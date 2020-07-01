#include "relay1250.h"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo>
#include <vector>
#include <ctime>
#include <omp.h>

#include "SuperLU_MT.h"
#include "qcustomplot.h"
#include "slu_mt_ddefs.h"


using namespace std;
using namespace arma;

#define PI 3.14159265358979323846
#define r 2

Relay1250::Relay1250():
    pmeshnode(nullptr),
    pmeshele(nullptr),
    pmeshnodelast(nullptr),
    pmeshelelast(nullptr),
    model(nullptr)
{

}

Relay1250::~Relay1250()
{
    gmsh::finalize();

    if(pmeshele) delete pmeshele;
    if(pmeshnode) delete pmeshnode;
    if(pmeshelelast) delete pmeshelelast;
    if(pmeshnodelast) delete pmeshnodelast;

    for(auto m : materialList){
        if(m) delete m;
    }
    materialList.clear();
    materialMap.clear();
}

/*!
 \brief 参数初始化

*/
void Relay1250::init()
{
    /** 文件名，需不带扩展名 **/
    sprintf(fileName,"%s","JRS1250/JRS1250bgm");
    /** 初始化时间步长 **/
    current_step = 0;
    min_time = 0;
    max_time = 4e-3;
    for(double t = 1e-5; t < 3e-3; t+=1e-5){
        timesteps.emplace_back(1e-5);
    }
    for(double t = 3e-3; t <= 3e-3; t+=5e-5){
        timesteps.emplace_back(5e-5);
    }

    /** 初始化变量 **/
    current_step = 0;
    Ddisplacements.resize(timesteps.size()+1);
    Ddisplacements.at(current_step) = 0;
    displacements.resize(timesteps.size()+1);
    displacements.at(current_step) = 0;
    velocities.resize(timesteps.size()+1);
    velocities.at(current_step) = 0;
    accelerations.resize(timesteps.size()+1);
    accelerations.at(current_step) = 0;
    PhiCoil.resize(timesteps.size()+1);
    PhiCoil.at(current_step) = 0;
    ICoil.resize(timesteps.size()+1);
    ICoil.at(current_step) = 0;
    magforces.resize(timesteps.size()+1);
    magforces.at(current_step) = 0;
    UCoil.resize(timesteps.size()+1);
    UCoil.at(current_step) = 420;

    MAX_NONLINEARSTEPS = 20;
    /** 设置材料参数 **/
    CMaterial* air_material = new CMaterial;
    materialList.emplace_back(air_material);

    CMaterial* pm_material = new CMaterial;
    materialList.emplace_back(pm_material);

    CMaterial* dt4e_material = new CMaterial;
    materialList.emplace_back(dt4e_material);

    CMaterial* coil_material = new CMaterial;
    materialList.emplace_back(coil_material);
    coil_material->tau = 60/(3e-3 * 17e-3);

    /** 设置映射关系 **/
    materialMap[1] = air_material;/** 衔铁 **/
    materialMap[2] = air_material;/**  **/
    materialMap[3] = air_material;/**  **/
    materialMap[4] = air_material;/**  **/
    materialMap[5] = air_material;/**  **/
    materialMap[6] = air_material;/**  **/
    materialMap[7] = air_material;/**  **/
    materialMap[8] = air_material;/**  **/
    materialMap[9] = air_material;/**  **/
    materialMap[10] = air_material;/**  **/
    /** 设置形变区域 **/
    tag_xiantie = 1;
    tag_air = 10;

    Precision = 1e-6;
}

/*!
 \brief 打开geo文件。注意，由于瞬态问题需要使用到上一次的
 数据，所以，最好在geo文件中，将发生形变的区域移到最后。

*/
void Relay1250::openGeo()
{
    int myargn = 6;
    char geoName[256];
    sprintf(geoName,"%s.geo",fileName);
    char *myargv[] = {(char*)"gmsh",(char*)"-setnumber",(char*)"disp",(char*)"1",(char*)"-format",(char*)"msh2",(char*)"-v",(char*)"1000"};
    gmsh::initialize(myargn,myargv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(geoName);
    //    gmsh::model::mesh::generate(2);
    //    GModel::current()->mesh(2);
    model = GModel::current();

//    /** 初次分网 **/
//    char nameMsh[256];
//    GenerateMesh(model, 2);
//    sprintf(nameMsh,"%s_%02d.msh",fileName,0);
//    gmsh::write(nameMsh);
}
int next_int(char **start)
{
    int i;
    char *end;

    i = strtol(*start,&end,10);
    *start = end;
    return(i);
}
/*!
 \brief 读取当前时间步的mesh，如果当前含有分网数据，先删除。
 gmsh网格的版本为22。

*/
void Relay1250::loadMesh()
{
    /** 先删除上一次的网格 **/
    if(pmeshelelast) delete pmeshelelast;
    if(pmeshnodelast) delete pmeshnodelast;
    /** 保存上一次的网格 **/
    if(pmeshele) pmeshelelast = pmeshele;
    if(pmeshnode) pmeshnodelast = pmeshnode;

    char nameMsh[256];
    sprintf(nameMsh,"%s_%02d.msh",fileName,current_step);
    char *ch = (char *)calloc(256,sizeof (char));
    //------------open file----------------------------------
    FILE * fp = nullptr;
    fp = fopen(nameMsh, "r");
    if (fp == nullptr) {
        printf("Error: openning file!");
        return ;
    }
    while(!feof(fp)){
        fgets(ch, 256, fp);

        if(strstr(ch,"$MeshFormat")){
            double version;
            int file_type;
            int data_size;
            if(fscanf(fp,"%lf %d %d\n",&version,&file_type,&data_size) != 3){
                printf("error reading format");
                return ;
            }else{
                if(version > 2.2){
                    printf("Can only open gmsh version 2.2 format");
                    return ;
                }
            }
            fgets(ch, 256, fp);
            if(!strstr(ch,"$EndMeshFormat")) {
                printf("$MeshFormat section should end to string $EndMeshFormat:\n%s\n",ch);
            }
        }else if(strstr(ch,"$Nodes")){
            if(fscanf(fp,"%d\n",&num_pts) != 1)
            {
                return ;
            }else{
                /** 读取节点坐标 **/
                pmeshnode = new CNode[num_pts];
                int index;
                for(int i = 0;i < num_pts;++i){
                    fscanf(fp,"%d %lf %lf %lf\n",&index,&pmeshnode[i].x,&pmeshnode[i].y,&pmeshnode[i].z);
                    //qDebug()<<index<<pmeshnode[i].x<<pmeshnode[i].y<<pmeshnode[i].z;
                }
            }
            fgets(ch, 256, fp);
            if(!strstr(ch,"$EndNodes")) {
                printf("$Node section should end to string $EndNodes:\n%s\n",ch);
            }
        }else if(strstr(ch,"$Elements")){
            int ele_number;
            //    int elm_type;
            int number_of_tags;
            char * chtmp;
            if(fscanf(fp,"%d\n",&num_ele) != 1){
                return ;
            }else{
                pmeshele = new CElement[num_ele];
                for(int i = 0;i < num_ele;++i){
                    chtmp = fgets(ch, 256, fp);
                    ele_number = next_int(&chtmp);
                    pmeshele[i].ele_type = next_int(&chtmp);
                    number_of_tags = next_int(&chtmp);
                    pmeshele[i].physic_tag = next_int(&chtmp);
                    pmeshele[i].geometry_tag = next_int(&chtmp);

                    int element_pmeshnode = 0;
                    switch (pmeshele[i].ele_type) {
                    case 15:/** point **/
                        element_pmeshnode = 1;
                        break;
                    case 1:/** line **/
                        element_pmeshnode = 2;
                        break;
                    case 2:/** triangle **/
                        element_pmeshnode = 3;
                        break;
                    default:
                        element_pmeshnode = 0;
                        break;
                    }

                    for(int j = 0; j < element_pmeshnode;++j)
                        pmeshele[i].n[j] = next_int(&chtmp)-1;
                }
            }
            fgets(ch, 256, fp);
            if(!strstr(ch,"$EndElements")) {
                printf("$Element section should end to string $EndElements:\n%s\n",ch);
            }
        }
    }
    delete ch;
    fclose(fp);
}

/*!
 \brief 默认查找外部边界，确定dof。

*/
void Relay1250::findBoundaryPoints(int index)
{
    /** 根据边界生成所有的点 **/
    int num_bndpoints = boundaryEdges.size() * 2;
    boundaryPoints.resize(num_bndpoints);
    for(int i = 0; i < boundaryEdges.size();i++){
        boundaryPoints.at(i*2 + 0) = boundaryEdges.at(i).start;
        boundaryPoints.at(i*2 + 1) = boundaryEdges.at(i).end;
    }
    /** 去重 **/
    std::vector<int>::iterator it_1 = boundaryPoints.begin();
    std::vector<int>::iterator it_2 = boundaryPoints.end();
    std::vector<int>::iterator new_end;

    sort(it_1,it_2);
    new_end = Myunique(it_1,it_2);
    boundaryPoints.erase(new_end,it_2);
}


/*!
 \brief 查找外边界上的边

 \param index 查找区域index，如果未指定，则是全局。
*/
void Relay1250::findBoundaryEdges(int index)
{
    /** 生成所有的棱 **/
    int numEdges = num_ele * 3;

    allEdges.resize(numEdges);
    allPoints.resize(num_ele * 3);

    for(int i=0; i < num_ele;i++){
        int a,b,c,t;
        a = pmeshele[i].n[0];
        b = pmeshele[i].n[1];
        c = pmeshele[i].n[2];
        /** 按照升序排序 **/
        if(a>b)  {t=a;a=b;b=t;}
        if(a>c)  {t=a;a=c;c=t;}
        if(b>c)  {t=b;b=c;c=t;}

        allEdges.at(i*3 + 0).start = a;
        allEdges.at(i*3 + 0).end   = b;
        allEdges.at(i*3 + 1).start = a;
        allEdges.at(i*3 + 1).end   = c;
        allEdges.at(i*3 + 2).start = b;
        allEdges.at(i*3 + 2).end   = c;

        allPoints.at(0) = a;
        allPoints.at(1) = b;
        allPoints.at(2) = c;
    }
    /** 去重找到边界 **/
    std::vector<FEMedge>::iterator it_1 = allEdges.begin();
    std::vector<FEMedge>::iterator it_2 = allEdges.end();
    std::vector<FEMedge>::iterator new_end;

    std::sort(it_1,it_2);

    new_end = Myunique(it_1,it_2);
    boundaryEdges = allEdges;
    boundaryEdges.erase(new_end,it_2);
    /** 去重找到所有的节点 **/
    std::vector<int>::iterator it_3 = allPoints.begin();
    std::vector<int>::iterator it_4 = allPoints.end();
    std::vector<int>::iterator new_end1;

    std::sort(it_3,it_4);

    new_end1 = std::unique(it_3,it_4);
    allPoints.erase(new_end1,it_4);
}

/*!
 \brief

 \param f
*/
void Relay1250::deleteFaceMesh(GFace *f)
{
    deMeshGFace dem;
    dem(f);
}

/*!
 \brief

 \param f
 \param dx
 \param dy
 \param dz
*/
void Relay1250::moveFace(GFace *f, double dx, double dy, double dz)
{
    /** 修改内部分网节点坐标 **/
    for(std::size_t j = 0; j < f->getNumMeshVertices(); j++) {
        f->getMeshVertex(j)->setXYZ(f->getMeshVertex(j)->x()-dx,
                                    f->getMeshVertex(j)->y()-dy,
                                    f->getMeshVertex(j)->z()-dz);
    }
    /** 修改边界棱上分网节点坐标 **/
    std::set<MVertex *, MVertexLessThanNum> all_vertices;
    for(auto e : f->edges()){
        for(auto line : e->lines){
            MVertex *v1 = line->getVertex(0);
            MVertex *v2 = line->getVertex(1);

            all_vertices.insert(v1);
            all_vertices.insert(v2);
        }
    }
    /** all_vertices比e->getMeshVertex多了顶点 **/
    for(std::set<MVertex *, MVertexLessThanNum>::iterator ite = all_vertices.begin();ite != all_vertices.end(); ite++) {
        (*ite)->setXYZ((*ite)->x()-dx,(*ite)->y()-dy,(*ite)->z()-dz);
    }
    /** 修改几何顶点坐标 **/
    for(auto v : f->vertices()){
        for(std::size_t j = 0; j < v->getNumMeshVertices(); j++){
            GPoint p(v->x()-dx,v->y()-dy,v->z()-dz);
            v->setPosition(p);
        }
    }
    /** 更新distance field **/
    FieldManager* fM = GModel::current()->getFields();
    std::map<int, Field *>::iterator iter;
    for(iter = fM->begin(); iter != fM->end(); iter++){
        iter->second->update_needed = true;
        iter->second->update();
    }
}

/*!
 \brief 对模型进行重分网，对衔铁进行移动，对可压缩
 区域进行重分网。

 \param dx
 \param dy
 \param dz
*/
void Relay1250::remesh(double dx, double dy)
{
    char nameMsh[256];
    /** 删掉空气的分网 **/
    GFace* f_xiantie = nullptr;
    GFace* f_air = nullptr;
    for(GModel::fiter it = model->firstFace(); it != model->lastFace(); ++it){
        if((*it)->tag() == tag_air){
            f_air = (*it);
        }
        if((*it)->tag() == tag_xiantie){
            f_xiantie = (*it);
        }
    }

    printf("-------------mesh step %d-------------\n",current_step);
    /** 删除空气的分网 **/
    printf("-------------deleting air surface %d-------------\n",tag_air);
    deleteFaceMesh(f_air);
    /** 移动衔铁的分网 **/
    printf("-------------moving armature mesh...-------------\n");
    moveFace(f_xiantie,dx,dy,0);

    /** 对空气进行重分网 **/
    printf("-------------remesh air domain...-------------\n");
    f_air->mesh(true);

    sprintf(nameMsh,"%s_%02d.msh",fileName,current_step);
    gmsh::write(nameMsh);

    printf("-------------Finish-------------\n");
}

/*!
 \brief 使用传统的牛顿法进行求解

*/
void Relay1250::AxisNRsolve()
{
    /** 计时 **/
    clock_t time[10];
    int tt = 0;
    time[tt++] = SuperLU_timer_();
    /** 变量 **/

    /** 时间步循环 **/
    current_step = 0;
    double totalTime = 0;
    for(int i = 0; i < timesteps.size();i++){
        totalTime += timesteps.at(i);
        printf("time step %d, current time is %lf\n",i,totalTime);
        /** remesh，要使用增量位移，向下为负，向上为正 **/
        remesh(0,Ddisplacements.at(i));
        /** 读取第i步的分网 **/
        loadMesh();
        /** 查找边界，默认一次边界 **/
        findBoundaryEdges(-1);
        findBoundaryPoints(-1);

        /** 非线性迭代求解 **/
        for(int non_iter = 0;non_iter < MAX_NONLINEARSTEPS;non_iter++){

            /** 有限元装配 **/
            for(int i_ele = 0; i_ele < num_ele; i_ele++){

            }

        }


        /** 更新结果 **/
        current_step += 1;
        Ddisplacements.at(current_step) = 0;
        displacements.at(current_step) = 0;
        velocities.at(current_step) = 0;
        accelerations.at(current_step) = 0;
        PhiCoil.at(current_step) = 0;
        ICoil.at(current_step) = 0;
        magforces.at(current_step) = 0;
        UCoil.at(current_step) = 0;
    }


}

/*!
 \brief 使用传输线法进行求解

*/
void Relay1250::AxisTLMsolve()
{

}

/*!
 \brief 运行测试案例

*/
void Relay1250::run()
{
    init();
    openGeo();

    AxisNRsolve();

    /** 输出信息 **/
    outputResults();
}

/*!
 \brief 将一些结果进行输出。

*/
void Relay1250::outputResults()
{

}
