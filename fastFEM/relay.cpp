#include "relay.h"

Relay::Relay():
    current_step(0),
    pmeshnode(nullptr),
    pmeshele(nullptr),
    pmeshnodelast(nullptr),
    pmeshelelast(nullptr),
    model(nullptr)
{

}

Relay::~Relay()
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

void Relay::init()
{

}

void Relay::setFileName(char fn[])
{
    sprintf(fileName,"%s",fn);
}

/*!
 \brief 打开geo文件。注意，由于瞬态问题需要使用到上一次的
 数据，所以，最好在geo文件中，将发生形变的区域移到最后。

*/
void Relay::openGeo()
{
    int myargn = 4;
    char geoName[256];
    sprintf(geoName,"%s.geo",fileName);
    char *myargv[] = {(char*)"gmsh",(char*)"-format",(char*)"msh2",(char*)"-v",(char*)"1000"};
    gmsh::initialize(myargn,myargv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(geoName);
    //    gmsh::model::mesh::generate(2);
    //    GModel::current()->mesh(2);
    model = GModel::current();

    /** 初次分网，重分网必须先有一个网格 **/
    char nameMsh[256];
    GenerateMesh(model, 3);
    sprintf(nameMsh,"%s_%02d.msh",fileName,0);
    gmsh::write(nameMsh);
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
void Relay::loadMesh()
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
                    case 3:/** quadrangle **/
                        element_pmeshnode = 4;
                        break;
                    case 4:/** tetrahedron **/
                        element_pmeshnode = 4;
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
void Relay::findBoundaryPoints(int index)
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
    new_end = unique(it_1,it_2);
    boundaryPoints.erase(new_end,it_2);
    printf("Total %llu boundary points.\n",boundaryPoints.size());
}

/*!
 \brief 查找外边界上的边

 \param index 查找区域index，如果未指定，则是全局。
*/
void Relay::findBoundaryEdges(int index)
{
    /** 计算出三角单元的数目 **/
    num_triangle = 0;
    for(int i=0; i < num_ele;i++){
        if(pmeshele[i].ele_type == 2){
            num_triangle++;
        }
    }
    /** 生成所有的棱 **/
    boundaryEdges.resize(num_triangle * 3);
    allPoints.resize(num_triangle * 3);

    for(int i=0; i < num_triangle;i++){
        int a,b,c,t;
        int i1 = num_ele-num_triangle+i;
        a = pmeshele[i1].n[0];
        b = pmeshele[i1].n[1];
        c = pmeshele[i1].n[2];
        /** 按照升序排序 **/
        if(a>b)  {t=a;a=b;b=t;}
        if(a>c)  {t=a;a=c;c=t;}
        if(b>c)  {t=b;b=c;c=t;}

        boundaryEdges.at(i*3 + 0).start = a;
        boundaryEdges.at(i*3 + 0).end   = b;
        boundaryEdges.at(i*3 + 1).start = a;
        boundaryEdges.at(i*3 + 1).end   = c;
        boundaryEdges.at(i*3 + 2).start = b;
        boundaryEdges.at(i*3 + 2).end   = c;

        allPoints.at(i*3 + 0) = a;
        allPoints.at(i*3 + 1) = b;
        allPoints.at(i*3 + 2) = c;
    }

    /** 去重找到边界 **/
    std::vector<FEMedge>::iterator it_1 = boundaryEdges.begin();
    std::vector<FEMedge>::iterator it_2 = boundaryEdges.end();
    std::vector<FEMedge>::iterator new_end;

    std::sort(it_1,it_2);

    new_end = Myunique(it_1,it_2);
    allEdges = boundaryEdges;
    boundaryEdges.erase(new_end,it_2);
//    for(int k = 0;k < boundaryEdges.size();k++){
//        printf("%d: start: %d,end:%d\n",k,boundaryEdges.at(k).start,boundaryEdges.at(k).end);
//    }
    /** 去重找到所有的节点 **/
    std::vector<int>::iterator it_3 = allPoints.begin();
    std::vector<int>::iterator it_4 = allPoints.end();
    std::vector<int>::iterator new_end1;

    std::sort(it_3,it_4);

    new_end1 = std::unique(it_3,it_4);
    allPoints.erase(new_end1,it_4);

    printf("Total %d triangles.\n",num_triangle);
    printf("Total %llu points.\n",allPoints.size());
    printf("Total %llu boundary edges.\n",boundaryEdges.size());
}

/*!
 \brief 三维版的查找边界棱。计算出三维分网棱的连接关系。
 三维的话，要计算边界棱的话，还得先计算出边界面，感觉还是
 应该计算出连接关系。

 \param index
*/
void Relay::findBoundaryEdges3D(int index)
{
    /** 计算出四面体单元的数目 **/
    num_tet = 0;
    for(int i=0; i < num_ele;i++){
        if(pmeshele[i].ele_type == 4){
            num_tet++;
        }
    }
    /** 生成所有的棱 **/
    boundaryEdges.resize(num_tet * 4);
    allPoints.resize(num_tet * 4);

    for(int i=0; i < num_tet;i++){
        int a,b,c,d,t;
        int i1 = num_ele-num_tet+i;
        a = pmeshele[i1].n[0];
        b = pmeshele[i1].n[1];
        c = pmeshele[i1].n[2];
        d = pmeshele[i1].n[3];
        /** 按照升序排序 **/
        if(a>b)  {t=a;a=b;b=t;}
        if(a>c)  {t=a;a=c;c=t;}
        if(a>d)  {t=a;a=d;d=t;}
        if(b>c)  {t=b;b=c;c=t;}
        if(b>d)  {t=b;b=d;d=t;}
        if(c>d)  {t=c;c=d;d=t;}

        boundaryEdges.at(i*6 + 0).start = a;
        boundaryEdges.at(i*6 + 0).end   = b;
        boundaryEdges.at(i*6 + 1).start = a;
        boundaryEdges.at(i*6 + 1).end   = c;
        boundaryEdges.at(i*6 + 2).start = a;
        boundaryEdges.at(i*6 + 2).end   = d;
        boundaryEdges.at(i*6 + 3).start = b;
        boundaryEdges.at(i*6 + 3).end   = c;
        boundaryEdges.at(i*6 + 4).start = b;
        boundaryEdges.at(i*6 + 4).end   = d;
        boundaryEdges.at(i*6 + 5).start = c;
        boundaryEdges.at(i*6 + 5).end   = d;

        allPoints.at(i*4 + 0) = a;
        allPoints.at(i*4 + 1) = b;
        allPoints.at(i*4 + 2) = c;
        allPoints.at(i*4 + 3) = d;
    }

    /** 去重找到边界 **/
    std::vector<FEMedge>::iterator it_1 = boundaryEdges.begin();
    std::vector<FEMedge>::iterator it_2 = boundaryEdges.end();
    std::vector<FEMedge>::iterator new_end;

    std::sort(it_1,it_2);

    new_end = Myunique(it_1,it_2);
    allEdges = boundaryEdges;
    boundaryEdges.erase(new_end,it_2);
//    for(int k = 0;k < boundaryEdges.size();k++){
//        printf("%d: start: %d,end:%d\n",k,boundaryEdges.at(k).start,boundaryEdges.at(k).end);
//    }
    /** 去重找到所有的节点 **/
    std::vector<int>::iterator it_3 = allPoints.begin();
    std::vector<int>::iterator it_4 = allPoints.end();
    std::vector<int>::iterator new_end1;

    std::sort(it_3,it_4);

    new_end1 = std::unique(it_3,it_4);
    allPoints.erase(new_end1,it_4);

    printf("Total %d triangles.\n",num_triangle);
    printf("Total %llu points.\n",allPoints.size());
    printf("Total %llu boundary edges.\n",boundaryEdges.size());
}

/*!
 \brief 查找三维体的边界上的面单元。

 \param index 为-1时，计算整体的边界。
*/
void Relay::findBoundaryFaces(int index)
{
    /** 计算出四面体单元的数目 **/
    num_tet = 0;
    for(int i=0; i < num_ele;i++){
        if(pmeshele[i].ele_type == 4){
            num_tet++;
        }
    }
    /** 生成所有的面 **/
    boundaryFaces.resize(num_tet * 6);
    allPoints.resize(num_tet * 4);

    for(unsigned int i=0; i < num_tet;i++){
        int a,b,c,d,t;
        int i1 = num_ele-num_tet+i;
        a = pmeshele[i1].n[0];
        b = pmeshele[i1].n[1];
        c = pmeshele[i1].n[2];
        d = pmeshele[i1].n[3];
        /** 按照升序排序 **/
        if(a>b)  {t=a;a=b;b=t;}
        if(a>c)  {t=a;a=c;c=t;}
        if(a>d)  {t=a;a=d;d=t;}
        if(b>c)  {t=b;b=c;c=t;}
        if(b>d)  {t=b;b=d;d=t;}
        if(c>d)  {t=c;c=d;d=t;}

        boundaryEdges.at(i*6 + 0).start = a;
        boundaryEdges.at(i*6 + 0).end   = b;
        boundaryEdges.at(i*6 + 1).start = a;
        boundaryEdges.at(i*6 + 1).end   = c;
        boundaryEdges.at(i*6 + 2).start = a;
        boundaryEdges.at(i*6 + 2).end   = d;
        boundaryEdges.at(i*6 + 3).start = b;
        boundaryEdges.at(i*6 + 3).end   = c;
        boundaryEdges.at(i*6 + 4).start = b;
        boundaryEdges.at(i*6 + 4).end   = d;
        boundaryEdges.at(i*6 + 5).start = c;
        boundaryEdges.at(i*6 + 5).end   = d;

        allPoints.at(i*4 + 0) = a;
        allPoints.at(i*4 + 1) = b;
        allPoints.at(i*4 + 2) = c;
        allPoints.at(i*4 + 3) = d;
    }

    /** 去重找到边界 **/
    std::vector<FEMedge>::iterator it_1 = boundaryEdges.begin();
    std::vector<FEMedge>::iterator it_2 = boundaryEdges.end();
    std::vector<FEMedge>::iterator new_end;

    std::sort(it_1,it_2);

    new_end = Myunique(it_1,it_2);
    allEdges = boundaryEdges;
    boundaryEdges.erase(new_end,it_2);
//    for(int k = 0;k < boundaryEdges.size();k++){
//        printf("%d: start: %d,end:%d\n",k,boundaryEdges.at(k).start,boundaryEdges.at(k).end);
//    }
    /** 去重找到所有的节点 **/
    std::vector<int>::iterator it_3 = allPoints.begin();
    std::vector<int>::iterator it_4 = allPoints.end();
    std::vector<int>::iterator new_end1;

    std::sort(it_3,it_4);

    new_end1 = std::unique(it_3,it_4);
    allPoints.erase(new_end1,it_4);

    printf("Total %d triangles.\n",num_triangle);
    printf("Total %llu points.\n",allPoints.size());
    printf("Total %llu boundary edges.\n",boundaryEdges.size());
}

void Relay::setXiantieTag(int xiantie)
{
    tag_xiantie = xiantie;
}

void Relay::setAirTag(int air)
{
    tag_air = air;
}

/*!
 \brief

 \param f
*/
void Relay::deleteFaceMesh(GFace *f)
{
    deMeshGFace dem;
    dem(f);
}

/**
 * @brief 删除区域内的网格。应该是内部的网格。
 *
 * @param r
 */
void Relay::deleteRegionMesh(GRegion* r)
{
    deMeshGRegion dem;
    dem(r);
}

/*!
 \brief 位移为正值，验证轴的正方向运动，反之，为负。

 \param f
 \param dx
 \param dy
 \param dz
*/
void Relay::moveFace(GFace *f, double dx, double dy, double dz)
{
    /** 修改内部分网节点坐标 **/
    for(std::size_t j = 0; j < f->getNumMeshVertices(); j++) {
        f->getMeshVertex(j)->setXYZ(f->getMeshVertex(j)->x()+dx,
                                    f->getMeshVertex(j)->y()+dy,
                                    f->getMeshVertex(j)->z()+dz);
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
        (*ite)->setXYZ((*ite)->x()+dx,(*ite)->y()+dy,(*ite)->z()+dz);
    }
    /** 修改几何顶点坐标 **/
    for(auto v : f->vertices()){
        for(std::size_t j = 0; j < v->getNumMeshVertices(); j++){
            GPoint p(v->x()+dx,v->y()+dy,v->z()+dz);
            v->setPosition(p);
        }
    }
    /** 更新distance field **/
    updateField();
}

/**
 * @brief 移动某个区域内的网格。
 *
 * @param r
 * @param dx
 * @param dy
 * @param dz
 */
void Relay::moveRegion(GRegion* r, double dx, double dy, double dz)
{
//    printf("r->getNumMeshVertices():%llu\n",r->getNumMeshVertices());
//    printf("r->embeddedVertices():%llu\n",r->embeddedVertices().size());
//    printf("r->vertices():%llu\n",r->vertices().size());

    /** 移动内部的节点 **/
    for(std::size_t j = 0;j <r->getNumMeshVertices();j++){
        r->getMeshVertex(j)->setXYZ(r->getMeshVertex(j)->x() + dx,
                                    r->getMeshVertex(j)->y() + dy,
                                    r->getMeshVertex(j)->z() + dz);
    }

    /** 修改边界面上分网节点坐标 **/
    std::set<MVertex *, MVertexLessThanNum> all_vertices;
    for(auto f : r->faces()){
        /** 修改内部分网节点坐标 **/
        for(std::size_t j = 0; j < f->getNumMeshVertices(); j++) {
            f->getMeshVertex(j)->setXYZ(f->getMeshVertex(j)->x()+dx,
                                        f->getMeshVertex(j)->y()+dy,
                                        f->getMeshVertex(j)->z()+dz);
        }
        /** 由于两个面相邻的话，分别移动边界上的点会造成重复移动 **/
        for(auto e : f->edges()){
            for(auto line : e->lines){
                MVertex *v1 = line->getVertex(0);
                MVertex *v2 = line->getVertex(1);

                all_vertices.insert(v1);
                all_vertices.insert(v2);
            }
        }
    }
    /** 移动面的边上的分网点**/
    for(std::set<MVertex *, MVertexLessThanNum>::iterator ite = all_vertices.begin();ite != all_vertices.end(); ite++) {
        (*ite)->setXYZ((*ite)->x()+dx,(*ite)->y()+dy,(*ite)->z()+dz);
    }
    /** 修改几何顶点坐标 **/
    for(auto v : r->vertices()){
        for(std::size_t j = 0; j < v->getNumMeshVertices(); j++){
            GPoint p(v->x()+dx,v->y()+dy,v->z()+dz);
            v->setPosition(p);
        }
    }
    /** 更新distance field **/
    updateField();
}

/**
 * @brief 旋转某个区域内的网格。正方向为右手握住轴的方向。
 *
 * @param r
 * @param dangle -180~180
 * @param x1 向量的位置
 * @param y1
 * @param z1
 * @param v1 旋转单位向量，过原点。
 * @param v2
 * @param v3
 */
void Relay::rotateRegion(GRegion* r, double dangle, double x1, double y1, double z1, double v1, double v2, double v3)
{
    Point3f p;
    /** 移动内部的节点 **/
    for(std::size_t j = 0;j <r->getNumMeshVertices();j++){
        p = RotateByVector(dangle,r->getMeshVertex(j)->x()-x1,r->getMeshVertex(j)->y()-y1,r->getMeshVertex(j)->z()-z1,v1,v2,v3);
        r->getMeshVertex(j)->setXYZ(p.x+x1,p.y+y1,p.z+z1);
    }

    /** 修改边界面上分网节点坐标 **/
    std::set<MVertex *, MVertexLessThanNum> all_vertices;
    for(auto f : r->faces()){
        /** 修改内部分网节点坐标 **/
        for(std::size_t j = 0; j < f->getNumMeshVertices(); j++) {
            p = RotateByVector(dangle,f->getMeshVertex(j)->x()-x1,f->getMeshVertex(j)->y()-y1,f->getMeshVertex(j)->z()-z1,v1,v2,v3);
            f->getMeshVertex(j)->setXYZ(p.x+x1,p.y+y1,p.z+z1);
        }
        /** 由于两个面相邻的话，分别移动边界上的点会造成重复移动 **/
        for(auto e : f->edges()){
            for(auto line : e->lines){
                MVertex *v1 = line->getVertex(0);
                MVertex *v2 = line->getVertex(1);

                all_vertices.insert(v1);
                all_vertices.insert(v2);
            }
        }
    }
    /** 移动面的边上的分网点**/
    for(std::set<MVertex *, MVertexLessThanNum>::iterator ite = all_vertices.begin();ite != all_vertices.end(); ite++) {
        p = RotateByVector(dangle,(*ite)->x()-x1,(*ite)->y()-y1,(*ite)->z()-z1,v1,v2,v3);
        (*ite)->setXYZ(p.x+x1,p.y+y1,p.z+z1);
    }
    /** 修改几何顶点坐标 **/
    for(auto v : r->vertices()){
        for(std::size_t j = 0; j < v->getNumMeshVertices(); j++){
            p = RotateByVector(dangle,v->x()-x1,v->y()-y1,v->z()-z1,v1,v2,v3);
            GPoint p1(p.x+x1,p.y+y1,p.z+z1);
            v->setPosition(p1);
        }
    }
    /** 更新distance field **/
    updateField();
}

/**
 * @brief 空间点绕任意轴的旋转公式。
 *
 * @param theta 角度
 * @param old_x
 * @param old_y
 * @param old_z
 * @param vx 单位向量。
 * @param vy
 * @param vz
 * @return Point3f
 */
Point3f Relay::RotateByVector(double theta, double old_x, double old_y, double old_z, double vx, double vy, double vz)
{
    double r = theta * M_PI / 180;
    double c = cos(r);
    double s = sin(r);

    double new_x = (vx*vx*(1 - c) + c) * old_x + (vx*vy*(1 - c) - vz*s) * old_y + (vx*vz*(1 - c) + vy*s) * old_z;
    double new_y = (vy*vx*(1 - c) + vz*s) * old_x + (vy*vy*(1 - c) + c) * old_y + (vy*vz*(1 - c) - vx*s) * old_z;
    double new_z = (vx*vz*(1 - c) - vy*s) * old_x + (vy*vz*(1 - c) + vx*s) * old_y + (vz*vz*(1 - c) + c) * old_z;
    return Point3f(new_x, new_y, new_z);
}

/*!
 \brief 对模型进行重分网，对衔铁进行移动，对可压缩
 区域进行重分网。

 \param dx
 \param dy
 \param dz
*/
void Relay::remesh(double dx, double dy)
{
    char nameMsh[256];
    sprintf(nameMsh,"%s_%02d.msh",fileName,current_step);
    /** 未发生位移就不要分了 **/
    if(dx == 0 && dy == 0){
        gmsh::write(nameMsh);
        return;
    }

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

    gmsh::write(nameMsh);

    printf("-------------Finish-------------\n");
}

/*!
 \brief 对3D分网进行平行移动。

 \param dx 移动的距离，默认大小与坐标轴一致。
 \param dy
 \param dz
*/
void Relay::remesh3DParallel(double dx, double dy, double dz)
{
    char nameMsh[256];
    sprintf(nameMsh,"%s_%02d.msh",fileName,current_step);
    /** 未发生位移就不要分了 **/
    if(dx == 0 && dy == 0 && dz == 0){
        gmsh::write(nameMsh);
        return;
    }

    GRegion* r_xiantie = nullptr;
    GRegion* r_air = nullptr;
    for(GModel::riter it = model->firstRegion(); it != model->lastRegion(); ++it){
        printf("rigion tag:%d\n",(*it)->tag());
        if((*it)->tag() == tag_air){
            r_air = (*it);
        }
        if((*it)->tag() == tag_xiantie){
            r_xiantie = (*it);
        }
    }

    printf("-------------mesh step %d-------------\n",current_step);
    /** 删除空气的分网 **/
    printf("-------------deleting air surface %d-------------\n",tag_air);
    deleteRegionMesh(r_air);
    /** 移动衔铁的分网 **/
    printf("-------------moving armature mesh...-------------\n");
    moveRegion(r_xiantie,dx,dy,dz);

    /** 对空气进行重分网 **/
    printf("-------------remesh air domain...-------------\n");
    r_air->mesh(true);

    gmsh::write(nameMsh);

    printf("-------------Finish-------------\n");
}

/*!
 \brief 对分网进行旋转移动。

 \param dangle 移动的角度，方向是右手握住轴的方向，即逆时针为正。
 \param x1 轴上的点x
 \param y1 轴上的点y
 \param z1 轴上的点z
 \param v1 轴向量x
 \param v2 轴向量y
 \param v3 轴向量z
*/
void Relay::remesh3DRotate(double dangle, double x1, double y1, double z1, double v1, double v2, double v3)
{
    char nameMsh[256];
    sprintf(nameMsh,"%s_%02d.msh",fileName,current_step);
    /** 未发生位移就不要分了 **/
    if(dangle == 0){
        gmsh::write(nameMsh);
        return;
    }

    GRegion* r_xiantie = nullptr;
    GRegion* r_air = nullptr;
    for(GModel::riter it = model->firstRegion(); it != model->lastRegion(); ++it){
        printf("region tag:%d\n",(*it)->tag());
        if((*it)->tag() == tag_air){
            r_air = (*it);
        }
        if((*it)->tag() == tag_xiantie){
            r_xiantie = (*it);
        }
    }


    printf("-------------mesh step %d-------------\n",current_step);
    /** 删除空气的分网 **/
    printf("-------------deleting air surface %d-------------\n",tag_air);
    deleteRegionMesh(r_air);
    /** 移动衔铁的分网 **/
    printf("-------------moving armature mesh...-------------\n");
    rotateRegion(r_xiantie,dangle,x1,y1,z1,v1,v2,v3);

    /** 对空气进行重分网 **/
    printf("-------------remesh air domain...-------------\n");
    r_air->mesh(true);

    gmsh::write(nameMsh);

    printf("-------------Finish-------------\n");
}

/**
 * @brief 更新gmsh的分网大小控制。
 *
 */
void Relay::updateField()
{
    /** 更新distance field **/
    FieldManager* fM = GModel::current()->getFields();
    std::map<int, Field *>::iterator iter;
    for(iter = fM->begin(); iter != fM->end(); iter++){
        iter->second->update_needed = true;
        iter->second->update();
    }
}

void Relay::run()
{

}

void Relay::stepIncrement()
{
    current_step++;
}

void Relay::calcMagForce(int index)
{

}

void Relay::calcMagTorque(int index)
{

}

void Relay::outputResults()
{

}

