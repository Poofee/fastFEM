//#include "plot.h"

//#include <QApplication>

//#include "FastFEMcore.h"

#include "relay1250.h"

#include "slu_mt_ddefs.h"

#include <armadillo> 

#include <stdio.h>

#include "meshGFace.h"
#include "Generator.h"
#include "GModel.h"
#include "gmsh.h"
#include "MLine.h"
#include "Field.h"

#define PI 3.14159265358979323846
using namespace arma;
//void triangletest();
//void triangletest1();
//void triangletestgroup();
//void quadtest();
//void armatest();
//void quadtlmtest();
//void triangletestvtm();
//void triangletestvtm2();
//void triangletestvtm3();
//void quadvtmtest();
//void T3NRtest();
//void T3NRTLMtest();
//void triangletestvtmsingle();
void test2DRemesh();

void test1250time(){
    Relay1250 relay;
    relay.run();
}

void armatest(){
	mat A (3, 3);
	colvec b (3);
	A(0, 0) = 3; A(0, 1) = -1; A(0, 2) = -2;
	A(1, 0) = -1; A(1, 1) = 4; A(1, 2) = -3;
	A(2, 0) = -2; A(2, 1) = -3; A(2, 2) = 5;

	b(0) = 4; b(1) = -5; b(2) = 1;
	colvec x2;//5 3 4
	bool status = solve(x2, A, b);
//	qDebug() << x2(0) << x2(1) << x2(2);
}
//void quadtest(){
//    CFastFEMcore fem;
//    //读取工程文件
//    fem.openProject("..\\model\\project1.mag");
//    //读取分网
//    if(fem.LoadQ4MeshCOMSOL("..\\model\\reg1.mphtxt") == 0){
////        qDebug()<<"OK";
////        qDebug()<<"number of elements:"<<fem.num_ele;
////        qDebug()<<"number of points:"<<fem.num_pts;

//		fem.StaticAxisQ4NR();
//    }
//}
//void quadtlmtest(){
//	CFastFEMcore fem;
//	//读取工程文件
//	fem.openProject("..\\model\\project1.mag");
//	//读取分网
//	if (fem.LoadQ4MeshCOMSOL("..\\model\\reg1.mphtxt") == 0){
////		qDebug() << "OK";
////		qDebug() << "number of elements:" << fem.num_ele;
////		qDebug() << "number of points:" << fem.num_pts;
//		//先进行一次牛顿求解
//		fem.StaticAxisQ4NR();
//		//设置一个猜测值
//		for (int i = 0; i < fem.num_pts; i++){
//			fem.pmeshnode[i].A *= fem.pmeshnode[i].x;
//		}
//		fem.StaticAxisQ4TLM();
//	}
//}

//void quadvtmtest(){
//	CFastFEMcore fem;
//	//读取工程文件
//	fem.openProject("..\\model\\project1.mag");
//	//读取分网
//	if (fem.LoadQ4MeshCOMSOL("..\\model\\reg1.mphtxt") == 0){
////		qDebug() << "OK";
////		qDebug() << "number of elements:" << fem.num_ele;
////		qDebug() << "number of points:" << fem.num_pts;
//		//先进行一次牛顿求解
//		fem.StaticAxisQ4NR();
//		//设置一个猜测值
//		for (int i = 0; i < fem.num_pts; i++){
//			fem.pmeshnode[i].A *= 0.5;
//		}
//		fem.StaticAxisQ4VTM();
//	}
//}

//void T3NRtest() {
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\mesh04.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miut = 0.9*fem.pmeshele[i].miut;
//			}
//		}
//		//t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		//t1 = SuperLU_timer_() - t1;
//		//qDebug() << "NR:" << t1;
//		//fem.CalcForce();
//	}
//}
//void triangletestgroup() {
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\mesh46k.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miu = 1*fem.pmeshele[i].miut;
//			}
//		}
//		t1 = SuperLU_timer_();
//		fem.StaticAxisT3TLMgroup();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "TLM:" << t1;
//		//fem.CalcForce();
//	}
//}
//void triangletestvtm() {
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\ruijiao.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miu = 1*fem.pmeshele[i].miut;
//			}
//		}
//		t1 = SuperLU_timer_();
//		fem.StaticAxisT3VTM();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "TLM:" << t1;
//		//fem.CalcForce();
//	}
//}
//void triangletestvtmsingle() {
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\ruijiao.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miu = 1 * fem.pmeshele[i].miut;
//			}
//		}
//		t1 = SuperLU_timer_();
//		fem.StaticAxisT3VTMsingle();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "TLM:" << t1;
//		//fem.CalcForce();
//	}
//}
//void triangletestvtm2() {
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miu = 0.9 * fem.pmeshele[i].miut;
//			}
//		}
//		t1 = SuperLU_timer_();
//		fem.StaticAxisT3VTM2();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "TLM:" << t1;
//		//fem.CalcForce();
//	}
//}

//void triangletestvtm3() {
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miu = 1 * fem.pmeshele[i].miut;
//			}
//		}
//		t1 = SuperLU_timer_();
//		fem.StaticAxisT3VTM3();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "TLM:" << t1;
//		//fem.CalcForce();
//	}
//}
//void triangletest1() {
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\mesh26k.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miu = 0.9*fem.pmeshele[i].miut;
//			}
//		}
//		t1 = SuperLU_timer_();
//		//fem.StaticAxisTLMNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "TLM:" << t1;
//		//fem.CalcForce();
//	}
//}
//void triangletest(){
//    CFastFEMcore fem;
//    double t1;
////    qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////    qDebug() << "current currentPath: " << QDir::currentPath();

//    fem.openProject("..\\model\\project1.mag");

//    if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
//        fem.preCalculation();
//        t1 = SuperLU_timer_();
//        fem.StaticAxisymmetricNR();
//        t1 = SuperLU_timer_() - t1;
////        qDebug() << "NR:" << t1;
//        for (int i = 0; i < fem.num_ele; i++) {
//            if (!fem.pmeshele[i].LinearFlag) {
//                fem.pmeshele[i].miu = 0.9*fem.pmeshele[i].miut;
//            }
//        }
//        t1 = SuperLU_timer_();
//        fem.StaticAxisymmetricTLM();
//        t1 = SuperLU_timer_() - t1;
////        qDebug() << "TLM:" << t1;
//        //fem.CalcForce();
//    }
//}
//void T3NRTLMtest(){
//	CFastFEMcore fem;
//	double t1;
////	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
////	qDebug() << "current currentPath: " << QDir::currentPath();

//	fem.openProject("..\\model\\project1.mag");

//	if (fem.Load2DMeshCOMSOL("..\\model\\mesh04.mphtxt") == 0) {
//		fem.preCalculation();
//		t1 = SuperLU_timer_();
//		fem.StaticAxisymmetricNR();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		for (int i = 0; i < fem.num_ele; i++) {
//			if (!fem.pmeshele[i].LinearFlag) {
//				fem.pmeshele[i].miut = 1*fem.pmeshele[i].miut;
//			}
//			fem.pmeshele[i].B = 0;
//		}
//		t1 = SuperLU_timer_();
//		fem.StaticAxisT3NRTLM();
//		t1 = SuperLU_timer_() - t1;
////		qDebug() << "NR:" << t1;
//		//fem.CalcForce();
//	}
//}



/*!
 \brief

 \param f
 \param dx
 \param dy
 \param dz
*/
void moveFace(GFace *f,double dx,double dy,double dz){
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
 \brief

 \param f
*/
void deleteFaceMesh(GFace* f){
    deMeshGFace dem;
    dem(f);
}

void test2DRemesh(){
    int myargn = 6;
    char fileName[] = "JRS1250/JRS1250bgm";
    char geoName[100];
    char nameMsh[100];
    sprintf(geoName,"%s.geo",fileName);
    char *myargv[] = {(char*)"gmsh",(char*)"-setnumber",(char*)"disp",(char*)"1",(char*)"-format",(char*)"msh2",(char*)"-v",(char*)"1000"};
    gmsh::initialize(myargn,myargv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(geoName);
    //    gmsh::model::mesh::generate(2);
    //    GModel::current()->mesh(2);
    GModel* model = GModel::current();
    GenerateMesh(model, 2);
    sprintf(nameMsh,"%s_%02d.msh",fileName,0);
    gmsh::write(nameMsh);
    /** 删掉空气的分网 **/
    int tag_xiantie = 1;
    int tag_air = 10;
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
    double xdisp = 0;
    double ydisp = 0.1;
    double xingcheng = 2.2;

    for(int i = 1; i <= 20;i++){
        printf("-------------mesh step %d-------------\n",i);
        /** 删除空气的分网 **/
        printf("-------------deleting air surface %d-------------\n",tag_air);
        deleteFaceMesh(f_air);
        /** 移动衔铁的分网 **/
        printf("-------------moving armature mesh...-------------\n");
        moveFace(f_xiantie,xdisp,ydisp,0);

        /** 对空气进行重分网 **/
        printf("-------------remesh air domain...-------------\n");
        f_air->mesh(true);

        sprintf(nameMsh,"%s_%02d.msh",fileName,i);
        gmsh::write(nameMsh);
    }

    gmsh::finalize();
    printf("-------------Finish-------------\n");
}

int main(int argc, char *argv[])
{
//    QApplication a(argc, argv);
//	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
//	qDebug() << "current currentPath: " << QDir::currentPath();

    //armatest();
    //triangletest1();
    //quadtlmtest();
    //triangletestvtm2();
    //triangletestvtm();
    //triangletestvtmsingle();
//    T3NRtest();
    //quadvtmtest();
    //T3NRTLMtest();
    //triangletestgroup();
    //quadtest();
//    test2DRemesh();

//	Plot myplot;
//	myplot.show();


    test1250time();

    return 0;
}
