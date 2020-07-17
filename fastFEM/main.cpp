// Copyright 2020 Poofee (https://github.com/Poofee)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------
/*****************************************************************************
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Poofee                                                          *
 *  Email:   poofee@qq.com                                                   *
 *  Address:                                                                 *
 *  Original Date: 2020-07-15                                                *
 *                                                                           *
 *****************************************************************************/
//#include "plot.h"

//#include <QApplication>

//#include "FastFEMcore.h"

#include "relay1250.h"
#include "testremesh.h"

#include "slu_mt_ddefs.h"
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
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


//    test1250time();

    testRemesh tr;
//    tr.testRemesh2D();
//    tr.testRemesh3DRotate();
//    tr.testRemesh3DParallel();
//    tr.testMove3DMeshParallel();
    tr.testRotate3DMesh();

    return 0;
}
