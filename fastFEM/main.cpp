#include "plot.h"
#include <QApplication>
#include "FastFEMcore.h"
#include "slu_mt_ddefs.h"
#include <armadillo> 

#define PI 3.14159265358979323846
using namespace arma;
void triangletest();
void triangletest1();
void triangletestgroup();
void quadtest();
void armatest();
void quadtlmtest();
void triangletestvtm();
void triangletestvtm2();
void triangletestvtm3();
void quadvtmtest();
void T3NRtest();

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	//armatest();
	//triangletest1();
	//quadtlmtest();
	//triangletestvtm2();
	triangletestvtm();
	T3NRtest();
	//quadvtmtest();
	Plot myplot;
	myplot.show();
    return a.exec();
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
	qDebug() << x2(0) << x2(1) << x2(2);
}
void quadtest(){
    CFastFEMcore fem;
    //读取工程文件
    fem.openProject("..\\model\\project1.mag");
    //读取分网
    if(fem.LoadQ4MeshCOMSOL("..\\model\\reg1.mphtxt") == 0){
        qDebug()<<"OK";
        qDebug()<<"number of elements:"<<fem.num_ele;
        qDebug()<<"number of points:"<<fem.num_pts;

		fem.StaticAxisQ4NR();
    }
}
void quadtlmtest(){
	CFastFEMcore fem;
	//读取工程文件
	fem.openProject("..\\model\\project1.mag");
	//读取分网
	if (fem.LoadQ4MeshCOMSOL("..\\model\\reg1.mphtxt") == 0){
		qDebug() << "OK";
		qDebug() << "number of elements:" << fem.num_ele;
		qDebug() << "number of points:" << fem.num_pts;
		//先进行一次牛顿求解
		fem.StaticAxisQ4NR();
		//设置一个猜测值
		for (int i = 0; i < fem.num_pts; i++){
			fem.pmeshnode[i].A *= 0.9;
		}
		fem.StaticAxisQ4TLM();
	}
}

void quadvtmtest(){
	CFastFEMcore fem;
	//读取工程文件
	fem.openProject("..\\model\\project1.mag");
	//读取分网
	if (fem.LoadQ4MeshCOMSOL("..\\model\\reg1.mphtxt") == 0){
		qDebug() << "OK";
		qDebug() << "number of elements:" << fem.num_ele;
		qDebug() << "number of points:" << fem.num_pts;
		//先进行一次牛顿求解
		fem.StaticAxisQ4NR();
		//设置一个猜测值
		for (int i = 0; i < fem.num_pts; i++){
			fem.pmeshnode[i].A *= 1;
		}
		fem.StaticAxisQ4VTM();
	}
}

void T3NRtest() {
	CFastFEMcore fem;
	double t1;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	fem.openProject("..\\model\\project1.mag");

	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
		fem.preCalculation();
		t1 = SuperLU_timer_();
		fem.StaticAxisymmetricNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "NR:" << t1;
		for (int i = 0; i < fem.num_ele; i++) {
			if (!fem.pmeshele[i].LinearFlag) {
				fem.pmeshele[i].miu = 0.9*fem.pmeshele[i].miut;
			}
		}
		t1 = SuperLU_timer_();
		fem.StaticAxisymmetricNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "NR:" << t1;
		//fem.CalcForce();
	}
}
void triangletestgroup() {
	CFastFEMcore fem;
	double t1;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	fem.openProject("..\\model\\project1.mag");

	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
		fem.preCalculation();
		t1 = SuperLU_timer_();
		fem.StaticAxisymmetricNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "NR:" << t1;
		for (int i = 0; i < fem.num_ele; i++) {
			if (!fem.pmeshele[i].LinearFlag) {
				fem.pmeshele[i].miu = 0.9*fem.pmeshele[i].miut;
			}
		}
		t1 = SuperLU_timer_();
		fem.StaticAxisT3TLMgroup();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "TLM:" << t1;
		//fem.CalcForce();
	}
}
void triangletestvtm() {
	CFastFEMcore fem;
	double t1;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	fem.openProject("..\\model\\project1.mag");

	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
		fem.preCalculation();
		t1 = SuperLU_timer_();
		fem.StaticAxisymmetricNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "NR:" << t1;
		for (int i = 0; i < fem.num_ele; i++) {
			if (!fem.pmeshele[i].LinearFlag) {
				fem.pmeshele[i].miu = 1*fem.pmeshele[i].miut;
			}
		}
		t1 = SuperLU_timer_();
		fem.StaticAxisT3VTM();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "TLM:" << t1;
		//fem.CalcForce();
	}
}

void triangletestvtm2() {
	CFastFEMcore fem;
	double t1;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	fem.openProject("..\\model\\project1.mag");

	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
		fem.preCalculation();
		t1 = SuperLU_timer_();
		fem.StaticAxisymmetricNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "NR:" << t1;
		for (int i = 0; i < fem.num_ele; i++) {
			if (!fem.pmeshele[i].LinearFlag) {
				fem.pmeshele[i].miu = 0.9 * fem.pmeshele[i].miut;
			}
		}
		t1 = SuperLU_timer_();
		fem.StaticAxisT3VTM2();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "TLM:" << t1;
		//fem.CalcForce();
	}
}

void triangletestvtm3() {
	CFastFEMcore fem;
	double t1;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	fem.openProject("..\\model\\project1.mag");

	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
		fem.preCalculation();
		t1 = SuperLU_timer_();
		fem.StaticAxisymmetricNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "NR:" << t1;
		for (int i = 0; i < fem.num_ele; i++) {
			if (!fem.pmeshele[i].LinearFlag) {
				fem.pmeshele[i].miu = 1 * fem.pmeshele[i].miut;
			}
		}
		t1 = SuperLU_timer_();
		fem.StaticAxisT3VTM3();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "TLM:" << t1;
		//fem.CalcForce();
	}
}
void triangletest1() {
	CFastFEMcore fem;
	double t1;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	fem.openProject("..\\model\\project1.mag");

	if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
		fem.preCalculation();
		t1 = SuperLU_timer_();
		//fem.StaticAxisymmetricNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "NR:" << t1;
		for (int i = 0; i < fem.num_ele; i++) {
			if (!fem.pmeshele[i].LinearFlag) {
				fem.pmeshele[i].miu = 0.9*fem.pmeshele[i].miut;
			}
		}
		t1 = SuperLU_timer_();
		fem.StaticAxisTLMNR();
		t1 = SuperLU_timer_() - t1;
		qDebug() << "TLM:" << t1;
		//fem.CalcForce();
	}
}
void triangletest(){
    CFastFEMcore fem;
    double t1;
    qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
    qDebug() << "current currentPath: " << QDir::currentPath();

    fem.openProject("..\\model\\project1.mag");

    if (fem.Load2DMeshCOMSOL("..\\model\\mesh24.mphtxt") == 0) {
        fem.preCalculation();
        t1 = SuperLU_timer_();
        fem.StaticAxisymmetricNR();
        t1 = SuperLU_timer_() - t1;
        qDebug() << "NR:" << t1;
        for (int i = 0; i < fem.num_ele; i++) {
            if (!fem.pmeshele[i].LinearFlag) {
                fem.pmeshele[i].miu = 0.9*fem.pmeshele[i].miut;
            }
        }
        t1 = SuperLU_timer_();
        fem.StaticAxisymmetricTLM();
        t1 = SuperLU_timer_() - t1;
        qDebug() << "TLM:" << t1;
        //fem.CalcForce();
    }
}
