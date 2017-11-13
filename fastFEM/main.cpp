#include "plot.h"
#include <QApplication>
#include "FastFEMcore.h"
#include "slu_mt_ddefs.h"

#define PI 3.14159265358979323846

void triangletest();
void quadtest();

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

	triangletest();
	Plot myplot;
	myplot.show();
    return a.exec();
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
