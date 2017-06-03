#include "plot.h"
#include <QApplication>
#include "FastFEMcore.h"
#include "slu_mt_ddefs.h"

#define PI 3.14159265358979323846

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	CFastFEMcore fem;
	double t1;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();

    fem.openProject("..\\..\\model\\project1.mag");

    if (fem.Load2DMeshCOMSOL("..\\..\\model\\mesh00.mphtxt") == 0) {
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
    return a.exec();
}
