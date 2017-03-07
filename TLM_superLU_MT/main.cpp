#include "plot.h"
#include <QApplication>
#include "FastFEMcore.h"

#define PI 3.14159265358979323846

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	CFastFEMcore fem;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();
	fem.openProject("..\\model\\project1.mag");
	if (fem.Load2DMeshCOMSOL("..\\model\\mesh.mphtxt") == 0) {
		fem.preCalculation();
		/*fem.StaticAxisymmetricTLM();
		for (int i = 0; i < fem.num_ele; i++) {
			if (!fem.pmeshele[i].LinearFlag) {
				fem.pmeshele[i].miu = 0.9*fem.pmeshele[i].miut;
			}			
		}
		fem.StaticAxisymmetricTLM();*/
		fem.StaticAxisymmetricNR();
	}	
    return a.exec();
}
