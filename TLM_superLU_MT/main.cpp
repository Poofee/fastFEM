#include "plot.h"
#include <QApplication>
#include "FastFEMcore.h"

#define PI 3.14159265358979323846

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	Plot w;
    w.show();
	CFastFEMcore fem;
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();
	fem.openProject("..\\model\\project1.mag");
	fem.LoadMeshCOMSOL("..\\model\\mesh00.mphtxt");
	fem.preCalculation();
	fem.StaticAxisymmetricNR();
    return a.exec();
}
