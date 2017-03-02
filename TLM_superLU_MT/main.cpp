#include "plot.h"
#include <QApplication>
#include "FastFEMcore.h"

#define PI 3.14159265358979323846

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	//Plot w;
    //w.show();
	CFastFEMcore fem;
	/*QMessageBox msgBox;
	msgBox.setText("Press Ok.");
	msgBox.exec();*/
	//QProgressDialog dialog;
	//dialog.exec();
	qDebug() << "current applicationDirPath: " << QCoreApplication::applicationDirPath();
	qDebug() << "current currentPath: " << QDir::currentPath();
	fem.openProject("..\\model\\project2.mag");
	fem.Load2DMeshCOMSOL("..\\model\\model_normal.mphtxt");
	fem.preCalculation();
	fem.StaticAxisymmetricTLM();
	//fem.StaticAxisymmetricNR();
    return a.exec();
}
