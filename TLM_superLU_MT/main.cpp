#include "mainwindow.h"
#include <QApplication>





#define PI 3.14159265358979323846

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
	//w.TLMcalculation();
    return a.exec();
}
