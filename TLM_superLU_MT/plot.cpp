#include "plot.h"
#include "ui_mainwindow.h"


Plot::Plot(QWidget *parent) :
QMainWindow(parent),
ui(new Ui::MainWindow) {
	ui->setupUi(this);

}
Plot::~Plot() {
	delete ui;
}
QCustomPlot* Plot::getQcustomPlot() {
	return ui->qplot;
}



