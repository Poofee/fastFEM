#include "plot.h"
#include "ui_mainwindow.h"


Plot::Plot(QWidget *parent) :
QMainWindow(parent),
ui(new Ui::MainWindow) {
	ui->setupUi(this);

	//fileMenu = menuBar()->addMenu(tr("&File"));

	//// File -> Open file
	//openAct = new QAction(QIcon(":/icons/document-open.png"), tr("&Open..."), this);
	//openAct->setShortcut(tr("Ctrl+O"));
	//openAct->setStatusTip(tr("Open geometry input file"));
	//connect(openAct, SIGNAL(triggered()), this, SLOT(triangletest()));

	//// File -> Load mesh...
	//loadAct = new QAction(QIcon(":/icons/document-open-folder.png"), tr("&Load mesh..."), this);
	//loadAct->setStatusTip(tr("Load Elmer mesh files"));
	//connect(loadAct, SIGNAL(triggered()), this, SLOT(quadtest()));

	//fileMenu->addAction(openAct);
	//fileMenu->addAction(loadAct);

	//connect(menuBar(), SIGNAL(triggered(QAction*)), this, SLOT(menuBarTriggeredSlot(QAction*)));
}
Plot::~Plot() {
	delete ui;
}
QCustomPlot* Plot::getQcustomPlot() {
	return ui->qplot;
}



