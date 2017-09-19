#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"


namespace Ui {
	class MainWindow;
}

class Plot : public QMainWindow
{
    Q_OBJECT

public:
	explicit Plot(QWidget *parent = 0);
	~Plot();
	QCustomPlot* getQcustomPlot();
public slots:

private:
	Ui::MainWindow *ui;
	QMenu *fileMenu;

	QAction *openAct;               // File -> Open...
	QAction *loadAct;               // File -> Load...

};

#endif // MAINWINDOW_H
