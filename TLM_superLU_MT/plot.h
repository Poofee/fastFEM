#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include "datatype.h"

namespace Ui {
	class Plot;
}

class Plot : public QMainWindow
{
    Q_OBJECT

public:
	explicit Plot(QWidget *parent = 0);
	~Plot();
	
public slots:
	void TLMcalculation();

private:
	Ui::Plot *ui;

};

#endif // MAINWINDOW_H
