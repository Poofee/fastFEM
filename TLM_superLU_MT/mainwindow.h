#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include "datatype.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
	
	double HB(double B);
	void readDataFile(char *fileName, int & num_pts, int & num_ele, CNode **ppmeshnode, CElement ** ppmeshele);
    void TLM(int & num_pts, int & num_ele, CNode *pmeshnode, CElement * m_l);
	void CalcForce(int & num_pts, int & num_ele, CNode *pmeshnode, CElement * m_l);
public slots:
	void TLMcalculation();

private:
    Ui::MainWindow *ui;

};

#endif // MAINWINDOW_H
