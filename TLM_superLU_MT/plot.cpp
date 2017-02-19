//2016-02-26 by poofee
//add the readfile function to simplify the code.
//add the sourcecode in the git using version control.
//2016-01-20 by poofee
//cheng the Y0 as the diffrent numbers
//as in the papers say
//by poofee
//add the LU fact directly using superLU
//and we don't process the LU fact in every iteration.
#include "plot.h"
#include "ui_mainwindow.h"

#include <stdio.h>

//#include <stdlib.h>
//#include <omp.h>
#include "datatype.h"

#include <cmath>



Plot::Plot(QWidget *parent) :
QMainWindow(parent),
ui(new Ui::MainWindow) {
	ui->setupUi(this);
	QAction *start = new QAction("TLM", this);
	addAction(start);
	connect(start, SIGNAL(triggered()), this, SLOT(TLMcalculation()));
	setContextMenuPolicy(Qt::ActionsContextMenu);
}
Plot::~Plot() {
	delete ui;
}
typedef struct{
	int index;
	double x;
	double y;
} ele;
bool compareele(const ele ele1, const ele ele2) {
	if (ele1.x < ele2.x) {
		return true;
	} else {
		return false;
	}
}
double mymax(double a, double b) {
	return ((a > b) ? a : b);
}
double mymin(double a, double b) {
	return ((a < b) ? a : b);
}


void Plot::TLMcalculation() {
	///*read coarse mesh data from file*/
	//double U0 = 1;// PI*4e-7;
	//char  coarseFN[] = "E:\\Projects\\axispmmodel\\mesh04.mphtxt";
	//CNode * cmeshnode = NULL;	int cnum_pts = 0;
	//CElement * cmeshele = NULL;	int cnum_ele = 0;
	//readDataFile(coarseFN, cnum_pts, cnum_ele, &cmeshnode, &cmeshele);

	///*using TLM get a coarse miu value*/
	//TLM(cnum_pts, cnum_ele, cmeshnode, cmeshele);
	///*read fine mesh data from file*/
	//char fineFN[] = "E:\\Projects\\axispmmodel\\mesh04.mphtxt";
	//CNode * fmeshnode = NULL;	int fnum_pts = 0;
	//CElement * fmeshele = NULL;	int fnum_ele = 0;
	//readDataFile(fineFN, fnum_pts, fnum_ele, &fmeshnode, &fmeshele);
	///*set the miu initial value */
	//double minx, maxx, miny, maxy;
	//QVector <ele> d34;
	//ele eletmp;
	//for (int i = 0; i < fnum_ele; i++) {
	//	if (fmeshele[i].domain == 3 || fmeshele[i].domain == 4) {
	//		eletmp.index = i;
	//		eletmp.x = fmeshele[i].rc;
	//		eletmp.y = fmeshele[i].zc;
	//		d34.push_back(eletmp);
	//	}
	//}
	//qSort(d34.begin(), d34.end(), compareele);//ascend

	//for (int i = 0; i < cnum_ele; i++) {
	//	if (cmeshele[i].domain == 3 || cmeshele[i].domain == 4) {
	//		/*find the minx max miny maxy in the ele*/
	//		double a, b, c;
	//		a = cmeshnode[cmeshele[i].n[0]].x;
	//		b = cmeshnode[cmeshele[i].n[1]].x;
	//		c = cmeshnode[cmeshele[i].n[2]].x;
	//		maxx = mymax(mymax(a, b), c);
	//		minx = mymin(mymin(a, b), c);
	//		a = cmeshnode[cmeshele[i].n[0]].y;
	//		b = cmeshnode[cmeshele[i].n[1]].y;
	//		c = cmeshnode[cmeshele[i].n[2]].y;
	//		maxy = mymax(mymax(a, b), c);
	//		miny = mymin(mymin(a, b), c);
	//		/*binary search to find the index */
	//		int bottom = 0, top = d34.size() - 1;
	//		if (top == -1) {
	//			break;
	//		}
	//		int position;
	//		while (bottom <= top) {
	//			position = (bottom + top) >> 1;
	//			if (d34[position].x >= minx && position == 0) {
	//				break;
	//			}
	//			if (d34[position].x >= minx && d34[position - 1].x < minx) {
	//				break;
	//			}
	//			if (d34[position].x < minx) {
	//				bottom = position + 1;
	//			} else {
	//				top = position - 1;
	//			}
	//		}
	//		/*start from the position*/
	//		for (int j = position; j < d34.size(); j++) {
	//			/*find the first x between minx~maxx*/
	//			if (d34[j].x > maxx) {
	//				break;
	//			} else if (d34[j].x > minx) {
	//				if (d34[j].y > miny && d34[j].y < maxy) {
	//					fmeshele[d34[j].index].miut = cmeshele[i].miu * 1;
	//					d34.remove(j);
	//					j -= 1;
	//				}
	//			}
	//		}
	//	}
	//}
	///*recall the TLM*/
	//TLM(fnum_pts, fnum_ele, fmeshnode, fmeshele);
	//CalcForce(fnum_pts, fnum_ele, fmeshnode, fmeshele);
	///*free the space*/
	//free(cmeshele); free(cmeshnode);
	//free(fmeshele); free(fmeshnode);
}



