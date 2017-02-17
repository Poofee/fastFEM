#include "FastFEMcore.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo> 
#include <vector>
#include <ctime>
#include <omp.h>
#include "SuperLU_MT.h"
#include <QtAlgorithms>
#include <QVector>
#include "spline.h"

using namespace std;
using namespace arma;

#define PI 3.14159265358979323846
#define r 2
#define INF 1
#define AIR 2
#define FIX 3
#define MOV 4
#define PM 5
#define DOWNCOIL 6
#define UPCOIL 7
const double ste = 4;
const double miu0 = PI*4e-7;

CFastFEMcore::CFastFEMcore() {
	Precision = 1e-6;;
	LengthUnits = 0;
	pmeshnode = NULL;
	pmeshele = NULL;
	materialList = NULL;
}


CFastFEMcore::~CFastFEMcore() {
	//should free the space allocated
	if (pmeshnode != NULL) free(pmeshnode);
	if (pmeshele != NULL) free(pmeshele);
	if (materialList != NULL) free(materialList);
}


// load mesh
int CFastFEMcore::LoadMesh() {
	char ch[256];
	//------------open file----------------------------------
	FILE * fp = NULL;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		//QMessageBox::warning(NULL, "Error:", "Error: opening file!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	//--------------Read the head-----------------------------
	for (int i = 0; i < 18; i++) {
		fgets(ch, 256, fp);
	}
	//-----------------mesh point-----------------------------	

	if (fscanf(fp, "%d # number of mesh points\n", &num_pts)) {
		pmeshnode = (CNode*)calloc(num_pts, sizeof(CNode));

		for (int i = 0; i < num_pts; i++) {
			pmeshnode[i].I = 0;
			pmeshnode[i].pm = 0;
		}
	} else {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_pts!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	int pts_ind;//the beginning of the points index

	if (fscanf(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading pts_ind!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);
	for (int i = pts_ind; i < num_pts; i++) {
		if (fscanf(fp, "%lf %lf \n", &(pmeshnode[i].x), &(pmeshnode[i].y)) != 2) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading mesh point!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	//---------------vertex-------------------------------
	for (int i = 0; i < 7; i++)
		fgets(ch, 256, fp);
	int num_vtx_ns, num_vtx_ele;
	if (fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns) != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ns!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}

	if (fscanf(fp, "%d # number of elements\n", &num_vtx_ele) != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);

	int *vtx;
	vtx = (int*)calloc(num_vtx_ele, sizeof(int));
	for (int i = 0; i < num_vtx_ele; i++) {
		if (fscanf(fp, "%d \n", vtx + i) != 1) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading vertex condition!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	if (vtx != NULL) free(vtx); vtx = NULL;
	//---------------vertex-------------------------------
	int num_vtx_ele2;
	fscanf(fp, "%d # number of geometric entity indices\n", &num_vtx_ele2);
	fgets(ch, 256, fp);
	int *vtx2;
	vtx2 = (int*)calloc(num_vtx_ele2, sizeof(int));
	for (int i = 0; i < num_vtx_ele2; i++) {
		if (fscanf(fp, "%d \n", vtx2 + i) != 1) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading vertex condition!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	if (vtx2 != NULL) free(vtx2); vtx2 = NULL;
	//--------------boundary--------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int num_bdr_ns, num_bdr_ele;//number of nodes per element;number of elements
	if (fscanf(fp, "%d # number of nodes per element\n", &num_bdr_ns) != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_bdr_ns!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}

	if (fscanf(fp, "%d # number of elements\n", &num_bdr_ele) != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_bdr_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);

	int *p1, *p2;
	p1 = (int*)calloc(num_bdr_ele, sizeof(int));
	p2 = (int*)calloc(num_bdr_ele, sizeof(int));
	for (int i = 0; i < num_bdr_ele; i++) {
		if (fscanf(fp, "%d %d\n", p1 + i, p2 + i) == 2) {
			//-----process the A=0 boundary--------------------------
			/*if (abs(pmeshnode[p1[i]].length() - 0.05) < 5e-3 || abs(pmeshnode[p1[i]].x) < 5e-5) {
			pmeshnode[p1[i]].bdr = 1;
			}*/
		} else {
			//QMessageBox::warning(NULL, "Error:", "Error: reading boundary condition!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	if (p1 != NULL) free(p1); p1 = NULL;
	if (p2 != NULL) free(p2); p2 = NULL;
	//---------------entity----------------------------------
	int num_entity;
	fscanf(fp, "%d # number of geometric entity indices\n", &num_entity);
	fgets(ch, 256, fp);
	int * entity;
	entity = (int*)calloc(num_entity, sizeof(int));
	for (int i = 0; i < num_entity; i++) {
		if (fscanf(fp, "%d \n", entity + i) != 1) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading boundary condition!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	if (entity != NULL) free(entity); entity = NULL;
	//----------------elements------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int ns_per_ele;// num_ele;//number of nodes per element;number of elements
	if (fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele) != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading ns_per_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}

	if (fscanf(fp, "%d # number of elements\n", &num_ele) == 1) {
		pmeshele = (CElement*)calloc(num_ele, sizeof(CElement));
	} else {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);

	for (int i = 0; i < num_ele; i++) {
		if (fscanf(fp, "%d %d %d \n", &pmeshele[i].n[0], &pmeshele[i].n[1], &pmeshele[i].n[2]) != 3) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading elements points!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	//---------------Domain----------------------------------
	int num_domain;
	fscanf(fp, "%d # number of geometric entity indices\n", &num_domain);
	fgets(ch, 256, fp);

	for (int i = 0; i < num_domain; i++) {
		if (fscanf(fp, "%d \n", &pmeshele[i].domain) != 1) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading domain points!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	fclose(fp);
	return 0;
}


bool CFastFEMcore::StaticAxisymmetric() {
	double U0 = 1;// PI*4e-7;
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele[i].LinearFlag) {
			D34.push_back(i);
		}
	}

	//------------build C Matrix-----------------------------
	umat locs(2, 9 * num_ele);
	locs.zeros();
	vec vals = zeros<vec>(9 * num_ele);
	double ce[3][3];
	ResistMarix *rm = (ResistMarix*)malloc(num_ele * sizeof(ResistMarix));
	vec bbJz = zeros<vec>(num_pts);
	vec b = zeros<vec>(num_pts);//b = bbJz + INL;
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	vec INL = zeros<vec>(num_pts);
	vec bb = zeros<vec>(num_pts);
	vec rpm = zeros<vec>(num_pts);
	double * ydot = (double*)malloc(num_ele);
	//轴对称：A'=rA,v'=v/r,
	for (int i = 0; i < num_ele; i++) {
		//确定单元的近似半径
		int flag = 0;
		for (int f = 0; f < 3; f++)
			if (pmeshnode[pmeshele[i].n[f]].x < 1e-6)
				flag++;

		if (flag == 2) {
			ydot = pmeshele[i].rc;
		} else {
			ydot[i] = 1 / (pmeshnode[pmeshele[i].n[0]].x + pmeshnode[pmeshele[i].n[1]].x);
			ydot[i] += 1 / (pmeshnode[pmeshele[i].n[0]].x + pmeshnode[pmeshele[i].n[2]].x);
			ydot[i] += 1 / (pmeshnode[pmeshele[i].n[1]].x + pmeshnode[pmeshele[i].n[2]].x);
			ydot[i] = 1.5 / ydot[i];
		}
		//计算单元导纳
		rm[i].Y11 = pmeshele[i].Q[0] * pmeshele[i].Q[0] + pmeshele[i].P[0] * pmeshele[i].P[0];
		rm[i].Y12 = pmeshele[i].Q[0] * pmeshele[i].Q[1] + pmeshele[i].P[0] * pmeshele[i].P[1];
		rm[i].Y13 = pmeshele[i].Q[0] * pmeshele[i].Q[2] + pmeshele[i].P[0] * pmeshele[i].P[2];
		rm[i].Y22 = pmeshele[i].Q[1] * pmeshele[i].Q[1] + pmeshele[i].P[1] * pmeshele[i].P[1];
		rm[i].Y23 = pmeshele[i].Q[1] * pmeshele[i].Q[2] + pmeshele[i].P[1] * pmeshele[i].P[2];
		rm[i].Y33 = pmeshele[i].Q[2] * pmeshele[i].Q[2] + pmeshele[i].P[2] * pmeshele[i].P[2];

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i];
		//生成单元矩阵，线性与非线性
		if (pmeshele[i].LinearFlag) {
			ce[0][0] = abs(rm[i].Y11) / pmeshele[i].miu;
			ce[1][1] = abs(rm[i].Y22) / pmeshele[i].miu;
			ce[2][2] = abs(rm[i].Y33) / pmeshele[i].miu;

			ce[0][1] = -abs(rm[i].Y12) / pmeshele[i].miu;
			ce[0][2] = -abs(rm[i].Y13) / pmeshele[i].miu;
			ce[1][2] = -abs(rm[i].Y23) / pmeshele[i].miu;

			ce[1][0] = ce[0][1];
			ce[2][0] = ce[0][2];
			ce[2][1] = ce[1][2];
		} else {
			ce[0][0] = abs(rm[i].Y11) / pmeshele[i].miut;
			ce[1][1] = abs(rm[i].Y22) / pmeshele[i].miut;
			ce[2][2] = abs(rm[i].Y33) / pmeshele[i].miut;

			ce[0][1] = -abs(rm[i].Y12) / pmeshele[i].miut;
			ce[0][2] = -abs(rm[i].Y13) / pmeshele[i].miut;
			ce[1][2] = -abs(rm[i].Y23) / pmeshele[i].miut;

			ce[1][0] = ce[0][1];
			ce[2][0] = ce[0][2];
			ce[2][1] = ce[1][2];
		}
		//将单元矩阵进行存储
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				locs(0, i * 9 + row * 3 + col) = pmeshele[i].n[row];
				locs(1, i * 9 + row * 3 + col) = pmeshele[i].n[col];
				vals(i * 9 + row * 3 + col) = ce[row][col];
			}
		}
		//计算电流密度//要注意domain会不会越界
		double jr = pmeshele[i].AREA*materialList[pmeshele[i].domain].Jr / 3;
		pmeshnode[pmeshele[i].n[0]].I += jr;
		pmeshnode[pmeshele[i].n[1]].I += jr;
		pmeshnode[pmeshele[i].n[2]].I += jr;
		//计算永磁部分
		pmeshnode[pmeshele[i].n[0]].pm += materialList[pmeshele[i].domain].H_c / 2.*pmeshele[i].Q[0];
		pmeshnode[pmeshele[i].n[1]].pm += materialList[pmeshele[i].domain].H_c / 2.*pmeshele[i].Q[1];
		pmeshnode[pmeshele[i].n[2]].pm += materialList[pmeshele[i].domain].H_c / 2.*pmeshele[i].Q[2];
	}
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts, num_pts, true, true);

	for (int i = 0; i < num_pts; i++) {
		bbJz(i) = pmeshnode[i].I;//set the current
		rpm(i) = pmeshnode[i].pm;
	}
	b = bbJz + INL + rpm;
	//---------------------superLU_MT---------------------------------------
	CSuperLU_MT superlumt(num_pts, X, b);
	superlumt.solve();
	double *sol = NULL;
	sol = superlumt.getResult();
	for (int i = 0; i < num_pts; i++) {
		pmeshnode[i].A = sol[i] * miu0;//the A is r*A_real
		A(i) = sol[i] * miu0;
	}
	//---------------------superLU--end----------------------------------
	QVector<double> x1(D34.size()), y1(D34.size());
	QVector<double> x2(num_pts), y2(num_pts);
	for (int i = 0; i < D34.size(); ++i) {
		x1[i] = i;
		x2[i] = i;
	}
	//---------the main loop---------------------------------
	double Precision = 1e-6;
	int steps = 250;

	int count;
	//int i,j;
	Voltage *Vr = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	Voltage *Vi = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	for (count = 0; count < steps; count++) {
		//------update miu----------------
		for (int i = 0; i < D34.size(); i++) {
			int d34 = D34[i];
			pmeshele[d34].Bx = (pmeshele[d34].Q[0] * pmeshnode[pmeshele[d34].n[0]].A
				+ pmeshele[d34].Q[1] * pmeshnode[pmeshele[d34].n[1]].A
				+ pmeshele[d34].Q[2] * pmeshnode[pmeshele[d34].n[2]].A) / 2. / pmeshele[d34].AREA / ydot[d34];
			pmeshele[d34].By = (pmeshele[d34].P[0] * pmeshnode[pmeshele[d34].n[0]].A
				+ pmeshele[d34].P[1] * pmeshnode[pmeshele[d34].n[1]].A
				+ pmeshele[d34].P[2] * pmeshnode[pmeshele[d34].n[2]].A) / 2. / pmeshele[d34].AREA / ydot[d34];

			pmeshele[d34].B = sqrt(pmeshele[d34].By*pmeshele[d34].By +
				pmeshele[d34].Bx*pmeshele[d34].Bx);
			pmeshele[d34].miu = materialList[pmeshele[d34].domain].GetMiu(pmeshele[d34].B);
			y1[i] = pmeshele[d34].B;
		}
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *meshelement = &pmeshele[i];
			double rtmp;//do this to mark it as private
			rtmp = (meshelement->miu - meshelement->miut) / (meshelement->miu + meshelement->miut);

			Vr[j].V12 = (pmeshnode[meshelement->n[0]].A - pmeshnode[meshelement->n[1]].A) - Vi[j].V12;
			Vr[j].V23 = (pmeshnode[meshelement->n[1]].A - pmeshnode[meshelement->n[2]].A) - Vi[j].V23;
			Vr[j].V13 = (pmeshnode[meshelement->n[2]].A - pmeshnode[meshelement->n[0]].A) - Vi[j].V13;


			Vi[j].V12 = rtmp*Vr[j].V12;
			INL(pmeshele[i].n[1]) += -2. *Vi[j].V12*abs(rm[i].Y12) / pmeshele[i].miut;
			INL(pmeshele[i].n[0]) += 2. * Vi[j].V12 *abs(rm[i].Y12) / pmeshele[i].miut;


			Vi[j].V23 = rtmp*Vr[j].V23;
			INL(pmeshele[i].n[1]) += 2. * Vi[j].V23*abs(rm[i].Y23) / pmeshele[i].miut;
			INL(pmeshele[i].n[2]) += -2. *Vi[j].V23*abs(rm[i].Y23) / pmeshele[i].miut;

			Vi[j].V13 = rtmp*Vr[j].V13;
			INL(pmeshele[i].n[2]) += 2. * Vi[j].V13*abs(rm[i].Y13) / pmeshele[i].miut;
			INL(pmeshele[i].n[0]) += -2.0 *Vi[j].V13*abs(rm[i].Y13) / pmeshele[i].miut;
		}
		b = bbJz + INL + rpm;
		A_old = A;

		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		superlumt.LUsolve();
		for (int i = 0; i < num_pts; i++) {
			pmeshnode[i].A = sol[i] * miu0;
			A(i) = sol[i] * miu0;
		}
		double error = norm((A_old - A), 2) / norm(A, 2);
		if (error < Precision) {
			break;
		}
		INL.zeros();
	}
	for (int i = 0; i < num_ele; i++) {
		pmeshele[i].Bx = (pmeshele[i].Q[0] * pmeshnode[pmeshele[i].n[0]].A
			+ pmeshele[i].Q[1] * pmeshnode[pmeshele[i].n[1]].A
			+ pmeshele[i].Q[2] * pmeshnode[pmeshele[i].n[2]].A) / 2. / pmeshele[i].AREA / ydot[i];
		pmeshele[i].By = (pmeshele[i].P[0] * pmeshnode[pmeshele[i].n[0]].A
			+ pmeshele[i].P[1] * pmeshnode[pmeshele[i].n[1]].A
			+ pmeshele[i].P[2] * pmeshnode[pmeshele[i].n[2]].A) / 2. / pmeshele[i].AREA / ydot[i];

		pmeshele[i].B = sqrt(pmeshele[i].By*pmeshele[i].By +
			pmeshele[i].Bx*pmeshele[i].Bx);
	}

	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vi != NULL) free(Vi);
	if (Vr != NULL) free(Vr);
	return true;
}

double CFastFEMcore::HB(double B) {
	double miu;
	double hint;
	double h[11] = { 0, 200, 500, 1000, 2500, 5000, 10000, 15000, 100000, 500000, 550000 };
	double b[11] = { 0, 1.2, 1.4, 1.5, 1.62, 1.71, 1.85, 1.851, 1.8511, 5, 5 };
	std::vector<double> BB(b, b + 11);
	std::vector<double> HH(h, h + 11);
	tk::spline s;
	s.set_points(BB, HH, false);
	hint = s(B);
	if (B > 4) {
		B = 2;
		miu = B / (1e4 + 5e6*(B - 1.85)) / miu0;
	} else {
		miu = B / hint / miu0;
	}
	return miu;
}


double CFastFEMcore::CalcForce() {
	/*******2016-12-28 by Poofee*************/
	/********Find the first layer***********/
	double xLeft = 1.25e-3;
	double xRight = 5.8e-3;
	double yDown = -4.7e-3;
	double yUp = 5.3e-3;
	double delta = 1e-10;
	double xForce = 0;
	double yForce = 0;
	yUp += ste * 1e-4;
	yDown += ste * 1e-4;

	/*QCustomPlot *customPlot;
	customPlot = ui->widget;
	customPlot->xAxis->setLabel("x");
	customPlot->xAxis->setRange(0, 0.09);
	customPlot->xAxis->setAutoTickStep(false);
	customPlot->xAxis->setTicks(false);
	customPlot->yAxis->setLabel("y");
	customPlot->yAxis->setRange(-0.09, 0.09);
	customPlot->xAxis2->setTicks(false);
	customPlot->yAxis->setScaleRatio(ui->widget->xAxis, 1.0);*/


	/*customPlot->yAxis->setAutoTickStep(false);
	ui->widget->yAxis->setAutoTickLabels(false);
	customPlot->yAxis->setTicks(false);
	customPlot->yAxis->grid()->setVisible(false);
	customPlot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
	customPlot->yAxis2->setTicks(false);
	customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);*/

	/*QCPCurve *newCurve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);

	QCPGraph *graph1 = customPlot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 5), QBrush(Qt::black), 5));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	QCPGraph *graph2 = customPlot->addGraph();
	graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, QPen(Qt::red, 5), QBrush(Qt::red), 5));
	graph2->setPen(QPen(QColor(120, 120, 120), 2));
	graph2->setLineStyle(QCPGraph::lsNone);*/
	QVector <double> gpx1, gpx2, gpy1, gpy2;
	//QColor cc[7];
	///*cc[0] = QColor(0,0,0);
	//cc[1] = QColor(0, 120, 0);
	//cc[2] = QColor(0, 0, 120);
	//cc[3] = QColor(120, 0, 0);
	//cc[4] = QColor(120, 120, 0);
	//cc[5] = QColor(0, 120, 120);
	//cc[6] = QColor(120, 0, 120);*/
	//cc[0] = QColor(0, 0, 0);
	//cc[1] = QColor(0, 120, 0);
	//cc[2] = QColor(68, 49, 242);
	//cc[3] = QColor(254, 31, 87);
	//cc[4] = QColor(255, 128, 0);
	//cc[5] = QColor(232, 223, 60);
	//cc[6] = QColor(232, 223, 60);
	//customPlot->xAxis->grid()->setSubGridVisible(false);
	//customPlot->xAxis->grid()->setSubGridPen(Qt::NoPen);
	//customPlot->yAxis->grid()->setSubGridVisible(false);
	//customPlot->yAxis->grid()->setSubGridPen(Qt::NoPen);
	FILE * fp = NULL;
	fp = fopen("E:\\index.txt", "w+");//delete exist, read and write
	for (int i = 0; i < num_ele; i++) {
		//QCPCurve *newCurve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
		QVector <double> x1(4);
		QVector <double> y1(4);
		int BCx = 0;
		int BCy = 0;
		for (int j = 0; j < 3; j++) {
			x1[j] = pmeshnode[pmeshele[i].n[j]].x;
			y1[j] = pmeshnode[pmeshele[i].n[j]].y;
		}
		x1[3] = x1[0];
		y1[3] = y1[0];
		if (pmeshele[i].domain == AIR) {
			for (int k = 0; k < 3; k++) {
				if (fabs(x1[k] - xLeft) < delta) {
					if (y1[k] >= yDown - delta && y1[k] <= yUp + delta) {
						BCx += 1;
					}
				}
				if (fabs(x1[k] - xRight) < delta) {
					if (y1[k] >= yDown - delta && y1[k] <= yUp + delta) {
						BCx += 3;
					}
				}
				if (fabs(y1[k] - yDown) < delta) {
					if (x1[k] >= xLeft - delta && x1[k] <= xRight + delta) {
						BCy += 1;
					}
				}
				if (fabs(y1[k] - yUp) < delta) {
					if (x1[k] >= xLeft - delta && x1[k] <= xRight + delta) {
						BCy += 3;
					}
				}
			}
		}

		int ind = 10;
		int s = 0;
		int num = 0;
		int index[3];
		if (BCx == 0 && BCy == 0) {//outside
			//newCurve->setBrush(Qt::NoBrush);
		} else if (BCy == 2 || BCy == 6) {//Y,2,edge
			//newCurve->setBrush(QColor(0, 0, 255));//blue
			/****Find which point is moved,1?2?3?*****/
			if (fabs(y1[0] - y1[1]) < delta) {
				ind = 2;
				num = 2;
				index[0] = 0;
				index[1] = 1;
			} else if (fabs(y1[1] - y1[2]) < delta) {
				ind = 0;
				num = 2;
				index[0] = 1;
				index[1] = 2;
			} else if (fabs(y1[0] - y1[2]) < delta) {
				ind = 1;
				num = 2;
				index[0] = 0;
				index[1] = 2;
			}
			s = 1;

			gpx1.push_back(x1[ind]);
			gpy1.push_back(y1[ind]);
		} else if (BCx == 0 && BCy % 2 == 1) {//Y,1,point
			//newCurve->setBrush(QColor(255, 0, 0));
			/****Find which point is moved,1?2?3?*****/
			for (int k = 0; k < 3; k++) {
				if (fabs(y1[k] - yUp) < delta) {
					ind = k;
					num = 1;
					index[0] = k;
					break;
				}
				if (fabs(y1[k] - yDown) < delta) {
					ind = k;
					num = 1;
					index[0] = k;
					break;
				}
			}
			s = 1;
			gpx2.push_back(x1[ind]);
			gpy2.push_back(y1[ind]);
		} else if (BCx == 2 || BCx == 6) {//X,2,edge
			//newCurve->setBrush(QColor(0, 0, 255));//blue
			/****Find which point is moved,1?2?3?*****/
			if (fabs(x1[0] - x1[1]) < delta) {
				ind = 2;
				num = 2;
				index[0] = 0;
				index[1] = 1;
			} else if (fabs(x1[1] - x1[2]) < delta) {
				ind = 0;
				num = 2;
				index[0] = 1;
				index[1] = 2;
			} else if (fabs(x1[0] - x1[2]) < delta) {
				ind = 1;
				num = 2;
				index[0] = 0;
				index[1] = 2;
			}
			s = 1;
			gpx1.push_back(x1[ind]);
			gpy1.push_back(y1[ind]);
		} else if (BCx % 2 == 1 && BCy == 0) {//X,1,point
			//newCurve->setBrush(QColor(255, 0, 0));
			/****Find which point is moved,1?2?3?*****/
			for (int k = 0; k < 3; k++) {
				if (fabs(x1[k] - xLeft) < delta) {
					ind = k;
					num = 1;
					index[0] = k;
					break;
				}
				if (fabs(x1[k] - xRight) < delta) {
					ind = k;
					num = 1;
					index[0] = k;
					break;
				}
			}
			s = 1;
			gpx2.push_back(x1[ind]);
			gpy2.push_back(y1[ind]);
		} else {//corner,point
			//newCurve->setBrush(QColor(0, 255, 0));//green	
			/****Find which point is moved,1?2?3?*****/
			for (int k = 0; k < 3; k++) {
				if (fabs(y1[k] - yUp) < delta) {
					ind = k;
					num = 1;
					index[0] = k;
					break;
				}
				if (fabs(y1[k] - yDown) < delta) {
					ind = k;
					num = 1;
					index[0] = k;
					break;
				}
			}
			s = 1;
			gpx2.push_back(x1[ind]);
			gpy2.push_back(y1[ind]);
		}

		/*****Cacl the Force*******/
		for (int n = 0; n < num; n++) {
			double xf1 = 0;
			double xf2 = 0;
			double xf3 = 0;
			double beta2 = 0;
			double beta3;
			double Ac = 0;
			double A1, A2, A3;
			double tmp = PI / 2 * PI * pmeshele[i].rc * pmeshele[i].AREA / (4 * PI*1e-7);
			A1 = pmeshnode[pmeshele[i].n[0]].A;
			A2 = pmeshnode[pmeshele[i].n[1]].A;
			A3 = pmeshnode[pmeshele[i].n[2]].A;
			Ac = (A1 + A2 + A3) / 3;
			if (ind != 10) {
				ind = index[n];
				xf1 = pmeshele[i].B * pmeshele[i].B / pmeshele[i].AREA;
				xf1 *= pmeshele[i].Q[ind] / 2.0;

				beta2 = pmeshele[i].Q[0] * A1 +
					pmeshele[i].Q[1] * A2 +
					pmeshele[i].Q[2] * A3;
				beta2 /= (2. * pmeshele[i].AREA);
				xf2 = -(2. * beta2 + 2.0 * Ac / pmeshele[i].rc)*beta2 / (2 * pmeshele[i].AREA)*(pmeshele[i].Q[ind]);

				beta3 = pmeshele[i].P[0] * A1 +
					pmeshele[i].P[1] * A2 +
					pmeshele[i].P[2] * A3;
				beta3 /= (2. * pmeshele[i].AREA);
				if (ind == 0) {
					xf3 = 2. * beta3*((A3 - A2)*(2. * pmeshele[i].AREA) - beta3*(2. * pmeshele[i].AREA)*(pmeshele[i].Q[0]));
					xf3 /= (2. * pmeshele[i].AREA)*(2. * pmeshele[i].AREA);
				} else if (ind == 1) {
					xf3 = 2. * beta3*((A1 - A3)*(2. * pmeshele[i].AREA) - beta3*(2. * pmeshele[i].AREA)*(pmeshele[i].Q[1]));
					xf3 /= (2. * pmeshele[i].AREA)*(2. * pmeshele[i].AREA);
				} else if (ind == 2) {
					xf3 = 2. * beta3*((A2 - A1)*(2. * pmeshele[i].AREA) - beta3*(2. * pmeshele[i].AREA)*(pmeshele[i].Q[2]));
					xf3 /= (2. * pmeshele[i].AREA)*(2. * pmeshele[i].AREA);
				}
				yForce -= (tmp*(xf1 + xf2 + xf3));
				//fprintf(fp, "%d\t%lf\t%lf\n", i, pmeshele[i].B, tmp*(xf1 + xf2 + xf3));
			}
		}

		//if (pmeshele[i].domain != AIR && pmeshele[i].domain != INF) {
		//newCurve->setData(x1, y1);
		//newCurve->setPen(QPen(cc[pmeshele[i].domain - 1]));
		//newCurve->setBrush(cc[pmeshele[i].domain - 1]);
		//}		
		//if (ind != 10) {
		//	fprintf(fp, "%d\n", ind);			
		//} 

	}
	//graph1->setData(gpx1, gpy1);
	//graph2->setData(gpx2, gpy2);
	fclose(fp);
	//this->setWindowTitle(QString::number(yForce));
	//customPlot->replot();
	return 0;
}


int CFastFEMcore::openProject() {

	return 0;
}


int CFastFEMcore::preCalculation() {
	for (int i = 0; i < num_ele; i++) {
		pmeshele[i].P[0] = pmeshnode[pmeshele[i].n[1]].y - pmeshnode[pmeshele[i].n[2]].y;
		pmeshele[i].P[1] = pmeshnode[pmeshele[i].n[2]].y - pmeshnode[pmeshele[i].n[0]].y;
		pmeshele[i].P[2] = pmeshnode[pmeshele[i].n[0]].y - pmeshnode[pmeshele[i].n[1]].y;

		pmeshele[i].Q[0] = pmeshnode[pmeshele[i].n[2]].x - pmeshnode[pmeshele[i].n[1]].x;
		pmeshele[i].Q[1] = pmeshnode[pmeshele[i].n[0]].x - pmeshnode[pmeshele[i].n[2]].x;
		pmeshele[i].Q[2] = pmeshnode[pmeshele[i].n[1]].x - pmeshnode[pmeshele[i].n[0]].x;

		pmeshele[i].AREA = 0.5*abs(pmeshele[i].P[1] * pmeshele[i].Q[2] - pmeshele[i].Q[1] * pmeshele[i].P[2]);
		pmeshele[i].rc = (pmeshnode[pmeshele[i].n[0]].x +
			pmeshnode[pmeshele[i].n[1]].x +
			pmeshnode[pmeshele[i].n[2]].x) / 3;
		pmeshele[i].zc = (pmeshnode[pmeshele[i].n[0]].y +
			pmeshnode[pmeshele[i].n[1]].y +
			pmeshnode[pmeshele[i].n[2]].y) / 3;

		//主要根据材料属性完成单元当中miu,miut,的赋值；
		//由于I,pm与形函数有关系，为实现分离，不在此计算
		CMaterial  mat;
		mat = materialList[pmeshele[i].domain];

		if (mat.BHpoints == 0) {
			pmeshele[i].miu = 1;
			pmeshele[i].miut = 1;
			pmeshele[i].LinearFlag = true;
		} else {
			pmeshele[i].miu = 1;
			pmeshele[i].miut = 100;
			pmeshele[i].LinearFlag = false;
		}
	}

	return 0;
}

//调用其他的子函数完成求解任务，这个求解可以有多个选项，
//使用NR或者TLM迭代算法，或者选择不同的形函数。
int CFastFEMcore::solve() {
	LoadMesh();
	preCalculation();

	return 0;
}
