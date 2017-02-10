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
	Precision = 0.0;
	LengthUnits = 0;
}


CFastFEMcore::~CFastFEMcore() {
}


// load mesh
int CFastFEMcore::LoadMesh() {
	char ch[256];
	int re;
	double U0 = 1;// PI*4e-7;
	double nonlinear_miu = 1;
	double Scoil = 14.4e-3*12.7e-3;//11.687e-3 * 14.7e-3;
	double I = 1;// 12 / 3.2;
	double N = 1420;
	double J = I * N / Scoil;
	double Jdown = 1 * 1250 / (12.4*18.4*1e-6);
	double Hc = 883310;
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
	re = fscanf(fp, "%d # number of mesh points\n", &num_pts);
	if (re == 1) {
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
	re = fscanf(fp, "%d # lowest mesh point index\n", &pts_ind);
	if (re != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading pts_ind!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);
	for (int i = pts_ind; i < num_pts; i++) {
		re = fscanf(fp, "%lf %lf \n", &(pmeshnode[i].x), &(pmeshnode[i].y));
		if (re != 2) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading mesh point!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	//---------------vertex-------------------------------
	for (int i = 0; i < 7; i++)
		fgets(ch, 256, fp);
	int num_vtx_ns, num_vtx_ele;
	re = fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns);
	if (re != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ns!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}

	re = fscanf(fp, "%d # number of elements\n", &num_vtx_ele);
	if (re != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);

	int *vtx;
	vtx = (int*)calloc(num_vtx_ele, sizeof(int));
	for (int i = 0; i < num_vtx_ele; i++) {
		re = fscanf(fp, "%d \n", vtx + i);
		if (re != 1) {
			//QMessageBox::warning(NULL, "Error:", "Error: reading vertex condition!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	if (vtx != NULL) free(vtx); vtx = NULL;
	//---------------vertex-------------------------------
	int num_vtx_ele2;
	re = fscanf(fp, "%d # number of geometric entity indices\n", &num_vtx_ele2);
	fgets(ch, 256, fp);
	int *vtx2;
	vtx2 = (int*)calloc(num_vtx_ele2, sizeof(int));
	for (int i = 0; i < num_vtx_ele2; i++) {
		re = fscanf(fp, "%d \n", vtx2 + i);
		if (re != 1) {
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
	re = fscanf(fp, "%d # number of nodes per element\n", &num_bdr_ns);
	if (re != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_bdr_ns!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}

	re = fscanf(fp, "%d # number of elements\n", &num_bdr_ele);
	if (re != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_bdr_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);

	int *p1, *p2;
	p1 = (int*)calloc(num_bdr_ele, sizeof(int));
	p2 = (int*)calloc(num_bdr_ele, sizeof(int));
	for (int i = 0; i < num_bdr_ele; i++) {
		re = fscanf(fp, "%d %d\n", p1 + i, p2 + i);
		if (re == 2) {
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
	re = fscanf(fp, "%d # number of geometric entity indices\n", &num_entity);
	fgets(ch, 256, fp);
	int * entity;
	entity = (int*)calloc(num_entity, sizeof(int));
	for (int i = 0; i < num_entity; i++) {
		re = fscanf(fp, "%d \n", entity + i);
		if (re != 1) {
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
	re = fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele);
	if (re != 1) {
		//QMessageBox::warning(NULL, "Error:", "Error: reading ns_per_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}

	
	re = fscanf(fp, "%d # number of elements\n", &num_ele);
	if (re == 1) {
		pmeshele = (CElement*)calloc(num_ele, sizeof(CElement));
	} else {
		//QMessageBox::warning(NULL, "Error:", "Error: reading num_ele!",
		//	QMessageBox::Ok, QMessageBox::Ok);
		return 1;
	}
	fgets(ch, 256, fp);

	for (int i = 0; i < num_ele; i++) {
		re = fscanf(fp, "%d %d %d \n", &pmeshele[i].n[0], &pmeshele[i].n[1], &pmeshele[i].n[2]);
		if (re == 3) {
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

			pmeshele[i].y12 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[0] * pmeshele[i].P[1] +
				pmeshele[i].Q[0] * pmeshele[i].Q[1]);
			pmeshele[i].y12 += (pmeshele[i].P[0] + pmeshele[i].P[1]) / 6;
			pmeshele[i].y12 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y12 *= 2 * PI;

			pmeshele[i].y23 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[1] * pmeshele[i].P[2] +
				pmeshele[i].Q[1] * pmeshele[i].Q[2]);
			pmeshele[i].y23 += (pmeshele[i].P[1] + pmeshele[i].P[2]) / 6;
			pmeshele[i].y23 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y23 *= 2 * PI;

			pmeshele[i].y31 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[2] * pmeshele[i].P[0] +
				pmeshele[i].Q[2] * pmeshele[i].Q[0]);
			pmeshele[i].y31 += (pmeshele[i].P[2] + pmeshele[i].P[0]) / 6;
			pmeshele[i].y31 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y31 *= 2 * PI;

			pmeshele[i].y10 = 3 * pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y10 += (4 * pmeshele[i].P[0] + pmeshele[i].P[1] + pmeshele[i].P[2]) / 6;
			pmeshele[i].y10 *= 2 * PI;
			pmeshele[i].y20 = 3 * pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y20 += (pmeshele[i].P[0] + 4 * pmeshele[i].P[1] + pmeshele[i].P[2]) / 6;
			pmeshele[i].y20 *= 2 * PI;
			pmeshele[i].y30 = 3 * pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y30 += (pmeshele[i].P[0] + pmeshele[i].P[1] + 4 * pmeshele[i].P[2]) / 6;
			pmeshele[i].y30 *= 2 * PI;


			pmeshele[i].y11 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[0] * pmeshele[i].P[0] +
				pmeshele[i].Q[0] * pmeshele[i].Q[0]);
			pmeshele[i].y11 += (pmeshele[i].P[0] + pmeshele[i].P[0]) / 6;
			pmeshele[i].y11 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y11 *= 2 * PI;

			pmeshele[i].y22 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[1] * pmeshele[i].P[1] +
				pmeshele[i].Q[1] * pmeshele[i].Q[1]);
			pmeshele[i].y22 += (pmeshele[i].P[1] + pmeshele[i].P[1]) / 6;
			pmeshele[i].y22 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y22 *= 2 * PI;

			pmeshele[i].y33 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[2] * pmeshele[i].P[2] +
				pmeshele[i].Q[2] * pmeshele[i].Q[2]);
			pmeshele[i].y33 += (pmeshele[i].P[2] + pmeshele[i].P[2]) / 6;
			pmeshele[i].y33 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
			pmeshele[i].y33 *= 2 * PI;

			pmeshele[i].vi10 = 0; pmeshele[i].vi20 = 0; pmeshele[i].vi30 = 0;
			//pmeshele[i].vi40 = 0; pmeshele[i].vi50 = 0; pmeshele[i].vi60 = 0;
			pmeshele[i].vi12 = 0; pmeshele[i].vi23 = 0; pmeshele[i].vi31 = 0;
		} else {
			//QMessageBox::warning(NULL, "Error:", "Error: reading elements points!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	//---------------Domain----------------------------------
	int num_domain;
	re = fscanf(fp, "%d # number of geometric entity indices\n", &num_domain);
	fgets(ch, 256, fp);

	for (int i = 0; i < num_domain; i++) {
		re = fscanf(fp, "%d \n", &pmeshele[i].domain);
		if (re == 1) {
			switch (pmeshele[i].domain) {
			case INF:
			case AIR:
				pmeshele[i].miu = U0;
				pmeshele[i].miut = U0;
				pmeshnode[pmeshele[i].n[0]].I += 0;
				pmeshnode[pmeshele[i].n[1]].I += 0;
				pmeshnode[pmeshele[i].n[2]].I += 0;
				pmeshnode[pmeshele[i].n[0]].pm += 0;
				pmeshnode[pmeshele[i].n[1]].pm += 0;
				pmeshnode[pmeshele[i].n[2]].pm += 0;
				break;
			case DOWNCOIL:
				pmeshele[i].miu = U0;
				pmeshele[i].miut = U0;
				pmeshnode[pmeshele[i].n[0]].I += 2. / 3 * PI*pmeshele[i].rc*pmeshele[i].AREA*Jdown;
				pmeshnode[pmeshele[i].n[1]].I += 2. / 3 * PI*pmeshele[i].rc*pmeshele[i].AREA*Jdown;
				pmeshnode[pmeshele[i].n[2]].I += 2. / 3 * PI*pmeshele[i].rc*pmeshele[i].AREA*Jdown;
				pmeshnode[pmeshele[i].n[0]].pm += 0;
				pmeshnode[pmeshele[i].n[1]].pm += 0;
				pmeshnode[pmeshele[i].n[2]].pm += 0;
				break;
			case UPCOIL:
				pmeshele[i].miu = U0;
				pmeshele[i].miut = U0;
				pmeshnode[pmeshele[i].n[0]].I += 2. / 3 * PI*pmeshele[i].rc*pmeshele[i].AREA*J;
				pmeshnode[pmeshele[i].n[1]].I += 2. / 3 * PI*pmeshele[i].rc*pmeshele[i].AREA*J;
				pmeshnode[pmeshele[i].n[2]].I += 2. / 3 * PI*pmeshele[i].rc*pmeshele[i].AREA*J;
				pmeshnode[pmeshele[i].n[0]].pm += 0;
				pmeshnode[pmeshele[i].n[1]].pm += 0;
				pmeshnode[pmeshele[i].n[2]].pm += 0;
				break;
			case FIX:
			case MOV:
				pmeshele[i].miu = 1.0 * U0;
				pmeshele[i].miut = (nonlinear_miu)* U0;
				pmeshnode[pmeshele[i].n[0]].I += 0;
				pmeshnode[pmeshele[i].n[1]].I += 0;
				pmeshnode[pmeshele[i].n[2]].I += 0;
				pmeshnode[pmeshele[i].n[0]].pm += 0;
				pmeshnode[pmeshele[i].n[1]].pm += 0;
				pmeshnode[pmeshele[i].n[2]].pm += 0;
				break;
			case PM:
				int k;
				double r0, K, be[3];
				be[0] = 0; be[1] = 0; be[2] = 0;
				pmeshele[i].miu = U0; pmeshele[i].miut = U0;
				for (int j = 0; j < 3; j++) {
					k = j + 1;
					if (k == 3) {
						k = 0;
					}
					r0 = (pmeshnode[pmeshele[i].n[j]].x + pmeshnode[pmeshele[i].n[k]].x) / 2;
					K = 2 * PI*r0*Hc / 2 * (pmeshnode[pmeshele[i].n[j]].x - pmeshnode[pmeshele[i].n[k]].x);
					be[j] += K;
					be[k] += K;
				}
				pmeshnode[pmeshele[i].n[0]].pm += be[0];
				pmeshnode[pmeshele[i].n[1]].pm += be[1];
				pmeshnode[pmeshele[i].n[2]].pm += be[2];
				pmeshnode[pmeshele[i].n[0]].I += 0;
				pmeshnode[pmeshele[i].n[1]].I += 0;
				pmeshnode[pmeshele[i].n[2]].I += 0;
				break;
			default:
				break;
			}
		} else {
			//QMessageBox::warning(NULL, "Error:", "Error: reading domain points!",
			//	QMessageBox::Ok, QMessageBox::Ok);
			return 1;
		}
	}
	fclose(fp);
	return 0;
}


bool CFastFEMcore::StaticAxisymmetric() {
	//QCustomPlot *customPlot;
	//customPlot = ui->widget;
	double U0 = 1;// PI*4e-7;
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (pmeshele[i].domain == FIX || pmeshele[i].domain == MOV) {
			D34.push_back(i);
		}
	}
	
	//long long T1, T2;
	//timeval tv1;
	//T1 = tv1.tv_sec * 1000 * 1000 + tv1.tv_usec;
	//------------build C Matrix-----------------------------
	umat locs(2, 9 * num_ele); locs.zeros();
	vec vals = zeros<vec>(9 * num_ele);
	double ce[3][3];
	//double cc = 20;
	for (int i = 0; i < num_ele; i++) {
		if (pmeshele[i].domain == 3 || pmeshele[i].domain == 4) {
			if (pmeshele[i].y10 < 0) {
				ce[0][0] = 0;// -abs(pmeshele[i].y10) / pmeshele[i].miut;
			} else {
				ce[0][0] = abs(pmeshele[i].y10) / pmeshele[i].miut;
			}
			if (pmeshele[i].y20 < 0) {
				ce[1][1] = 0;//-abs(pmeshele[i].y20) / pmeshele[i].miut;
			} else {
				ce[1][1] = abs(pmeshele[i].y20) / pmeshele[i].miut;
			}
			if (pmeshele[i].y30 < 0) {
				ce[2][2] = 0;//-abs(pmeshele[i].y30) / pmeshele[i].miut;
			} else {
				ce[2][2] = abs(pmeshele[i].y30) / pmeshele[i].miut;
			}
			if (pmeshele[i].y12 > 0) {
				ce[0][0] += 0;//-abs(pmeshele[i].y12) / pmeshele[i].miut;
				ce[1][1] += 0;//-abs(pmeshele[i].y12) / pmeshele[i].miut;
				ce[0][1] = 0;//abs(pmeshele[i].y12) / pmeshele[i].miut;
			} else {
				ce[0][0] += abs(pmeshele[i].y12) / pmeshele[i].miut;
				ce[1][1] += abs(pmeshele[i].y12) / pmeshele[i].miut;
				ce[0][1] = -abs(pmeshele[i].y12) / pmeshele[i].miut;
			}
			if (pmeshele[i].y31 > 0) {
				ce[0][0] += 0;//-abs(pmeshele[i].y31) / pmeshele[i].miut;
				ce[2][2] += 0;//-abs(pmeshele[i].y31) / pmeshele[i].miut;
				ce[0][2] = 0;//abs(pmeshele[i].y31) / pmeshele[i].miut;
			} else {
				ce[0][0] += abs(pmeshele[i].y31) / pmeshele[i].miut;
				ce[2][2] += abs(pmeshele[i].y31) / pmeshele[i].miut;
				ce[0][2] = -abs(pmeshele[i].y31) / pmeshele[i].miut;
			}
			if (pmeshele[i].y23 > 0) {
				ce[1][1] += 0;//-abs(pmeshele[i].y23) / pmeshele[i].miut;
				ce[2][2] += 0;//-abs(pmeshele[i].y23) / pmeshele[i].miut;
				ce[1][2] = 0;//abs(pmeshele[i].y23) / pmeshele[i].miut;
			} else {
				ce[1][1] += abs(pmeshele[i].y23) / pmeshele[i].miut;
				ce[2][2] += abs(pmeshele[i].y23) / pmeshele[i].miut;
				ce[1][2] = -abs(pmeshele[i].y23) / pmeshele[i].miut;
			}
			ce[1][0] = ce[0][1];
			ce[2][0] = ce[0][2];
			ce[2][1] = ce[1][2];
		} else {
			ce[0][0] = pmeshele[i].y11 / pmeshele[i].miu;
			ce[0][1] = pmeshele[i].y12 / pmeshele[i].miu;
			ce[0][2] = pmeshele[i].y31 / pmeshele[i].miu;

			ce[1][0] = ce[0][1];
			ce[1][1] = pmeshele[i].y22 / pmeshele[i].miu;
			ce[1][2] = pmeshele[i].y23 / pmeshele[i].miu;

			ce[2][0] = ce[0][2];
			ce[2][1] = ce[1][2];
			ce[2][2] = pmeshele[i].y33 / pmeshele[i].miu;
		}

		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				locs(0, i * 9 + row * 3 + col) = pmeshele[i].n[row];
				locs(1, i * 9 + row * 3 + col) = pmeshele[i].n[col];
				vals(i * 9 + row * 3 + col) = ce[row][col];
			}
		}
	}
	
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts, num_pts, true, true);

	vec bbJz = zeros<vec>(num_pts);	vec b = zeros<vec>(num_pts);//b = bbJz + INL;
	vec A = zeros<vec>(num_pts);	vec A_old = A;
	vec INL = zeros<vec>(num_pts); vec bb = zeros<vec>(num_pts);
	vec rpm = zeros<vec>(num_pts);
	for (int i = 0; i < num_pts; i++) {
		bbJz(i) = pmeshnode[i].I;//set the current
		rpm(i) = pmeshnode[i].pm;
	}
	b = bbJz + INL + rpm;
	//---------------------superLU_MT---------------------------------------
	
	

	
	

	

	
	//---------------------superLU--end----------------------------------
	QVector<double> x1(D34.size()), y1(D34.size());
	QVector<double> x2(num_pts), y2(num_pts);
	for (int i = 0; i < D34.size(); ++i) {
		x1[i] = i;
		x2[i] = i;
	}
	////----------QCustomPlot setting---------------------------
	//QCPGraph *graph1 = customPlot->addGraph();
	//graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 3));
	//graph1->setPen(QPen(QColor(120, 120, 120), 2));
	//graph1->setLineStyle(QCPGraph::lsNone);
	//QCPGraph *graph2 = customPlot->addGraph();
	//graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 1.5), QBrush(Qt::red), 3));
	//graph2->setPen(QPen(QColor(120, 120, 120), 2));
	//graph2->setLineStyle(QCPGraph::lsNone);
	//customPlot->xAxis->setLabel("x");
	//customPlot->yAxis->setLabel("y");
	//customPlot->rescaleAxes(true);
	//customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	//---------openMP setting--------------------------------
	//int num_threads = omp_get_num_procs();
	//omp_set_num_threads(nprocs);
	//---------the main loop---------------------------------
	double Precision = 1e-6;
	int steps = 250;
	//clock_t t2 = clock();

	//int rowequ, colequ, notran;
	//int ldb = Bstore->lda;
	//Gstat_t   Gstat;
	int count;
	//int i,j;
	for (count = 0; count < steps; count++) {
		//#pragma omp parallel for
		//------update miu----------------
		for (int i = 0; i < D34.size(); i++) {
			//meshelement = &pmeshele[D34[i]];
			int d34 = D34[i];
			pmeshele[d34].Bx = (pmeshele[d34].Q[0] * pmeshnode[pmeshele[d34].n[0]].A
				+ pmeshele[d34].Q[1] * pmeshnode[pmeshele[d34].n[1]].A
				+ pmeshele[d34].Q[2] * pmeshnode[pmeshele[d34].n[2]].A) / 2. / pmeshele[d34].AREA;
			pmeshele[d34].By = (pmeshele[d34].P[0] * pmeshnode[pmeshele[d34].n[0]].A
				+ pmeshele[d34].P[1] * pmeshnode[pmeshele[d34].n[1]].A
				+ pmeshele[d34].P[2] * pmeshnode[pmeshele[d34].n[2]].A) / 2. / pmeshele[d34].AREA;
			pmeshele[d34].By += (pmeshnode[pmeshele[d34].n[0]].A
				+ pmeshnode[pmeshele[d34].n[1]].A
				+ pmeshnode[pmeshele[d34].n[2]].A) / 3. / pmeshele[d34].rc;
			pmeshele[d34].B = sqrt(pmeshele[d34].By*pmeshele[d34].By +
				pmeshele[d34].Bx*pmeshele[d34].Bx);
			pmeshele[d34].miu = HB(pmeshele[d34].B);// *miu0;
			y1[i] = pmeshele[d34].B;
			//y2[i] = pmeshele[D34[i]].miut / r;
		}
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *meshelement = &pmeshele[i];
			double rtmp;//do this to mark it as private
			rtmp = (meshelement->miu - meshelement->miut) / (meshelement->miu + meshelement->miut);

			meshelement->vr12 = (pmeshnode[meshelement->n[0]].A - pmeshnode[meshelement->n[1]].A) - meshelement->vi12;
			meshelement->vr23 = (pmeshnode[meshelement->n[1]].A - pmeshnode[meshelement->n[2]].A) - meshelement->vi23;
			meshelement->vr31 = (pmeshnode[meshelement->n[2]].A - pmeshnode[meshelement->n[0]].A) - meshelement->vi31;

			meshelement->vr10 = (pmeshnode[meshelement->n[0]].A - 0) - meshelement->vi10;
			meshelement->vr20 = (pmeshnode[meshelement->n[1]].A - 0) - meshelement->vi20;
			meshelement->vr30 = (pmeshnode[meshelement->n[2]].A - 0) - meshelement->vi30;

			if (meshelement->y10 > 0) {//the conductor
				meshelement->vi10 = rtmp*meshelement->vr10;
				INL(pmeshele[i].n[0]) += 2. / pmeshele[i].miut*pmeshele[i].vi10*abs(pmeshele[i].y10) / miu0;
			} else {//the controlled current source
				meshelement->vi10 = (pmeshnode[meshelement->n[0]].A - 0);//rtmp*meshelement->vr10;//
				INL(pmeshele[i].n[0]) += 1. / pmeshele[i].miu*pmeshele[i].vi10*abs(pmeshele[i].y10) / miu0;
			}
			if (meshelement->y20 > 0) {
				meshelement->vi20 = rtmp*meshelement->vr20;
				INL(pmeshele[i].n[1]) += 2. / pmeshele[i].miut*pmeshele[i].vi20*abs(pmeshele[i].y20) / miu0;
			} else {
				meshelement->vi20 = (pmeshnode[meshelement->n[1]].A - 0);//rtmp*meshelement->vr20;
				INL(pmeshele[i].n[1]) += 1. / pmeshele[i].miu*pmeshele[i].vi20*abs(pmeshele[i].y20) / miu0;
			}
			if (meshelement->y30 > 0) {
				meshelement->vi30 = rtmp*meshelement->vr30;
				INL(pmeshele[i].n[2]) += 2. / pmeshele[i].miut*pmeshele[i].vi30*abs(pmeshele[i].y30) / miu0;
			} else {
				meshelement->vi30 = (pmeshnode[meshelement->n[2]].A - 0);//rtmp*meshelement->vr30;
				INL(pmeshele[i].n[2]) += 1. / pmeshele[i].miu*pmeshele[i].vi30*abs(pmeshele[i].y30) / miu0;
			}
			if (meshelement->y12 < 0) {
				meshelement->vi12 = rtmp*meshelement->vr12;
				INL(pmeshele[i].n[1]) += -2. / pmeshele[i].miut*pmeshele[i].vi12*abs(pmeshele[i].y12) / miu0;
				INL(pmeshele[i].n[0]) += 2. / pmeshele[i].miut* pmeshele[i].vi12 *abs(pmeshele[i].y12) / miu0;
			} else {
				meshelement->vi12 = (pmeshnode[meshelement->n[0]].A - pmeshnode[meshelement->n[1]].A);//rtmp*meshelement->vr12;
				INL(pmeshele[i].n[1]) += -1. / pmeshele[i].miu*pmeshele[i].vi12*abs(pmeshele[i].y12) / miu0;
				INL(pmeshele[i].n[0]) += 1. / pmeshele[i].miu* pmeshele[i].vi12 *abs(pmeshele[i].y12) / miu0;
			}
			if (meshelement->y23 < 0) {
				meshelement->vi23 = rtmp*meshelement->vr23;
				INL(pmeshele[i].n[1]) += 2. / pmeshele[i].miut* pmeshele[i].vi23*abs(pmeshele[i].y23) / miu0;
				INL(pmeshele[i].n[2]) += -2. / pmeshele[i].miut*pmeshele[i].vi23*abs(pmeshele[i].y23) / miu0;
			} else {
				meshelement->vi23 = (pmeshnode[meshelement->n[1]].A - pmeshnode[meshelement->n[2]].A);//rtmp*meshelement->vr23;
				INL(pmeshele[i].n[1]) += 1. / pmeshele[i].miu* pmeshele[i].vi23*abs(pmeshele[i].y23) / miu0;
				INL(pmeshele[i].n[2]) += -1. / pmeshele[i].miu*pmeshele[i].vi23*abs(pmeshele[i].y23) / miu0;
			}
			if (meshelement->y31 < 0) {
				meshelement->vi31 = rtmp*meshelement->vr31;
				INL(pmeshele[i].n[2]) += 2. / pmeshele[i].miut* pmeshele[i].vi31*abs(pmeshele[i].y31) / miu0;
				INL(pmeshele[i].n[0]) += -2.0 / pmeshele[i].miut*pmeshele[i].vi31*abs(pmeshele[i].y31) / miu0;
			} else {
				meshelement->vi31 = (pmeshnode[meshelement->n[2]].A - pmeshnode[meshelement->n[0]].A);//rtmp*meshelement->vr31;
				INL(pmeshele[i].n[2]) += 1. / pmeshele[i].miu* pmeshele[i].vi31*abs(pmeshele[i].y31) / miu0;
				INL(pmeshele[i].n[0]) += -1.0 / pmeshele[i].miu*pmeshele[i].vi31*abs(pmeshele[i].y31) / miu0;
			}


		}
		b = bbJz + INL + rpm; 
		A_old = A;

		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		


		/*double error = 0;
		for (int i = 0; i < num_pts; i++) {
			pmeshnode[i].A = sol[i] * miu0;
			A(i) = sol[i] * miu0;
		}*/

		double error = norm((A_old - A), 2) / norm(A, 2);
		if (error < Precision) {
			break;
		}
		INL.zeros();

		//this->setWindowTitle(QString::number(t1));
		//graph1->setData(x1, y1);

		//customPlot->rescaleAxes(true);
		//customPlot->replot();//necessary
		//if (QString::number(error) == "nan") {
		//	return;
		//} else {
		//	this->setWindowTitle(QString::number(count));
		//}
	}
	for (int i = 0; i < num_ele; i++) {
		pmeshele[i].Bx = (pmeshele[i].Q[0] * pmeshnode[pmeshele[i].n[0]].A
			+ pmeshele[i].Q[1] * pmeshnode[pmeshele[i].n[1]].A
			+ pmeshele[i].Q[2] * pmeshnode[pmeshele[i].n[2]].A) / 2. / pmeshele[i].AREA;
		pmeshele[i].By = (pmeshele[i].P[0] * pmeshnode[pmeshele[i].n[0]].A
			+ pmeshele[i].P[1] * pmeshnode[pmeshele[i].n[1]].A
			+ pmeshele[i].P[2] * pmeshnode[pmeshele[i].n[2]].A) / 2. / pmeshele[i].AREA;
		pmeshele[i].By += (pmeshnode[pmeshele[i].n[0]].A
			+ pmeshnode[pmeshele[i].n[1]].A
			+ pmeshnode[pmeshele[i].n[2]].A) / 3. / pmeshele[i].rc;
		pmeshele[i].B = sqrt(pmeshele[i].By*pmeshele[i].By +
			pmeshele[i].Bx*pmeshele[i].Bx);
	}
	//A.save("A.txt", arma_ascii);
	/*using another loop*/
	//gettimeofday(&tv1,NULL);
	//T2 = tv1.tv_sec * 1000 * 1000 + tv1.tv_usec;
	//this->setWindowTitle(QString::number(T2 - T1));
	//this->setWindowTitle(QString::number(count) + "," + QString::number(T2 - T1));
	//graph1->setData(x1, y1);
	//graph2->setData(x2, y2);
	//customPlot->rescaleAxes(true);
	//customPlot->replot();//necessary
	//customPlot->removeGraph(graph1);
	//customPlot->removeGraph(graph2);
	
	//customPlot->removeGraph(graph1);
	return true;
}

double CFastFEMcore::HB(double B) {
	double miu;
	double hint;
	double h[11] = { 0, 200, 500, 1000, 2500, 5000, 10000, 15000, 100000, 500000, 500000 };
	double b[11] = { 0, 1.2, 1.4, 1.5, 1.62, 1.71, 1.85, 1.851, 1.8511, 5 , 5.5};
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
