#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "datatype.h"
#include "spline.h"
#include <cmath>
#include<QtAlgorithms>
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo> 
#include <vector>
#include <time.h>
#include <omp.h>
#include <QMessageBox>
#include "slu_mt_ddefs.h"

using namespace std;
using namespace arma;


#define r 0.1
#define PI 3.14159265358979323846

const double ste = 20;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    QAction *start = new QAction("NR", this);
    addAction(start);
    connect(start, SIGNAL(triggered()), this, SLOT(NRFEM()));
    setContextMenuPolicy(Qt::ActionsContextMenu);

}

MainWindow::~MainWindow() {
    delete ui;
}
typedef struct{
    int index;
    double x;
    double y;
} ele;
bool compareele(const ele ele1,const ele ele2){
    if(ele1.x < ele2.x){
        return true;
    }else{
        return false;
    }
}
double mymax(double a,double b){
    return ((a > b) ? a : b);
}
double mymin(double a, double b){
    return ((a < b) ? a : b);
}
void MainWindow::NRFEM() {
    /*read coarse mesh data from file*/
    double U0 = PI*4e-7;
    char  coarseFN[] = "E:\\Projects\\axispmmodel\\mesh20.mphtxt";
    CNode * cmeshnode = NULL;	int cnum_pts = 0;
    CElement * cmeshele = NULL;	int cnum_ele = 0;
    readDataFile(coarseFN, cnum_pts, cnum_ele, &cmeshnode, &cmeshele);
    /*using TLM get a coarse miu value*/
    NRCalculate(cnum_pts, cnum_ele, cmeshnode, cmeshele);
	CalcForce(cnum_pts, cnum_ele, cmeshnode, cmeshele);
	
	///*read fine mesh data from file*/
 //   char fineFN[] = "E:\\Projects\\axispmmodel\\fine.mphtxt";
 //   CNode * fmeshnode = NULL;	int fnum_pts = 0;
 //   CElement * fmeshele = NULL;	int fnum_ele = 0;
 //   readDataFile(fineFN, fnum_pts, fnum_ele, &fmeshnode, &fmeshele);
 //   /*set the miu initial value */
 //   double minx,maxx,miny,maxy;
 //   QVector <ele> d34;
 //   ele eletmp;
 //   for(int i = 0; i < fnum_ele;i++){
 //       if(fmeshele[i].domain == 3 || fmeshele[i].domain == 4){
 //           eletmp.index = i;
 //           eletmp.x = fmeshele[i].rc;
 //           eletmp.y = fmeshele[i].zc;
 //           d34.push_back(eletmp);
 //       }
 //   }
 //   qSort(d34.begin(),d34.end(),compareele);//ascend

 //   for (int i = 0; i < cnum_ele; i++) {
 //       if (cmeshele[i].domain == 3 || cmeshele[i].domain == 4) {
 //           /*find the minx max miny maxy in the ele*/
 //           double a,b,c;
 //           a = cmeshnode[cmeshele[i].n[0]].x;
 //           b = cmeshnode[cmeshele[i].n[1]].x;
 //           c = cmeshnode[cmeshele[i].n[2]].x;
 //           maxx = mymax(mymax(a,b),c);
 //           minx = mymin(mymin(a,b),c);
 //           a = cmeshnode[cmeshele[i].n[0]].y;
 //           b = cmeshnode[cmeshele[i].n[1]].y;
 //           c = cmeshnode[cmeshele[i].n[2]].y;
 //           maxy = mymax(mymax(a,b),c);
 //           miny = mymin(mymin(a,b),c);
 //           for(int j = 0; j < d34.size(); j++){
 //               /*find the first x between minx~maxx*/
 //               if(d34[j].x > maxx){
 //                   break;
 //               }else if(d34[j].x > minx){
 //                   if(d34[j].y > miny && d34[j].y < maxy){
 //                       fmeshele[d34[j].index].miu = cmeshele[i].miu * r;
 //                       d34.remove(j);
 //                       j -= 1;
 //                   }
 //               }
 //           }
 //       }
 //   }
 //   /*recall the TLM*/
 //   NRCalculate(fnum_pts, fnum_ele, fmeshnode, fmeshele);
 //   /*free the space*/
 //   free(cmeshele); free(cmeshnode);
 //   free(fmeshele); free(fmeshnode);
	
}

void MainWindow::readDataFile(char *fileName, int & num_pts, int & num_ele, CNode **ppmeshnode, CElement ** ppmeshele) {
    char ch[256];
    int re;
    double U0 = PI*4e-7;
    double nonlinear_miu = 1;
	double Scoil = 14.4e-3*12.7e-3;// 11.687e-3 * 14.7e-3;
	double I = 1;// 12 / 3.2;
	double N = 1420;// 518;
    double J = I * N / Scoil;
	double Jdown = 1 * 1250 / (12.4*18.4*1e-6);
	double Hc = 883310;
    //------------open file----------------------------------
    FILE * fp = NULL;
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        //printf("Error: opening file!\n");
        QMessageBox::warning(NULL, "Error:", "Error: opening file!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }
    //--------------Read the head-----------------------------
    for (int i = 0; i < 18; i++) {
        fgets(ch, 256, fp);
        //printf("%s", ch);
    }
    //-----------------mesh point-----------------------------
    CNode *meshnode = NULL;

    re = fscanf(fp, "%d # number of mesh points\n", &num_pts);
    if (re == 1) {
        //printf("num_pts:%d\n", num_pts);
        *ppmeshnode = (CNode*)calloc(num_pts, sizeof(CNode));
        meshnode = *ppmeshnode;
		for (int i = 0; i < num_pts; i++) {
			meshnode[i].I = 0;
			meshnode[i].pm = 0;
		}
    } else {//printf("Error: reading num_pts!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading num_pts!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }

    int pts_ind;//the beginning of the points index
    re = fscanf(fp, "%d # lowest mesh point index\n", &pts_ind);
    if (re == 1) {//printf("pts_ind:%d\n", pts_ind);
        ;
    } else {//printf("Error: reading pts_ind!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading pts_ind!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }

    fgets(ch, 256, fp);
    //x = (double*)calloc(num_pts, sizeof(double));
    //y = (double*)calloc(num_pts, sizeof(double));
    for (int i = pts_ind; i < num_pts; i++) {
        re = fscanf(fp, "%lf %lf \n", &meshnode[i].x, &meshnode[i].y);
        if (re == 2) {//printf("%lf %lf\n", x[i], y[i]);
            ;
        } else {//printf("Error:loading mesh point!\n");
            QMessageBox::warning(NULL, "Error:", "Error: reading mesh point!",
                                 QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        //getchar();
    }
    //---------------vertex-------------------------------
    for (int i = 0; i < 7; i++)
        fgets(ch, 256, fp);
    int num_vtx_ns, num_vtx_ele;
    re = fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns);
    if (re == 1) {//printf("num_bdr_ns:%d\n", num_vtx_ns);
        ;
    } else {	//printf("Error: reading num_vtx_ns!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ns!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }

    re = fscanf(fp, "%d # number of elements\n", &num_vtx_ele);
    if (re == 1) {
        //printf("num_vtx_ele:%d\n", num_vtx_ele);
        ;
    } else {//printf("Error: reading num_vtx_ele!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ele!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }
    fgets(ch, 256, fp);

    int *vtx;
    vtx = (int*)calloc(num_vtx_ele, sizeof(int));
    for (int i = 0; i < num_vtx_ele; i++) {
        re = fscanf(fp, "%d \n", vtx + i);
        if (re == 1) {//printf("%d %d\n", p1[i], p2[i]);
            ;
        } else {	//printf("Error: loading vertex condition!\n");
            QMessageBox::warning(NULL, "Error:", "Error: reading vertex condition!",
                                 QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        //getchar();
    }
    //---------------vertex-------------------------------
    int num_vtx_ele2;
    re = fscanf(fp, "%d # number of geometric entity indices\n", &num_vtx_ele2);
    fgets(ch, 256, fp);
    int *vtx2;
    vtx2 = (int*)calloc(num_vtx_ele2, sizeof(int));
    for (int i = 0; i < num_vtx_ele2; i++) {
        re = fscanf(fp, "%d \n", vtx2 + i);
        if (re == 1) {//printf("%d %d\n", p1[i], p2[i]);
            ;
        } else {	//printf("Error: loading vertex condition!\n");
            QMessageBox::warning(NULL, "Error:", "Error: reading vertex condition!",
                                 QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        //getchar();
    }
    //--------------boundary--------------------------------
    for (int i = 0; i < 5; i++)
        fgets(ch, 256, fp);
    int num_bdr_ns, num_bdr_ele;//number of nodes per element;number of elements
    re = fscanf(fp, "%d # number of nodes per element\n", &num_bdr_ns);
    if (re == 1) {
        ;
    }//printf("num_bdr_ns:%d\n", num_bdr_ns);
    else {//printf("Error: reading num_bdr_ns!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading num_bdr_ns!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }

    re = fscanf(fp, "%d # number of elements\n", &num_bdr_ele);
    if (re == 1) {
        ;
    }	//printf("num_bdr_ele:%d\n", num_bdr_ele);
    else {	//printf("Error: reading num_bdr_ele!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading num_bdr_ele!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }
    fgets(ch, 256, fp);

    int *p1, *p2;
    p1 = (int*)calloc(num_bdr_ele, sizeof(int));
    p2 = (int*)calloc(num_bdr_ele, sizeof(int));
    for (int i = 0; i < num_bdr_ele; i++) {
        re = fscanf(fp, "%d %d\n", p1 + i, p2 + i);
        if (re == 2) {
            //printf("%d %d\n", p1[i], p2[i]);
            //-----process the A=0 boundary--------------------------
            /*if (abs(meshnode[p1[i]].length() - 0.05) < 5e-3 || abs(meshnode[p1[i]].x) < 5e-5) {
                meshnode[p1[i]].bdr = 1;
            }*/
        } else {	//printf("Error: loading boundary condition!\n");
            QMessageBox::warning(NULL, "Error:", "Error: reading boundary condition!",
                                 QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        //getchar();
    }
    //---------------entity----------------------------------
    int num_entity;
    re = fscanf(fp, "%d # number of geometric entity indices\n", &num_entity);
    fgets(ch, 256, fp);
    int * entity;
    entity = (int*)calloc(num_entity, sizeof(int));
    for (int i = 0; i < num_entity; i++) {
        re = fscanf(fp, "%d \n", entity + i);
        if (re == 1) {//printf("%d %d\n", p1[i], p2[i]);
            ;
        } else {	//printf("Error: loading boundary condition!\n");
            QMessageBox::warning(NULL, "Error:", "Error: reading boundary condition!",
                                 QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        //getchar();
    }
    //----------------elements------------------------------
    for (int i = 0; i < 5; i++)
        fgets(ch, 256, fp);
    int ns_per_ele;//number of nodes per element;number of elements
    re = fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele);
    if (re == 1) {	//printf("ns_per_ele:%d\n", ns_per_ele);
        ;
    } else {	//printf("Error: reading ns_per_ele!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading ns_per_ele!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }
    CElement * meshele = NULL;
    re = fscanf(fp, "%d # number of elements\n", &num_ele);
    if (re == 1) {
        //printf("num_ele:%d\n", num_ele);
        *ppmeshele = (CElement*)calloc(num_ele, sizeof(CElement));
        meshele = *ppmeshele;
    } else {	//printf("Error: reading num_ele!\n");
        QMessageBox::warning(NULL, "Error:", "Error: reading num_ele!",
                             QMessageBox::Ok, QMessageBox::Ok);
        return;
    }
    fgets(ch, 256, fp);

    for (int i = 0; i < num_ele; i++) {
        re = fscanf(fp, "%d %d %d \n", &meshele[i].n[0], &meshele[i].n[1], &meshele[i].n[2]);
        if (re == 3) {	//printf("%d %d %d \n", pi[i], pj[i], pk[i]);
            meshele[i].P[0] = meshnode[meshele[i].n[1]].y - meshnode[meshele[i].n[2]].y;
            meshele[i].P[1] = meshnode[meshele[i].n[2]].y - meshnode[meshele[i].n[0]].y;
            meshele[i].P[2] = meshnode[meshele[i].n[0]].y - meshnode[meshele[i].n[1]].y;

            meshele[i].Q[0] = meshnode[meshele[i].n[2]].x - meshnode[meshele[i].n[1]].x;
            meshele[i].Q[1] = meshnode[meshele[i].n[0]].x - meshnode[meshele[i].n[2]].x;
            meshele[i].Q[2] = meshnode[meshele[i].n[1]].x - meshnode[meshele[i].n[0]].x;

            meshele[i].AREA = 0.5*abs(meshele[i].P[1] * meshele[i].Q[2] - meshele[i].Q[1] * meshele[i].P[2]);
            meshele[i].rc = (meshnode[meshele[i].n[0]].x +
                    meshnode[meshele[i].n[1]].x +
                    meshnode[meshele[i].n[2]].x) / 3;
            meshele[i].y12 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[0] * meshele[i].P[1] +
                    meshele[i].Q[0] * meshele[i].Q[1]);
            meshele[i].y12 += (meshele[i].P[0] + meshele[i].P[1]) / 6;
            meshele[i].y12 += meshele[i].AREA / 9 / meshele[i].rc;
            meshele[i].y12 *= 2 * PI; //meshele[i].y12 = abs(meshele[i].y12);

            meshele[i].y23 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[1] * meshele[i].P[2] +
                    meshele[i].Q[1] * meshele[i].Q[2]);
            meshele[i].y23 += (meshele[i].P[1] + meshele[i].P[2]) / 6;
            meshele[i].y23 += meshele[i].AREA / 9 / meshele[i].rc;
            meshele[i].y23 *= 2 * PI; //meshele[i].y23 = abs(meshele[i].y23);

            meshele[i].y13 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[2] * meshele[i].P[0] +
                    meshele[i].Q[2] * meshele[i].Q[0]);
            meshele[i].y13 += (meshele[i].P[2] + meshele[i].P[0]) / 6;
            meshele[i].y13 += meshele[i].AREA / 9 / meshele[i].rc;
            meshele[i].y13 *= 2 * PI; //meshele[i].y31 = abs(meshele[i].y31);

            meshele[i].y11 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[0] * meshele[i].P[0] +
                    meshele[i].Q[0] * meshele[i].Q[0]);
            meshele[i].y11 += (meshele[i].P[0] + meshele[i].P[0]) / 6;
            meshele[i].y11 += meshele[i].AREA / 9 / meshele[i].rc;
            meshele[i].y11 *= 2 * PI;

            meshele[i].y22 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[1] * meshele[i].P[1] +
                    meshele[i].Q[1] * meshele[i].Q[1]);
            meshele[i].y22 += (meshele[i].P[1] + meshele[i].P[1]) / 6;
            meshele[i].y22 += meshele[i].AREA / 9 / meshele[i].rc;
            meshele[i].y22 *= 2 * PI;

            meshele[i].y33 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[2] * meshele[i].P[2] +
                    meshele[i].Q[2] * meshele[i].Q[2]);
            meshele[i].y33 += (meshele[i].P[2] + meshele[i].P[2]) / 6;
            meshele[i].y33 += meshele[i].AREA / 9 / meshele[i].rc;
            meshele[i].y33 *= 2 * PI;
        } else {	//printf("Error: loading elements points!\n");
            QMessageBox::warning(NULL, "Error:", "Error: reading elements points!",
                                 QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        //getchar();
    }
    //---------------Domain----------------------------------
    int num_domain;
    re = fscanf(fp, "%d # number of geometric entity indices\n", &num_domain);
    fgets(ch, 256, fp);
    //int *domain;
    //domain = (int*)calloc(num_domain, sizeof(int));

    for (int i = 0; i < num_domain; i++) {
        re = fscanf(fp, "%d \n", &meshele[i].domain);
        if (re == 1) {	//printf("%d %d %d \n", pi[i], pj[i], pk[i]);
            switch (meshele[i].domain) {
            case 1:
            case 2:
                meshele[i].miu = U0;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
				meshnode[meshele[i].n[0]].I += 0;
				meshnode[meshele[i].n[1]].I += 0;
				meshnode[meshele[i].n[2]].I += 0;
                break;
            case 6:
				meshele[i].miu = U0;
				meshnode[meshele[i].n[0]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*Jdown;
				meshnode[meshele[i].n[1]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*Jdown;
				meshnode[meshele[i].n[2]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*Jdown;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
				break;
			case 7:
                meshele[i].miu = U0;
                meshnode[meshele[i].n[0]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*J;
                meshnode[meshele[i].n[1]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*J;
                meshnode[meshele[i].n[2]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*J;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
                break;
            case 3:
            case 4:
                meshele[i].miu = 1 * U0;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
				meshnode[meshele[i].n[0]].I += 0;
				meshnode[meshele[i].n[1]].I += 0;
				meshnode[meshele[i].n[2]].I += 0;
                break;
			case 5:
				int k;
				double r0, K, be[3];
				be[0] = 0; be[1] = 0; be[2] = 0;
				meshele[i].miu = U0; 
				for (int j = 0; j < 3; j++) {
					k = j + 1;
					if (k == 3) {
						k = 0;
					}
					r0 = (meshnode[meshele[i].n[j]].x + meshnode[meshele[i].n[k]].x) / 2;
					K = 2 * PI*r0*Hc / 2 * (meshnode[meshele[i].n[j]].x - meshnode[meshele[i].n[k]].x);
					be[j] += K;
					be[k] += K;
				}
				meshnode[meshele[i].n[0]].pm += be[0];
				meshnode[meshele[i].n[1]].pm += be[1];
				meshnode[meshele[i].n[2]].pm += be[2];
				meshnode[meshele[i].n[0]].I += 0;
				meshnode[meshele[i].n[1]].I += 0;
				meshnode[meshele[i].n[2]].I += 0;
				break;
            default:
                break;
            }
        } else {	//printf("Error: loading domain points!\n");
            QMessageBox::warning(NULL, "Error:", "Error: reading domain points!",
                                 QMessageBox::Ok, QMessageBox::Ok);
            return;
        }
        //getchar();

    }
    fclose(fp);
}

void MainWindow::NRCalculate(int & num_pts, int & num_ele, CNode *m_n, CElement * m_l){
    QCustomPlot *customPlot;
    double U0 = 4 * PI*1e-7;
    customPlot = ui->qcustomplot;
    std::vector <int> D34;
    for (int i = 0; i < num_ele; i++) {
        if (m_l[i].domain == 3 || m_l[i].domain == 4) {
            D34.push_back(i);
        }
    }
    long long T1,T2;
    timeval tv1;
    //gettimeofday(&tv1,NULL);
    T1 = tv1.tv_sec * 1000 * 1000 + tv1.tv_usec;

    //------------get C Matrix-----------------------------
    umat locs(2, 9 * num_ele);
    vec vals(9 * num_ele);

	QVector<double> x1(D34.size()), y1(D34.size()), y2(D34.size());
    for (int i = 0; i < D34.size(); ++i) {
        x1[i] = i;
    }
    double ce[3][3];
    double err = 1e-6;
    double maxur = 1;
    vec b = zeros<vec>(num_pts);
    vec bbjz = zeros<vec>(num_pts);
	vec rpm = zeros<vec>(num_pts);
    vec A = ones<vec>(num_pts);
    vec A_old = zeros<vec>(num_pts);
    vec bnr = zeros<vec>(D34.size());
    vec miunr = zeros<vec>(D34.size());
    for (int i = 0; i < num_pts; i++) {
        bbjz(i) = m_n[i].I;
		rpm(i) = m_n[i].pm;
    }

    QCPGraph *graph1 = customPlot->addGraph();
	QCPGraph *graph2 = customPlot->addGraph();	
	graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 1.5), QBrush(Qt::red), 3));
	graph2->setPen(QPen(QColor(120, 120, 120), 2));
	graph2->setLineStyle(QCPGraph::lsNone);
    graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 3));
    graph1->setPen(QPen(QColor(120, 120, 120), 2));
    graph1->setLineStyle(QCPGraph::lsNone);
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");
    // set axes ranges, so we see all data:
    //customPlot->xAxis->setRange(0, D34.size());
    //customPlot->yAxis->setRange(-1e-4, 1e-4);
    customPlot->rescaleAxes(true);
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    //while (maxur > err)

    //---------------------superLU_MT---------------------------------------
    SuperMatrix sluA, L, U;
    SuperMatrix sluB, sluX;
    NCformat    *Astore;
    SCPformat   *Lstore;DNformat	   *Bstore;
    NCPformat   *Ustore;
    int_t         nprocs;
    fact_t      fact;
    trans_t     trans;
    yes_no_t    refact, usepr;
    equed_t     equed;
    double      *a;double t1;
    int_t         *asub, *xa;
    int_t         *perm_c; /* column permutation vector */
    int_t         *perm_r; /* row permutations from partial pivoting */
    void        *work;
    superlumt_options_t superlumt_options;
    int_t         info, lwork, nrhs, ldx, panel_size, relax;
    int_t         m, n, nnz, permc_spec;
    double      *rhsb, *rhsx, *xact;
    double      *R, *C;
    double      *ferr, *berr;
    double      u, drop_tol, rpg, rcond;
    superlu_memusage_t superlu_memusage;
    void parse_command_line();

    /* Default parameters to control factorization. */
    nprocs = 1;
    fact  = EQUILIBRATE;
    trans = NOTRANS;
    equed = NOEQUIL;
    refact= NO;
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    u     = 1.0;
    usepr = NO;
    drop_tol = 0.0;
    lwork = 0;
    nrhs  = 1;

    if ( lwork > 0 ) {
        work = SUPERLU_MALLOC(lwork);
        printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
        if ( !work ) {
            SUPERLU_ABORT("DLINSOLX: cannot allocate work[]");
        }
    }


    /*   * Solve the system and compute the condition number
         * and error bounds using pdgssvx.
         */

    double error = 5;
    double precision = 1e-6;
    double ratio = 0.5;
    clock_t t2 = clock();
    //timeval tv1;
    //gettimeofday(&tv1,NULL);
    T1 = tv1.tv_sec * 1000 * 1000 + tv1.tv_usec;
    /*the loop*/
    int count;
    for ( count = 0; count < 100; count++) {
        for (int i = 0; i < num_ele; i++) {
            ce[0][0] = m_l[i].y11 / m_l[i].miu;
            ce[0][1] = m_l[i].y12 / m_l[i].miu;
            ce[0][2] = m_l[i].y13 / m_l[i].miu;

            ce[1][0] = m_l[i].y12 / m_l[i].miu;
            ce[1][1] = m_l[i].y22 / m_l[i].miu;
            ce[1][2] = m_l[i].y23 / m_l[i].miu;

            ce[2][0] = m_l[i].y13 / m_l[i].miu;
            ce[2][1] = m_l[i].y23 / m_l[i].miu;
            ce[2][2] = m_l[i].y33 / m_l[i].miu;

            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    locs(0, i * 9 + row * 3 + col) = m_l[i].n[row];
                    locs(1, i * 9 + row * 3 + col) = m_l[i].n[col];
                    vals(i * 9 + row * 3 + col) = ce[row][col];
                }
            }
        }
        sp_mat X;
        X = sp_mat(true, locs, vals, num_pts, num_pts, true, true);
        /* Read the matrix in Harwell-Boeing format. */
        /*create A*/
        /* create matrix A in Harwell-Boeing format.*/
        m = num_pts; n = num_pts; nnz = X.n_nonzero;
        a = const_cast<double *>(X.values);
        asub = (int*)const_cast<u32 *>(X.row_indices);
        xa = (int*)const_cast<u32 *>(X.col_ptrs);
        dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

        //------create B and X-------------------
        if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
        dCreate_Dense_Matrix(&sluX, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
        b = bbjz + rpm;
        rhsb = const_cast<double*>(b.mem);
        dCreate_Dense_Matrix(&sluB, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
        Bstore = (DNformat*)sluB.Store;
        double *Bmat = (double*)Bstore->nzval;
        double *Xmat = (double*)((DNformat*)sluX.Store)->nzval;
        ldx = m;
        if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
        if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
        if (!(R = (double *) SUPERLU_MALLOC(sluA.nrow * sizeof(double))))
            SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
        if ( !(C = (double *) SUPERLU_MALLOC(sluA.ncol * sizeof(double))) )
            SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
        if ( !(ferr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
            SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
        if ( !(berr = (double *) SUPERLU_MALLOC(nrhs * sizeof(double))) )
            SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");

        /*   * Get column permutation vector perm_c[], according to permc_spec:
             *   permc_spec = 0: natural ordering
             *   permc_spec = 1: minimum degree ordering on structure of A'*A
             *   permc_spec = 2: minimum degree ordering on structure of A'+A
             *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
             */
        permc_spec = 1;
        get_perm_c(permc_spec, &sluA, perm_c);

        superlumt_options.nprocs = nprocs;
        superlumt_options.fact = EQUILIBRATE;
        superlumt_options.trans = NOTRANS;
        superlumt_options.refact = NO;
        superlumt_options.panel_size = panel_size;
        superlumt_options.relax = relax;
        superlumt_options.usepr = NO;
        superlumt_options.drop_tol = drop_tol;
        superlumt_options.diag_pivot_thresh = u;
        superlumt_options.SymmetricMode = NO;
        superlumt_options.PrintStat = NO;
        superlumt_options.perm_c = perm_c;
        superlumt_options.perm_r = perm_r;
        superlumt_options.work = work;
        superlumt_options.lwork = 0;
        if ( !(superlumt_options.etree = intMalloc(n)) )
            SUPERLU_ABORT("Malloc fails for etree[].");
        if ( !(superlumt_options.colcnt_h = intMalloc(n)) )
            SUPERLU_ABORT("Malloc fails for colcnt_h[].");
        if ( !(superlumt_options.part_super_h = intMalloc(n)) )
            SUPERLU_ABORT("Malloc fails for colcnt_h[].");

        /*   * Solve the system and compute the condition number
             * and error bounds using pdgssvx.
             */
        t1 = SuperLU_timer_();
        pdgssvx(nprocs, &superlumt_options, &sluA, perm_c, perm_r,
                &equed, R, C, &L, &U, &sluB, &sluX, &rpg, &rcond,
                ferr, berr, &superlu_memusage, &info);
        t1 = SuperLU_timer_() - t1;
        //printf("%lf\n",t1);
        /* This is how you could access the solution matrix. */
        double *sol = NULL;
        A_old = A;
        if (info == 0 || info == n + 1) {
            sol = (double*)((DNformat*)sluX.Store)->nzval;
            for (int i = 0; i < num_pts; i++) {
                m_n[i].A = sol[i];
                A(i) = sol[i];
            }
        } else if (info > 0 && lwork == -1) {

        }

        error = norm((A_old - A),2) / norm(A,2);
        if (error < precision) {
            break;
        }
        SUPERLU_FREE(rhsx);
        //        SUPERLU_FREE(etree);
        SUPERLU_FREE(perm_r);
        SUPERLU_FREE(perm_c);
        SUPERLU_FREE(R);
        SUPERLU_FREE(C);
        SUPERLU_FREE(ferr);
        SUPERLU_FREE(berr);
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
        /*-------------superLU end-----------*/
        for (int i = 0; i < D34.size(); i++) {
            m_l[D34[i]].Bx = (m_l[D34[i]].Q[0] * m_n[m_l[D34[i]].n[0]].A
                    + m_l[D34[i]].Q[1] * m_n[m_l[D34[i]].n[1]].A
                    + m_l[D34[i]].Q[2] * m_n[m_l[D34[i]].n[2]].A) / 2 / m_l[D34[i]].AREA;
            m_l[D34[i]].By = (m_l[D34[i]].P[0] * m_n[m_l[D34[i]].n[0]].A
                    + m_l[D34[i]].P[1] * m_n[m_l[D34[i]].n[1]].A
                    + m_l[D34[i]].P[2] * m_n[m_l[D34[i]].n[2]].A) / 2 / m_l[D34[i]].AREA;
            m_l[D34[i]].By += (m_n[m_l[D34[i]].n[0]].A
                    + m_n[m_l[D34[i]].n[1]].A
                    + m_n[m_l[D34[i]].n[2]].A) / 3 / m_l[D34[i]].rc;
            m_l[D34[i]].B = sqrt(m_l[D34[i]].By*m_l[D34[i]].By +
                    m_l[D34[i]].Bx*m_l[D34[i]].Bx);
            m_l[D34[i]].miu_t = HB(m_l[D34[i]].B)*U0;
            m_l[D34[i]].miu += ratio*(m_l[D34[i]].miu_t - m_l[D34[i]].miu);

			y1[i] = m_l[D34[i]].B;
			/*if (count == 0) {
				y2[i] = y1[i];
			}*/
            bnr(i) = m_l[D34[i]].B;
            miunr(i) = m_l[D34[i]].miu;
        }
		
        graph1->setData(x1, y1);
		//graph2->setData(x1, y2);
        customPlot->rescaleAxes(true);
        customPlot->replot();//necessary
        this->setWindowTitle(QString::number(count));
    }
//    A.save("NR_A.txt",arma_ascii);
//    bnr.save("NR_B.txt",arma_ascii);
//    miunr.save("NR_miu.txt",arma_ascii);
    //gettimeofday(&tv1,NULL);
    T2 = tv1.tv_sec * 1000 * 1000 + tv1.tv_usec;
    clock_t t3 = clock();
    graph1->setData(x1, y1);
    customPlot->rescaleAxes(true);
    customPlot->replot();//necessary
    customPlot->removeGraph(graph1);
	customPlot->removeGraph(graph2);
    //this->setWindowTitle(QString::number(T2 - T1));
    this->setWindowTitle(QString::number(count)+","+QString::number(T2 - T1));
}

double MainWindow::HB(double B) {
    /*double miu;
    if (B <= 0.6) {
        miu = 1 / 1.7e-4;
    } else {
        miu = 1 / (1e-4*(B - 0.6)*(B - 0.6) + 1.7e-4);
    }*/
	double miu;
	double hint;
	double h[10] = { 0, 200, 500, 1000, 2500, 5000, 10000, 15000, 100000,500000 };
	double b[10] = { 0, 1.2, 1.4, 1.5, 1.62, 1.71, 1.85, 1.851, 1.8511 ,5};
	std::vector<double> BB(b, b + 9);
	std::vector<double> HH(h, h + 9);
	tk::spline s;
	s.set_points(BB, HH,false);
	hint = s(B);
	if (B > 4) {
		B = 3;
		miu = B / (1e4 + 5e6*(B - 1.85)) / (4*PI*1e-7);
	} else {
		miu = B / hint / (4 * PI*1e-7);
	}	
	
    return miu;
}
void MainWindow::CalcForce(int & num_pts, int & num_ele, CNode *m_n, CElement * m_l) {
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

	QCustomPlot *customPlot;
	customPlot = ui->qcustomplot;
	customPlot->xAxis->setLabel("x");
	customPlot->xAxis->setRange(0, 0.09);
	customPlot->yAxis->setLabel("y");
	customPlot->yAxis->setRange(-0.09, 0.09);
	customPlot->yAxis->setScaleRatio(ui->qcustomplot->xAxis, 1.0);
	customPlot->yAxis->setAutoTickStep(false);
	customPlot->yAxis->setTickStep(ui->qcustomplot->xAxis->tickStep());
	ui->qcustomplot->yAxis->setAutoTickLabels(true);
	customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

	QCPCurve *newCurve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);

	QCPGraph *graph1 = customPlot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 5), QBrush(Qt::black), 5));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	QCPGraph *graph2 = customPlot->addGraph();
	graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, QPen(Qt::red, 5), QBrush(Qt::red), 5));
	graph2->setPen(QPen(QColor(120, 120, 120), 2));
	graph2->setLineStyle(QCPGraph::lsNone);
	QVector <double> gpx1, gpx2, gpy1, gpy2;
	QColor cc[7];
	cc[0] = QColor(0, 0, 0);
	cc[1] = QColor(0, 120, 0);
	cc[2] = QColor(0, 0, 120);
	cc[3] = QColor(120, 0, 0);
	cc[4] = QColor(120, 120, 0);
	cc[5] = QColor(0, 120, 120);
	cc[6] = QColor(120, 0, 120);

	//FILE * fp = NULL;
	//fp = fopen("E:\\index.txt", "w+");//delete exist, read and write
	for (int i = 0; i < num_ele; i++) {
		QCPCurve *newCurve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
		QVector <double> x1(4);
		QVector <double> y1(4);
		int BCx = 0;
		int BCy = 0;
		for (int j = 0; j < 3; j++) {
			x1[j] = m_n[m_l[i].n[j]].x;
			y1[j] = m_n[m_l[i].n[j]].y;
		}
		x1[3] = x1[0];
		y1[3] = y1[0];
		if (m_l[i].domain == 2) {
			for (int k = 0; k < 3; k++) {
				if (abs(x1[k] - xLeft) < delta) {
					if (y1[k] >= yDown - delta && y1[k] <= yUp + delta) {
						BCx += 1;
					}
				}
				if (abs(x1[k] - xRight) < delta) {
					if (y1[k] >= yDown - delta && y1[k] <= yUp + delta) {
						BCx += 3;
					}
				}
				if (abs(y1[k] - yDown) < delta) {
					if (x1[k] >= xLeft - delta && x1[k] <= xRight + delta) {
						BCy += 1;
					}
				}
				if (abs(y1[k] - yUp) < delta) {
					if (x1[k] >= xLeft - delta && x1[k] <= xRight + delta) {
						BCy += 3;
					}
				}
			}
		}

		int ind = 10;
		int s = 0;
		if (BCx == 0 && BCy == 0) {//outside
			newCurve->setBrush(Qt::NoBrush);
		} else if (BCy == 2 || BCy == 6) {//Y,2,edge
			newCurve->setBrush(QColor(0, 0, 255));//blue
			/****Find which point is moved,1?2?3?*****/
			if (abs(y1[0] - y1[1]) < delta) {
				ind = 2;
			} else if (abs(y1[1] - y1[2]) < delta) {
				ind = 0;
			} else if (abs(y1[0] - y1[2]) < delta) {
				ind = 1;
			}
			s = -1;
			gpx1.push_back(x1[ind]);
			gpy1.push_back(y1[ind]);
		} else if (BCx == 0 && BCy % 2 == 1) {//Y,1,point
			newCurve->setBrush(QColor(255, 0, 0));
			/****Find which point is moved,1?2?3?*****/
			for (int k = 0; k < 3; k++) {
				if (abs(y1[k] - yUp) < delta) {
					ind = k;
					break;
				}
				if (abs(y1[k] - yDown) < delta) {
					ind = k;
					break;
				}
			}
			s = 1;
			gpx2.push_back(x1[ind]);
			gpy2.push_back(y1[ind]);
		} else if (BCx == 2 || BCx == 6) {//X,2,edge
			newCurve->setBrush(QColor(0, 0, 255));//blue
			/****Find which point is moved,1?2?3?*****/
			if (abs(x1[0] - x1[1]) < delta) {
				ind = 2;
			} else if (abs(x1[1] - x1[2]) < delta) {
				ind = 0;
			} else if (abs(x1[0] - x1[2]) < delta) {
				ind = 1;
			}
			s = -1;
			gpx1.push_back(x1[ind]);
			gpy1.push_back(y1[ind]);
		} else if (BCx % 2 == 1 && BCy == 0) {//X,1,point
			newCurve->setBrush(QColor(255, 0, 0));
			/****Find which point is moved,1?2?3?*****/
			for (int k = 0; k < 3; k++) {
				if (abs(x1[k] - xLeft) < delta) {
					ind = k;
					break;
				}
				if (abs(x1[k] - xRight) < delta) {
					ind = k;
					break;
				}
			}
			s = 1;
			gpx2.push_back(x1[ind]);
			gpy2.push_back(y1[ind]);
		} else {//corner,point
			newCurve->setBrush(QColor(0, 255, 0));//green	
			/****Find which point is moved,1?2?3?*****/
			for (int k = 0; k < 3; k++) {
				if (abs(y1[k] - yUp) < delta) {
					ind = k;
					break;
				}
				if (abs(y1[k] - yDown) < delta) {
					ind = k;
					break;
				}
			}
			s = 1;
			gpx2.push_back(x1[ind]);
			gpy2.push_back(y1[ind]);
		}

		/*****Cacl the Force*******/
		double xf1 = 0;
		double xf2 = 0;
		double xf3 = 0;
		double beta2 = 0;
		double beta3;
		double Ac = 0;
		double A1, A2, A3;
		double tmp = PI / 2 * PI * m_l[i].rc * m_l[i].AREA / (4 * PI*1e-7);
		A1 = m_n[m_l[i].n[0]].A;
		A2 = m_n[m_l[i].n[1]].A;
		A3 = m_n[m_l[i].n[2]].A;
		Ac = (A1 + A2 + A3) / 3;
		if (ind != 10) {
			xf1 += m_l[i].B * m_l[i].B / m_l[i].AREA;
			xf1 *= m_l[i].Q[ind] / 2.0*s;

			beta2 = m_l[i].Q[0] * A1 +
				m_l[i].Q[1] * A2 +
				m_l[i].Q[2] * A3;
			beta2 /= (2. * m_l[i].AREA);
			xf2 = (2. * beta2 + 2.0 * Ac / m_l[i].rc)*beta2 / (2.0 * m_l[i].AREA)*(-m_l[i].Q[ind] * s);

			beta3 = m_l[i].P[0] * A1 +
				m_l[i].P[1] * A2 +
				m_l[i].P[2] * A3;
			beta3 /= (2. * m_l[i].AREA);
			if (ind == 0) {
				xf3 = 2. * beta3*(s*(A3 - A2)*(2. * m_l[i].AREA) - beta3*(2. * m_l[i].AREA)*(m_l[i].Q[0])*s);
				xf3 /= (2. * m_l[i].AREA)*(2. * m_l[i].AREA);
			} else if (ind == 1) {
				xf3 = 2. * beta3*(s*(A1 - A3)*(2. * m_l[i].AREA) - beta3*(2. * m_l[i].AREA)*(m_l[i].Q[1])*s);
				xf3 /= (2. * m_l[i].AREA)*(2. * m_l[i].AREA);
			} else {
				xf3 = 2. * beta3*(s*(A2 - A1)*(2. * m_l[i].AREA) - beta3*(2. * m_l[i].AREA)*(m_l[i].Q[2])*s);
				xf3 /= (2. * m_l[i].AREA)*(2. * m_l[i].AREA);
			}
			yForce -= (tmp*(xf1 + xf2 + xf3));
		}
		newCurve->setData(x1, y1);
		newCurve->setPen(QPen(cc[m_l[i].domain - 1]));
		//if (ind != 10) {
		//	fprintf(fp, "%d\n", ind);			
		//} 

	}
	graph1->setData(gpx1, gpy1);
	graph2->setData(gpx2, gpy2);
	//fclose(fp);
	this->setWindowTitle(QString::number(yForce));
	customPlot->replot();
}