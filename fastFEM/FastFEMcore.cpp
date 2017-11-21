#include <QFile>
#include <QXmlStreamWriter>
#include <QXmlStreamReader>
#include <QDebug>
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
#include <QtAlgorithms>
#include <QVector>
#include "FastFEMcore.h"
#include "SuperLU_MT.h"
#include "qcustomplot.h"
#include "slu_mt_ddefs.h"

using namespace std;
using namespace arma;

#define PI 3.14159265358979323846
#define r 2

double GaussPoint3[] = { -0.774596669241483, 0, 0.774596669241483 };
double GaussWeight3[] = { 0.555555555555556, 0.888888888888889, 0.555555555555556 };

double GaussPoint5[] = { -0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664 };
double GaussWeight5[] = { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 };

const double ste = 0;
const double miu0 = PI*4e-7;
typedef struct {
	int level;
	int index;
} ele;

bool compareele(const ele ele1, const ele ele2) {
	if (ele1.level < ele2.level) {
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
CFastFEMcore::CFastFEMcore() {
	Precision = 1e-6;;
	LengthUnits = 0;
	pmeshnode = NULL;
	pmeshele = NULL;
	materialList = NULL;

	thePlot = new Plot();
	thePlot->show();
	qApp->processEvents();//强制刷新界面
}


CFastFEMcore::~CFastFEMcore() {
	//should free the space allocated
	if (pmeshnode != NULL) free(pmeshnode);
	if (pmeshele != NULL) free(pmeshele);
	if (materialList != NULL) delete[]materialList;
	if (thePlot != NULL) delete thePlot;
}


// load mesh
int CFastFEMcore::Load2DMeshCOMSOL(const char fn[]) {
	char ch[256];
	//------------open file----------------------------------
	FILE * fp = NULL;
	fp = fopen(fn, "r");
	if (fp == NULL) {
		qDebug() << "Error: openning file!";
		return 1;
	}
	//--------------Read the head-----------------------------
	for (int i = 0; i < 18; i++) {
		fgets(ch, 256, fp);
	}
	//-----------------mesh point-----------------------------
	//读取节点数目
	if (fscanf(fp, "%d # number of mesh points\n", &num_pts)) {
		pmeshnode = (CNode*)calloc(num_pts, sizeof(CNode));

		for (int i = 0; i < num_pts; i++) {
			pmeshnode[i].I = 0;
			pmeshnode[i].pm = 0;
		}
	} else {
		qDebug() << "Error: reading num_pts!";
		return 1;
	}
	int pts_ind;//the beginning of the points index
	//读取节点索引，默认从0开始
	if (fscanf(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
		qDebug() << "Error: reading pts_ind!";
		return 1;
	}
	fgets(ch, 256, fp);

	for (int i = pts_ind; i < num_pts; i++) {
		//读取x,y坐标
		if (fscanf(fp, "%lf %lf \n", &(pmeshnode[i].x), &(pmeshnode[i].y)) != 2) {
			qDebug() << "Error: reading mesh point!";
			return 1;
		}
	}
	//---------------vertex-------------------------------
	for (int i = 0; i < 7; i++)
		fgets(ch, 256, fp);
	int num_vtx_ns, num_vtx_ele;
	//
	if (fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns) != 1) {
		qDebug() << "Error: reading num_vtx_ns!";
		return 1;
	}

	if (fscanf(fp, "%d # number of elements\n", &num_vtx_ele) != 1) {
		qDebug() << "Error: reading num_vtx_ele!";
		return 1;
	}
	fgets(ch, 256, fp);

	int *vtx;
	vtx = (int*)calloc(num_vtx_ele, sizeof(int));
	for (int i = 0; i < num_vtx_ele; i++) {
		//好象是每一个域的顶点编号
		if (fscanf(fp, "%d \n", vtx + i) != 1) {
			qDebug() << "Error: reading vertex condition!";
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
			qDebug() << "Error: reading vertex condition!";
			return 1;
		}
	}
	if (vtx2 != NULL) free(vtx2); vtx2 = NULL;
	//--------------boundary--------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int num_bdr_ns, num_bdr_ele;//number of nodes per element;number of elements
	//读取一个边界单元中的数目，2D的话为2，表示线段
	if (fscanf(fp, "%d # number of nodes per element\n", &num_bdr_ns) != 1) {
		qDebug() << "Error: reading num_bdr_ns!";
		return 1;
	}
	//读取线段边界数目
	if (fscanf(fp, "%d # number of elements\n", &num_bdr_ele) != 1) {
		qDebug() << "Error: reading num_bdr_ele!";
		return 1;
	}
	fgets(ch, 256, fp);

	int *p1, *p2;
	p1 = (int*)calloc(num_bdr_ele, sizeof(int));
	p2 = (int*)calloc(num_bdr_ele, sizeof(int));
	for (int i = 0; i < num_bdr_ele; i++) {
		//读取线段边界的起点和终点
		if (fscanf(fp, "%d %d\n", p1 + i, p2 + i) == 2) {
			pmeshnode[p1[i]].bdr = 1;
		} else {
			qDebug() << "Error: reading boundary condition!";
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
			qDebug() << "Error: reading boundary condition!";
			return 1;
		}
	}
	if (entity != NULL) free(entity); entity = NULL;
	//----------------elements------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int ns_per_ele;// num_ele;//number of nodes per element;number of elements
	if (fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele) != 1) {
		qDebug() << "Error: reading ns_per_ele!";
		return 1;
	}
	//读取分网单元数目
	if (fscanf(fp, "%d # number of elements\n", &num_ele) == 1) {
		pmeshele = (CElement*)calloc(num_ele, sizeof(CElement));
	} else {
		qDebug() << "Error: reading num_ele!";
		return 1;
	}
	fgets(ch, 256, fp);
	//读取分网三角单元的三个节点索引
	for (int i = 0; i < num_ele; i++) {
		if (fscanf(fp, "%d %d %d \n", &pmeshele[i].n[0], &pmeshele[i].n[1], &pmeshele[i].n[2]) != 3) {
			qDebug() << "Error: reading elements points!";
			return 1;
		}
	}
	//---------------Domain----------------------------------
	int num_domain;
	//读取domain数目
	fscanf(fp, "%d # number of geometric entity indices\n", &num_domain);
	fgets(ch, 256, fp);

	for (int i = 0; i < num_domain; i++) {
		//读取每个单元所在的domain
		if (fscanf(fp, "%d \n", &pmeshele[i].domain) != 1) {
			qDebug() << "Error: reading domain points!";
			return 1;
		}
	}
	fclose(fp);
	return 0;
}

void CFastFEMcore::myTriSolve(int ncore, SuperMatrix *L, SuperMatrix *U,
	int_t *perm_r, int_t *perm_c, SuperMatrix *B, int_t *info){
	//--------------Plot start-------------------
	QCustomPlot * customplot;
	customplot = thePlot->getQcustomPlot();
	customplot->yAxis->setScaleRatio(customplot->xAxis, 1);

	SCPformat *Lstore = (SCPformat *)L->Store;
	NCPformat *Ustore = (NCPformat *)U->Store;
	//-----------------------------------------
	//绘制U矩阵
	QCPGraph *graphU = customplot->addGraph();
	graphU->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 1), QBrush(Qt::red), 1));
	graphU->setPen(QPen(QColor(120, 120, 120), 2));
	graphU->setLineStyle(QCPGraph::lsNone);

	//绘制边框
	int m = L->ncol;
	int nrhs = 1;
	QVector<double> xU(4), yU(4);
	xU[0] = 0; xU[1] = m - 1; xU[2] = m - 1; xU[3] = 0;
	yU[0] = 0; yU[1] = 0; yU[2] = m - 1; yU[3] = m - 1;
	QCPCurve *newCurveU = new QCPCurve(customplot->xAxis, customplot->yAxis);
	newCurveU->setBrush(QColor(255, 0, 0, 100));


	//绘制L矩阵
	QCPGraph *graphL = customplot->addGraph();
	graphL->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 1));
	graphL->setPen(QPen(QColor(120, 120, 120), 2));
	graphL->setLineStyle(QCPGraph::lsNone);

	QVector<double> xL(2), yL(2);
	xL[0] = 0 + m; xL[1] = m - 1 + m;
	yL[0] = m - 1; yL[1] = 0;
	QCPCurve *newCurveL = new QCPCurve(customplot->xAxis, customplot->yAxis);
	newCurveL->setBrush(QColor(255, 0, 0, 100));
	newCurveL->setData(xL, yL);

	int nnzL = Lstore->nnz;
	//asub = Lstore->rowind;

	//qDebug()<<"L:"<<nnzL;
	//qDebug()<<"nsuper:"<<Lstore->nsuper;//start from zero
	double t1 = SuperLU_timer_();
	QVector<int> xnL(nnzL), ynL(nnzL);
	QVector<double> valueL(nnzL);
	//    int *xnL = (int*)malloc(nnzL*sizeof(int));
	//    int *ynL = (int*)malloc(nnzL*sizeof(int));
	//    double *valueL = (double*)malloc(nnzL*sizeof(double));
	t1 = SuperLU_timer_() - t1; qDebug() << t1;
	int nnzU = Ustore->nnz;
	//qDebug()<<"U:"<<nnzU;
	//asubU = Ustore->rowind;
	t1 = SuperLU_timer_();
	QVector<int> xnU(nnzU), ynU(nnzU);
	QVector<double> valueU(nnzU);
	//    int *xnU = (int*)malloc(nnzU*sizeof(int));
	//    int *ynU = (int*)malloc(nnzU*sizeof(int));
	//    double *valueU = (double*)malloc(nnzU*sizeof(double));
	t1 = SuperLU_timer_() - t1; qDebug() << t1;
	int fsupc, istart, nsupr, nsupc, nrow;//
	int count = 0; int countU = 0;
	double value;
	t1 = SuperLU_timer_();
	QVector<int> level(m), sortlevel(m), sortlevelU(m);
	t1 = SuperLU_timer_() - t1; qDebug() << t1;
	for (int i = 0; i < m; i++){
		level[i] = 0;
	}
	QVector<double> Udiag(m);
	//对Lstore进行读取,nsuper,超级节点个数，从0开始
	for (int ksupno = 0; ksupno <= Lstore->nsuper; ++ksupno) {
		fsupc = Lstore->sup_to_colbeg[ksupno];//第ksupno个超级节点长方形区域的起始列号
		istart = Lstore->rowind_colbeg[fsupc];//第fsupc列的起始行号，这个列必须是超级节点第一列，
		//超级节点是长方形结构，行号一样,但是某些位置是被零填充的
		nsupr = Lstore->rowind_colend[fsupc] - istart;//第fsupc列的结束行号，算出来是超级节点的列宽
		nsupc = Lstore->sup_to_colend[ksupno] - fsupc;//算出来是超级节点的行高度
		nrow = nsupr - nsupc;//算出来是超级节点不在对角区域的行高度

		if (nsupc == 1) {//超级节点只有一列
			for (int j = 0; j < nrhs; j++) {//对B是多列的情况
				int luptr = Lstore->nzval_colbeg[fsupc];//非零元素的序号
				for (int iptr = istart; iptr < Lstore->rowind_colend[fsupc]; iptr++){
					int irow = Lstore->rowind[iptr];//行号
					value = ((double*)Lstore->nzval)[luptr];//值

					if (fabs(value) > 1e-9){//非零元素,如果不考虑这个，那个level就没法算了
						if (irow >= fsupc){//下三角
							xnL[count] = fsupc;
							ynL[count] = irow;

							if (irow == fsupc){//对角线上元素
								level[irow] += 1;//计算level=max(level)+1;
								valueL[count] = 1;
							} else if (irow > fsupc){//非对角线元素
								//计算level，取该行的最大值
								level[irow] = std::max(level[fsupc], level[irow]);
								valueL[count] = value;
							}
							count++;
						}
					}
					++luptr;
				}
			}
		} else {
			for (int j = 0; j < nrhs; j++){
				for (int i = 0; i < nsupc; i++){
					int luptr = Lstore->nzval_colbeg[fsupc + i];
					for (int iptr = istart; iptr < Lstore->rowind_colend[fsupc]; iptr++){
						int irow = Lstore->rowind[iptr];
						value = ((double*)Lstore->nzval)[luptr];

						if (fabs(value) > 1e-9){
							if (irow >= fsupc + i){
								xnL[count] = fsupc + i;
								ynL[count] = irow;

								if (irow == fsupc + i){
									level[irow] += 1;
									valueL[count] = 1;
								} else if (irow > fsupc + i){
									level[irow] = std::max(level[fsupc + i], level[irow]);
									valueL[count] = value;
								}
								count++;
							}
						}
						++luptr;
					}
				}
			}
		} /* if-else: nsupc == 1 ... */
	} /* for L-solve */

	//对Ustore进行读取
	QVector <int> levelU(m); levelU.fill(0);
	for (int ksupno = Lstore->nsuper; ksupno >= 0; --ksupno) {
		fsupc = Lstore->sup_to_colbeg[ksupno];//第ksupno个超级节点长方形区域的起始列号
		istart = Lstore->rowind_colbeg[fsupc];//第fsupc列的起始行号，这个列必须是超级节点第一列，
		//超级节点是长方形结构，行号一样,但是某些位置是被零填充的
		nsupr = Lstore->rowind_colend[fsupc] - istart;//第fsupc列的结束行号，算出来是超级节点的列宽
		nsupc = Lstore->sup_to_colend[ksupno] - fsupc;//算出来是超级节点的行高度
		nrow = nsupr - nsupc;//算出来是超级节点不在对角区域的行高度
		int iend = Lstore->rowind_colend[fsupc] - 1;

		if (nsupc == 1) {//超级节点只有一列
			for (int j = 0; j < nrhs; j++) {//对B是多列的情况
				int luptr = Lstore->nzval_colend[fsupc] - 1;//非零元素的序号
				for (int iptr = iend; iptr >= Lstore->rowind_colbeg[fsupc]; iptr--){
					int irow = Lstore->rowind[iptr];//行号
					value = ((double*)Lstore->nzval)[luptr];//值

					if (fabs(value) > 1e-9){//非零元素,如果不考虑这个，那个level就没法算了
						if (irow < fsupc){//上三角
							xnU[countU] = fsupc;
							ynU[countU] = irow;
							valueU[countU] = value;
							levelU[irow] = std::max(levelU[fsupc], levelU[irow]);

							countU++;
						} else if (irow == fsupc){
							Udiag[irow] = value;
							levelU[irow] += 1;//计算level=max(level)+1;
						}
					}
					--luptr;
				}
			}
			//读取Ustore中数据
			for (int u = Ustore->colend[fsupc] - 1; u >= Ustore->colbeg[fsupc]; u--){
				double value1 = ((double*)Ustore->nzval)[u];
				int irow = Ustore->rowind[u];//row
				if (fabs(value1) > 1e-9){
					xnU[countU] = fsupc;//col
					ynU[countU] = irow;
					valueU[countU] = value1;
					levelU[irow] = std::max(levelU[fsupc], levelU[irow]);

					countU++;
				}
			}
		} else {
			for (int j = 0; j < nrhs; j++){
				for (int i = nsupc - 1; i >= 0; i--){
					int luptr = Lstore->nzval_colend[fsupc + i] - 1;
					for (int iptr = iend; iptr >= Lstore->rowind_colbeg[fsupc]; iptr--){
						int irow = Lstore->rowind[iptr];
						value = ((double*)Lstore->nzval)[luptr];

						if (fabs(value) > 1e-9){
							if (irow < fsupc + i){
								xnU[countU] = fsupc + i;
								ynU[countU] = irow;
								valueU[countU] = value;
								levelU[irow] = std::max(levelU[fsupc + i], levelU[irow]);

								countU++;
							} else if (irow == fsupc + i){
								Udiag[irow] = value;
								levelU[irow] += 1;
							}
						}
						--luptr;
					}
					//读取该列的Ustore数据
					for (int u = Ustore->colend[fsupc + i] - 1; u >= Ustore->colbeg[fsupc + i]; u--){
						double value1 = ((double*)Ustore->nzval)[u];
						int irow = Ustore->rowind[u];//row
						if (fabs(value1) > 1e-9){
							xnU[countU] = fsupc + i;//col
							ynU[countU] = irow;
							valueU[countU] = value1;
							levelU[irow] = std::max(levelU[fsupc + i], levelU[irow]);

							countU++;
						}
					}
				}
			}
		} /* if-else: nsupc == 1 ... */
	} /* for U-solve */


	qDebug() << "countL :" << count;
	qDebug() << "countU :" << countU;
	int maxLevel = 0;
	int maxLevelU = 0;
	QVector <ele> slevel(m);
	QVector <ele> slevelU(m);
	for (int i = 0; i < m; i++){
		//计算最大的level
		maxLevel = maxLevel > level[i] ? maxLevel : level[i];
		maxLevelU = maxLevelU > levelU[i] ? maxLevelU : levelU[i];
		slevel[i].index = i;
		slevel[i].level = level[i];
		slevelU[i].index = i;
		slevelU[i].level = levelU[i];
	}
	qDebug() << "maxLevel:" << maxLevel;
	qDebug() << "maxLevelU:" << maxLevelU;

	qSort(slevel.begin(), slevel.end(), compareele);//ascend
	qSort(slevelU.begin(), slevelU.end(), compareele);//ascend
	QVector <int> Lcol(maxLevel);
	QVector <int> Ucol(maxLevelU);
	Lcol[0] = 0;
	int lcount = 0;
	int ucount = 0;
	for (int i = 1; i < m; i++){
		if (slevel[i].level != slevel[i - 1].level){
			Lcol[++lcount] = i;
		}
		if (slevelU[i].level != slevelU[i - 1].level){
			Ucol[++ucount] = i;
		}
	}
	qDebug() << "lcount:" << lcount;
	qDebug() << "ucount:" << ucount;

	for (int i = 0; i < m; i++){
		sortlevel[slevel[i].index] = i;
		sortlevelU[slevelU[i].index] = i;
	}

	QVector <int> xxnL(count), yynL(count);
	QVector <int> xxnU(countU), yynU(countU);

	for (int i = 0; i < count; i++){
		xxnL[i] = xnL[i];
		yynL[i] = ynL[i];
	}
	QVector <ele> Lloc(count), Uloc(countU);
	//转换坐标
	for (int i = 0; i < count; i++){
		xxnL[i] = sortlevel[xnL[i]];
		yynL[i] = sortlevel[ynL[i]];
		Lloc[i].level = yynL[i];//row
		Lloc[i].index = xxnL[i];//col
	}
	for (int i = 0; i < countU; i++){
		xxnU[i] = sortlevelU[xnU[i]];
		yynU[i] = sortlevelU[ynU[i]];
		Uloc[i].level = yynL[i];//row
		Uloc[i].index = xxnL[i];//col
	}
	qSort(Lloc.begin(), Lloc.end(), compareele);//ascend
	qSort(Uloc.begin(), Uloc.end(), compareele);//ascend

	QVector <int> Lrowindex(m), Urowindex(m);
	Lrowindex[0] = 0;
	for (int i = 1; i < count; i++){
		if (Lloc[i].level != Lloc[i - 1].level){
			Lrowindex[Lloc[i].level] = i;
		}
	}
	for (int i = 1; i < countU; i++){
		if (Uloc[i].level != Uloc[i - 1].level){
			Urowindex[Uloc[i].level] = i;
		}
	}

	qDebug() << m;
	t1 = SuperLU_timer_();

	DNformat *Bstore = (DNformat*)B->Store;
	double *Bmat = (double*)Bstore->nzval;
	//L solve
	int cc = 0;
	omp_set_num_threads(8);
	qDebug() << omp_get_num_procs();
	for (int i = 0; i < maxLevel - 1; i++){
#pragma omp parallel for
		for (int j = Lcol[i + 1]; j < m; j++){//row
			for (int k = Lrowindex[j]; k < Lrowindex[j + 1]; k++){
				//(Lloc[k].index >=Lcol[i] && Lloc[k].index <= Lcol[i+1])
				Bmat[j] -= Bmat[Lloc[k].index] * valueL[k];
				//cc++;//because of this
			}
		}
	}

	//U solve
	for (int i = 0; i < maxLevelU; i++){

	}
	t1 = SuperLU_timer_() - t1; qDebug() << t1 << "\t" << cc;
	//绘制level分割线等
	double lastline = 0;
	QPen pen;
	pen.setStyle(Qt::DashLine);
	pen.setWidth(1);
	pen.setColor(QColor(0, 0, 0));
	QVector <int> levelLrow(maxLevel + 1); levelLrow[0] = 0;
	for (int i = 0; i < m - 1; i++){
		if (slevel[i].level != slevel[i + 1].level){
			levelLrow[slevel[i + 1].level] = i + 1;
			//qDebug()<<i+1;
			QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
			QVector <double> xline(2), yline(2);
			xline[0] = 0; xline[1] = i + 0.5;
			yline[0] = m - 1 - (i + 0.5); yline[1] = m - 1 - (i + 0.5);
			//newCurve->setData(xline, yline);
			//newCurve->setPen(pen);
			//newCurve->pen().setStyle(Qt::DotLine);
			//newCurve->setBrush(QColor(0, 120, 0));
			//--绘制并行区域
			QCPCurve *rec1 = new QCPCurve(customplot->xAxis, customplot->yAxis);
			QCPCurve *rec2 = new QCPCurve(customplot->xAxis, customplot->yAxis);
			QVector <double> recdatax(5), recdatay(5);
			recdatax[0] = lastline; recdatax[1] = i;
			recdatax[2] = i; recdatax[3] = lastline; recdatax[4] = lastline;
			recdatay[0] = m - 1 - i; recdatay[1] = m - 1 - i;
			recdatay[2] = m - 1 - lastline; recdatay[3] = m - 1 - lastline; recdatay[4] = m - 1 - i;
			rec1->setData(recdatax, recdatay);
			rec1->setPen(Qt::NoPen);
			rec1->setBrush(QColor(255, 0, 0, 100));

			recdatax[0] = lastline; recdatax[1] = i;
			recdatax[2] = i; recdatax[3] = lastline; recdatax[4] = lastline;
			recdatay[0] = 0; recdatay[1] = 0;
			recdatay[2] = m - 1 - i; recdatay[3] = m - 1 - i; recdatay[4] = 0;
			if (slevel[i].level % 2 == 1){
				rec2->setData(recdatax, recdatay);
				rec2->setPen(Qt::NoPen);
				rec2->setBrush(QColor(0, 255, 0, 100));
			}

			lastline = i + 1;
		}
	}
	//绘制U矩阵的分割线
	lastline = 0;
	for (int i = 0; i < m - 1; i++){
		if (slevelU[i].level != slevelU[i + 1].level){
			QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
			QVector <double> xline(2), yline(2);
			xline[0] = m - 1; xline[1] = m - 1 - i - 0.5;
			yline[0] = (i + 0.5); yline[1] = (i + 0.5);
			//newCurve->setData(xline, yline);
			//newCurve->setPen(QPen(QColor(0, 120, 255),1));
			//newCurve->setBrush(QColor(0, 120, 255));
			//--绘制并行区域
			QCPCurve *rec1 = new QCPCurve(customplot->xAxis, customplot->yAxis);
			QCPCurve *rec2 = new QCPCurve(customplot->xAxis, customplot->yAxis);
			QVector <double> recdatax(5), recdatay(5);
			recdatax[0] = m - 1 - lastline; recdatax[1] = m - 1 - i;
			recdatax[2] = m - 1 - i; recdatax[3] = m - 1 - lastline; recdatax[4] = m - 1 - lastline;
			recdatay[0] = lastline; recdatay[1] = lastline;
			recdatay[2] = i; recdatay[3] = i; recdatay[4] = lastline;
			rec1->setData(recdatax, recdatay);
			rec1->setPen(Qt::NoPen);
			rec1->setBrush(QColor(255, 0, 20, 100));

			recdatax[0] = m - 1 - lastline; recdatax[1] = m - 1 - i;
			recdatax[2] = m - 1 - i; recdatax[3] = m - 1 - lastline; recdatax[4] = m - 1 - lastline;
			recdatay[0] = i; recdatay[1] = i;
			recdatay[2] = m - 1; recdatay[3] = m - 1; recdatay[4] = i;
			if (slevelU[i].level % 2 == 1){
				rec2->setData(recdatax, recdatay);
				rec2->setPen(Qt::NoPen);
				rec2->setBrush(QColor(0, 255, 0, 100));
			}
			lastline = i + 1;
		}
	}
	//graphU->setData(xxnU, yynU);
	//graphL->setData(xxnL, yynL);
	//--------------Plot End---------------------
}

// 这部分使用TLM进行求解
bool CFastFEMcore::StaticAxisymmetricTLM() {
	double time[10];
	int tt = 0;
	time[tt++] = SuperLU_timer_();
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele[i].LinearFlag) {
			D34.push_back(i);
		}
	}
	uvec node_reorder = zeros<uvec>(num_pts);
	uvec node_pos = zeros<uvec>(num_pts);
	//------------build C Matrix-----------------------------
	umat locs(2, 9 * num_ele);
	locs.zeros();
	mat vals(1, 9 * num_ele);
	double ce[3][3] = { 0 };
	ResistMarix *rm = (ResistMarix*)malloc(num_ele * sizeof(ResistMarix));
	vec bbJz = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	vec INL = zeros<vec>(num_pts);
	double * ydot = (double*)malloc(num_ele*sizeof(double));
	//重新对节点进行编号，将边界点分离
	int node_bdr = 0;
	for (int i = 0; i < num_pts; i++) {
		if (pmeshnode[i].bdr == 3) {
			node_bdr++;
			node_reorder(num_pts - node_bdr) = i;
			node_pos(i) = num_pts - node_bdr;
			pmeshnode[i].A = 0;
			A(i) = 0;
		} else {
			node_reorder(i - node_bdr) = i;
			node_pos(i) = i - node_bdr;
		}
	}
	double* unknown_b = (double*)calloc(num_pts - node_bdr, sizeof(double));
	int pos = 0;
	//轴对称：A'=rA,v'=v/r,
	for (int i = 0; i < num_ele; i++) {
		//确定单元的近似半径
		int flag = 0;
		for (int f = 0; f < 3; f++)
			if (pmeshnode[pmeshele[i].n[f]].x < 1e-7)
				flag++;

		if (flag == 2) {
			ydot[i] = pmeshele[i].rc;
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

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;//猜测值
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;

		//生成单元矩阵，线性与非线性
		// 因为线性与非线性的差不多，所以不再分开讨论了
		ce[0][0] = rm[i].Y11;
		ce[1][1] = rm[i].Y22;
		ce[2][2] = rm[i].Y33;

		if (pmeshele[i].LinearFlag) {//线性区，不用计算
			ce[0][1] = rm[i].Y12;
			ce[0][2] = rm[i].Y13;
			ce[1][2] = rm[i].Y23;
		} else {
			if (rm[i].Y12 < 0){//电阻
				ce[0][1] = rm[i].Y12;
			} else{
				ce[0][0] += rm[i].Y12;//在对角项上减去受控源
				ce[1][1] += rm[i].Y12;
				ce[0][1] = 0;//受控源在右侧，所以为0
			}
			if (rm[i].Y13 < 0){
				ce[0][2] = rm[i].Y13;
			} else{
				ce[0][0] += rm[i].Y13;
				ce[2][2] += rm[i].Y13;
				ce[0][2] = 0;
			}
			if (rm[i].Y23 < 0){
				ce[1][2] = rm[i].Y23;
			} else{
				ce[1][1] += rm[i].Y23;
				ce[2][2] += rm[i].Y23;
				ce[1][2] = 0;
			}
		}
		ce[1][0] = ce[0][1];
		ce[2][0] = ce[0][2];
		ce[2][1] = ce[1][2];

		//将单元矩阵进行存储
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				//判断节点是否在未知节点内
				//得到排序之后的编号
				int n_row = node_pos(pmeshele[i].n[row]);
				int n_col = node_pos(pmeshele[i].n[col]);
				if (n_row < num_pts - node_bdr && n_col < num_pts - node_bdr) {
					locs(0, pos) = n_row;
					locs(1, pos) = n_col;
					vals(0, pos) = ce[row][col];
					pos++;
				}
			}
		}
		//计算电流密度//要注意domain会不会越界
		double jr = pmeshele[i].AREA*materialList[pmeshele[i].domain - 1].Jr / 3;
		for (int j = 0; j < 3; j++) {
			bbJz(pmeshele[i].n[j]) += jr;
			// 计算永磁部分
			bbJz(pmeshele[i].n[j]) -= materialList[pmeshele[i].domain - 1].H_c / 2.*pmeshele[i].Q[j];
		}
	}//end for
	time[tt++] = SuperLU_timer_();
	locs.reshape(2, pos);//重新调整大小
	vals.reshape(1, pos);
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

	INL += bbJz;
	for (int i = 0; i < num_pts - node_bdr; i++) {
		unknown_b[i] = INL(node_reorder(i));
	}
	time[tt++] = SuperLU_timer_();
	//---------------------superLU_MT---------------------------------------
	//CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
	SuperMatrix   sluA; SuperMatrix sluB, sluX;
	//NCformat *Astore;
	double   *a;
	int_t      *asub, *xa;
	int_t      *perm_r; /* row permutations from partial pivoting */
	int_t      *perm_c; /* column permutation vector */
	SuperMatrix   L;       /* factor L */
	SCPformat *Lstore;
	SuperMatrix   U;       /* factor U */
	NCPformat *Ustore;
	int_t      nrhs, info, m, n, nnz;
	int_t      nprocs; /* maximum number of processors to use. */
	int_t      panel_size, relax, maxsup;
	int_t      permc_spec;
	trans_t  trans;
	//double   *rhs;
	superlu_memusage_t   superlu_memusage;
	DNformat	   *Bstore;
	double      *rhsb, *rhsx;
	Gstat_t  Gstat1;


	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	maxsup = sp_ienv(3);

	nprocs = 1;
	nrhs = 1;
	trans = NOTRANS;
	/* create matrix A in Harwell-Boeing format.*/
	m = num_pts - node_bdr; n = num_pts - node_bdr; nnz = X.n_nonzero;
	a = const_cast<double *>(X.values);

	StatAlloc(n, nprocs, panel_size, relax, &Gstat1);
	StatInit(n, nprocs, &Gstat1);

	asub = (int*)const_cast<unsigned int*>(X.row_indices);
	xa = (int*)const_cast<unsigned int*>(X.col_ptrs);
	dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

	//------create B and X-------------------
	if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&sluX, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);

	rhsb = unknown_b;
	dCreate_Dense_Matrix(&sluB, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
	Bstore = (DNformat*)sluB.Store;

	if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");

	/*
	* Get column permutation vector perm_c[], according to permc_spec:
	*   permc_spec = 0: natural ordering
	*   permc_spec = 1: minimum degree ordering on structure of A'*A
	*   permc_spec = 2: minimum degree ordering on structure of A'+A
	*   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	*/
	permc_spec = 1;
	get_perm_c(permc_spec, &sluA, perm_c);

	pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);

	if (info != 0) {
		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = (double*)((DNformat*)sluB.Store)->nzval;
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
			A(node_reorder(i)) = sol[i];
		}
	}
	time[tt++] = SuperLU_timer_();
	//---------------------superLU--end----------------------------------
	//-----------绘图----------------------------------------
	QVector<double> x(num_ele), y(num_ele);
	for (int i = 0; i < num_ele; ++i) {
		x[i] = i;
	}
	QCustomPlot * customplot;
	customplot = thePlot->getQcustomPlot();
	QCPGraph *graph1 = customplot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 3));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	customplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	customplot->xAxis->setLabel("x");
	customplot->xAxis->setRange(0, num_ele);
	//customplot->xAxis->setAutoTickStep(false);
	//customplot->xAxis->setTicks(false);
	customplot->yAxis->setLabel("y");
	//customplot->yAxis->setRange(0, 4);
	//customplot->yAxis->setTicks(false);
	//customplot->xAxis2->setTicks(false);
	//customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);
	//---------the main loop---------------------------------
	int steps = 300;
	int count;
	double alpha = 1;
	Voltage *Vr = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	Voltage *Vi = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	time[tt++] = SuperLU_timer_();
	for (count = 0; count < steps; count++) {
		//------update miu----------------
		for (int i = 0; i < num_ele; i++) {
			double bx = 0;
			double by = 0;
			for (int j = 0; j < 3; j++) {
				bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
				by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
			}
			pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
			pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
			pmeshele[i].miut = y[i] + (count>90 ? 0.02 : (0.9 - count*0.001))*(pmeshele[i].miut - y[i]);
			//pmeshele[i].miut = y[i] + (count>30 ? 0.02 : (0.9-count*0.03))*(pmeshele[i].miut - y[i]);
			//            if(fabs(y[i]-pmeshele[i].B)/y[i] > 0.2 && count > 180){
			//            //if(pmeshele[i].domain == 4 || pmeshele[i].domain == 3){
			//                QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
			//                QVector <double> x11(4);
			//                QVector <double> y11(4);
			//                x11[0] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[0] = pmeshnode[pmeshele[i].n[0]].y;
			//                x11[1] = pmeshnode[pmeshele[i].n[1]].x;
			//                y11[1] = pmeshnode[pmeshele[i].n[1]].y;
			//                x11[2] = pmeshnode[pmeshele[i].n[2]].x;
			//                y11[2] = pmeshnode[pmeshele[i].n[2]].y;
			//                x11[3] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[3] = pmeshnode[pmeshele[i].n[0]].y;
			//                newCurve->setData(x11, y11);
			//                newCurve->setPen(QPen(Qt::black, 1));
			//                newCurve->setBrush(QColor(255, 0, 0,100));
			//                qDebug()<<pmeshele[i].B;
			//            }
			//            if(count == 1){
			//            //if(pmeshele[i].domain == 4 || pmeshele[i].domain == 3){
			//                QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
			//                QVector <double> x11(4);
			//                QVector <double> y11(4);
			//                x11[0] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[0] = pmeshnode[pmeshele[i].n[0]].y;
			//                x11[1] = pmeshnode[pmeshele[i].n[1]].x;
			//                y11[1] = pmeshnode[pmeshele[i].n[1]].y;
			//                x11[2] = pmeshnode[pmeshele[i].n[2]].x;
			//                y11[2] = pmeshnode[pmeshele[i].n[2]].y;
			//                x11[3] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[3] = pmeshnode[pmeshele[i].n[0]].y;
			//                newCurve->setData(x11, y11);
			//                newCurve->setPen(QPen(Qt::black, 1));
			//                newCurve->setBrush(QColor(255, 255, 0,100));
			//            }


			//y[i] = (A(pmeshele[i].n[0]) + A(pmeshele[i].n[1]) + A(pmeshele[i].n[2]))/3;
			y[i] = pmeshele[i].miut;
		}
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			double rtmp;//do this to mark it as private
			rtmp = (m_e->miut - m_e->miu) / (m_e->miu + m_e->miut);
			double hehe = -(m_e->miu - m_e->miut) / m_e->miut;
			//qDebug()<<m_e->miut<<"\t"<<m_e->miu<<"\t"<<rtmp;

			Vr[j].V12 = (pmeshnode[k].A - pmeshnode[m].A) - Vi[j].V12;
			Vr[j].V23 = (pmeshnode[m].A - pmeshnode[n].A) - Vi[j].V23;
			Vr[j].V13 = (pmeshnode[n].A - pmeshnode[k].A) - Vi[j].V13;

			if (rm[i].Y12 < 0) {
				Vi[j].V12 = Vr[j].V12 * rtmp;
			} else {
				Vi[j].V12 = 0;// 0.5*(pmeshnode[k].A - pmeshnode[m].A) * 1;
			}
			INL(k) += 2. *Vi[j].V12*abs(rm[i].Y12);
			INL(m) += -2. * Vi[j].V12 *abs(rm[i].Y12);
			if (rm[i].Y23 < 0) {
				Vi[j].V23 = Vr[j].V23*rtmp;
			} else {
				Vi[j].V23 = 0;// 0.5*(pmeshnode[m].A - pmeshnode[n].A) * 1;
			}
			INL(m) += 2. * Vi[j].V23*abs(rm[i].Y23);
			INL(n) += -2. *Vi[j].V23*abs(rm[i].Y23);
			if (rm[i].Y13 < 0) {
				Vi[j].V13 = Vr[j].V13 * rtmp;
			} else {
				Vi[j].V13 = 0;// 0.5*(pmeshnode[n].A - pmeshnode[k].A) * 1;
			}
			INL(n) += 2. * Vi[j].V13*abs(rm[i].Y13);
			INL(k) += -2.0 *Vi[j].V13*abs(rm[i].Y13);
		}
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		if (info != 0) {
			qDebug() << "Error: superlumt.slove";
			//qDebug() << "info: " << superlumt.info;
			break;
		} else {
			double *sol = NULL;
			A_old = A;
			sol = (double*)((DNformat*)sluB.Store)->nzval;

			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
				A(node_reorder(i)) = sol[i];
			}
		}
		double error = norm((A_old - A), 2) / norm(A, 2);
		//qDebug() << "iter: " << count;
		//qDebug() << "error: " << error;

		graph1->setData(x, y);
		customplot->rescaleAxes(true);
		//customplot->xAxis->setRange(0, 0.09);
		//customplot->yAxis->setRange(-0.04, 0.04);
		//customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);
		customplot->replot();


		if (error < Precision) {
			break;
		}
		INL.zeros();
	}
	time[tt++] = SuperLU_timer_();
	// 求解结束，更新B
	for (int i = 0; i < num_ele; i++) {
		double bx = 0;
		double by = 0;
		for (int j = 0; j < 3; j++) {
			bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
			by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
		}
		pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
		pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
	}
	for (int i = 0; i < num_pts - node_bdr; i++) {
		pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;// / pmeshnode[i].x;//the A is r*A_real
	}
	//output the time
	for (int i = 1; i < tt; i++){
		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vi != NULL) free(Vi);
	if (Vr != NULL) free(Vr);
	if (unknown_b != NULL) free(unknown_b);

	//SUPERLU_FREE(rhs);
	//SUPERLU_FREE(xact);
	//SUPERLU_FREE(perm_r);
	//SUPERLU_FREE(perm_c);
	//Destroy_CompCol_Matrix(&sluA);
	//Destroy_SuperMatrix_Store(&sluB);
	//Destroy_SuperNode_SCP(&L);
	//Destroy_CompCol_NCP(&U);
	//StatFree(&Gstat1);
	return true;
}
//尝试将三角形单元的电路看作为一个整体
bool CFastFEMcore::StaticAxisQ3TLMgroup() {
	double time[10];
	int tt = 0;
	time[tt++] = SuperLU_timer_();
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele[i].LinearFlag) {
			D34.push_back(i);
		}
	}
	uvec node_reorder = zeros<uvec>(num_pts);
	uvec node_pos = zeros<uvec>(num_pts);
	//------------build C Matrix-----------------------------
	umat locs(2, 9 * num_ele);
	locs.zeros();
	mat vals(1, 9 * num_ele);
	double ce[3][3] = { 0 };
	ResistMarix *rm = (ResistMarix*)malloc(num_ele * sizeof(ResistMarix));
	vec bbJz = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	vec INL = zeros<vec>(num_pts);
	double * ydot = (double*)malloc(num_ele*sizeof(double));
	//重新对节点进行编号，将边界点分离
	int node_bdr = 0;
	for (int i = 0; i < num_pts; i++) {
		if (pmeshnode[i].bdr == 3) {
			node_bdr++;
			node_reorder(num_pts - node_bdr) = i;
			node_pos(i) = num_pts - node_bdr;
			pmeshnode[i].A = 0;
			A(i) = 0;
		} else {
			node_reorder(i - node_bdr) = i;
			node_pos(i) = i - node_bdr;
		}
	}
	double* unknown_b = (double*)calloc(num_pts - node_bdr, sizeof(double));
	int pos = 0;
	//轴对称：A'=rA,v'=v/r,
	for (int i = 0; i < num_ele; i++) {
		//确定单元的近似半径
		int flag = 0;
		for (int f = 0; f < 3; f++)
			if (pmeshnode[pmeshele[i].n[f]].x < 1e-7)
				flag++;

		if (flag == 2) {
			ydot[i] = pmeshele[i].rc;
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

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;//猜测值
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;

		//生成单元矩阵，线性与非线性
		// 因为线性与非线性的差不多，所以不再分开讨论了
		ce[0][0] = rm[i].Y11;
		ce[1][1] = rm[i].Y22;
		ce[2][2] = rm[i].Y33;

		if (pmeshele[i].LinearFlag) {//线性区，不用计算
			ce[0][1] = rm[i].Y12;
			ce[0][2] = rm[i].Y13;
			ce[1][2] = rm[i].Y23;
		} else {
			if (rm[i].Y12 < 0) {//电阻
				ce[0][1] = rm[i].Y12;
			} else {
				ce[0][0] += rm[i].Y12;//在对角项上减去受控源
				ce[1][1] += rm[i].Y12;
				ce[0][1] = 0;//受控源在右侧，所以为0
			}
			if (rm[i].Y13 < 0) {
				ce[0][2] = rm[i].Y13;
			} else {
				ce[0][0] += rm[i].Y13;
				ce[2][2] += rm[i].Y13;
				ce[0][2] = 0;
			}
			if (rm[i].Y23 < 0) {
				ce[1][2] = rm[i].Y23;
			} else {
				ce[1][1] += rm[i].Y23;
				ce[2][2] += rm[i].Y23;
				ce[1][2] = 0;
			}
		}
		ce[1][0] = ce[0][1];
		ce[2][0] = ce[0][2];
		ce[2][1] = ce[1][2];

		//将单元矩阵进行存储
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				//判断节点是否在未知节点内
				//得到排序之后的编号
				int n_row = node_pos(pmeshele[i].n[row]);
				int n_col = node_pos(pmeshele[i].n[col]);
				if (n_row < num_pts - node_bdr && n_col < num_pts - node_bdr) {
					locs(0, pos) = n_row;
					locs(1, pos) = n_col;
					vals(0, pos) = ce[row][col];
					pos++;
				}
			}
		}
		//计算电流密度//要注意domain会不会越界
		double jr = pmeshele[i].AREA*materialList[pmeshele[i].domain - 1].Jr / 3;
		for (int j = 0; j < 3; j++) {
			bbJz(pmeshele[i].n[j]) += jr;
			// 计算永磁部分
			bbJz(pmeshele[i].n[j]) -= materialList[pmeshele[i].domain - 1].H_c / 2.*pmeshele[i].Q[j];
		}
	}//end for
	time[tt++] = SuperLU_timer_();
	locs.reshape(2, pos);//重新调整大小
	vals.reshape(1, pos);
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

	INL += bbJz;
	for (int i = 0; i < num_pts - node_bdr; i++) {
		unknown_b[i] = INL(node_reorder(i));
	}
	time[tt++] = SuperLU_timer_();
	//---------------------superLU_MT---------------------------------------
	//CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
	SuperMatrix   sluA; SuperMatrix sluB, sluX;
	//NCformat *Astore;
	double   *a;
	int_t      *asub, *xa;
	int_t      *perm_r; /* row permutations from partial pivoting */
	int_t      *perm_c; /* column permutation vector */
	SuperMatrix   L;       /* factor L */
	SCPformat *Lstore;
	SuperMatrix   U;       /* factor U */
	NCPformat *Ustore;
	int_t      nrhs, info, m, n, nnz;
	int_t      nprocs; /* maximum number of processors to use. */
	int_t      panel_size, relax, maxsup;
	int_t      permc_spec;
	trans_t  trans;
	//double   *rhs;
	superlu_memusage_t   superlu_memusage;
	DNformat	   *Bstore;
	double      *rhsb, *rhsx;
	Gstat_t  Gstat1;


	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	maxsup = sp_ienv(3);

	nprocs = 1;
	nrhs = 1;
	trans = NOTRANS;
	/* create matrix A in Harwell-Boeing format.*/
	m = num_pts - node_bdr; n = num_pts - node_bdr; nnz = X.n_nonzero;
	a = const_cast<double *>(X.values);

	StatAlloc(n, nprocs, panel_size, relax, &Gstat1);
	StatInit(n, nprocs, &Gstat1);

	asub = (int*)const_cast<unsigned int*>(X.row_indices);
	xa = (int*)const_cast<unsigned int*>(X.col_ptrs);
	dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

	//------create B and X-------------------
	if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&sluX, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);

	rhsb = unknown_b;
	dCreate_Dense_Matrix(&sluB, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
	Bstore = (DNformat*)sluB.Store;

	if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");

	/*
	* Get column permutation vector perm_c[], according to permc_spec:
	*   permc_spec = 0: natural ordering
	*   permc_spec = 1: minimum degree ordering on structure of A'*A
	*   permc_spec = 2: minimum degree ordering on structure of A'+A
	*   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	*/
	permc_spec = 1;
	get_perm_c(permc_spec, &sluA, perm_c);

	pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);

	if (info != 0) {
		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = (double*)((DNformat*)sluB.Store)->nzval;
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
			A(node_reorder(i)) = sol[i];
		}
	}
	time[tt++] = SuperLU_timer_();
	//---------------------superLU--end----------------------------------
	//-----------绘图----------------------------------------
	QVector<double> x(num_ele), y(num_ele);
	for (int i = 0; i < num_ele; ++i) {
		x[i] = i;
	}
	QCustomPlot * customplot;
	customplot = thePlot->getQcustomPlot();
	QCPGraph *graph1 = customplot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 3));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	customplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	customplot->xAxis->setLabel("x");
	customplot->xAxis->setRange(0, num_ele);
	//customplot->xAxis->setAutoTickStep(false);
	//customplot->xAxis->setTicks(false);
	customplot->yAxis->setLabel("y");
	//customplot->yAxis->setRange(0, 4);
	//customplot->yAxis->setTicks(false);
	//customplot->xAxis2->setTicks(false);
	//customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);
	//---------the main loop---------------------------------
	int steps = 300;
	int count;
	double alpha = 1;
	Voltage *Vr = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	Voltage *Vi = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	time[tt++] = SuperLU_timer_();
	for (count = 0; count < steps; count++) {
		//------update miu----------------
		for (int i = 0; i < num_ele; i++) {
			double bx = 0;
			double by = 0;
			for (int j = 0; j < 3; j++) {
				bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
				by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
			}
			pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
			pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
			pmeshele[i].miut = y[i] + (count>90 ? 0.02 : (0.9 - count*0.001))*(pmeshele[i].miut - y[i]);
			//pmeshele[i].miut = y[i] + (count>30 ? 0.02 : (0.9-count*0.03))*(pmeshele[i].miut - y[i]);
			//            if(fabs(y[i]-pmeshele[i].B)/y[i] > 0.2 && count > 180){
			//            //if(pmeshele[i].domain == 4 || pmeshele[i].domain == 3){
			//                QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
			//                QVector <double> x11(4);
			//                QVector <double> y11(4);
			//                x11[0] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[0] = pmeshnode[pmeshele[i].n[0]].y;
			//                x11[1] = pmeshnode[pmeshele[i].n[1]].x;
			//                y11[1] = pmeshnode[pmeshele[i].n[1]].y;
			//                x11[2] = pmeshnode[pmeshele[i].n[2]].x;
			//                y11[2] = pmeshnode[pmeshele[i].n[2]].y;
			//                x11[3] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[3] = pmeshnode[pmeshele[i].n[0]].y;
			//                newCurve->setData(x11, y11);
			//                newCurve->setPen(QPen(Qt::black, 1));
			//                newCurve->setBrush(QColor(255, 0, 0,100));
			//                qDebug()<<pmeshele[i].B;
			//            }
			//            if(count == 1){
			//            //if(pmeshele[i].domain == 4 || pmeshele[i].domain == 3){
			//                QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
			//                QVector <double> x11(4);
			//                QVector <double> y11(4);
			//                x11[0] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[0] = pmeshnode[pmeshele[i].n[0]].y;
			//                x11[1] = pmeshnode[pmeshele[i].n[1]].x;
			//                y11[1] = pmeshnode[pmeshele[i].n[1]].y;
			//                x11[2] = pmeshnode[pmeshele[i].n[2]].x;
			//                y11[2] = pmeshnode[pmeshele[i].n[2]].y;
			//                x11[3] = pmeshnode[pmeshele[i].n[0]].x;
			//                y11[3] = pmeshnode[pmeshele[i].n[0]].y;
			//                newCurve->setData(x11, y11);
			//                newCurve->setPen(QPen(Qt::black, 1));
			//                newCurve->setBrush(QColor(255, 255, 0,100));
			//            }


			//y[i] = (A(pmeshele[i].n[0]) + A(pmeshele[i].n[1]) + A(pmeshele[i].n[2]))/3;
			y[i] = pmeshele[i].miut;
		}
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			double rtmp;//do this to mark it as private
			rtmp = (m_e->miut - m_e->miu) / (m_e->miu + m_e->miut);
			double hehe = -(m_e->miu - m_e->miut) / m_e->miut;
			//qDebug()<<m_e->miut<<"\t"<<m_e->miu<<"\t"<<rtmp;

			Vr[j].V12 = (pmeshnode[k].A - pmeshnode[m].A) - Vi[j].V12;
			Vr[j].V23 = (pmeshnode[m].A - pmeshnode[n].A) - Vi[j].V23;
			Vr[j].V13 = (pmeshnode[n].A - pmeshnode[k].A) - Vi[j].V13;

			if (rm[i].Y12 < 0) {
				Vi[j].V12 = Vr[j].V12 * rtmp;
			} else {
				Vi[j].V12 = 0;// 0.5*(pmeshnode[k].A - pmeshnode[m].A) * 1;
			}
			INL(k) += 2. *Vi[j].V12*abs(rm[i].Y12);
			INL(m) += -2. * Vi[j].V12 *abs(rm[i].Y12);
			if (rm[i].Y23 < 0) {
				Vi[j].V23 = Vr[j].V23*rtmp;
			} else {
				Vi[j].V23 = 0;// 0.5*(pmeshnode[m].A - pmeshnode[n].A) * 1;
			}
			INL(m) += 2. * Vi[j].V23*abs(rm[i].Y23);
			INL(n) += -2. *Vi[j].V23*abs(rm[i].Y23);
			if (rm[i].Y13 < 0) {
				Vi[j].V13 = Vr[j].V13 * rtmp;
			} else {
				Vi[j].V13 = 0;// 0.5*(pmeshnode[n].A - pmeshnode[k].A) * 1;
			}
			INL(n) += 2. * Vi[j].V13*abs(rm[i].Y13);
			INL(k) += -2.0 *Vi[j].V13*abs(rm[i].Y13);
		}
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		if (info != 0) {
			qDebug() << "Error: superlumt.slove";
			//qDebug() << "info: " << superlumt.info;
			break;
		} else {
			double *sol = NULL;
			A_old = A;
			sol = (double*)((DNformat*)sluB.Store)->nzval;

			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
				A(node_reorder(i)) = sol[i];
			}
		}
		double error = norm((A_old - A), 2) / norm(A, 2);
		//qDebug() << "iter: " << count;
		//qDebug() << "error: " << error;

		graph1->setData(x, y);
		customplot->rescaleAxes(true);
		//customplot->xAxis->setRange(0, 0.09);
		//customplot->yAxis->setRange(-0.04, 0.04);
		//customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);
		customplot->replot();


		if (error < Precision) {
			break;
		}
		INL.zeros();
	}
	time[tt++] = SuperLU_timer_();
	// 求解结束，更新B
	for (int i = 0; i < num_ele; i++) {
		double bx = 0;
		double by = 0;
		for (int j = 0; j < 3; j++) {
			bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
			by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
		}
		pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
		pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
	}
	for (int i = 0; i < num_pts - node_bdr; i++) {
		pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;// / pmeshnode[i].x;//the A is r*A_real
	}
	//output the time
	for (int i = 1; i < tt; i++) {
		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vi != NULL) free(Vi);
	if (Vr != NULL) free(Vr);
	if (unknown_b != NULL) free(unknown_b);

	//SUPERLU_FREE(rhs);
	//SUPERLU_FREE(xact);
	//SUPERLU_FREE(perm_r);
	//SUPERLU_FREE(perm_c);
	//Destroy_CompCol_Matrix(&sluA);
	//Destroy_SuperMatrix_Store(&sluB);
	//Destroy_SuperNode_SCP(&L);
	//Destroy_CompCol_NCP(&U);
	//StatFree(&Gstat1);
	return true;
}
//该方法采用牛顿的迭代思想将非线性元件线性化为电阻与电流源的并联
bool CFastFEMcore::StaticAxisTLMNR() {
	double time[10];
	int tt = 0;
	time[tt++] = SuperLU_timer_();
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele[i].LinearFlag) {
			D34.push_back(i);
		}
	}
	uvec node_reorder = zeros<uvec>(num_pts);
	uvec node_pos = zeros<uvec>(num_pts);
	//------------build C Matrix-----------------------------
	umat locs(2, 9 * num_ele);
	locs.zeros();
	mat vals(1, 9 * num_ele);
	double ce[3][3] = { 0 };
	ResistMarix *rm = (ResistMarix*)malloc(num_ele * sizeof(ResistMarix));
	vec bbJz = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	vec INL = zeros<vec>(num_pts);
	double * ydot = (double*)malloc(num_ele*sizeof(double));
	//重新对节点进行编号，将边界点分离
	int node_bdr = 0;
	for (int i = 0; i < num_pts; i++) {
		if (pmeshnode[i].bdr == 3) {
			node_bdr++;
			node_reorder(num_pts - node_bdr) = i;
			node_pos(i) = num_pts - node_bdr;
			pmeshnode[i].A = 0;
			A(i) = 0;
		} else {
			node_reorder(i - node_bdr) = i;
			node_pos(i) = i - node_bdr;
		}
	}
	double* unknown_b = (double*)calloc(num_pts - node_bdr, sizeof(double));
	int pos = 0;
	//轴对称：A'=rA,v'=v/r,
	for (int i = 0; i < num_ele; i++) {
		//确定单元的近似半径
		int flag = 0;
		for (int f = 0; f < 3; f++)
			if (pmeshnode[pmeshele[i].n[f]].x < 1e-7)
				flag++;

		if (flag == 2) {
			ydot[i] = pmeshele[i].rc;
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

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;//猜测值
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i] * pmeshele[i].miu;

		//生成单元矩阵，线性与非线性
		// 因为线性与非线性的差不多，所以不再分开讨论了
		

		if (pmeshele[i].LinearFlag) {//线性区，不用计算
			ce[0][1] = rm[i].Y12;
			ce[0][2] = rm[i].Y13;
			ce[1][2] = rm[i].Y23;
			ce[0][0] = rm[i].Y11;
			ce[1][1] = rm[i].Y22;
			ce[2][2] = rm[i].Y33;
		} else {
			ce[0][1] = -abs(rm[i].Y12);
			ce[0][2] = -abs(rm[i].Y13);
			ce[1][2] = -abs(rm[i].Y23);
			ce[0][0] = -ce[0][1] - ce[0][2];
			ce[1][1] = -ce[0][1] - ce[1][2];
			ce[2][2] = -ce[0][2]- ce[1][2];
		}
		ce[1][0] = ce[0][1];
		ce[2][0] = ce[0][2];
		ce[2][1] = ce[1][2];

		//将单元矩阵进行存储
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				//判断节点是否在未知节点内
				//得到排序之后的编号
				int n_row = node_pos(pmeshele[i].n[row]);
				int n_col = node_pos(pmeshele[i].n[col]);
				if (n_row < num_pts - node_bdr && n_col < num_pts - node_bdr) {
					locs(0, pos) = n_row;
					locs(1, pos) = n_col;
					vals(0, pos) = ce[row][col];
					pos++;
				}
			}
		}
		//计算电流密度//要注意domain会不会越界
		double jr = pmeshele[i].AREA*materialList[pmeshele[i].domain - 1].Jr / 3;
		for (int j = 0; j < 3; j++) {
			bbJz(pmeshele[i].n[j]) += jr;
			// 计算永磁部分
			bbJz(pmeshele[i].n[j]) -= materialList[pmeshele[i].domain - 1].H_c / 2.*pmeshele[i].Q[j];
		}
	}//end for
	time[tt++] = SuperLU_timer_();
	locs.reshape(2, pos);//重新调整大小
	vals.reshape(1, pos);
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

	INL += bbJz;
	for (int i = 0; i < num_pts - node_bdr; i++) {
		unknown_b[i] = INL(node_reorder(i));
	}
	time[tt++] = SuperLU_timer_();
	//---------------------superLU_MT---------------------------------------
	//CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
	SuperMatrix   sluA; SuperMatrix sluB, sluX;
	//NCformat *Astore;
	double   *a;
	int_t      *asub, *xa;
	int_t      *perm_r; /* row permutations from partial pivoting */
	int_t      *perm_c; /* column permutation vector */
	SuperMatrix   L;       /* factor L */
	SCPformat *Lstore;
	SuperMatrix   U;       /* factor U */
	NCPformat *Ustore;
	int_t      nrhs, info, m, n, nnz;
	int_t      nprocs; /* maximum number of processors to use. */
	int_t      panel_size, relax, maxsup;
	int_t      permc_spec;
	trans_t  trans;
	//double   *rhs;
	superlu_memusage_t   superlu_memusage;
	DNformat	   *Bstore;
	double      *rhsb, *rhsx;
	Gstat_t  Gstat1;


	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	maxsup = sp_ienv(3);

	nprocs = 1;
	nrhs = 1;
	trans = NOTRANS;
	/* create matrix A in Harwell-Boeing format.*/
	m = num_pts - node_bdr; n = num_pts - node_bdr; nnz = X.n_nonzero;
	a = const_cast<double *>(X.values);

	StatAlloc(n, nprocs, panel_size, relax, &Gstat1);
	StatInit(n, nprocs, &Gstat1);

	asub = (int*)const_cast<unsigned int*>(X.row_indices);
	xa = (int*)const_cast<unsigned int*>(X.col_ptrs);
	dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

	//------create B and X-------------------
	if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&sluX, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);

	rhsb = unknown_b;
	dCreate_Dense_Matrix(&sluB, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
	Bstore = (DNformat*)sluB.Store;

	if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");

	/*
	* Get column permutation vector perm_c[], according to permc_spec:
	*   permc_spec = 0: natural ordering
	*   permc_spec = 1: minimum degree ordering on structure of A'*A
	*   permc_spec = 2: minimum degree ordering on structure of A'+A
	*   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	*/
	permc_spec = 1;
	get_perm_c(permc_spec, &sluA, perm_c);

	pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);

	if (info != 0) {
		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = (double*)((DNformat*)sluB.Store)->nzval;
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
			A(node_reorder(i)) = sol[i];
		}
	}
	time[tt++] = SuperLU_timer_();
	//---------------------superLU--end----------------------------------
	//-----------绘图----------------------------------------
	QVector<double> x(num_ele), y(num_ele);
	for (int i = 0; i < num_ele; ++i) {
		x[i] = i;
	}
	QCustomPlot * customplot;
	customplot = thePlot->getQcustomPlot();
	QCPGraph *graph1 = customplot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 3));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	customplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	customplot->xAxis->setLabel("x");
	customplot->xAxis->setRange(0, num_ele);
	//customplot->xAxis->setAutoTickStep(false);
	//customplot->xAxis->setTicks(false);
	customplot->yAxis->setLabel("y");
	//customplot->yAxis->setRange(0, 4);
	//customplot->yAxis->setTicks(false);
	//customplot->xAxis2->setTicks(false);
	//customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);
	//---------the main loop---------------------------------
	int steps = 200;
	int count;
	double alpha = 1;
	Voltage *Vr = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	Voltage *Vi = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	time[tt++] = SuperLU_timer_();
	for (count = 0; count < steps; count++) {
		//------update miu----------------
		for (int i = 0; i < num_ele; i++) {
			double bx = 0;
			double by = 0;
			for (int j = 0; j < 3; j++) {
				bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
				by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
			}
			pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
			pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
			y[i] = pmeshele[i].miut;
			if (std::isnan(y[i])) {
				qDebug() << A(pmeshele[i].n[0]);
				qDebug() << A(pmeshele[i].n[1]);
				qDebug() << A(pmeshele[i].n[2]);
			}
		}
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			//计算从线性网络反射回的电压
			Vr[j].V12 = (pmeshnode[k].A - pmeshnode[m].A) - Vi[j].V12;
			Vr[j].V23 = (pmeshnode[m].A - pmeshnode[n].A) - Vi[j].V23;
			Vr[j].V13 = (pmeshnode[n].A - pmeshnode[k].A) - Vi[j].V13;
			//求解小电路，计算入射向线性网络电压
			double J12, J23, J13, D12, D23, D13;
			
			double dvdB = materialList[m_e->domain - 1].getdvdB(m_e->B);
			dvdB /= ydot[i] * ydot[i];
			double c12, c23, c13,c11,c22,c33;
			c12 = rm[i].Y12*ydot[i] * m_e->miu;
			c23 = rm[i].Y23*ydot[i] * m_e->miu;
			c13 = rm[i].Y13*ydot[i] * m_e->miu;
			c11 = -(c12 + c13);
			c22 = -(c12 + c23);
			c33 = -(c13 + c23);
			double b = m_e->B * ydot[i];
			double miu = m_e->miut*ydot[i];
			double tmp, tmp1;
			tmp = c11*A(k)*A(k) + c22*A(m)*A(m) + c33*A(n)*A(n);
			tmp += 2 * c12*A(k)*A(m) + 2 * c13*A(k)*A(n) + 2 * c23*A(m)*A(n);
			tmp *= dvdB / b / m_e->AREA;
			tmp1 = tmp;
			tmp += 1 / miu;

			J12 = c12*tmp;
			J23 = c23*tmp;
			J13 = c13*tmp;
			
			Vi[j].V12 = Vr[j].V12*(fabs(rm[i].Y12) + J12) + c12*tmp1*(A(m) - A(k)) - c13*tmp1*(A(k) - A(n));
			Vi[j].V12 /= (fabs(rm[i].Y12) - J12);
			Vi[j].V23 = Vr[j].V23*(fabs(rm[i].Y23) + J23) + c23*tmp1*(A(n) - A(m)) - c12*tmp1*(A(m) - A(k));
			Vi[j].V23 /= (fabs(rm[i].Y23) - J23);
			Vi[j].V13 = Vr[j].V13*(fabs(rm[i].Y13) + J13) + c13*tmp1*(A(k) - A(n)) - c23*tmp1*(A(n) - A(m));
			Vi[j].V13 /= (fabs(rm[i].Y13) - J13);
			if (fabs(Vi[j].V12 / Vr[j].V12) > 1) {
				y[i] = 0;
			}
			if (std::isnan(Vi[j].V12)) {
				qDebug() << Vi[j].V12;
			}
			if (std::isnan(Vi[j].V23)) {
				qDebug() << Vi[j].V23;
			}
			if (std::isnan(Vi[j].V13)) {
				qDebug() << Vi[j].V13;
			}
			INL(k) += 2.*Vi[j].V12*fabs(rm[i].Y12) ;
			INL(m) -= 2.*Vi[j].V12*fabs(rm[i].Y12) ;
			
			INL(m) += 2.*Vi[j].V23*fabs(rm[i].Y23);
			INL(n) -= 2.*Vi[j].V23*fabs(rm[i].Y23) ;
			
			INL(n) += 2.*Vi[j].V13*fabs(rm[i].Y13);
			INL(k) -= 2.*Vi[j].V13*fabs(rm[i].Y13) ;
		}
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		if (info != 0) {
			qDebug() << "Error: superlumt.slove";
			//qDebug() << "info: " << superlumt.info;
			break;
		} else {
			double *sol = NULL;
			A_old = A;
			sol = (double*)((DNformat*)sluB.Store)->nzval;

			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
				A(node_reorder(i)) = sol[i];
			}
		}
		double error = norm((A_old - A), 2) / norm(A, 2);
		//qDebug() << "iter: " << count;
		//qDebug() << "error: " << error;

		graph1->setData(x, y);
		customplot->rescaleAxes(true);
		//customplot->xAxis->setRange(0, 0.09);
		//customplot->yAxis->setRange(-0.04, 0.04);
		//customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);
		customplot->replot();


		if (error < Precision) {
			break;
		}
		INL.zeros();
	}
	time[tt++] = SuperLU_timer_();
	// 求解结束，更新B
	for (int i = 0; i < num_ele; i++) {
		double bx = 0;
		double by = 0;
		for (int j = 0; j < 3; j++) {
			bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
			by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
		}
		pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
		pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
	}
	for (int i = 0; i < num_pts - node_bdr; i++) {
		pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;// / pmeshnode[i].x;//the A is r*A_real
	}
	//output the time
	for (int i = 1; i < tt; i++) {
		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vi != NULL) free(Vi);
	if (Vr != NULL) free(Vr);
	if (unknown_b != NULL) free(unknown_b);

	//SUPERLU_FREE(rhs);
	//SUPERLU_FREE(xact);
	//SUPERLU_FREE(perm_r);
	//SUPERLU_FREE(perm_c);
	//Destroy_CompCol_Matrix(&sluA);
	//Destroy_SuperMatrix_Store(&sluB);
	//Destroy_SuperNode_SCP(&L);
	//Destroy_CompCol_NCP(&U);
	//StatFree(&Gstat1);
	return true;
}

double CFastFEMcore::CalcForce() {
	/*******2016-12-28 by Poofee*************/
	/********Find the first layer***********/
	double xLeft = 1.25e-3;
	double xRight = 5.8e-3;
	double yDown = -4.7e-3;
	double yUp = 5.3e-3;
	double delta = 1e-10;
	//	double xForce = 0;
	double yForce = 0;
	yUp += ste * 1e-4;
	yDown += ste * 1e-4;

	QCustomPlot *customPlot;
	customPlot = thePlot->getQcustomPlot();
	customPlot->xAxis->setLabel("x");
	customPlot->xAxis->setRange(0, 0.09);
	customPlot->xAxis->setAutoTickStep(false);
	customPlot->xAxis->setTicks(false);
	customPlot->yAxis->setLabel("y");
	customPlot->yAxis->setRange(-0.09, 0.09);
	customPlot->xAxis2->setTicks(false);
	customPlot->yAxis->setScaleRatio(customPlot->xAxis, 1.0);


	customPlot->yAxis->setAutoTickStep(false);
	customPlot->yAxis->setAutoTickLabels(false);
	customPlot->yAxis->setTicks(false);
	customPlot->yAxis->grid()->setVisible(false);
	customPlot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
	customPlot->yAxis2->setTicks(false);
	customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

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
	cc[0] = QColor(0, 0, 0);
	cc[1] = QColor(0, 120, 0);
	cc[2] = QColor(68, 49, 242);
	cc[3] = QColor(254, 31, 87);
	cc[4] = QColor(255, 128, 0);
	cc[5] = QColor(232, 223, 60);
	cc[6] = QColor(232, 223, 60);
	customPlot->xAxis->grid()->setSubGridVisible(false);
	customPlot->xAxis->grid()->setSubGridPen(Qt::NoPen);
	customPlot->yAxis->grid()->setSubGridVisible(false);
	customPlot->yAxis->grid()->setSubGridPen(Qt::NoPen);
	FILE * fp = NULL;
	fp = fopen("E:\\index.txt", "w+");//delete exist, read and write
	for (int i = 0; i < num_ele; i++) {
		QCPCurve *newCurve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
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
		if (pmeshele[i].domain == 2) {
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
			newCurve->setBrush(Qt::NoBrush);
		} else if (BCy == 2 || BCy == 6) {//Y,2,edge
			newCurve->setBrush(QColor(0, 0, 255));//blue
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
			newCurve->setBrush(QColor(255, 0, 0));
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
			newCurve->setBrush(QColor(0, 0, 255));//blue
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
			newCurve->setBrush(QColor(255, 0, 0));
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
			newCurve->setBrush(QColor(0, 255, 0));//green
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
			double tmp = pmeshele[i].rc * pmeshele[i].AREA / (4 * PI*1e-7);
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

		//if (pmeshele[i].domain != 2 && pmeshele[i].domain != 1) {
		newCurve->setData(x1, y1);
		// newCurve->setPen(QPen(cc[pmeshele[i].domain - 1]));
		//newCurve->setBrush(cc[pmeshele[i].domain - 1]);
		//}
		//if (ind != 10) {
		//	fprintf(fp, "%d\n", ind);
		//}

	}
	graph1->setData(gpx1, gpy1);
	graph2->setData(gpx2, gpy2);
	fclose(fp);
	qDebug() << "yForce: " << yForce;
	customPlot->replot();
	return 0;
}

//打开工程文件，读取
int CFastFEMcore::openProject(QString proFile) {
	QFile ifile(proFile);
	if (ifile.open(QIODevice::ReadOnly | QIODevice::Text)) {

		QXmlStreamReader reader(&ifile);

		while (!reader.atEnd()) {
			if (reader.isStartElement()) {
				if (reader.name() == "Project") {
					readProjectElement(reader);
				} else {
					reader.raiseError("Not a valid Project file");
				}
			} else {
				reader.readNext();
			}
		}
	} else {
		qDebug() << "read inbox file error...";
	}
	ifile.close();
	return 0;
}


int CFastFEMcore::preCalculation() {
	int k, m, n;

	for (int i = 0; i < num_ele; i++) {
		k = pmeshele[i].n[0];
		m = pmeshele[i].n[1];
		n = pmeshele[i].n[2];
		pmeshele[i].P[0] = pmeshnode[m].y - pmeshnode[n].y;
		pmeshele[i].P[1] = pmeshnode[n].y - pmeshnode[k].y;
		pmeshele[i].P[2] = pmeshnode[k].y - pmeshnode[m].y;

		pmeshele[i].Q[0] = pmeshnode[n].x - pmeshnode[m].x;
		pmeshele[i].Q[1] = pmeshnode[k].x - pmeshnode[n].x;
		pmeshele[i].Q[2] = pmeshnode[m].x - pmeshnode[k].x;

		pmeshele[i].AREA = 0.5*abs(pmeshele[i].P[1] * pmeshele[i].Q[2] - pmeshele[i].Q[1] * pmeshele[i].P[2]);
		pmeshele[i].rc = (pmeshnode[k].x +
			pmeshnode[m].x +
			pmeshnode[n].x) / 3;
		pmeshele[i].zc = (pmeshnode[k].y +
			pmeshnode[m].y +
			pmeshnode[n].y) / 3;

		//主要根据材料属性完成单元当中miu,miut,的赋值；

		//由于I,pm与形函数有关系，为实现分离，不在此计算

		//由于I,pm与形函数有关系，为实现分离，不在此计算

		if (materialList[pmeshele[i].domain - 1].BHpoints == 0) {
			pmeshele[i].miu = materialList[pmeshele[i].domain - 1].miu;
			pmeshele[i].miut = materialList[pmeshele[i].domain - 1].miu;//must be 1
			pmeshele[i].LinearFlag = true;
		} else {
			pmeshele[i].miu = 1 * miu0;
			pmeshele[i].miut = 100 * miu0;
			pmeshele[i].LinearFlag = false;
		}
		//查找物理边界
		if (pmeshele[i].domain == 1) {
			pmeshnode[k].bdr = 3;

			pmeshnode[m].bdr = 3;

			pmeshnode[n].bdr = 3;

		} else {
			//对称轴上的部分边界点
			if (pmeshnode[k].bdr == 1 && pmeshnode[k].x < 1e-8) {
				pmeshnode[k].bdr = 3;
			}
			if (pmeshnode[m].bdr == 1 && pmeshnode[m].x < 1e-8) {
				pmeshnode[m].bdr = 3;
			}
			if (pmeshnode[n].bdr == 1 && pmeshnode[n].x < 1e-8) {
				pmeshnode[n].bdr = 3;
			}
		}
	}

	return 0;
}

//调用其他的子函数完成求解任务，这个求解可以有多个选项，
//使用NR或者TLM迭代算法，或者选择不同的形函数。
int CFastFEMcore::solve() {
	//openProject();
	//LoadMesh();
	preCalculation();
	StaticAxisymmetricTLM();
	return 0;
}


void CFastFEMcore::readProjectElement(QXmlStreamReader &reader) {
	Q_ASSERT(reader.isStartElement() && reader.name() == "Project");
	reader.readNext();
	while (!reader.atEnd()) {
		if (reader.isEndElement()) {
			//qDebug()<<reader.name();
			//break;
		}

		if (reader.isStartElement()) {
			if (reader.name() == "name") {
				reader.readElementText();
				//qDebug() << "name = " << reader.readElementText();
			} else if (reader.name() == "version") {
				reader.readElementText();
				//qDebug() << "version = " << reader.readElementText();
			} else if (reader.name() == "precision") {
				Precision = reader.readElementText().toDouble();
				//qDebug() << "precision = " << Precision;
			} else if (reader.name() == "unit") {
				LengthUnits = reader.readElementText().toInt();
				//qDebug() << "unit = " << LengthUnits;
			} else if (reader.name() == "proType") {
				reader.readElementText();
				//qDebug() << "proType = " << reader.readElementText();
			} else if (reader.name() == "coordinate") {
				reader.readElementText();
				//qDebug() << "coordinate = " << reader.readElementText();
			} else if (reader.name() == "Domains") {
				reader.readNextStartElement();
				if (reader.name() == "domainNum") {
					numDomain = reader.readElementText().toInt();
					materialList = new CMaterial[numDomain];
					qDebug() << "domainNum = " << numDomain;
				}
				for (int i = 0; i < numDomain; i++) {
					readDomainElement(reader, i);
				}

			}
		} else {
			reader.readNext();
		}
	}
}


void CFastFEMcore::readDomainElement(QXmlStreamReader &reader, int i) {
	reader.readNext();
	reader.readNext();
	//qDebug()<<reader.name();
	while (!(reader.isEndElement() && reader.name() == "Domain")) {
		reader.readNextStartElement();
		if (reader.name() == "domainName") {
			reader.readElementText();
			qDebug() << "domainName = " << reader.readElementText();
		} else if (reader.name() == "miu") {
			materialList[i].miu = reader.readElementText().toDouble() * 4 * PI*1e-7;
			qDebug() << "miu = " << materialList[i].miu;
		} else if (reader.name() == "BH") {
			readBHElement(reader, i);
		} else if (reader.name() == "Jr") {
			materialList[i].Jr = reader.readElementText().toDouble();
			qDebug() << "Jr = " << materialList[i].Jr;
		} else if (reader.name() == "H_c") {
			materialList[i].H_c = reader.readElementText().toDouble();
			qDebug() << "H_c = " << materialList[i].H_c;
		}
	}
}


void CFastFEMcore::readBHElement(QXmlStreamReader &reader, int i) {
	reader.readNextStartElement();
	//qDebug()<<reader.name();
	if (reader.name() == "BHpoints") {
		materialList[i].BHpoints = reader.readElementText().toInt();
		qDebug() << "BHpoints = " << materialList[i].BHpoints;
		if (materialList[i].BHpoints != 0) {
			materialList[i].Bdata = (double*)malloc(materialList[i].BHpoints*sizeof(double));
			materialList[i].Hdata = (double*)malloc(materialList[i].BHpoints*sizeof(double));
			reader.readNextStartElement();
			if (reader.name() == "Bdata") {
				for (int j = 0; j < materialList[i].BHpoints; j++) {
					reader.readNextStartElement();
					if (reader.name() == "B") {
						materialList[i].Bdata[j] = reader.readElementText().toDouble();
						//qDebug() << materialList[i].Bdata[j];
					}
				}
			}
			//qDebug()<<reader.name();
			reader.readNextStartElement();
			//qDebug()<<reader.name();
			reader.readNextStartElement();
			//qDebug()<<reader.name();
			if (reader.name() == "Hdata") {
				for (int j = 0; j < materialList[i].BHpoints; j++) {
					reader.readNextStartElement();
					if (reader.name() == "H") {
						materialList[i].Hdata[j] = reader.readElementText().toDouble();
						//qDebug() << materialList[i].Hdata[j];
					}

				}
			}
		} else {
			materialList[i].Bdata = NULL;
			materialList[i].Hdata = NULL;
		}

	}
	reader.skipCurrentElement();
}

//使用牛顿迭代实现非线性求解，是真的牛顿迭代而不是松弛迭代
//牛顿迭代的公式推导，参见颜威利数值分析教材P54
int CFastFEMcore::StaticAxisymmetricNR() {
	clock_t time[10];
	int tt = 0;
	time[tt++] = clock();
	//所需要的变量
	umat locs(2, 9 * num_ele);
	locs.zeros();
	mat vals(1, 9 * num_ele);
	double ce[3][3] = { 0 };
	double cn[3][3] = { 0 };//用于牛顿迭代
	vec bbJz = zeros<vec>(num_pts);
	uvec node_reorder = zeros<uvec>(num_pts);
	uvec node_pos = zeros<uvec>(num_pts);
	vec bn = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	double * ydot = (double*)malloc(num_ele*sizeof(double));
	ResistMarix *rm = (ResistMarix*)malloc(num_ele * sizeof(ResistMarix));

	QVector<double> x(num_ele);
	QVector<double> y(num_ele);
	QCustomPlot * customplot;
	customplot = thePlot->getQcustomPlot();
	QCPGraph *graph1 = customplot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 1.5), QBrush(Qt::black), 3));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	customplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	customplot->xAxis->setLabel("x");
	customplot->xAxis->setRange(0, num_ele);
	//customplot->xAxis->setAutoTickStep(false);
	//customplot->xAxis->setTicks(false);
	customplot->yAxis->setLabel("y");
	customplot->yAxis->setRange(0, 4);
	//customplot->xAxis2->setTicks(false);
	//customplot->yAxis->setScaleRatio(ui->widget->xAxis, 1.0);
	for (int i = 0; i < num_ele; i++) {
		x[i] = i;
	}
	//重新对节点进行编号，将边界点分离
	int node_bdr = 0;
	for (int i = 0; i < num_pts; i++) {
		if (pmeshnode[i].bdr == 3) {
			node_bdr++;
			node_reorder(num_pts - node_bdr) = i;
			node_pos(i) = num_pts - node_bdr;
			pmeshnode[i].A = 0;
			A(i) = 0;
		} else {
			node_reorder(i - node_bdr) = i;
			node_pos(i) = i - node_bdr;
		}
	}
	double* unknown_b = (double*)calloc(num_pts - node_bdr, sizeof(double));
	int iter = 0;//迭代步数
	int pos = 0;
	while (1) {
		//生成全局矩阵
		for (int i = 0; i < num_ele; i++) {
			//这部分只需计算一次即可
			if (iter == 0) {
				int flag = 0;
				for (int f = 0; f < 3; f++) {
					if (pmeshnode[pmeshele[i].n[f]].x < 1e-7) {
						flag++;
					}
				}
				//计算三角形重心半径
				if (flag == 2) {
					ydot[i] = pmeshele[i].rc;
				} else {
					ydot[i] = 1 / (pmeshnode[pmeshele[i].n[0]].x + pmeshnode[pmeshele[i].n[1]].x);
					ydot[i] += 1 / (pmeshnode[pmeshele[i].n[0]].x + pmeshnode[pmeshele[i].n[2]].x);
					ydot[i] += 1 / (pmeshnode[pmeshele[i].n[1]].x + pmeshnode[pmeshele[i].n[2]].x);
					ydot[i] = 1.5 / ydot[i];
				}
				//计算上三角矩阵，这些系数在求解过程当中是不变的
				rm[i].Y11 = pmeshele[i].Q[0] * pmeshele[i].Q[0] + pmeshele[i].P[0] * pmeshele[i].P[0];
				rm[i].Y12 = pmeshele[i].Q[0] * pmeshele[i].Q[1] + pmeshele[i].P[0] * pmeshele[i].P[1];
				rm[i].Y13 = pmeshele[i].Q[0] * pmeshele[i].Q[2] + pmeshele[i].P[0] * pmeshele[i].P[2];
				rm[i].Y22 = pmeshele[i].Q[1] * pmeshele[i].Q[1] + pmeshele[i].P[1] * pmeshele[i].P[1];
				rm[i].Y23 = pmeshele[i].Q[1] * pmeshele[i].Q[2] + pmeshele[i].P[1] * pmeshele[i].P[2];
				rm[i].Y33 = pmeshele[i].Q[2] * pmeshele[i].Q[2] + pmeshele[i].P[2] * pmeshele[i].P[2];

				rm[i].Y11 /= 4. * pmeshele[i].AREA;
				rm[i].Y12 /= 4. * pmeshele[i].AREA;
				rm[i].Y13 /= 4. * pmeshele[i].AREA;
				rm[i].Y22 /= 4. * pmeshele[i].AREA;
				rm[i].Y23 /= 4. * pmeshele[i].AREA;
				rm[i].Y33 /= 4. * pmeshele[i].AREA;

				//计算电流密度//要注意domain会不会越界
				double jr = pmeshele[i].AREA*materialList[pmeshele[i].domain - 1].Jr / 3;
				for (int j = 0; j < 3; j++) {
					bbJz(pmeshele[i].n[j]) += jr;
					// 计算永磁部分
					bbJz(pmeshele[i].n[j]) -= materialList[pmeshele[i].domain - 1].H_c / 2.*pmeshele[i].Q[j];
				}
			}//end of iter=0
			//miut对于线性就等于真值，对于非线性等于上一次的值
			//主要求解结果不要漏掉miu0

			ce[0][0] = rm[i].Y11 / pmeshele[i].miut / ydot[i];
			ce[1][1] = rm[i].Y22 / pmeshele[i].miut / ydot[i];
			ce[2][2] = rm[i].Y33 / pmeshele[i].miut / ydot[i];

			ce[0][1] = rm[i].Y12 / pmeshele[i].miut / ydot[i];
			ce[0][2] = rm[i].Y13 / pmeshele[i].miut / ydot[i];
			ce[1][2] = rm[i].Y23 / pmeshele[i].miut / ydot[i];

			//计算牛顿迭代部分的单元矩阵项,如果是第一次迭代的话，A=0，
			//所以就不计算了，参见颜威利书P56
			double v[3];

			v[0] = rm[i].Y11*A(pmeshele[i].n[0]) +
				rm[i].Y12*A(pmeshele[i].n[1]) +
				rm[i].Y13*A(pmeshele[i].n[2]);
			v[1] = rm[i].Y12*A(pmeshele[i].n[0]) +
				rm[i].Y22*A(pmeshele[i].n[1]) +
				rm[i].Y23*A(pmeshele[i].n[2]);
			v[2] = rm[i].Y13*A(pmeshele[i].n[0]) +
				rm[i].Y23*A(pmeshele[i].n[1]) +
				rm[i].Y33*A(pmeshele[i].n[2]);

			if (iter != 0) {
				double tmp;
				if (pmeshele[i].LinearFlag) {
					tmp = 0;
				} else {
					tmp = materialList[pmeshele[i].domain - 1].getdvdB(pmeshele[i].B);
					if (pmeshele[i].B > 1e-9){
						tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
						tmp /= ydot[i] * ydot[i] * ydot[i];
					}

				}
				cn[0][0] = v[0] * v[0] * tmp;
				cn[1][1] = v[1] * v[1] * tmp;
				cn[2][2] = v[2] * v[2] * tmp;

				cn[0][1] = v[0] * v[1] * tmp;
				cn[0][2] = v[0] * v[2] * tmp;
				cn[1][2] = v[1] * v[2] * tmp;

				cn[1][0] = cn[0][1];
				cn[2][0] = cn[0][2];
				cn[2][1] = cn[1][2];
			}
			ce[0][0] += cn[0][0];
			ce[1][1] += cn[1][1];
			ce[2][2] += cn[2][2];

			ce[0][1] += cn[0][1];
			ce[0][2] += cn[0][2];
			ce[1][2] += cn[1][2];

			ce[1][0] = ce[0][1];
			ce[2][0] = ce[0][2];
			ce[2][1] = ce[1][2];

			for (int row = 0; row < 3; row++) {
				for (int col = 0; col < 3; col++) {
					//判断节点是否在未知节点内
					//得到排序之后的编号
					int n_row = node_pos(pmeshele[i].n[row]);
					int n_col = node_pos(pmeshele[i].n[col]);
					if (n_row < num_pts - node_bdr && n_col < num_pts - node_bdr) {
						locs(0, pos) = n_row;
						locs(1, pos) = n_col;
						vals(0, pos) = ce[row][col];
						pos++;
					}
					//bn的顺序并没有改
					bn(pmeshele[i].n[row]) += cn[row][col] * A(pmeshele[i].n[col]);
				}
			}
		}//end of elememt iteration
		if (iter == 0) {
			locs.reshape(2, pos);
			vals.reshape(1, pos);
		}
		bn += bbJz;
		//使用构造函数来生成稀疏矩阵
		sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = bn(node_reorder(i));
		}
		//---------------------superLU_MT---------------------------------------
		CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
		if (superlumt.solve() == 1) {
			qDebug() << "Error: superlumt.slove";
			qDebug() << "info: " << superlumt.info;
			break;
		} else {
			double *sol = NULL;
			A_old = A;
			sol = superlumt.getResult();

			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
				A(node_reorder(i)) = sol[i];
			}
		}
		//有必要求出所有单元的B值
		for (int i = 0; i < num_ele; i++) {
			double bx = 0;
			double by = 0;
			for (int j = 0; j < 3; j++) {
				bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
				by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
			}
			pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
			pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);

			y[i] = pmeshele[i].miut;
		}
		double error = norm((A_old - A), 2) / norm(A, 2);
		iter++;
		//qDebug() << "iter: " << iter;
		//qDebug() << "error: " << error;
		if (error < Precision || iter > 20) {
			A.save("NRA.txt", arma::arma_ascii, false);
			break;
		}
		bn.zeros();
		pos = 0;

		graph1->setData(x, y);
		customplot->rescaleAxes(true);
		customplot->replot();
	}
	time[tt++] = clock();
	for (int i = 1; i < tt; i++){
		qDebug() << time[i] - time[i - 1];
	}
	qDebug() << "NR steps: " << iter;
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (unknown_b != NULL) free(unknown_b);
	return 0;
}

int CFastFEMcore::LoadQ4MeshCOMSOL(const char fn[]){
	char ch[256];
	//------------open file----------------------------------
	FILE * fp = NULL;
	fp = fopen(fn, "r");
	if (fp == NULL) {
		qDebug() << "Error: openning file!";
		return 1;
	}
	//--------------Read the head-----------------------------
	for (int i = 0; i < 18; i++) {
		fgets(ch, 256, fp);
	}
	//-----------------mesh point-----------------------------
	//读取节点数目
	if (fscanf(fp, "%d # number of mesh points\n", &num_pts)) {
		pmeshnode = (CNode*)calloc(num_pts, sizeof(CNode));

		for (int i = 0; i < num_pts; i++) {
			pmeshnode[i].I = 0;
			pmeshnode[i].pm = 0;
		}
	} else {
		qDebug() << "Error: reading num_pts!";
		return 1;
	}
	int pts_ind;//the beginning of the points index
	//读取节点索引，默认从0开始
	if (fscanf(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
		qDebug() << "Error: reading pts_ind!";
		return 1;
	}
	fgets(ch, 256, fp);

	for (int i = pts_ind; i < num_pts; i++) {
		//读取x,y坐标
		if (fscanf(fp, "%lf %lf \n", &(pmeshnode[i].x), &(pmeshnode[i].y)) != 2) {
			qDebug() << "Error: reading mesh point!";
			return 1;
		}
	}
	//---------------vertex-------------------------------
	for (int i = 0; i < 7; i++)
		fgets(ch, 256, fp);
	int num_vtx_ns, num_vtx_ele;
	//
	if (fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns) != 1) {
		qDebug() << "Error: reading num_vtx_ns!";
		return 1;
	}

	if (fscanf(fp, "%d # number of elements\n", &num_vtx_ele) != 1) {
		qDebug() << "Error: reading num_vtx_ele!";
		return 1;
	}
	fgets(ch, 256, fp);

	int *vtx;
	vtx = (int*)calloc(num_vtx_ele, sizeof(int));
	for (int i = 0; i < num_vtx_ele; i++) {
		//好象是每一个域的顶点编号
		if (fscanf(fp, "%d \n", vtx + i) != 1) {
			qDebug() << "Error: reading vertex condition!";
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
			qDebug() << "Error: reading vertex condition!";
			return 1;
		}
	}
	if (vtx2 != NULL) free(vtx2); vtx2 = NULL;
	//--------------boundary--------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int num_bdr_ns, num_bdr_ele;//number of nodes per element;number of elements
	//读取一个边界单元中的数目，2D的话为2，表示线段
	if (fscanf(fp, "%d # number of nodes per element\n", &num_bdr_ns) != 1) {
		qDebug() << "Error: reading num_bdr_ns!";
		return 1;
	}
	//读取线段边界数目
	if (fscanf(fp, "%d # number of elements\n", &num_bdr_ele) != 1) {
		qDebug() << "Error: reading num_bdr_ele!";
		return 1;
	}
	fgets(ch, 256, fp);

	int *p1, *p2;
	p1 = (int*)calloc(num_bdr_ele, sizeof(int));
	p2 = (int*)calloc(num_bdr_ele, sizeof(int));
	for (int i = 0; i < num_bdr_ele; i++) {
		//读取线段边界的起点和终点
		if (fscanf(fp, "%d %d\n", p1 + i, p2 + i) == 2) {
			pmeshnode[p1[i]].bdr = 1;
		} else {
			qDebug() << "Error: reading boundary condition!";
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
			qDebug() << "Error: reading boundary condition!";
			return 1;
		}
	}
	if (entity != NULL) free(entity); entity = NULL;
	//----------------elements------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int ns_per_ele;// num_ele;//number of nodes per element;number of elements
	if (fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele) != 1) {
		qDebug() << "Error: reading ns_per_ele!";
		return 1;
	}
	//读取分网单元数目
	if (fscanf(fp, "%d # number of elements\n", &num_ele) == 1) {
		pmeshele4 = (CElement4*)calloc(num_ele, sizeof(CElement4));
	} else {
		qDebug() << "Error: reading num_ele!";
		return 1;
	}
	fgets(ch, 256, fp);
	//读取分网四边形单元的四个节点索引
	//注意这里有点问题，comsol文件似乎不是按照逆时针存储的，而是1->2->4->3
	for (int i = 0; i < num_ele; i++) {
		if (fscanf(fp, "%d %d %d %d\n", &pmeshele4[i].n[0], &pmeshele4[i].n[1], &pmeshele4[i].n[3], &pmeshele4[i].n[2]) != 4) {
			qDebug() << "Error: reading elements points!";
			return 1;
		}
	}
	//---------------Domain----------------------------------
	int num_domain;
	//读取domain数目
	fscanf(fp, "%d # number of geometric entity indices\n", &num_domain);
	fgets(ch, 256, fp);

	for (int i = 0; i < num_domain; i++) {
		//读取每个单元所在的domain
		if (fscanf(fp, "%d \n", &pmeshele4[i].domain) != 1) {
			qDebug() << "Error: reading domain points!";
			return 1;
		}
	}
	fclose(fp);
	return 0;
}
//Ki,start from 0
double CFastFEMcore::getLocal4Matrix(int Ki, int Kj, int index){
	int gaussPoints = 3;
	double * gaussweight;
	double * gausspoint;
	double Kij = 0;
	if (gaussPoints == 3){
		gaussweight = GaussWeight3;
		gausspoint = GaussPoint3;
	} else if (gaussPoints == 5){
		gaussweight = GaussWeight5;
		gausspoint = GaussPoint5;
	}
	for (int i = 0; i < gaussPoints; i++){
		for (int j = 0; j < gaussPoints; j++){
			Kij += gaussweight[i] * gaussweight[j] * getP(Ki, Kj, gausspoint[i], gausspoint[j], index);
		}
	}
	//if (std::isnan(Kij))
	//	qDebug() << Kij;
	return Kij;
}
double CFastFEMcore::getTLMKij(int Ki, int Kj, int index) {
	int gaussPoints = 3;
	double * gaussweight;
	double * gausspoint;
	double Kij = 0;
	if (gaussPoints == 3) {
		gaussweight = GaussWeight3;
		gausspoint = GaussPoint3;
	} else if (gaussPoints == 5) {
		gaussweight = GaussWeight5;
		gausspoint = GaussPoint5;
	}
	for (int i = 0; i < gaussPoints; i++) {
		for (int j = 0; j < gaussPoints; j++) {
			Kij += gaussweight[i] * gaussweight[j] * getPx(Ki, Kj, gausspoint[i], gausspoint[j], index);
		}
	}
	return Kij;
}
double CFastFEMcore::getJi(int Ki, int index){
	int gaussPoints = 3;
	double * gaussweight;
	double * gausspoint;
	double Ji = 0;
	if (gaussPoints == 3){
		gaussweight = GaussWeight3;
		gausspoint = GaussPoint3;
	} else if (gaussPoints == 5){
		gaussweight = GaussWeight5;
		gausspoint = GaussPoint5;
	}
	for (int i = 0; i < gaussPoints; i++){
		for (int j = 0; j < gaussPoints; j++){
			Ji += gaussweight[i] * gaussweight[j] * Ne(gausspoint[i], gausspoint[j], Ki)* getJacobi(gausspoint[i], gausspoint[j], index);
		}
	}
	return Ji;
}
double CFastFEMcore::getD(int i, int j, double xi, double eta, int index){
	double D = 0;
	//计算x
	double x = getx(xi, eta, index);
	if (x == 0){
		x = pmeshele4[index].rc;
	}
	//计算B
	double bx = 0;
	double by = 0;
	double tmpA = 0;
	//by += 1 / getx(xi, eta, index);
	//by *= getA(xi, eta, index);//考虑到x可能为0，所以后乘A，
	double c1, c2;
	c1 = 0;
	c2 = 0;
	for (int iele = 0; iele < 4; iele++){
		bx += -pmeshnode[pmeshele4[index].n[iele]].A * getdNidy(iele, xi, eta, index); //qDebug() << bx;

		by += pmeshnode[pmeshele4[index].n[iele]].A * getdNidx(iele, xi, eta, index);

		c1 += getCij(i, iele, xi, eta, index)*pmeshnode[pmeshele4[index].n[iele]].A;
		c2 += getCij(j, iele, xi, eta, index)*pmeshnode[pmeshele4[index].n[iele]].A;
	}

	double B = sqrt(bx*bx + by*by);
	B /= getx(xi, eta, index);
	//计算dvdB
	double dvdB = materialList[pmeshele4[index].domain - 1].getdvdB(B);
	if (std::isnan(dvdB))
		qDebug() << "error in getD";
	//计算Cij
	//积分
	if (pmeshele4[index].LinearFlag){
		return 0;
	}

	if (std::isnan(D))
		qDebug() << "error in getD";
	D = dvdB*c1*c2 * getJacobi(xi, eta, index);
	if (B > 1e-9){
		D /= B * x * x * x;
	}
	if (std::isnan(D))
		qDebug() << "error in getD";
	return D;
}
double CFastFEMcore::getCij(int Ki, int Kj, double xi, double eta, int index){
	double Cij = 0;
	Cij = getdNidx(Ki, xi, eta, index)*getdNidx(Kj, xi, eta, index);
	Cij += getdNidy(Ki, xi, eta, index)*getdNidy(Kj, xi, eta, index);
	return Cij;
}
double CFastFEMcore::getDij(int Ki, int Kj, int index){
	int gaussPoints = 3;
	double * gaussweight;
	double * gausspoint;
	double Dij = 0;
	if (gaussPoints == 3){
		gaussweight = GaussWeight3;
		gausspoint = GaussPoint3;
	} else if (gaussPoints == 5){
		gaussweight = GaussWeight5;
		gausspoint = GaussPoint5;
	}
	for (int i = 0; i < gaussPoints; i++){
		for (int j = 0; j < gaussPoints; j++){
			Dij += gaussweight[i] * gaussweight[j] * getD(Ki, Kj, gausspoint[i], gausspoint[j], index);
		}
	}
	return Dij;
}
double CFastFEMcore::getPx(int Ki, int Kj, double xi, double eta, int index) {
	double p = 0;
	p = getdNidx(Ki, xi, eta, index)*getdNidx(Kj, xi, eta, index);
	p += getdNidy(Ki, xi, eta, index)*getdNidy(Kj, xi, eta, index);
	p *= getJacobi(xi, eta, index);	
	p /=  getx(xi, eta, index);//可能为0的bug
	return p;
}
double CFastFEMcore::getP(int Ki, int Kj, double xi, double eta, int index){
	double p = 0;
	p = getdNidx(Ki, xi, eta, index)*getdNidx(Kj, xi, eta, index);
	p += getdNidy(Ki, xi, eta, index)*getdNidy(Kj, xi, eta, index);
	p *= getJacobi(xi, eta, index);
	//应当考虑r*\nu部分，
	double bx = 0;
	double by = 0;
	double tmpA = 0;
	//by += 1 / getx(xi, eta, index);
	//by *= getA(xi,eta,index);//考虑到x可能为0，所以后乘A，
	for (int iele = 0; iele < 4; iele++){
		bx += -pmeshnode[pmeshele4[index].n[iele]].A * getdNidy(iele, xi, eta, index); //qDebug() << bx;

		by += pmeshnode[pmeshele4[index].n[iele]].A * getdNidx(iele, xi, eta, index);
	}

	double B = sqrt(bx*bx + by*by);
	B /= getx(xi, eta, index);
	double miu = materialList[pmeshele4[index].domain - 1].getMiu(B);
	p /= miu * getx(xi, eta, index);//可能为0的bug
	if (std::isnan(p))
		qDebug() << B << miu << getx(xi, eta, index) << by << p;
	return p;
}

double CFastFEMcore::getdNidx(int Ki, double xi, double eta, int index){
	double result = 0;
	result = getdydeta(xi, index)*getdNdxi(Ki, eta) - getdydxi(eta, index)*getdNdeta(Ki, xi);
	result /= getJacobi(xi, eta, index);
	return result;
}
double CFastFEMcore::getdNidy(int Ki, double xi, double eta, int index){
	double result = 0;
	result = -getdxdeta(xi, index)*getdNdxi(Ki, eta) + getdxdxi(eta, index)*getdNdeta(Ki, xi);
	result /= getJacobi(xi, eta, index);
	return result;
}
double CFastFEMcore::getJacobi(double xi, double eta, int index){
	double tmp = getdxdxi(eta, index)*getdydeta(xi, index) - getdydxi(eta, index)*getdxdeta(xi, index);
	//qDebug() << tmp;
	if ((tmp) == 0)
		qDebug() << tmp;
	return fabs(tmp);
}

double CFastFEMcore::getdxdeta(double xi, int index){
	double sum = 0;
	for (int i = 0; i < 4; i++){
		sum += getdNdeta(i, xi)*pmeshnode[pmeshele4[index].n[i]].x;
	}
	return sum;
}
double CFastFEMcore::getdxdxi(double eta, int index){
	double sum = 0;
	for (int i = 0; i < 4; i++){
		sum += getdNdxi(i, eta)*pmeshnode[pmeshele4[index].n[i]].x;
	}
	return sum;
}
double CFastFEMcore::getdydeta(double xi, int index){
	double sum = 0;
	for (int i = 0; i < 4; i++){
		sum += getdNdeta(i, xi)*pmeshnode[pmeshele4[index].n[i]].y;
	}
	return sum;
}
double CFastFEMcore::getdydxi(double eta, int index){
	double sum = 0;
	for (int i = 0; i < 4; i++){
		sum += getdNdxi(i, eta)*pmeshnode[pmeshele4[index].n[i]].y;
	}
	return sum;
}

double CFastFEMcore::getdNdeta(int i, double xi){
	double xxi[] = { -1, 1, 1, -1 };
	double yeta[] = { -1, -1, 1, 1 };
	return yeta[i] * (1 + xxi[i] * xi)*0.25;
}

double CFastFEMcore::getdNdxi(int i, double eta){
	double xxi[] = { -1, 1, 1, -1 };
	double yeta[] = { -1, -1, 1, 1 };
	return xxi[i] * (1 + yeta[i] * eta)*0.25;
}
//根据等参单元法求得xi,eta与A的对应关系
double CFastFEMcore::getA(double xi, double eta, int index){
	double A = 0;
	//教材里插值针对的是A‘，而不是A
	for (int i = 0; i < 4; i++){
		A += pmeshnode[pmeshele4[index].n[i]].A*Ne(xi, eta, i);
	}
	A /= getx(xi, eta, index);
	return A;
}
//根据等参单元法求得xi,eta与rA的对应关系
double CFastFEMcore::getrA(double xi, double eta, int index){
	double A = 0;
	//这里得到的是rA
	for (int i = 0; i < 4; i++){
		A += pmeshnode[pmeshele4[index].n[i]].A*Ne(xi, eta, i);
	}
	return A;
}
//根据等参单元法求得xi,eta与x的对应关系
double CFastFEMcore::getx(double xi, double eta, int index){
	double x = 0;
	for (int i = 0; i < 4; i++){
		x += pmeshnode[pmeshele4[index].n[i]].x*Ne(xi, eta, i);
	}
	return x;
}
//根据等参单元法求得xi,eta与y的对应关系
double CFastFEMcore::gety(double xi, double eta, int index){
	double y = 0;
	for (int i = 0; i < 4; i++){
		y += pmeshnode[pmeshele4[index].n[i]].y*Ne(xi, eta, i);
	}
	return y;
}
double CFastFEMcore::Ne(double xi, double eta, int index){
	double xxi[] = { -1, 1, 1, -1 };
	double yeta[] = { -1, -1, 1, 1 };
	return (1 + yeta[index] * eta)*(1 + xxi[index] * xi)*0.25;
}
bool CFastFEMcore::StaticAxisQ4NR(){
	//计算边界信息	
	for (int i = 0; i < num_ele; i++){
		int k, l, m, n;
		k = pmeshele4[i].n[0];
		l = pmeshele4[i].n[1];
		m = pmeshele4[i].n[2];
		n = pmeshele4[i].n[3];

		pmeshele4[i].rc = pmeshnode[k].x + pmeshnode[l].x + pmeshnode[m].x + pmeshnode[n].x;
		pmeshele4[i].rc /= 4;

		if (materialList[pmeshele4[i].domain - 1].BHpoints == 0) {
			pmeshele4[i].miu = materialList[pmeshele4[i].domain - 1].miu;
			pmeshele4[i].miut = materialList[pmeshele4[i].domain - 1].miu;//must be 1
			pmeshele4[i].LinearFlag = true;
		} else {
			pmeshele4[i].miu = 1 * miu0;
			pmeshele4[i].miut = 100 * miu0;
			pmeshele4[i].LinearFlag = false;
		}
		//查找物理边界
		if (pmeshele4[i].domain == 1) {
			pmeshnode[k].bdr = 3;

			pmeshnode[l].bdr = 3;

			pmeshnode[m].bdr = 3;

			pmeshnode[n].bdr = 3;
		} else {
			//对称轴上的部分边界点
			if (pmeshnode[k].bdr == 1 && pmeshnode[k].x < 1e-8) {
				pmeshnode[k].bdr = 3;
			}
			if (pmeshnode[l].bdr == 1 && pmeshnode[l].x < 1e-8) {
				pmeshnode[l].bdr = 3;
			}
			if (pmeshnode[m].bdr == 1 && pmeshnode[m].x < 1e-8) {
				pmeshnode[m].bdr = 3;
			}
			if (pmeshnode[n].bdr == 1 && pmeshnode[n].x < 1e-8) {
				pmeshnode[n].bdr = 3;
			}
		}
	}
	//定义一些变量
	umat locs(2, 16 * num_ele);    locs.zeros();
	mat vals(1, 16 * num_ele);  vals.zeros();
	vec bbJz = zeros<vec>(num_pts);
	uvec node_reorder = zeros<uvec>(num_pts);
	uvec node_pos = zeros<uvec>(num_pts);
	//vec bn = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	double ce[4][4] = { 0 };
	//处理边界点
	int node_bdr = 0;
	for (int i = 0; i < num_pts; i++) {
		if (pmeshnode[i].bdr == 3) {
			node_bdr++;
			node_reorder(num_pts - node_bdr) = i;
			node_pos(i) = num_pts - node_bdr;
			pmeshnode[i].A = 0;
			A(i) = 0;
		} else {
			node_reorder(i - node_bdr) = i;
			node_pos(i) = i - node_bdr;
		}
	}
	double* unknown_b = (double*)calloc(num_pts - node_bdr, sizeof(double));
	int iter = 0;//迭代步数
	int pos = 0;

	//主循环
	while (1){
		//生成全局矩阵
		for (int i = 0; i < num_ele; i++){
			//计算电流密度矩阵
			//计算电流密度//要注意domain会不会越界
			double jr = materialList[pmeshele4[i].domain - 1].Jr; 			
			double hc = materialList[pmeshele4[i].domain - 1].H_c;
			for (int row = 0; row < 4; row++) {
				bbJz(pmeshele4[i].n[row]) += jr*getJi(row, i);
				double ctmp = 0;
				for (int col = 0; col < 4; col++) {
					//计算系数
					ce[row][col] = getLocal4Matrix(row, col, i) + getDij(row, col, i);					
					ctmp += getDij(row, col, i)*pmeshnode[pmeshele4[i].n[col]].A;
					//判断节点是否在未知节点内
					//得到排序之后的编号
					int n_row = node_pos(pmeshele4[i].n[row]);
					int n_col = node_pos(pmeshele4[i].n[col]);
					if (n_row < num_pts - node_bdr && n_col < num_pts - node_bdr) {
						locs(0, pos) = n_row;
						locs(1, pos) = n_col;
						vals(0, pos) = ce[row][col];
						pos++;
					}				
				}
				bbJz(pmeshele4[i].n[row]) += ctmp;
				// 计算永磁部分
				int kk = row + 1; if (kk == 4) kk = 0;
				bbJz(pmeshele4[i].n[row]) += - hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
				bbJz(pmeshele4[i].n[kk]) += - hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
			}
		}
		//生成稀疏矩阵
		if (iter == 0) {
			locs.reshape(2, pos);
			vals.reshape(1, pos);
		}

		//使用构造函数来生成稀疏矩阵
		sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = bbJz(node_reorder(i));
		}
		//---------------------superLU_MT---------------------------------------
		CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
		if (superlumt.solve() == 1) {
			qDebug() << "Error: superlumt.slove";
			qDebug() << "info: " << superlumt.info;
			break;
		} else {
			double *sol = NULL;
			A_old = A;
			sol = superlumt.getResult();

			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A = sol[i];//the A is r*A_real
				A(node_reorder(i)) = sol[i];
			}
		}
		
		//double normA = norm(A, 2); //qDebug() << normA;
		double error = norm((A_old - A), 2);
		iter++;
		qDebug() << "iter: " << iter;
		qDebug() << "error: " << error;
		if (error < Precision && iter > 5) {
			//转换A
			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;//the A is r*A_real
				A(node_reorder(i)) /= pmeshnode[node_reorder(i)].x;
			}
			A.save("D:\\mypaper\\zhcore\\插图\\NRpmA.txt", arma::arma_ascii, false);
			bbJz.save("D:\\mypaper\\zhcore\\插图\\pmbn.txt", arma::arma_ascii, false);
			qDebug() << "solve over";
			break;
		}
		//重新初始化
		pos = 0;
		bbJz.zeros();
	}
	return true;
}
//采用传输线法求解轴对称静磁场，四边形分网
bool CFastFEMcore::StaticAxisQ4TLM(){
	//计算边界信息	
	for (int i = 0; i < num_ele; i++){
		int k, l, m, n;
		k = pmeshele4[i].n[0];
		l = pmeshele4[i].n[1];
		m = pmeshele4[i].n[2];
		n = pmeshele4[i].n[3];

		pmeshele4[i].rc = pmeshnode[k].x + pmeshnode[l].x + pmeshnode[m].x + pmeshnode[n].x;
		pmeshele4[i].rc /= 4;

		if (materialList[pmeshele4[i].domain - 1].BHpoints == 0) {
			pmeshele4[i].miu = materialList[pmeshele4[i].domain - 1].miu;
			pmeshele4[i].miut = materialList[pmeshele4[i].domain - 1].miu;//must be 1
			pmeshele4[i].LinearFlag = true;
		} else {
			pmeshele4[i].miu = 1 * miu0;
			pmeshele4[i].miut = 100 * miu0;
			pmeshele4[i].LinearFlag = false;
		}
		//查找物理边界
		if (pmeshele4[i].domain == 1) {
			pmeshnode[k].bdr = 3;

			pmeshnode[l].bdr = 3;

			pmeshnode[m].bdr = 3;

			pmeshnode[n].bdr = 3;
		} else {
			//对称轴上的部分边界点
			if (pmeshnode[k].bdr == 1 && pmeshnode[k].x < 1e-8) {
				pmeshnode[k].bdr = 3;
			}
			if (pmeshnode[l].bdr == 1 && pmeshnode[l].x < 1e-8) {
				pmeshnode[l].bdr = 3;
			}
			if (pmeshnode[m].bdr == 1 && pmeshnode[m].x < 1e-8) {
				pmeshnode[m].bdr = 3;
			}
			if (pmeshnode[n].bdr == 1 && pmeshnode[n].x < 1e-8) {
				pmeshnode[n].bdr = 3;
			}
		}
	}
	//寻找非线性单元
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele4[i].LinearFlag) {
			D34.push_back(i);
		}
	}
	umat locs(2, 16 * num_ele);
	locs.zeros();
	mat vals(1, 16 * num_ele);
	double ce[4][4] = { 0 };
	Resist4Matrix *rm = (Resist4Matrix*)malloc(D34.size() * sizeof(Resist4Matrix));
	vec bbJz = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	vec INL = zeros<vec>(num_pts);
	//根据边界条件对节点重新进行编号
	uvec node_reorder = zeros<uvec>(num_pts);
	uvec node_pos = zeros<uvec>(num_pts);
	int node_bdr = 0;
	for (int i = 0; i < num_pts; i++) {
		if (pmeshnode[i].bdr == 3) {
			node_bdr++;
			node_reorder(num_pts - node_bdr) = i;
			node_pos(i) = num_pts - node_bdr;
			pmeshnode[i].A = 0;
			A(i) = 0;
		} else {
			node_reorder(i - node_bdr) = i;
			node_pos(i) = i - node_bdr;
		}
	}
	double* unknown_b = (double*)calloc(num_pts - node_bdr, sizeof(double));
	int pos = 0;

	//计算导纳矩阵
	for (int i = 0; i < num_ele; i++) {
		if (pmeshele[i].LinearFlag) {//线性区，不用计算
			ce[0][1] = rm[i].Y12;
			ce[0][2] = rm[i].Y13;
			ce[1][2] = rm[i].Y23;
		} else {
			if (rm[i].Y12 < 0){//电阻
				ce[0][1] = rm[i].Y12;
			} else{
				ce[0][0] += rm[i].Y12;//在对角项上减去受控源
				ce[1][1] += rm[i].Y12;
				ce[0][1] = 0;//受控源在右侧，所以为0
			}
			if (rm[i].Y13 < 0){
				ce[0][2] = rm[i].Y13;
			} else{
				ce[0][0] += rm[i].Y13;
				ce[2][2] += rm[i].Y13;
				ce[0][2] = 0;
			}
			if (rm[i].Y23 < 0){
				ce[1][2] = rm[i].Y23;
			} else{
				ce[1][1] += rm[i].Y23;
				ce[2][2] += rm[i].Y23;
				ce[1][2] = 0;
			}
		}

		//将单元矩阵进行存储
		double jr = materialList[pmeshele4[i].domain - 1].Jr;
		double hc = materialList[pmeshele4[i].domain - 1].H_c;
		for (int row = 0; row < 4; row++) {			
			for (int col = 0; col < 4; col++) {	
				//线性部分，无需加传输线
				if (pmeshele4[i].LinearFlag) {
					ce[row][col] = getLocal4Matrix(row, col, i);
				} else {
					
				}
				//判断节点是否在未知节点内
				//得到排序之后的编号
				int n_row = node_pos(pmeshele[i].n[row]);
				int n_col = node_pos(pmeshele[i].n[col]);
				if (n_row < num_pts - node_bdr && n_col < num_pts - node_bdr) {
					locs(0, pos) = n_row;
					locs(1, pos) = n_col;
					vals(0, pos) = ce[row][col];
					pos++;
				}
			}
			//计算电流密度
			bbJz(pmeshele4[i].n[row]) += jr*getJi(row, i);
			// 计算永磁部分
			int kk = row + 1; if (kk == 4) kk = 0;
			bbJz(pmeshele4[i].n[row]) += -hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
			bbJz(pmeshele4[i].n[kk]) += -hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
		}		
	}//end for
	locs.reshape(2, pos);//重新调整大小
	vals.reshape(1, pos);
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

	INL += bbJz;
	for (int i = 0; i < num_pts - node_bdr; i++) {
		unknown_b[i] = INL(node_reorder(i));
	}
	//第一次求解
	//---------------------superLU_MT---------------------------------------
	//CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
	SuperMatrix   sluA; SuperMatrix sluB, sluX;
	//NCformat *Astore;
	double   *a;
	int_t      *asub, *xa;
	int_t      *perm_r; /* row permutations from partial pivoting */
	int_t      *perm_c; /* column permutation vector */
	SuperMatrix   L;       /* factor L */
	SCPformat *Lstore;
	SuperMatrix   U;       /* factor U */
	NCPformat *Ustore;
	int_t      nrhs, info, m, n, nnz;
	int_t      nprocs; /* maximum number of processors to use. */
	int_t      panel_size, relax, maxsup;
	int_t      permc_spec;
	trans_t  trans;
	//double   *rhs;
	superlu_memusage_t   superlu_memusage;
	DNformat	   *Bstore;
	double      *rhsb, *rhsx;
	Gstat_t  Gstat1;


	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	maxsup = sp_ienv(3);

	nprocs = 1;
	nrhs = 1;
	trans = NOTRANS;
	/* create matrix A in Harwell-Boeing format.*/
	m = num_pts - node_bdr; n = num_pts - node_bdr; nnz = X.n_nonzero;
	a = const_cast<double *>(X.values);

	StatAlloc(n, nprocs, panel_size, relax, &Gstat1);
	StatInit(n, nprocs, &Gstat1);

	asub = (int*)const_cast<unsigned int*>(X.row_indices);
	xa = (int*)const_cast<unsigned int*>(X.col_ptrs);
	dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

	//------create B and X-------------------
	if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&sluX, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);

	rhsb = unknown_b;
	dCreate_Dense_Matrix(&sluB, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
	Bstore = (DNformat*)sluB.Store;

	if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");

	/*
	* Get column permutation vector perm_c[], according to permc_spec:
	*   permc_spec = 0: natural ordering
	*   permc_spec = 1: minimum degree ordering on structure of A'*A
	*   permc_spec = 2: minimum degree ordering on structure of A'+A
	*   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	*/
	permc_spec = 1;
	get_perm_c(permc_spec, &sluA, perm_c);

	pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);

	if (info != 0) {
		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = (double*)((DNformat*)sluB.Store)->nzval;
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
			A(node_reorder(i)) = sol[i];
		}
	}
	//迭代

	//
	return true;
}
