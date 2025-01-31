#include <QFile>
#include <QXmlStreamWriter>
#include <QXmlStreamReader>
//#include <QDebug>
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
//		qDebug() << "Error: openning file!";
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
//		qDebug() << "Error: reading num_pts!";
		return 1;
	}
	int pts_ind;//the beginning of the points index
	//读取节点索引，默认从0开始
	if (fscanf(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
//		qDebug() << "Error: reading pts_ind!";
		return 1;
	}
	fgets(ch, 256, fp);

	for (int i = pts_ind; i < num_pts; i++) {
		//读取x,y坐标
		if (fscanf(fp, "%lf %lf \n", &(pmeshnode[i].x), &(pmeshnode[i].y)) != 2) {
//			qDebug() << "Error: reading mesh point!";
			return 1;
		}
	}
	//---------------vertex-------------------------------
	for (int i = 0; i < 7; i++)
		fgets(ch, 256, fp);
	int num_vtx_ns, num_vtx_ele;
	//
	if (fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns) != 1) {
//		qDebug() << "Error: reading num_vtx_ns!";
		return 1;
	}

	if (fscanf(fp, "%d # number of elements\n", &num_vtx_ele) != 1) {
//		qDebug() << "Error: reading num_vtx_ele!";
		return 1;
	}
	fgets(ch, 256, fp);

	int *vtx;
	vtx = (int*)calloc(num_vtx_ele, sizeof(int));
	for (int i = 0; i < num_vtx_ele; i++) {
		//好象是每一个域的顶点编号
		if (fscanf(fp, "%d \n", vtx + i) != 1) {
//			qDebug() << "Error: reading vertex condition!";
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
//			qDebug() << "Error: reading vertex condition!";
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
//		qDebug() << "Error: reading num_bdr_ns!";
		return 1;
	}
	//读取线段边界数目
	if (fscanf(fp, "%d # number of elements\n", &num_bdr_ele) != 1) {
//		qDebug() << "Error: reading num_bdr_ele!";
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
//			qDebug() << "Error: reading boundary condition!";
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
//			qDebug() << "Error: reading boundary condition!";
			return 1;
		}
	}
	if (entity != NULL) free(entity); entity = NULL;
	//----------------elements------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int ns_per_ele;// num_ele;//number of nodes per element;number of elements
	if (fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele) != 1) {
//		qDebug() << "Error: reading ns_per_ele!";
		return 1;
	}
	//读取分网单元数目
	if (fscanf(fp, "%d # number of elements\n", &num_ele) == 1) {
		pmeshele = (CElement*)calloc(num_ele, sizeof(CElement));
	} else {
//		qDebug() << "Error: reading num_ele!";
		return 1;
	}
	fgets(ch, 256, fp);
	//读取分网三角单元的三个节点索引
	for (int i = 0; i < num_ele; i++) {
		if (fscanf(fp, "%d %d %d \n", &pmeshele[i].n[0], &pmeshele[i].n[1], &pmeshele[i].n[2]) != 3) {
//			qDebug() << "Error: reading elements points!";
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
//			qDebug() << "Error: reading domain points!";
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
    t1 = SuperLU_timer_() - t1;
//    qDebug() << t1;
	int nnzU = Ustore->nnz;
	//qDebug()<<"U:"<<nnzU;
	//asubU = Ustore->rowind;
	t1 = SuperLU_timer_();
	QVector<int> xnU(nnzU), ynU(nnzU);
	QVector<double> valueU(nnzU);
	//    int *xnU = (int*)malloc(nnzU*sizeof(int));
	//    int *ynU = (int*)malloc(nnzU*sizeof(int));
	//    double *valueU = (double*)malloc(nnzU*sizeof(double));
    t1 = SuperLU_timer_() - t1;
//    qDebug() << t1;
	int fsupc, istart, nsupr, nsupc, nrow;//
	int count = 0; int countU = 0;
	double value;
	t1 = SuperLU_timer_();
	QVector<int> level(m), sortlevel(m), sortlevelU(m);
    t1 = SuperLU_timer_() - t1;
//    qDebug() << t1;
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


//	qDebug() << "countL :" << count;
//	qDebug() << "countU :" << countU;
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
//	qDebug() << "maxLevel:" << maxLevel;
//	qDebug() << "maxLevelU:" << maxLevelU;

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
//	qDebug() << "lcount:" << lcount;
//	qDebug() << "ucount:" << ucount;

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

//	qDebug() << m;
	t1 = SuperLU_timer_();

	DNformat *Bstore = (DNformat*)B->Store;
	double *Bmat = (double*)Bstore->nzval;
	//L solve
	int cc = 0;
	omp_set_num_threads(8);
//	qDebug() << omp_get_num_procs();
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
    t1 = SuperLU_timer_() - t1;
//    qDebug() << t1 << "\t" << cc;
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
//		qDebug() << "Error: superlumt.slove";
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
//			qDebug() << "Error: superlumt.slove";
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
//		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

//	qDebug() << "TLM steps:" << count;
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
//2018-03-04
//by Poofee
//该函数尝试对三角单元的地进行处理
bool CFastFEMcore::StaticAxisT3VTM2() {
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
	//传输线导纳矩阵，只保存上三角矩阵
	double * Ytl = (double*)malloc(D34.size() * 6 * sizeof(double));
	int nlin = 0;
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

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i];//猜测值
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i];

		//生成单元矩阵，线性与非线性
		//因为线性与非线性的差不多，所以不再分开讨论了
		ce[0][1] = rm[i].Y12 / pmeshele[i].miu;
		ce[0][2] = rm[i].Y13 / pmeshele[i].miu;
		ce[1][2] = rm[i].Y23 / pmeshele[i].miu;

		ce[0][0] = rm[i].Y11 / pmeshele[i].miu;
		ce[1][1] = rm[i].Y22 / pmeshele[i].miu;
		ce[2][2] = rm[i].Y33 / pmeshele[i].miu;

		//非线性区，计算传输线导纳矩阵
		if (!pmeshele[i].LinearFlag) {
			Ytl[6 * nlin + 0] = ce[0][0];
			Ytl[6 * nlin + 1] = ce[0][1];
			Ytl[6 * nlin + 2] = ce[0][2];
			Ytl[6 * nlin + 3] = ce[1][1];
			Ytl[6 * nlin + 4] = ce[1][2];
			Ytl[6 * nlin + 5] = ce[2][2];
			nlin++;
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
//		qDebug() << "Error: superlumt.slove";
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
	QVector<double> x(num_ele), y(num_ele),y1(num_ele);
	for (int i = 0; i < num_ele; ++i) {
		x[i] = i;
	}
	QCustomPlot * customplot;
	customplot = thePlot->getQcustomPlot();
	QCPGraph *graph1 = customplot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 3));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	QCPGraph *graph2 = customplot->addGraph();
	graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 1.5), QBrush(Qt::red), 3));
	graph2->setPen(QPen(QColor(120, 120, 120), 2));
	graph2->setLineStyle(QCPGraph::lsNone);
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
	//保存非线性单元节点电压，不是反射电压了
	Voltage *Vs = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of system
	Voltage *Ve = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of element
	//保存电流
	Voltage *Is = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into system
	Voltage *Ie = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into element
	time[tt++] = SuperLU_timer_();
	for (count = 0; count < steps; count++) {
		//PART A：计算非线性单元部分
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			//PART B：计算流入非线性单元侧的电流
			double Vlast[3];//保存上一步的值
			Vlast[0] = Ie[j].V12;
			Vlast[1] = Ie[j].V23;
			Vlast[2] = Ie[j].V13;
			

			Ie[j].V12 = 2 * (Ytl[6 * j + 0] * Vs[j].V12 + Ytl[6 * j + 1] * Vs[j].V23 + Ytl[6 * j + 2] * Vs[j].V13) - Is[j].V12;
			Ie[j].V23 = 2 * (Ytl[6 * j + 1] * Vs[j].V12 + Ytl[6 * j + 3] * Vs[j].V23 + Ytl[6 * j + 4] * Vs[j].V13) - Is[j].V23;
			Ie[j].V13 = 2 * (Ytl[6 * j + 2] * Vs[j].V12 + Ytl[6 * j + 4] * Vs[j].V23 + Ytl[6 * j + 5] * Vs[j].V13) - Is[j].V13;
			//PART C：计算流入线性系统侧的电流
			Is[j].V12 = 2 * (Ytl[6 * j + 0] * Ve[j].V12 + Ytl[6 * j + 1] * Ve[j].V23 + Ytl[6 * j + 2] * Ve[j].V13) - Vlast[0];
			Is[j].V23 = 2 * (Ytl[6 * j + 1] * Ve[j].V12 + Ytl[6 * j + 3] * Ve[j].V23 + Ytl[6 * j + 4] * Ve[j].V13) - Vlast[1];
			Is[j].V13 = 2 * (Ytl[6 * j + 2] * Ve[j].V12 + Ytl[6 * j + 4] * Ve[j].V23 + Ytl[6 * j + 5] * Ve[j].V13) - Vlast[2];

			INL(k) += Is[j].V12;
			INL(m) += Is[j].V23;
			INL(n) += Is[j].V13;
			//PART D：牛顿迭代求解小电路，计算电压
			mat C(2, 2);//单元系数矩阵，为了方便计算
			C(0, 0) = rm[i].Y11; C(0, 1) = rm[i].Y12;
			C(1, 0) = rm[i].Y12; C(1, 1) = rm[i].Y22; 
			mat AJ(2, 2);
			colvec b(2);
			colvec x2(3); x2.zeros();
			colvec x3(2); x3.zeros();
			double err1 = 0;

			for (int iter = 0; iter < 20; iter++){
				//1.初始化电流
				b(0) = Ie[j].V12;
				b(1) = Ie[j].V23;								
				//2.求解单元mu值
				double bx = 0;
				double by = 0;
				for (int j = 0; j < 3; j++) {
					bx += pmeshele[i].Q[j] * x2(j);
					by += pmeshele[i].P[j] * x2(j);
				}
				pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
				pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
				y[i] = Ve[j].V23 - Vs[j].V23;
				//y1[i] = Vs[j].V13;
				//3.计算雅可比矩阵
				double tmp = materialList[pmeshele[i].domain - 1].getdvdB(pmeshele[i].B);
				if (pmeshele[i].B > 1e-9){
					tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
					tmp /= ydot[i];
				}
				AJ.zeros();
				double ca[2];//理论上来说，这里的A应当是三个节点的A
				for (int row = 0; row < 2; row++){
					ca[row] = C(row, 0)*x3(0) + C(row, 1)*x3(1);
				}
				for (int row = 0; row < 2; row++){
					for (int col = 0; col < 2; col++){
						//注意C已经被除了一次ydot了
						AJ(row, col) = ca[row];
						AJ(row, col) *= ca[col];
						AJ(row, col) *= tmp;
						b(row) += AJ(row, col)*x3(col);
						AJ(row, col) += C(row, col) / m_e->miut;
					}
				}
				//4.加上传输线导纳
				AJ(0, 0) += Ytl[6 * j + 0];
				AJ(0, 1) += Ytl[6 * j + 1];
				AJ(1, 0) += Ytl[6 * j + 1];
				AJ(1, 1) += Ytl[6 * j + 3];
				//5.求解电压V
				bool status = arma::solve(x3, AJ, b);
				if (!status){
//					qDebug() << "error: solve !";
					return false;
				}
				x2(0) = x3(0) + Vs[j].V13;
				x2(1) = x3(1) + Vs[j].V13;
				x2(2) = Vs[j].V13;
			}

			//qDebug() << x2(0) << x2(1) << x2(2);
			Ve[j].V12 = x2(0);
			Ve[j].V23 = x2(1);
			Ve[j].V13 = x2(2);

			//计算电压
			Vs[j].V12 = pmeshnode[k].A - 0;
			Vs[j].V23 = pmeshnode[m].A - 0;
			Vs[j].V13 = pmeshnode[n].A - 0;
		}
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		//调用superLU三角求解
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		//读取求解结果
		if (info != 0) {
//			qDebug() << "Error: superlumt.slove";
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
		//误差判断
		double error = norm((A_old - A), 2) / norm(A, 2);
//		qDebug() << "iter: " << count;
//		qDebug() << "error: " << error;

		//绘制计算结果
		graph1->setData(x, y);
		//graph2->setData(x, y1);
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
//		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

//	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vs != NULL) free(Vs);
	if (Ve != NULL) free(Ve);
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
//2018-03-05
//by Poofee
//该函数不采用rA的变换方式进行求解
bool CFastFEMcore::StaticAxisT3VTM3() {
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
	//传输线导纳矩阵，只保存上三角矩阵
	double * Ytl = (double*)malloc(D34.size() * 6 * sizeof(double));
	int nlin = 0;
	int pos = 0;
	//
	for (int i = 0; i < num_ele; i++) {		
		//计算单元导纳
		rm[i].Y12 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[0] * pmeshele[i].P[1] +
			pmeshele[i].Q[0] * pmeshele[i].Q[1]);
		rm[i].Y12 += (pmeshele[i].P[0] + pmeshele[i].P[1]) / 6;
		rm[i].Y12 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
		rm[i].Y12 *= 2 * PI;

		rm[i].Y23 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[1] * pmeshele[i].P[2] +
			pmeshele[i].Q[1] * pmeshele[i].Q[2]);
		rm[i].Y23 += (pmeshele[i].P[1] + pmeshele[i].P[2]) / 6;
		rm[i].Y23 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
		rm[i].Y23 *= 2 * PI;

		rm[i].Y13 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[2] * pmeshele[i].P[0] +
			pmeshele[i].Q[2] * pmeshele[i].Q[0]);
		rm[i].Y13 += (pmeshele[i].P[2] + pmeshele[i].P[0]) / 6;
		rm[i].Y13 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
		rm[i].Y13 *= 2 * PI;

		rm[i].Y11 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[0] * pmeshele[i].P[0] +
			pmeshele[i].Q[0] * pmeshele[i].Q[0]);
		rm[i].Y11 += (pmeshele[i].P[0] + pmeshele[i].P[0]) / 6;
		rm[i].Y11 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
		rm[i].Y11 *= 2 * PI;

		rm[i].Y22 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[1] * pmeshele[i].P[1] +
			pmeshele[i].Q[1] * pmeshele[i].Q[1]);
		rm[i].Y22 += (pmeshele[i].P[1] + pmeshele[i].P[1]) / 6;
		rm[i].Y22 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
		rm[i].Y22 *= 2 * PI;

		rm[i].Y33 = pmeshele[i].rc / 4 / pmeshele[i].AREA*(pmeshele[i].P[2] * pmeshele[i].P[2] +
			pmeshele[i].Q[2] * pmeshele[i].Q[2]);
		rm[i].Y33 += (pmeshele[i].P[2] + pmeshele[i].P[2]) / 6;
		rm[i].Y33 += pmeshele[i].AREA / 9 / pmeshele[i].rc;
		rm[i].Y33 *= 2 * PI;

		//生成单元矩阵，线性与非线性
		//因为线性与非线性的差不多，所以不再分开讨论了
		ce[0][1] = rm[i].Y12 / pmeshele[i].miut;
		ce[0][2] = rm[i].Y13 / pmeshele[i].miut;
		ce[1][2] = rm[i].Y23 / pmeshele[i].miut;

		ce[0][0] = rm[i].Y11 / pmeshele[i].miut;
		ce[1][1] = rm[i].Y22 / pmeshele[i].miut;
		ce[2][2] = rm[i].Y33 / pmeshele[i].miut;

		//非线性区，计算传输线导纳矩阵
		if (!pmeshele[i].LinearFlag) {
			Ytl[6 * nlin + 0] = ce[0][0];
			Ytl[6 * nlin + 1] = ce[0][1];
			Ytl[6 * nlin + 2] = ce[0][2];
			Ytl[6 * nlin + 3] = ce[1][1];
			Ytl[6 * nlin + 4] = ce[1][2];
			Ytl[6 * nlin + 5] = ce[2][2];
			nlin++;
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
		double jr = 2 * PI*pmeshele[i].rc*pmeshele[i].AREA*materialList[pmeshele[i].domain - 1].Jr / 3;
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
//		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = (double*)((DNformat*)sluB.Store)->nzval;
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is A_real
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
	//保存非线性单元节点电压，不是反射电压了
	Voltage *Vs = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of system
	Voltage *Ve = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of element
	//保存电流
	Voltage *Is = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into system
	Voltage *Ie = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into element
	time[tt++] = SuperLU_timer_();
	for (count = 0; count < steps; count++) {
		//PART A：计算非线性单元部分
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			//PART B：计算流入非线性单元侧的电流
			double Vlast[3];//保存上一步的值
			Vlast[0] = Ie[j].V12;
			Vlast[1] = Ie[j].V23;
			Vlast[2] = Ie[j].V13;
			//计算电压
			Vs[j].V12 = pmeshnode[k].A - 0;
			Vs[j].V23 = pmeshnode[m].A - 0;
			Vs[j].V13 = pmeshnode[n].A - 0;

			Ie[j].V12 = 2 * (Ytl[6 * j + 0] * Vs[j].V12 + Ytl[6 * j + 1] * Vs[j].V23 + Ytl[6 * j + 2] * Vs[j].V13) - Is[j].V12;
			Ie[j].V23 = 2 * (Ytl[6 * j + 1] * Vs[j].V12 + Ytl[6 * j + 3] * Vs[j].V23 + Ytl[6 * j + 4] * Vs[j].V13) - Is[j].V23;
			Ie[j].V13 = 2 * (Ytl[6 * j + 2] * Vs[j].V12 + Ytl[6 * j + 4] * Vs[j].V23 + Ytl[6 * j + 5] * Vs[j].V13) - Is[j].V13;
			//PART C：计算流入线性系统侧的电流
			Is[j].V12 = 2 * (Ytl[6 * j + 0] * Ve[j].V12 + Ytl[6 * j + 1] * Ve[j].V23 + Ytl[6 * j + 2] * Ve[j].V13) - Vlast[0];
			Is[j].V23 = 2 * (Ytl[6 * j + 1] * Ve[j].V12 + Ytl[6 * j + 3] * Ve[j].V23 + Ytl[6 * j + 4] * Ve[j].V13) - Vlast[1];
			Is[j].V13 = 2 * (Ytl[6 * j + 2] * Ve[j].V12 + Ytl[6 * j + 4] * Ve[j].V23 + Ytl[6 * j + 5] * Ve[j].V13) - Vlast[2];

			INL(k) += Is[j].V12;
			INL(m) += Is[j].V23;
			INL(n) += Is[j].V13;
			//PART D：牛顿迭代求解小电路，计算电压
			mat C(3, 3);//单元系数矩阵，为了方便计算
			C(0, 0) = rm[i].Y11; C(0, 1) = rm[i].Y12; C(0, 2) = rm[i].Y13;
			C(1, 0) = rm[i].Y12; C(1, 1) = rm[i].Y22; C(1, 2) = rm[i].Y23;
			C(2, 0) = rm[i].Y13; C(2, 1) = rm[i].Y23; C(2, 2) = rm[i].Y33;
			mat AJ(3, 3);
			colvec b(3);
			colvec x2(3); x2.zeros();
			double err1 = 0;

			for (int iter = 0; iter < 20; iter++) {
				//1.初始化电流
				b(0) = Ie[j].V12;
				b(1) = Ie[j].V23;
				b(2) = Ie[j].V13;
				//
				AJ.zeros();
				double ca[3];//理论上来说，这里的A应当是三个节点的A
				for (int row = 0; row < 3; row++) {
					ca[row] = C(row, 0)*x2(0) + C(row, 1)*x2(1) + C(row, 2)*x2(2);
				}
				//2.求解单元mu值
				double bx = 0;
				double by = 0;
				for (int j = 0; j < 3; j++) {
					bx += pmeshele[i].Q[j] * x2(j);
					by += pmeshele[i].P[j] * x2(j);
				}
				pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
				pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
				y[i] = pmeshele[i].miut;
				//3.计算雅可比矩阵
				double tmp = materialList[pmeshele[i].domain - 1].getdvdB(pmeshele[i].B);
				if (pmeshele[i].B > 1e-9) {
					tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
					tmp /= ydot[i];
				}
				for (int row = 0; row < 3; row++) {
					for (int col = 0; col < 3; col++) {
						//注意C已经被除了一次ydot了
						AJ(row, col) = ca[row];
						AJ(row, col) *= ca[col];
						AJ(row, col) *= tmp;
						b(row) += AJ(row, col)*x2(col);
						AJ(row, col) += C(row, col) / m_e->miut;
					}
				}
				//4.加上传输线导纳
				AJ(0, 0) += Ytl[6 * j + 0];
				AJ(0, 1) += Ytl[6 * j + 1];
				AJ(0, 2) += Ytl[6 * j + 2];
				AJ(1, 0) += Ytl[6 * j + 1];
				AJ(1, 1) += Ytl[6 * j + 3];
				AJ(1, 2) += Ytl[6 * j + 4];
				AJ(2, 0) += Ytl[6 * j + 2];
				AJ(2, 1) += Ytl[6 * j + 4];
				AJ(2, 2) += Ytl[6 * j + 5];
				//5.求解电压V
				bool status = arma::solve(x2, AJ, b);
				if (!status) {
//					qDebug() << "error: solve !";
					return false;
				}
			}

			//qDebug() << x2(0) << x2(1) << x2(2);
			Ve[j].V12 = x2(0);
			Ve[j].V23 = x2(1);
			Ve[j].V13 = x2(2);
		}
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		//调用superLU三角求解
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		//读取求解结果
		if (info != 0) {
//			qDebug() << "Error: superlumt.slove";
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
		//误差判断
		double error = norm((A_old - A), 2) / norm(A, 2);
//		qDebug() << "iter: " << count;
//		qDebug() << "error: " << error;

		//绘制计算结果
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
//		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

//	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vs != NULL) free(Vs);
	if (Ve != NULL) free(Ve);
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

//2018-03-03
//by Poofee
//采用VTM来进行静磁场的求解
//每个单元都被隔离
bool CFastFEMcore::StaticAxisT3VTM() {
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
	//传输线导纳矩阵，只保存上三角矩阵
	double * Ytl = (double*)malloc(D34.size()*6*sizeof(double));
	int nlin = 0;
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

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i];//猜测值
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i];

		//生成单元矩阵，线性与非线性
		//因为线性与非线性的差不多，所以不再分开讨论了	

		//非线性区，计算传输线导纳矩阵
		double delta = 1e-3;//接地导纳
		if (!pmeshele[i].LinearFlag) {
			ce[0][1] = rm[i].Y12 / pmeshele[i].miu;
			ce[0][2] = rm[i].Y13 / pmeshele[i].miu;
			ce[1][2] = rm[i].Y23 / pmeshele[i].miu;

			ce[0][0] = rm[i].Y11 / pmeshele[i].miu;
			ce[1][1] = rm[i].Y22 / pmeshele[i].miu;
			ce[2][2] = rm[i].Y33 / pmeshele[i].miu*(1 + delta);

			if (ce[0][1] > 0){
				Ytl[6 * nlin + 1] = -ce[0][1];
			} else{
				Ytl[6 * nlin + 1] = ce[0][1];
			}
			if (ce[0][2] > 0){
				Ytl[6 * nlin + 2] = -ce[0][2];
			} else{
				Ytl[6 * nlin + 2] = ce[0][2];
			}
			if (ce[1][2] > 0){
				Ytl[6 * nlin + 4] = -ce[1][2];
			} else{
				Ytl[6 * nlin + 4] = ce[1][2];
			}

			Ytl[6 * nlin + 0] = ce[0][0];
			
			
			Ytl[6 * nlin + 3] = ce[1][1];
			
			Ytl[6 * nlin + 5] = ce[2][2];
			nlin++;
		} else{
			ce[0][1] = rm[i].Y12 / pmeshele[i].miu;
			ce[0][2] = rm[i].Y13 / pmeshele[i].miu;
			ce[1][2] = rm[i].Y23 / pmeshele[i].miu;

			ce[0][0] = rm[i].Y11 / pmeshele[i].miu;
			ce[1][1] = rm[i].Y22 / pmeshele[i].miu;
			ce[2][2] = rm[i].Y33 / pmeshele[i].miu;
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
//		qDebug() << "Error: superlumt.slove";
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
	//保存非线性单元节点电压，不是反射电压了
	Voltage *Vs = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of system
	Voltage *Ve = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of element
	//保存电流
	Voltage *Is = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into system
	Voltage *Ie = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into element
	time[tt++] = SuperLU_timer_();
	for (count = 0; count < steps; count++) {
		//PART A：计算非线性单元部分
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			//PART B：计算流入非线性单元侧的电流
			double Vlast[3];//保存上一步的值
			Vlast[0] = Ie[j].V12;
			Vlast[1] = Ie[j].V23;
			Vlast[2] = Ie[j].V13;
			//计算电压
			Vs[j].V12 = pmeshnode[k].A - 0;
			Vs[j].V23 = pmeshnode[m].A - 0;
			Vs[j].V13 = pmeshnode[n].A - 0;

			Ie[j].V12 = 2 * (Ytl[6 * j + 0] * Vs[j].V12 + Ytl[6 * j + 1] * Vs[j].V23 + Ytl[6 * j + 2] * Vs[j].V13) - Is[j].V12;
			Ie[j].V23 = 2 * (Ytl[6 * j + 1] * Vs[j].V12 + Ytl[6 * j + 3] * Vs[j].V23 + Ytl[6 * j + 4] * Vs[j].V13) - Is[j].V23;
			Ie[j].V13 = 2 * (Ytl[6 * j + 2] * Vs[j].V12 + Ytl[6 * j + 4] * Vs[j].V23 + Ytl[6 * j + 5] * Vs[j].V13) - Is[j].V13;
			//PART C：计算流入线性系统侧的电流
			Is[j].V12 = 2 * (Ytl[6 * j + 0] * Ve[j].V12 + Ytl[6 * j + 1] * Ve[j].V23 + Ytl[6 * j + 2] * Ve[j].V13) - Vlast[0];
			Is[j].V23 = 2 * (Ytl[6 * j + 1] * Ve[j].V12 + Ytl[6 * j + 3] * Ve[j].V23 + Ytl[6 * j + 4] * Ve[j].V13) - Vlast[1];
			Is[j].V13 = 2 * (Ytl[6 * j + 2] * Ve[j].V12 + Ytl[6 * j + 4] * Ve[j].V23 + Ytl[6 * j + 5] * Ve[j].V13) - Vlast[2];

			INL(k) += Is[j].V12;
			INL(m) += Is[j].V23;
			INL(n) += Is[j].V13;
			//PART D：牛顿迭代求解小电路，计算电压
			mat C(3, 3);//单元系数矩阵，为了方便计算
			C(0, 0) = rm[i].Y11; C(0, 1) = rm[i].Y12; C(0, 2) = rm[i].Y13;
			C(1, 0) = rm[i].Y12; C(1, 1) = rm[i].Y22; C(1, 2) = rm[i].Y23;
			C(2, 0) = rm[i].Y13; C(2, 1) = rm[i].Y23; C(2, 2) = rm[i].Y33;
			mat AJ(3, 3);
			colvec b(3);
			colvec x2(3); x2.zeros();
			double err1 = 0;

			for (int iter = 0; iter < 20; iter++){
				//1.初始化电流
				b(0) = Ie[j].V12;
				b(1) = Ie[j].V23;
				b(2) = Ie[j].V13;
				//
				AJ.zeros();
				double ca[3];//理论上来说，这里的A应当是三个节点的A
				for (int row = 0; row < 3; row++){
					ca[row] = C(row, 0)*x2(0) + C(row, 1)*x2(1) + C(row, 2)*x2(2);
				}
				//2.求解单元mu值
				double bx = 0;
				double by = 0;
				for (int j = 0; j < 3; j++) {
					bx += pmeshele[i].Q[j] * x2(j);
					by += pmeshele[i].P[j] * x2(j);
				}
				pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
				pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
				y[i] = pmeshele[i].miut;
				//3.计算雅可比矩阵
				double tmp = materialList[pmeshele[i].domain - 1].getdvdB(pmeshele[i].B);
				if (pmeshele[i].B > 1e-9){
					tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
					tmp /= ydot[i];
				}
				for (int row = 0; row < 3; row++){
					for (int col = 0; col < 3; col++){
						//注意C已经被除了一次ydot了
						AJ(row, col) = ca[row];
						AJ(row, col) *= ca[col];
						AJ(row, col) *= tmp;
						b(row) += AJ(row, col)*x2(col);
						AJ(row, col) += C(row, col) / m_e->miut;
					}
				}
				//4.加上传输线导纳
				AJ(0, 0) += Ytl[6 * j + 0];
				AJ(0, 1) += Ytl[6 * j + 1];
				AJ(0, 2) += Ytl[6 * j + 2];
				AJ(1, 0) += Ytl[6 * j + 1];
				AJ(1, 1) += Ytl[6 * j + 3];
				AJ(1, 2) += Ytl[6 * j + 4];
				AJ(2, 0) += Ytl[6 * j + 2];
				AJ(2, 1) += Ytl[6 * j + 4];
				AJ(2, 2) += Ytl[6 * j + 5];
				//5.求解电压V
				bool status = arma::solve(x2, AJ, b);
				if (!status){
//					qDebug() << "error: solve !";
					return false;
				}
			}

			//qDebug() << x2(0) << x2(1) << x2(2);
			Ve[j].V12 = x2(0);
			Ve[j].V23 = x2(1);
			Ve[j].V13 = x2(2);
		}
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		//调用superLU三角求解
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		//读取求解结果
		if (info != 0) {
//			qDebug() << "Error: superlumt.slove";
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
		//误差判断
		double error = norm((A_old - A), 2) / norm(A, 2);
//		qDebug() << "iter: " << count;
//		qDebug() << "error: " << error;

		//绘制计算结果
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
//		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

//	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vs != NULL) free(Vs);
	if (Ve != NULL) free(Ve);
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

//尝试将三角形单元的电路看作为一个整体，三角单元内部暂时还是使用线性来处理
bool CFastFEMcore::StaticAxisT3TLMgroup() {
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
//	qDebug() << "D34 size:" << D34.size();
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
	double * y10 = (double*)malloc(D34.size()*sizeof(double));
	double * y20 = (double*)malloc(D34.size()*sizeof(double));
	double * y30 = (double*)malloc(D34.size()*sizeof(double));
	int nlin = 0;
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

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i];//猜测值
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i];

		//生成单元矩阵，线性与非线性
		// 因为线性与非线性的差不多，所以不再分开讨论了


		if (pmeshele[i].LinearFlag) {//线性区，不用计算
			ce[0][1] = rm[i].Y12 / pmeshele[i].miu;
			ce[0][2] = rm[i].Y13 / pmeshele[i].miu;
			ce[1][2] = rm[i].Y23 / pmeshele[i].miu;

			ce[0][0] = rm[i].Y11 / pmeshele[i].miu;
			ce[1][1] = rm[i].Y22 / pmeshele[i].miu;
			ce[2][2] = rm[i].Y33 / pmeshele[i].miu;
		} else {
			ce[0][1] = 0;
			ce[0][2] = 0;
			ce[1][2] = 0;
			int a = 1000;
			//y10[nlin] = a * rm[i].Y11 / pmeshele[i].miu ;
			//y20[nlin] = a * rm[i].Y22 / pmeshele[i].miu;
			//y30[nlin] = a * rm[i].Y33 / pmeshele[i].miu;
			
			y10[nlin] = a * rm[i].Y11 / pmeshele[i].miu + a * rm[i].Y22 / pmeshele[i].miu + a * rm[i].Y33 / pmeshele[i].miu;
			y20[nlin] = a * rm[i].Y11 / pmeshele[i].miu + a * rm[i].Y22 / pmeshele[i].miu + a * rm[i].Y33 / pmeshele[i].miu;
			y30[nlin] = a * rm[i].Y11 / pmeshele[i].miu + a * rm[i].Y22 / pmeshele[i].miu + a * rm[i].Y33 / pmeshele[i].miu;

			/*y10[nlin] = rm[i].Y11*pmeshnode[pmeshele[i].n[0]].A +
				rm[i].Y12*pmeshnode[pmeshele[i].n[1]].A +
				rm[i].Y13*pmeshnode[pmeshele[i].n[2]].A;
			y10[nlin] /= pmeshele[i].miu;
			y10[nlin] /= pmeshnode[pmeshele[i].n[0]].A;
			if (std::isinf(y10[nlin])){
				y10[nlin] = abs(rm[i].Y12) / pmeshele[i].miu;
			}
			y20[nlin] = rm[i].Y12*pmeshnode[pmeshele[i].n[0]].A +
				rm[i].Y22*pmeshnode[pmeshele[i].n[1]].A +
				rm[i].Y23*pmeshnode[pmeshele[i].n[2]].A;
			y20[nlin] /= pmeshele[i].miu;
			y20[nlin] /= pmeshnode[pmeshele[i].n[1]].A;
			if (std::isinf(y20[nlin])){
				y20[nlin] = abs(rm[i].Y23) / pmeshele[i].miu;
			}
			y30[nlin] = rm[i].Y13*pmeshnode[pmeshele[i].n[0]].A +
				rm[i].Y23*pmeshnode[pmeshele[i].n[1]].A +
				rm[i].Y33*pmeshnode[pmeshele[i].n[2]].A;
			y30[nlin] /= pmeshele[i].miu;
		    y30[nlin] /= pmeshnode[pmeshele[i].n[2]].A;
			if (std::isinf(y30[nlin])){
				y30[nlin] = abs(rm[i].Y13) / pmeshele[i].miu;
			}*/
			
			ce[0][0] = abs(y10[nlin]);
			ce[1][1] = abs(y20[nlin]);
			ce[2][2] = abs(y30[nlin]);
			nlin++;
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
//		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = (double*)((DNformat*)sluB.Store)->nzval;
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			pmeshnode[node_reorder(i)].A *= 1;// / pmeshnode[i].x;//the A is r*A_real
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
	int steps = 3000;
	int count;
	double alpha = 1;
	Voltage *Vr = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	Voltage *Vi = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	//根据粗糙结果反过来计算出每一个单元当中的入射电压
	for (int j = 0; j < D34.size(); j++){
		int i = D34[j];
		CElement *m_e = pmeshele + i;
		int k, m, n;
		k = m_e->n[0];
		m = m_e->n[1];
		n = m_e->n[2];
		
		//计算磁导率
		double bx = 0;
		double by = 0;
		bx = pmeshele[i].Q[0] * pmeshnode[k].A+
			pmeshele[i].Q[1] * pmeshnode[m].A+
			pmeshele[i].Q[2] * pmeshnode[n].A;
		by = pmeshele[i].P[0] * pmeshnode[k].A +
			pmeshele[i].P[1] * pmeshnode[m].A +
			pmeshele[i].P[2] * pmeshnode[n].A;
		//pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
		//pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
		mat C(3, 3);//单元系数矩阵，为了方便计算
		C(0, 0) = rm[i].Y11 / pmeshele[i].miut; C(0, 1) = rm[i].Y12 / pmeshele[i].miut; C(0, 2) = rm[i].Y13 / pmeshele[i].miut;
		C(1, 0) = rm[i].Y12 / pmeshele[i].miut; C(1, 1) = rm[i].Y22 / pmeshele[i].miut; C(1, 2) = rm[i].Y23 / pmeshele[i].miut;
		C(2, 0) = rm[i].Y13 / pmeshele[i].miut; C(2, 1) = rm[i].Y23 / pmeshele[i].miut; C(2, 2) = rm[i].Y33 / pmeshele[i].miut;

		Vr[j].V12 = (abs(y10[j])+C(0, 0))*pmeshnode[k].A + C(0, 1)*pmeshnode[m].A + C(0, 2)*pmeshnode[n].A;
		Vr[j].V12 /=  2 * abs(y10[j]);
		Vr[j].V23 = C(1, 0)*pmeshnode[k].A + (abs(y20[j])+C(1, 1))*pmeshnode[m].A + C(1, 2)*pmeshnode[n].A;
		Vr[j].V23 /=  2 * abs(y20[j]);
		Vr[j].V13 = C(2, 0)*pmeshnode[k].A + C(2, 1)*pmeshnode[m].A + (abs(y30[j])+C(2, 2))*pmeshnode[n].A;
		Vr[j].V13 /=  2 * abs(y30[j]);
	}
	time[tt++] = SuperLU_timer_();
	double a1, a2;
	for (count = 0; count < steps; count++) {
		//------update miu----------------
		//for (int i = 0; i < num_ele; i++) {
		//	double bx = 0;
		//	double by = 0;
		//	for (int j = 0; j < 3; j++) {
		//		bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
		//		by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
		//	}
		//	pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
		//	pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
		//	//pmeshele[i].miut = y[i] + (count>90 ? 0.02 : (0.9 - count*0.001))*(pmeshele[i].miut - y[i]);
		//	y[i] = pmeshele[i].miut;
		//}
		//#pragma omp parallel for
		a1 = SuperLU_timer_();
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			//计算反射电压
			if (count != 0){
				Vr[j].V12 = (pmeshnode[k].A - 0) - Vi[j].V12;
				Vr[j].V23 = (pmeshnode[m].A - 0) - Vi[j].V23;
				Vr[j].V13 = (pmeshnode[n].A - 0) - Vi[j].V13;
			}
			
			//求解小电路，计算入射电压
			//采用线性代入不收敛，遂尝试牛顿迭代
			mat C(3, 3);//单元系数矩阵，为了方便计算
			C(0, 0) = rm[i].Y11; C(0, 1) = rm[i].Y12; C(0, 2) = rm[i].Y13;
			C(1, 0) = rm[i].Y12; C(1, 1) = rm[i].Y22; C(1, 2) = rm[i].Y23;
			C(2, 0) = rm[i].Y13; C(2, 1) = rm[i].Y23; C(2, 2) = rm[i].Y33;
			mat AJ(3, 3);
			colvec b(3);
			colvec x2(3); x2.zeros();
			colvec x2old(3); x2old.zeros();
			double err1 = 0;
			
			for (int iter = 0; iter < 20; iter++){
				//初始化电流
				b(0) = 2. * Vr[j].V12*abs(y10[j]);
				b(1) = 2. * Vr[j].V23*abs(y20[j]);
				b(2) = 2. * Vr[j].V13*abs(y30[j]);
				//
				AJ.zeros();
				double ca[3];//理论上来说，这里的A应当是三个节点的A
				for (int row = 0; row < 3; row++){
					ca[row] = C(row, 0)*x2(0) + C(row, 1)*x2(1) + C(row, 2)*x2(2);
				}
				//求解单元mu值
				double bx = 0;
				double by = 0;
				for (int j = 0; j < 3; j++) {
					bx += pmeshele[i].Q[j] * x2(j);
					by += pmeshele[i].P[j] * x2(j);
				}
				pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
				pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
				y[i] = pmeshele[i].miut;
				//计算雅可比矩阵
				double tmp = materialList[pmeshele[i].domain - 1].getdvdB(pmeshele[i].B);
				if (pmeshele[i].B > 1e-9){
					tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
					tmp /= ydot[i];
				}
				for (int row = 0; row < 3; row++){
					for (int col = 0; col < 3; col++){
						//注意C已经被除了一次ydot了
						AJ(row, col) = ca[row];
						AJ(row, col) *= ca[col];
						AJ(row, col) *= tmp;
						b(row) += AJ(row, col)*x2(col);
						AJ(row, col) += C(row, col) / m_e->miut;
					}
				}
				//加上对地导纳
				AJ(0, 0) += abs(y10[j]);
				AJ(1, 1) += abs(y20[j]);
				AJ(2, 2) += abs(y30[j]);
				//求解电压V
				x2old = x2;
				bool status = arma::solve(x2, AJ, b);
				if (!status){
//					qDebug() << "error: solve !";
					return false;
				}
				//判断误差
				if (iter > 2){
					double err = norm(x2old - x2) / norm(x2old);
					//qDebug() << err;
					if (err < 1e-6){
						break;
					}
				}
			}
			
			//qDebug() << x2(0) << x2(1) << x2(2);
			Vi[j].V12 = x2(0) - Vr[j].V12;
			Vi[j].V23 = x2(1) - Vr[j].V23;
			Vi[j].V13 = x2(2) - Vr[j].V13;
			INL(k) += 2. *Vi[j].V12*abs(y10[j]);

			INL(m) += 2. * Vi[j].V23*abs(y20[j]);

			INL(n) += 2. * Vi[j].V13*abs(y30[j]);
		}
		INL += bbJz;
		a2 = SuperLU_timer_();
//		qDebug() <<"t1: "<< a2 - a1;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		/*double a = SuperLU_timer_();
		if (count == 0){
			qDebug() << a - time[tt - 1];
		}*/
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		if (info != 0) {
//			qDebug() << "Error: superlumt.slove";
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
//		qDebug() << count;
//		qDebug() << error;

		//graph1->setData(x, y);
		//customplot->rescaleAxes(true);
		//customplot->replot();


		if (error < 1e-5) {
			break;
		}
		INL.zeros();
		a1 = SuperLU_timer_();
//		qDebug() <<"t2: "<< a1 - a2;
		/*double b = SuperLU_timer_();
		if (count == 0){
			qDebug() << b - a;
		}*/
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
		//A(node_reorder(i)) /= pmeshnode[node_reorder(i)].x;
	}
	//output the time
	for (int i = 1; i < tt; i++) {
//		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

//	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vi != NULL) free(Vi);
	if (Vr != NULL) free(Vr);
	if (unknown_b != NULL) free(unknown_b);

	A.save("TLM_T3_A.txt", arma::arma_ascii, false);
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
			ce[2][2] = -ce[0][2] - ce[1][2];
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
//		qDebug() << "Error: superlumt.slove";
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
//				qDebug() << A(pmeshele[i].n[0]);
//				qDebug() << A(pmeshele[i].n[1]);
//				qDebug() << A(pmeshele[i].n[2]);
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
			double c12, c23, c13, c11, c22, c33;
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
//				qDebug() << Vi[j].V12;
			}
			if (std::isnan(Vi[j].V23)) {
//				qDebug() << Vi[j].V23;
			}
			if (std::isnan(Vi[j].V13)) {
//				qDebug() << Vi[j].V13;
			}
			INL(k) += 2.*Vi[j].V12*fabs(rm[i].Y12);
			INL(m) -= 2.*Vi[j].V12*fabs(rm[i].Y12);

			INL(m) += 2.*Vi[j].V23*fabs(rm[i].Y23);
			INL(n) -= 2.*Vi[j].V23*fabs(rm[i].Y23);

			INL(n) += 2.*Vi[j].V13*fabs(rm[i].Y13);
			INL(k) -= 2.*Vi[j].V13*fabs(rm[i].Y13);
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
//			qDebug() << "Error: superlumt.slove";
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
//		qDebug() << i << "\t" << time[i] - time[i - 1];
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
//	qDebug() << "yForce: " << yForce;
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
//		qDebug() << "read inbox file error...";
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
//					qDebug() << "domainNum = " << numDomain;
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
//			qDebug() << "domainName = " << reader.readElementText();
		} else if (reader.name() == "miu") {
			materialList[i].miu = reader.readElementText().toDouble() * 4 * PI*1e-7;
//			qDebug() << "miu = " << materialList[i].miu;
		} else if (reader.name() == "BH") {
			readBHElement(reader, i);
		} else if (reader.name() == "Jr") {
			materialList[i].Jr = reader.readElementText().toDouble();
//			qDebug() << "Jr = " << materialList[i].Jr;
		} else if (reader.name() == "H_c") {
			materialList[i].H_c = reader.readElementText().toDouble();
//			qDebug() << "H_c = " << materialList[i].H_c;
		}
	}
}


void CFastFEMcore::readBHElement(QXmlStreamReader &reader, int i) {
	reader.readNextStartElement();
	//qDebug()<<reader.name();
	if (reader.name() == "BHpoints") {
		materialList[i].BHpoints = reader.readElementText().toInt();
//		qDebug() << "BHpoints = " << materialList[i].BHpoints;
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
//非线性迭代采用NR方法，每一步的线性方程组求解采用TLM迭代
bool CFastFEMcore::StaticAxisT3NRTLM(){
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
	vec INL = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	vec A_tlm = A;//用来保存tlm迭代中的值
	double * ydot = (double*)malloc(num_ele*sizeof(double));
	ResistMarix *rm = (ResistMarix*)malloc(num_ele * sizeof(ResistMarix));

	//绘图的一些变量
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
	//寻找非线性单元编号
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele[i].LinearFlag) {
			D34.push_back(i);
		}
	}
//	qDebug() << "num of points: " << num_pts;
//	qDebug() << "num of elements: " << num_ele;
//	qDebug() << "D34 size: " << D34.size();
	//装配传输线导纳稀疏矩阵，原则在于不添加传输线区域的，矩阵不变，
	//添加传输线的，导纳要大于零
	int pos = 0;
	for (int i = 0; i < num_ele; i++) {
		//这部分只需计算一次即可		
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

		// /=容易出错
		rm[i].Y11 /= 4. * pmeshele[i].AREA * pmeshele[i].miut * ydot[i]/1.5;
		rm[i].Y12 /= 4. * pmeshele[i].AREA * pmeshele[i].miut * ydot[i]/1.5;
		rm[i].Y13 /= 4. * pmeshele[i].AREA * pmeshele[i].miut * ydot[i]/1.5;
		rm[i].Y22 /= 4. * pmeshele[i].AREA * pmeshele[i].miut * ydot[i]/1.5;
		rm[i].Y23 /= 4. * pmeshele[i].AREA * pmeshele[i].miut * ydot[i]/1.5;
		rm[i].Y33 /= 4. * pmeshele[i].AREA * pmeshele[i].miut * ydot[i]/1.5;

		//计算电流密度//要注意domain会不会越界
		double jr = pmeshele[i].AREA*materialList[pmeshele[i].domain - 1].Jr / 3;
		for (int j = 0; j < 3; j++) {
			bbJz(pmeshele[i].n[j]) += jr;
			// 计算永磁部分
			bbJz(pmeshele[i].n[j]) -= materialList[pmeshele[i].domain - 1].H_c / 2.*pmeshele[i].Q[j];
		}
		//主要求解结果不要漏掉miu0
		if (pmeshele[i].LinearFlag){
			ce[0][0] = rm[i].Y11/1.5 ;
			ce[1][1] = rm[i].Y22/1.5 ;
			ce[2][2] = rm[i].Y33/1.5 ;

			ce[0][1] = rm[i].Y12/1.5 ;
			ce[0][2] = rm[i].Y13/1.5 ;
			ce[1][2] = rm[i].Y23/1.5 ;
		} else{
			ce[0][1] = -abs(rm[i].Y12);
			ce[0][2] = -abs(rm[i].Y13);
			ce[1][2] = -abs(rm[i].Y23);

			ce[0][0] = -ce[0][1] - ce[0][2];
			ce[1][1] = -ce[0][1] - ce[1][2];
			ce[2][2] = -ce[1][2] - ce[0][2];
		}
		

		//计算牛顿迭代部分的单元矩阵项,如果是第一次迭代的话，A=0，
		//所以就不计算了，参见颜威利书P56
		//double v[3];

		/*v[0] = rm[i].Y11*A(pmeshele[i].n[0]) +
			rm[i].Y12*A(pmeshele[i].n[1]) +
			rm[i].Y13*A(pmeshele[i].n[2]);
		v[1] = rm[i].Y12*A(pmeshele[i].n[0]) +
			rm[i].Y22*A(pmeshele[i].n[1]) +
			rm[i].Y23*A(pmeshele[i].n[2]);
		v[2] = rm[i].Y13*A(pmeshele[i].n[0]) +
			rm[i].Y23*A(pmeshele[i].n[1]) +
			rm[i].Y33*A(pmeshele[i].n[2]);*/

		//if (iter != 0) {
		//	double tmp;
		//	if (pmeshele[i].LinearFlag) {
		//		tmp = 0;
		//	} else {
		//		tmp = materialList[pmeshele[i].domain - 1].getdvdB(pmeshele[i].B);
		//		if (pmeshele[i].B > 1e-9){
		//			tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
		//			tmp /= ydot[i] * ydot[i] * ydot[i];
		//		}

		//	}
		//	cn[0][0] = v[0] * v[0] * tmp;
		//	cn[1][1] = v[1] * v[1] * tmp;
		//	cn[2][2] = v[2] * v[2] * tmp;

		//	cn[0][1] = v[0] * v[1] * tmp;
		//	cn[0][2] = v[0] * v[2] * tmp;
		//	cn[1][2] = v[1] * v[2] * tmp;

		//	cn[1][0] = cn[0][1];
		//	cn[2][0] = cn[0][2];
		//	cn[2][1] = cn[1][2];
		//}
		//ce[0][0] += cn[0][0];
		//ce[1][1] += cn[1][1];
		//ce[2][2] += cn[2][2];

		//ce[0][1] += cn[0][1];
		//ce[0][2] += cn[0][2];
		//ce[1][2] += cn[1][2];

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
				//bn(pmeshele[i].n[row]) += cn[row][col] * A(pmeshele[i].n[col]);
			}
		}
	}//end of elememt iteration

	locs.reshape(2, pos);
	vals.reshape(1, pos);
	//使用构造函数来生成稀疏矩阵
	sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);
	double* unknown_b = (double*)calloc(num_pts - node_bdr, sizeof(double));
	bn += bbJz;
	//bn.save("bn.txt", arma::arma_ascii, false);
	for (int i = 0; i < num_pts - node_bdr; i++) {
		unknown_b[i] = bn(node_reorder(i));
	}
	//第一次求解，主要是完成LU分解
	CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
	if (superlumt.solve1() == 1) {
//		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		sol = superlumt.get1Result();
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
			A(node_reorder(i)) = pmeshnode[node_reorder(i)].A;
		}
	}
	//牛顿迭代法	
	int iter = 0;//迭代步数	
	double tlm_tol = 1e-10;//TLM迭代精度
	Voltage *Vr = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	Voltage *Vi = (Voltage*)calloc(D34.size(), sizeof(Voltage));
	while (1) {
		//PART A:迭代开始，更新非线性区域的导纳矩阵，以及电流矩阵
		//保存上一步的A
		A_old = A;
		A.zeros();
		int negY = 0;
		for (int iter_tlm = 0; iter_tlm < 10000; iter_tlm++){
			//非线性区域计算
			for (int j = 0; j < D34.size(); j++){
				int i = D34[j];
				CElement * m_e = pmeshele + i;
				int k, m, n;
				k = m_e->n[0];
				m = m_e->n[1];
				n = m_e->n[2];
				
				//计算非线性导纳
				double Y11, Y12, Y22, Y23, Y13, Y33;
				Y11 = m_e->Q[0] * m_e->Q[0] + m_e->P[0] * m_e->P[0];
				Y12 = m_e->Q[0] * m_e->Q[1] + m_e->P[0] * m_e->P[1];
				Y13 = m_e->Q[0] * m_e->Q[2] + m_e->P[0] * m_e->P[2];
				Y22 = m_e->Q[1] * m_e->Q[1] + m_e->P[1] * m_e->P[1];
				Y23 = m_e->Q[1] * m_e->Q[2] + m_e->P[1] * m_e->P[2];
				Y33 = m_e->Q[2] * m_e->Q[2] + m_e->P[2] * m_e->P[2];

				Y11 /= 4. * m_e->AREA;
				Y12 /= 4. * m_e->AREA;
				Y13 /= 4. * m_e->AREA;
				Y22 /= 4. * m_e->AREA;
				Y23 /= 4. * m_e->AREA;
				Y33 /= 4. * m_e->AREA;

				double v[3];
				v[0] = Y11*A_old(k) + Y12*A_old(m) + Y13*A_old(n);
				v[1] = Y12*A_old(k) + Y22*A_old(m) + Y23*A_old(n);
				v[2] = Y13*A_old(k) + Y23*A_old(m) + Y33*A_old(n);

				Y11 /= m_e->miut * ydot[i];
				Y12 /= m_e->miut * ydot[i];
				Y13 /= m_e->miut * ydot[i];
				Y22 /= m_e->miut * ydot[i];
				Y23 /= m_e->miut * ydot[i];
				Y33 /= m_e->miut * ydot[i];

				//计算偏导数
				double tmp=0;
				tmp = materialList[m_e->domain - 1].getdvdB(m_e->B);
				if (m_e->B > 1e-9){
					tmp /= m_e->B * m_e->AREA;//B==0?
					tmp /= ydot[i] * ydot[i] * ydot[i];
				}
				Y12 += v[0] * v[1] * tmp;//注意使用的时候要取负号
				Y23 += v[1] * v[2] * tmp;
				Y13 += v[0] * v[2] * tmp;				
				
				//计算牛顿迭代产生的电流
				INL(k) += v[0] * v[0] * tmp*A_old(k) + v[0] * v[1] * tmp*A_old(m) + v[0] * v[2] * tmp*A_old(n);
				INL(m) += v[1] * v[0] * tmp*A_old(k) + v[1] * v[1] * tmp*A_old(m) + v[1] * v[2] * tmp*A_old(n);
				INL(n) += v[2] * v[0] * tmp*A_old(k) + v[2] * v[1] * tmp*A_old(m) + v[2] * v[2] * tmp*A_old(n);
				//计算从线性系统反射回电压
				Vr[j].V12 = (A(k) - A(m)) - Vi[j].V12;
				Vr[j].V23 = (A(m) - A(n)) - Vi[j].V23;
				Vr[j].V13 = (A(n) - A(k)) - Vi[j].V13;
				//计算入射向线性系统的电压，如果电阻为负，则可能出错
				if (Y12 == 0){//受控电流源
					Vi[j].V12 = Vr[j].V12 + (A(k) - A(m))*Y12 / abs(rm[i].Y12);
				} else{//纯电阻
					Vi[j].V12 = Vr[j].V12*(abs(rm[i].Y12) + (Y12)) / (abs(rm[i].Y12) - (Y12));					
				}
				if (Y23 == 0){//受控电流源
					Vi[j].V23 = Vr[j].V23 + (A(m) - A(n))*Y23 / abs(rm[i].Y23);
				} else{//纯电阻
					Vi[j].V23 = Vr[j].V23*(abs(rm[i].Y23) + (Y23)) / (abs(rm[i].Y23) - (Y23));					
				}
				if (Y13 == 0){//受控电流源
					Vi[j].V13 = Vr[j].V13 + (A(n) - A(k))*Y13 / abs(rm[i].Y13);
				} else{//纯电阻
					Vi[j].V13 = Vr[j].V13*(abs(rm[i].Y13) + (Y13)) / (abs(rm[i].Y13) - (Y13));					
				}	
				INL(k) += 2 * Vi[j].V12*abs(rm[i].Y12);
				INL(m) -= 2 * Vi[j].V12*abs(rm[i].Y12);

				INL(m) += 2 * Vi[j].V23*abs(rm[i].Y23);
				INL(n) -= 2 * Vi[j].V23*abs(rm[i].Y23);

				INL(n) += 2 * Vi[j].V13*abs(rm[i].Y13);
				INL(k) -= 2 * Vi[j].V13*abs(rm[i].Y13);
			}
			//更新电流
			INL += bbJz;
			//INL.save("INL.txt", arma::arma_ascii, false);
			for (int i = 0; i < num_pts - node_bdr; i++) {
				unknown_b[i] = INL(node_reorder(i));
			}			
			
			//使用superLU_MT进行三角求解线性系统
			if (superlumt.triangleSolve() == 1) {
//				qDebug() << "Error: superlumt.slove";
//				qDebug() << "info: " << superlumt.info;
				break;
			} else {
				double *sol = NULL;
				A_tlm = A;//保存上一步的求解值
				sol = superlumt.get1Result();

				for (int i = 0; i < num_pts - node_bdr; i++) {
					//pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
					A(node_reorder(i)) = sol[i];
				}
			}
			//A_tlm.save("A_tlmNRA.txt", arma::arma_ascii, false);
			//A.save("NRA.txt", arma::arma_ascii, false);
			//if (iter == 1){
			//	int a = 1;
			//}
			//判断收敛
			double inner_error = norm((A_tlm - A), 2) / norm(A, 2);
//			qDebug() << iter_tlm << inner_error;
			if (inner_error < tlm_tol && iter_tlm>5) {
				//A.save("NRA.txt", arma::arma_ascii, false);
				//A_tlm.save("A_tlmNRA.txt", arma::arma_ascii, false);
//				qDebug() << "TLM steps: " << iter_tlm << " in " << iter << "th NR steps";
				
				break;
			}
			//重置电流
			bn.zeros();
			INL.zeros();
		}		
		
		//有必要更新所有单元的B值
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
		//计算牛顿迭代误差
		double error = norm((A_old - A), 2) / norm(A, 2);
		//qDebug() << "negY: " << negY << " in " << iter << "th NR steps";
		iter++;
//		qDebug() << "iter: " << iter;
//		qDebug() << "error: " << error;
		if (error < Precision && iter > 3) {
			//A.save("NRA.txt", arma::arma_ascii, false);
			break;
		}
		bn.zeros();
		
		pos = 0;

		//graph1->setData(x, y);
		//customplot->rescaleAxes(true);
		//customplot->replot();
	}
	//信息输出
	time[tt++] = clock();
	for (int i = 1; i < tt; i++){
//		qDebug() << time[i] - time[i - 1];
	}
//	qDebug() << "NR steps: " << iter;

	//空间回收
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (unknown_b != NULL) free(unknown_b);

	return true;
}
//使用牛顿迭代实现非线性求解，是真的牛顿迭代而不是松弛迭代
//牛顿迭代的公式推导，参见颜威利数值分析教材P54
int CFastFEMcore::StaticAxisymmetricNR() {
	clock_t time[10];
	int tt = 0;
	time[tt++] = SuperLU_timer_();
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
	int negY = 0;
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

			if (ce[0][1] > 0){
				negY++;
			}
			if (ce[0][2] > 0){
				negY++;
			}
			if (ce[1][2] > 0){
				negY++;
			}
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
//			qDebug() << "Error: superlumt.slove";
//			qDebug() << "info: " << superlumt.info;
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
		FILE *fp1 = fopen("B_T3_NR.txt", "w");
		for (int i = 0; i < num_ele; i++) {
			double bx = 0;
			double by = 0;
			for (int j = 0; j < 3; j++) {
				bx += pmeshele[i].Q[j] * A(pmeshele[i].n[j]);
				by += pmeshele[i].P[j] * A(pmeshele[i].n[j]);
			}
			pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
			pmeshele[i].Bx = bx / 2. / pmeshele[i].AREA / ydot[i] ;
			pmeshele[i].By = by / 2. / pmeshele[i].AREA / ydot[i];
			pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
			fprintf(fp1, "%lf \t %lf \t %lf \t %lf \t %lf\n", pmeshele[i].rc, pmeshele[i].zc, pmeshele[i].Bx, pmeshele[i].By, pmeshele[i].B);
			y[i] = pmeshele[i].miut;
		}
		fclose(fp1);
		double error = norm((A_old - A), 2) / norm(A, 2);
		iter++;
//		qDebug() << "iter: " << iter;
//		qDebug() << "error: " << error;
		//qDebug() << negY << " in " << iter << "th NR steps";
		negY = 0;
		if (error < Precision || iter > 20) {
			
			break;
		}
		bn.zeros();
		pos = 0;

		//graph1->setData(x, y);
		//customplot->rescaleAxes(true);
		//customplot->replot();
	}
	//write B resluts
	time[tt++] = SuperLU_timer_();
	for (int i = 1; i < tt; i++){
//		qDebug() << time[i] - time[i - 1];
	}
	A.save("NR_T3_A.txt", arma::arma_ascii, false);
//	qDebug() << "NR steps: " << iter;
//	qDebug() << "NR time: " << time[1] - time[1 - 1];
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
//		qDebug() << "Error: openning file!";
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
//		qDebug() << "Error: reading num_pts!";
		return 1;
	}
	int pts_ind;//the beginning of the points index
	//读取节点索引，默认从0开始
	if (fscanf(fp, "%d # lowest mesh point index\n", &pts_ind) != 1) {
//		qDebug() << "Error: reading pts_ind!";
		return 1;
	}
	fgets(ch, 256, fp);

	for (int i = pts_ind; i < num_pts; i++) {
		//读取x,y坐标
		if (fscanf(fp, "%lf %lf \n", &(pmeshnode[i].x), &(pmeshnode[i].y)) != 2) {
//			qDebug() << "Error: reading mesh point!";
			return 1;
		}
	}
	//---------------vertex-------------------------------
	for (int i = 0; i < 7; i++)
		fgets(ch, 256, fp);
	int num_vtx_ns, num_vtx_ele;
	//
	if (fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns) != 1) {
//		qDebug() << "Error: reading num_vtx_ns!";
		return 1;
	}

	if (fscanf(fp, "%d # number of elements\n", &num_vtx_ele) != 1) {
//		qDebug() << "Error: reading num_vtx_ele!";
		return 1;
	}
	fgets(ch, 256, fp);

	int *vtx;
	vtx = (int*)calloc(num_vtx_ele, sizeof(int));
	for (int i = 0; i < num_vtx_ele; i++) {
		//好象是每一个域的顶点编号
		if (fscanf(fp, "%d \n", vtx + i) != 1) {
//			qDebug() << "Error: reading vertex condition!";
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
//			qDebug() << "Error: reading vertex condition!";
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
//		qDebug() << "Error: reading num_bdr_ns!";
		return 1;
	}
	//读取线段边界数目
	if (fscanf(fp, "%d # number of elements\n", &num_bdr_ele) != 1) {
//		qDebug() << "Error: reading num_bdr_ele!";
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
//			qDebug() << "Error: reading boundary condition!";
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
//			qDebug() << "Error: reading boundary condition!";
			return 1;
		}
	}
	if (entity != NULL) free(entity); entity = NULL;
	//----------------elements------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int ns_per_ele;// num_ele;//number of nodes per element;number of elements
	if (fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele) != 1) {
//		qDebug() << "Error: reading ns_per_ele!";
		return 1;
	}
	//读取分网单元数目
	if (fscanf(fp, "%d # number of elements\n", &num_ele) == 1) {
		pmeshele4 = (CElement4*)calloc(num_ele, sizeof(CElement4));
	} else {
//		qDebug() << "Error: reading num_ele!";
		return 1;
	}
	fgets(ch, 256, fp);
	//读取分网四边形单元的四个节点索引
	//注意这里有点问题，comsol文件似乎不是按照逆时针存储的，而是1->2->4->3
	for (int i = 0; i < num_ele; i++) {
		if (fscanf(fp, "%d %d %d %d\n", &pmeshele4[i].n[0], &pmeshele4[i].n[1], &pmeshele4[i].n[3], &pmeshele4[i].n[2]) != 4) {
//			qDebug() << "Error: reading elements points!";
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
//			qDebug() << "Error: reading domain points!";
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
	if (std::isnan(dvdB)){
//		qDebug() << "error in getD: dvdB" << getx(xi, eta, index);
		return 1;
	}
		
	//计算Cij
	//积分
	if (pmeshele4[index].LinearFlag){
		return 0;
	}

	if (std::isnan(D))
//		qDebug() << "error in getD: D";
	D = dvdB*c1*c2 * getJacobi(xi, eta, index);
	if (B > 1e-9){
		D /= B * x * x * x;
	}
	if (std::isnan(D))
//		qDebug() << "error in getD";
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
	p /= getx(xi, eta, index);//可能为0的bug
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
//		qDebug() << B << miu << getx(xi, eta, index) << by << p;
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
	if (tmp <= 0)
//		qDebug() << "error: Jacobi.";
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
	double time[10]; int tt = 10;
	time[tt++] = SuperLU_timer_();
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
	double error;
	int pos = 0;
	time[tt++] = SuperLU_timer_();
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
				bbJz(pmeshele4[i].n[row]) += -hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
				bbJz(pmeshele4[i].n[kk]) += -hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
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
//			qDebug() << "Error: superlumt.slove";
//			qDebug() << "info: " << superlumt.info;
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
		error = norm((A_old - A), 2) / norm(A, 2);
		iter++;
//		qDebug() << "iter: " << iter;
//		qDebug() << "error: " << error;
		if (error < Precision && iter > 5) {
			//转换A
			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;//the A is r*A_real
				A(node_reorder(i)) /= pmeshnode[node_reorder(i)].x;
			}
			//A.save("D:\\mypaper\\zhcore\\插图\\NRpmA.txt", arma::arma_ascii, false);
			//bbJz.save("D:\\mypaper\\zhcore\\插图\\pmbn.txt", arma::arma_ascii, false);
			//qDebug() << "solve over";
			break;
		}
		//重新初始化
		pos = 0;
		bbJz.zeros();
	}
	time[tt++] = SuperLU_timer_();
//	qDebug() << "iter: " << iter;
//	qDebug() << "error: " << error;
//	qDebug() <<"Single step time of NR: "<< (time[tt - 1] - time[tt - 2]) / iter;
//	qDebug() << "Total step time of NR: " << (time[tt - 1] - time[tt - 2]);
	A.save("D:\\NRpmA.txt", arma::arma_ascii, false);
	return true;
}
//采用传输线法求解轴对称静磁场，四边形分网
//该方法的主要思路是将每个单元的等效电路当作一个黑盒子，不去考虑具体内部是什么电路
bool CFastFEMcore::StaticAxisQ4TLM(){
	double time[10]; int  tt = 0;
	time[tt++] = SuperLU_timer_();
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
			pmeshele4[i].miut = 1000 * miu0;
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
//	qDebug() << "D34 size: " << D34.size();
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
	time[tt++] = SuperLU_timer_();
	//计算导纳矩阵
	int nonlinr = -1;
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele4[i].LinearFlag){
			nonlinr++;
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
					//非线性部分，采用传输线,仅含有对地支路
					if (row == col){
						rm[nonlinr].Y[row] = 500*getLocal4Matrix(0, 0, i);// *pmeshnode[pmeshele4[i].n[0]].A;
						rm[nonlinr].Y[row] += 500*getLocal4Matrix(1, 1, i);// *pmeshnode[pmeshele4[i].n[1]].A;
						rm[nonlinr].Y[row] += 500*getLocal4Matrix(2, 2, i);// *pmeshnode[pmeshele4[i].n[2]].A;
						rm[nonlinr].Y[row] += 500*getLocal4Matrix(3, 3, i);// *pmeshnode[pmeshele4[i].n[3]].A;
						//rm[nonlinr].Y[row] = abs(rm[nonlinr].Y[row] / pmeshnode[pmeshele4[i].n[row]].A );
						if (std::isinf(rm[nonlinr].Y[row])){
							//qDebug() << pmeshnode[pmeshele4[i].n[row]].A;
							rm[nonlinr].Y[row] = abs(getLocal4Matrix(row, col, i));
						}
							
						ce[row][col] = rm[nonlinr].Y[row];
					} else{
						ce[row][col] = 0;
					}
				}
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
	time[tt++] = SuperLU_timer_();
	INL += bbJz;
	for (int i = 0; i < num_pts - node_bdr; i++) {
		unknown_b[i] = INL(node_reorder(i));
	}
	//第一次求解
	//---------------------superLU_MT---------------------------------------
	CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
	if (superlumt.solve1() == 1) {
//		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = superlumt.get1Result();
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			//pmeshnode[node_reorder(i)].A *= pmeshnode[node_reorder(i)].x;// / pmeshnode[i].x;//the A is r*A_real
			A(node_reorder(i)) = pmeshnode[node_reorder(i)].A ;
		}
	}
	time[tt++] = SuperLU_timer_();
	double a1, a2;
	//迭代
	int count;//迭代步数
	double error;
	VoltageQ4 *Vr = (VoltageQ4*)calloc(D34.size(), sizeof(VoltageQ4));
	VoltageQ4 *Vi = (VoltageQ4*)calloc(D34.size(), sizeof(VoltageQ4));
	for (count = 0; count < 500; count++){
		a1 = SuperLU_timer_();
		//反射到单个非线性单元进行迭代
		for (int j = 0; j < D34.size(); j++){
			int i = D34[j];
			CElement4 *m_e = pmeshele4 + i;
			int n1, n2, n3, n4;
			n1 = m_e->n[0];
			n2 = m_e->n[1];
			n3 = m_e->n[2];
			n4 = m_e->n[3];
			//计算反射电压
			Vr[j].V[0] = A(n1) - Vi[j].V[0];
			Vr[j].V[1] = A(n2) - Vi[j].V[1];
			Vr[j].V[2] = A(n3) - Vi[j].V[2];
			Vr[j].V[3] = A(n4) - Vi[j].V[3];
			//使用牛顿迭代求解小电路
			mat AJ(4, 4);
			colvec b(4);
			colvec x2(4); x2.zeros();
			double err1 = 0;

			for (int iter = 0; iter < 5; iter++){
				//流向节点的电流源
				b(0) = 2 * Vr[j].V[0] * rm[j].Y[0];
				b(1) = 2 * Vr[j].V[1] * rm[j].Y[1];
				b(2) = 2 * Vr[j].V[2] * rm[j].Y[2];
				b(3) = 2 * Vr[j].V[3] * rm[j].Y[3];

				for (int row = 0; row < 4; row++){
					for (int col = 0; col < 4; col++){
						AJ(row, col) = getLocal4Matrix(row, col, i) + getDij(row, col, i);
						b(row) += getDij(row, col, i)*pmeshnode[m_e->n[col]].A;
					}
				}
				//加上对地导纳
				AJ(0, 0) += rm[j].Y[0];
				AJ(1, 1) += rm[j].Y[1];
				AJ(2, 2) += rm[j].Y[2];
				AJ(3, 3) += rm[j].Y[3];
				//更新电压
				pmeshnode[n1].A = x2(0);
				pmeshnode[n2].A = x2(1);
				pmeshnode[n3].A = x2(2);
				pmeshnode[n4].A = x2(3);
				//求解电压V
				bool status = arma::solve(x2, AJ, b);
				if (!status){
//					qDebug() << "error: solve !";
					return false;
				}			
				
				//判断收敛
				/*double error1 = 0;
				for (int ii = 0; ii < 4; ii++){
					error1 += (x2(ii) -)* x2(ii);
				}
				double error1 = abs(x2(0)) + abs(x2(1)) + abs(x2(2)) + abs(x2(3));
				error1 -= abs(pmeshnode[n1].A) + abs(pmeshnode[n2].A) + abs(pmeshnode[n3].A) + abs(pmeshnode[n4].A);
				error1 /= 4;
				if (abs(error1) < 1e-5) break;*/
			}
			//计算入射电压
			Vi[j].V[0] = x2(0) - Vr[j].V[0];
			Vi[j].V[1] = x2(1) - Vr[j].V[1];
			Vi[j].V[2] = x2(2) - Vr[j].V[2];
			Vi[j].V[3] = x2(3) - Vr[j].V[3];
			//入射电流源项
			INL(n1) += 2 * Vi[j].V[0] * rm[j].Y[0];
			INL(n2) += 2 * Vi[j].V[1] * rm[j].Y[1];
			INL(n3) += 2 * Vi[j].V[2] * rm[j].Y[2];
			INL(n4) += 2 * Vi[j].V[3] * rm[j].Y[3];
		}
		a2 = SuperLU_timer_();
//		qDebug() << "t1: " << a2 - a1;
		//入射到线性网络过程
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}		
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		if (superlumt.triangleSolve() == 1) {
//			qDebug() << "Error: superlumt.slove";
			//qDebug() << "info: " << superlumt.info;
			break;
		} else {
			double *sol = NULL;
			A_old = A;
			sol = superlumt.get1Result();

			for (int i = 0; i < num_pts - node_bdr; i++) {
				//pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
				A(node_reorder(i)) = sol[i];
			}
		}
		error = norm((A_old - A), 2) / norm(A, 2);
//		qDebug() << "steps: " << count;
//		qDebug() <<"error: "<< error;
		if (error < Precision && count > 10) {
			//转换A
			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;//the A is r*A_real
				A(node_reorder(i)) /= pmeshnode[node_reorder(i)].x;
			}
			//A.save("D:\\mypaper\\zhcore\\插图\\TLMpmA.txt", arma::arma_ascii, false);
			//bbJz.save("D:\\mypaper\\zhcore\\插图\\pmbn.txt", arma::arma_ascii, false);
//			qDebug() << "Solve Successfully!";
			break;
		}
		INL.zeros();
		a1 = SuperLU_timer_();
//		qDebug() << "t2: " << a1 - a2;
	}
	time[tt++] = SuperLU_timer_();
//	qDebug() << "steps: " << count;
//	qDebug() <<"error: "<< error;
//	qDebug() <<"Single step time of TLM: "<< (time[tt - 1] - time[tt - 2]) / count;
//	qDebug() << "Total time of TLM: " << (time[tt - 1] - time[tt - 2]);
	//
	for (int i = 0; i < num_pts - node_bdr; i++) {
		pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;//the A is r*A_real
		A(node_reorder(i)) /= pmeshnode[node_reorder(i)].x;
	}
	A.save("D:\\TLMpmA.txt", arma::arma_ascii, false);
	return true;
}
//2018-03-05
//by Poofee
//采用VTM方法对四边形问题进行求解
bool CFastFEMcore::StaticAxisQ4VTM(){
	double time[10]; int  tt = 0;
	time[tt++] = SuperLU_timer_();
	//PART A:计算边界信息	
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
			pmeshele4[i].miut = 1000 * miu0;
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
	//PART B: 寻找非线性单元
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (!pmeshele4[i].LinearFlag) {
			D34.push_back(i);
		}
	}
//	qDebug() << "D34 size: " << D34.size();
	umat locs(2, 16 * num_ele);
	locs.zeros();
	mat vals(1, 16 * num_ele);
	double ce[4][4] = { 0 };
	Resist4Matrix *rm = (Resist4Matrix*)malloc(D34.size() * sizeof(Resist4Matrix));
	vec bbJz = zeros<vec>(num_pts);
	vec A = zeros<vec>(num_pts);
	vec A_old = A;
	vec INL = zeros<vec>(num_pts);
	//PART C: 根据边界条件对节点重新进行编号
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
	time[tt++] = SuperLU_timer_();
	//PART D: 有限元装配
	int nonlinr = -1;
	for (int i = 0; i < num_ele; i++) {
		//将单元矩阵进行存储
		double jr = materialList[pmeshele4[i].domain - 1].Jr;
		double hc = materialList[pmeshele4[i].domain - 1].H_c;
		for (int row = 0; row < 4; row++) {
			for (int col = 0; col < 4; col++) {
				//利用对称性，只计算一半
				if (row > col){
					ce[row][col] = ce[col][row];
				} else{
					ce[row][col] = getLocal4Matrix(row, col, i);
				}							
				
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
			//计算电流密度
			bbJz(pmeshele4[i].n[row]) += jr*getJi(row, i);
			// 计算永磁部分
			int kk = row + 1; if (kk == 4) kk = 0;
			bbJz(pmeshele4[i].n[row]) += -hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
			bbJz(pmeshele4[i].n[kk]) += -hc / 2.*(pmeshnode[pmeshele4[i].n[kk]].x - pmeshnode[pmeshele4[i].n[row]].x);
		}
		//计算导纳矩阵
		if (!pmeshele4[i].LinearFlag){
			nonlinr++;
			rm[nonlinr].Y[0] = abs(ce[0][0]);
			rm[nonlinr].Y[1] = -abs(ce[0][1]);
			rm[nonlinr].Y[2] = -abs(ce[0][2]);
			rm[nonlinr].Y[3] = -abs(ce[0][3]);
			rm[nonlinr].Y[4] = abs(ce[1][1]);
			rm[nonlinr].Y[5] = -abs(ce[1][2]);
			rm[nonlinr].Y[6] = -abs(ce[1][3]);
			rm[nonlinr].Y[7] = abs(ce[2][2]);
			rm[nonlinr].Y[8] = -abs(ce[2][3]);
			rm[nonlinr].Y[9] = abs(ce[3][3]);
		}		
	}//end for
	locs.reshape(2, pos);//重新调整大小
	vals.reshape(1, pos);
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);
	time[tt++] = SuperLU_timer_();
	INL += bbJz;
	for (int i = 0; i < num_pts - node_bdr; i++) {
		unknown_b[i] = INL(node_reorder(i));
	}
	//第一次求解
	//---------------------superLU_MT---------------------------------------
	CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
	if (superlumt.solve1() == 1) {
//		qDebug() << "Error: superlumt.slove";
		//qDebug() << "info: " << superlumt.info;
	} else {
		double *sol = NULL;
		A_old = A;
		sol = superlumt.get1Result();
		//取得结果
		for (int i = 0; i < num_pts - node_bdr; i++) {
			//pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
			A(node_reorder(i)) = sol[i];
		}
	}
	time[tt++] = SuperLU_timer_();
	//迭代
	int count;//迭代步数
	double error;
	//边界上的节点电压
	VoltageQ4 *Ve = (VoltageQ4*)calloc(D34.size(), sizeof(VoltageQ4));
	VoltageQ4 *Vs = (VoltageQ4*)calloc(D34.size(), sizeof(VoltageQ4));
	//边界上的节点电流
	VoltageQ4 *Ie = (VoltageQ4*)calloc(D34.size(), sizeof(VoltageQ4));
	VoltageQ4 *Is = (VoltageQ4*)calloc(D34.size(), sizeof(VoltageQ4));
	//传输线迭代
	for (count = 0; count < 500; count++){
		////PART 1：计算非线性单元部分
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++){
			int i = D34[j];
			CElement4 *m_e = pmeshele4 + i;
			int n1, n2, n3, n4;
			n1 = m_e->n[0];
			n2 = m_e->n[1];
			n3 = m_e->n[2];
			n4 = m_e->n[3];
			//PART 2：计算流入非线性单元侧的电流
			double Ielast[4];//保存上一步的值
			Ielast[0] = Ie[j].V[0];
			Ielast[1] = Ie[j].V[1];
			Ielast[2] = Ie[j].V[2];
			Ielast[3] = Ie[j].V[3];
			//计算电压
			Vs[j].V[0] = A(n1);
			Vs[j].V[1] = A(n2);
			Vs[j].V[2] = A(n3);
			Vs[j].V[3] = A(n4);

			Ie[j].V[0] = 2 * (rm[j].Y[0] * Vs[j].V[0] + rm[j].Y[1] * Vs[j].V[1] + rm[j].Y[2] * Vs[j].V[2] + rm[j].Y[3] * Vs[j].V[3]) - Is[j].V[0];
			Ie[j].V[1] = 2 * (rm[j].Y[1] * Vs[j].V[0] + rm[j].Y[4] * Vs[j].V[1] + rm[j].Y[5] * Vs[j].V[2] + rm[j].Y[6] * Vs[j].V[3]) - Is[j].V[1];
			Ie[j].V[2] = 2 * (rm[j].Y[2] * Vs[j].V[0] + rm[j].Y[5] * Vs[j].V[1] + rm[j].Y[7] * Vs[j].V[2] + rm[j].Y[8] * Vs[j].V[3]) - Is[j].V[2];
			Ie[j].V[3] = 2 * (rm[j].Y[3] * Vs[j].V[0] + rm[j].Y[6] * Vs[j].V[1] + rm[j].Y[8] * Vs[j].V[2] + rm[j].Y[9] * Vs[j].V[3]) - Is[j].V[3];
			//PART 3：计算流入线性系统侧的电流
			Is[j].V[0] = 2 * (rm[j].Y[0] * Ve[j].V[0] + rm[j].Y[1] * Ve[j].V[1] + rm[j].Y[2] * Ve[j].V[2] + rm[j].Y[3] * Ve[j].V[3]) - Ielast[0];
			Is[j].V[1] = 2 * (rm[j].Y[1] * Ve[j].V[0] + rm[j].Y[4] * Ve[j].V[1] + rm[j].Y[5] * Ve[j].V[2] + rm[j].Y[6] * Ve[j].V[3]) - Ielast[1];
			Is[j].V[2] = 2 * (rm[j].Y[2] * Ve[j].V[0] + rm[j].Y[5] * Ve[j].V[1] + rm[j].Y[7] * Ve[j].V[2] + rm[j].Y[8] * Ve[j].V[3]) - Ielast[2];
			Is[j].V[3] = 2 * (rm[j].Y[3] * Ve[j].V[0] + rm[j].Y[6] * Ve[j].V[1] + rm[j].Y[8] * Ve[j].V[2] + rm[j].Y[9] * Ve[j].V[3]) - Ielast[3];

			//入射电流源项
			INL(n1) += Is[j].V[0];
			INL(n2) += Is[j].V[1];
			INL(n3) += Is[j].V[2];
			INL(n4) += Is[j].V[3];
			//PART 4: 使用牛顿迭代求解小电路
			mat AJ(4, 4);
			colvec b(4);
			colvec x2(4); x2.zeros();
			double err1 = 0;


			pmeshnode[n1].A = 0;// Ve[j].V[0];
			pmeshnode[n2].A = 0;// Ve[j].V[1];
			pmeshnode[n3].A = 0;// Ve[j].V[2];
			pmeshnode[n4].A = 0;// Ve[j].V[3];
			for (int iter = 0; iter < 30; iter++){
				//1.初始化电流
				b(0) = Ie[j].V[0];
				b(1) = Ie[j].V[1];
				b(2) = Ie[j].V[2];
				b(3) = Ie[j].V[3];
				//2.牛顿迭代
				for (int row = 0; row < 4; row++){
					for (int col = 0; col < 4; col++){
						if (row > col){
							AJ(row, col) = AJ(col, row);
						} else{
							AJ(row, col) = getLocal4Matrix(row, col, i) + getDij(row, col, i);
						}						
						b(row) += getDij(row, col, i)*pmeshnode[m_e->n[col]].A;
					}
				}
				//3.加上传输线导纳
				AJ(0, 0) += rm[j].Y[0];
				AJ(0, 1) += rm[j].Y[1];
				AJ(0, 2) += rm[j].Y[2];
				AJ(0, 3) += rm[j].Y[3];

				AJ(1, 0) += rm[j].Y[1];
				AJ(1, 1) += rm[j].Y[4];
				AJ(1, 2) += rm[j].Y[5];
				AJ(1, 3) += rm[j].Y[6];

				AJ(2, 0) += rm[j].Y[2];
				AJ(2, 1) += rm[j].Y[5];
				AJ(2, 2) += rm[j].Y[7];
				AJ(2, 3) += rm[j].Y[8];

				AJ(3, 0) += rm[j].Y[3];
				AJ(3, 1) += rm[j].Y[6];
				AJ(3, 2) += rm[j].Y[8];
				AJ(3, 3) += rm[j].Y[9];
				//更新电压
				pmeshnode[n1].A = x2(0);
				pmeshnode[n2].A = x2(1);
				pmeshnode[n3].A = x2(2);
				pmeshnode[n4].A = x2(3);
				//求解电压V
				bool status = arma::solve(x2, AJ, b);
				if (!status){
//					qDebug() << "error: solve !";
					return false;
				}				
				//判断收敛，x2比较小，要用相对误差比较！
				double error1 = x2(0) + x2(1) + x2(2) + x2(3);
				error1 -= pmeshnode[n1].A + pmeshnode[n2].A + pmeshnode[n3].A + pmeshnode[n4].A;
				error1 /= 4;
				//if (abs(error1) < 1e-5) break;
			}
			//计算电压
			Ve[j].V[0] = x2(0);
			Ve[j].V[1] = x2(1);
			Ve[j].V[2] = x2(2);
			Ve[j].V[3] = x2(3);			
		}
		//入射到线性网络过程
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		if (superlumt.triangleSolve() == 1) {
//			qDebug() << "Error: superlumt.slove";
			//qDebug() << "info: " << superlumt.info;
			break;
		} else {
			double *sol = NULL;
			A_old = A;
			sol = superlumt.get1Result();

			for (int i = 0; i < num_pts - node_bdr; i++) {
				//pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
				A(node_reorder(i)) = sol[i];
			}
		}
		error = norm((A_old - A), 2) / norm(A, 2);
//		qDebug() << "steps: " << count;
//		qDebug() << "error: " << error;
		if (error < Precision && count > 10) {
			//转换A
			for (int i = 0; i < num_pts - node_bdr; i++) {
				pmeshnode[node_reorder(i)].A /= pmeshnode[node_reorder(i)].x;//the A is r*A_real
				A(node_reorder(i)) /= pmeshnode[node_reorder(i)].x;
			}
			//A.save("D:\\mypaper\\zhcore\\插图\\TLMpmA.txt", arma::arma_ascii, false);
			//bbJz.save("D:\\mypaper\\zhcore\\插图\\pmbn.txt", arma::arma_ascii, false);
//			qDebug() << "Solve Successfully!";
			break;
		}
		INL.zeros();
	}
	time[tt++] = SuperLU_timer_();
//	qDebug() << "steps: " << count;
//	qDebug() << "error: " << error;
//	qDebug() << "Single step time of TLM: " << (time[tt - 1] - time[tt - 2]) / count;
//	qDebug() << "Total time of TLM: " << (time[tt - 1] - time[tt - 2]);
	//
	return true;
}
//2018-03-12
//by Poofee
//本函数实现用VTM方法来求解NR每一迭代步的线性方程组
bool CFastFEMcore::StaticAxisT3NRVTM(){
	
	return true;
}
//2018-03-29
//by Poofee
//本函数采用VTM模型实现非线性静磁场的求解
//三角单元的每一个电阻被添加传输线
bool CFastFEMcore::StaticAxisT3VTMsingle(){
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
	//传输线导纳矩阵，只保存上三角矩阵
	double * Ytl = (double*)malloc(D34.size() * 6 * sizeof(double));
	int nlin = 0;
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

		rm[i].Y11 /= 4. * pmeshele[i].AREA*ydot[i];//猜测值
		rm[i].Y12 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y13 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y22 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y23 /= 4. * pmeshele[i].AREA*ydot[i];
		rm[i].Y33 /= 4. * pmeshele[i].AREA*ydot[i];

		//生成单元矩阵，线性与非线性
		//因为线性与非线性的差不多，所以不再分开讨论了	

		//非线性区，计算传输线导纳矩阵
		double delta = 1e-3;//接地导纳
		if (!pmeshele[i].LinearFlag) {
			ce[0][1] = rm[i].Y12 / pmeshele[i].miu;
			ce[0][2] = rm[i].Y13 / pmeshele[i].miu;
			ce[1][2] = rm[i].Y23 / pmeshele[i].miu;

			ce[0][0] = rm[i].Y11 / pmeshele[i].miu;
			ce[1][1] = rm[i].Y22 / pmeshele[i].miu;
			ce[2][2] = rm[i].Y33 / pmeshele[i].miu;

			if (ce[0][1] > 0){
				int a = 1;
			}
			if (ce[0][2] > 0){
				int a = 0;
			}
			if (ce[1][2] > 0){
				int a = 1;
			}

			Ytl[6 * nlin + 0] = ce[0][0];
			Ytl[6 * nlin + 1] = ce[0][1];
			Ytl[6 * nlin + 2] = ce[0][2];
			Ytl[6 * nlin + 3] = ce[1][1];
			Ytl[6 * nlin + 4] = ce[1][2];
			Ytl[6 * nlin + 5] = ce[2][2];
			nlin++;
		} else{
			ce[0][1] = rm[i].Y12 / pmeshele[i].miu;
			ce[0][2] = rm[i].Y13 / pmeshele[i].miu;
			ce[1][2] = rm[i].Y23 / pmeshele[i].miu;

			ce[0][0] = rm[i].Y11 / pmeshele[i].miu;
			ce[1][1] = rm[i].Y22 / pmeshele[i].miu;
			ce[2][2] = rm[i].Y33 / pmeshele[i].miu;
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
//		qDebug() << "Error: superlumt.slove";
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
	//保存非线性单元节点电压，不是反射电压了
	Voltage *Vs = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of system
	Voltage *Ve = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Voltage of element
	//保存电流
	Voltage *Is = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into system
	Voltage *Ie = (Voltage*)calloc(D34.size(), sizeof(Voltage));//Current flowing into element
	time[tt++] = SuperLU_timer_();
	for (count = 0; count < steps; count++) {
		//PART A：计算非线性单元部分
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *m_e = pmeshele + i;
			int k, m, n;
			k = m_e->n[0];
			m = m_e->n[1];
			n = m_e->n[2];
			//PART B：计算流入非线性单元侧的电流
			double Vlast[3];//保存上一步的值
			Vlast[0] = Ie[j].V12;
			Vlast[1] = Ie[j].V23;
			Vlast[2] = Ie[j].V13;
			//计算电压
			Vs[j].V12 = pmeshnode[k].A - 0;
			Vs[j].V23 = pmeshnode[m].A - 0;
			Vs[j].V13 = pmeshnode[n].A - 0;

			Ie[j].V12 = 2 * (Ytl[6 * j + 0] * Vs[j].V12 + Ytl[6 * j + 1] * Vs[j].V23 + Ytl[6 * j + 2] * Vs[j].V13) - Is[j].V12;
			Ie[j].V23 = 2 * (Ytl[6 * j + 1] * Vs[j].V12 + Ytl[6 * j + 3] * Vs[j].V23 + Ytl[6 * j + 4] * Vs[j].V13) - Is[j].V23;
			Ie[j].V13 = 2 * (Ytl[6 * j + 2] * Vs[j].V12 + Ytl[6 * j + 4] * Vs[j].V23 + Ytl[6 * j + 5] * Vs[j].V13) - Is[j].V13;
			//PART C：计算流入线性系统侧的电流
			Is[j].V12 = 2 * (Ytl[6 * j + 0] * Ve[j].V12 + Ytl[6 * j + 1] * Ve[j].V23 + Ytl[6 * j + 2] * Ve[j].V13) - Vlast[0];
			Is[j].V23 = 2 * (Ytl[6 * j + 1] * Ve[j].V12 + Ytl[6 * j + 3] * Ve[j].V23 + Ytl[6 * j + 4] * Ve[j].V13) - Vlast[1];
			Is[j].V13 = 2 * (Ytl[6 * j + 2] * Ve[j].V12 + Ytl[6 * j + 4] * Ve[j].V23 + Ytl[6 * j + 5] * Ve[j].V13) - Vlast[2];

			INL(k) += Is[j].V12;
			INL(m) += Is[j].V23;
			INL(n) += Is[j].V13;
			//PART D：牛顿迭代求解小电路，计算电压
			mat C(3, 3);//单元系数矩阵，为了方便计算
			C(0, 0) = -rm[i].Y12; C(0, 1) = rm[i].Y12; C(0, 2) = 0;
			C(1, 0) = 0; C(1, 1) = -rm[i].Y23; C(1, 2) = rm[i].Y23;
			C(2, 0) = rm[i].Y13; C(2, 1) = 0; C(2, 2) = -rm[i].Y13;
			mat D(3, 3);//单元系数矩阵，为了方便计算
			D(0, 0) = rm[i].Y11; D(0, 1) = rm[i].Y12; D(0, 2) = rm[i].Y13;
			D(1, 0) = rm[i].Y12; D(1, 1) = rm[i].Y22; D(1, 2) = rm[i].Y23;
			D(2, 0) = rm[i].Y13; D(2, 1) = rm[i].Y23; D(2, 2) = rm[i].Y33;
			mat AJ(3, 3);
			colvec b(3);
			colvec x2(3); x2.zeros();
			double err1 = 0;

			for (int iter = 0; iter < 20; iter++){
				//1.初始化电流
				b(0) = Ie[j].V12;
				b(1) = Ie[j].V23;
				b(2) = Ie[j].V13;
				//
				AJ.zeros();
				double ca[3];//理论上来说，这里的A应当是三个节点的A
				double cb[3];
				for (int row = 0; row < 3; row++){
					ca[row] = C(row, 0)*x2(0) + C(row, 1)*x2(1) + C(row, 2)*x2(2);
					cb[row] = D(row, 0)*x2(0) + D(row, 1)*x2(1) + D(row, 2)*x2(2);
				}
				//2.求解单元mu值
				double bx = 0;
				double by = 0;
				for (int j = 0; j < 3; j++) {
					bx += pmeshele[i].Q[j] * x2(j);
					by += pmeshele[i].P[j] * x2(j);
				}
				pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / ydot[i];
				pmeshele[i].miut = materialList[pmeshele[i].domain - 1].getMiu(pmeshele[i].B);
				y[i] = pmeshele[i].miut;
				//3.计算雅可比矩阵
				double tmp = materialList[pmeshele[i].domain - 1].getdvdB(pmeshele[i].B);
				if (pmeshele[i].B > 1e-9){
					tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
					tmp /= ydot[i];
				}
				for (int row = 0; row < 3; row++){
					for (int col = 0; col < 3; col++){
						//注意C已经被除了一次ydot了
						AJ(row, col) = ca[row];
						AJ(row, col) *= cb[col];
						AJ(row, col) *= tmp;
						b(row) += AJ(row, col)*x2(col);
						AJ(row, col) += C(row, col) / m_e->miut;
					}
				}
				//4.加上传输线导纳
				AJ(0, 0) += -Ytl[6 * j + 1];
				AJ(0, 1) += Ytl[6 * j + 1];
				//AJ(0, 2) += Ytl[6 * j + 2];
				//AJ(1, 0) += Ytl[6 * j + 1];
				AJ(1, 1) += -Ytl[6 * j + 4];
				AJ(1, 2) += Ytl[6 * j + 4];
				AJ(2, 0) += Ytl[6 * j + 2];
				//AJ(2, 1) += Ytl[6 * j + 4];
				AJ(2, 2) += -Ytl[6 * j + 2];
				//5.求解电压V
				bool status = arma::solve(x2, AJ, b);
				if (!status){
//					qDebug() << "error: solve !";
					return false;
				}
			}

			//qDebug() << x2(0) << x2(1) << x2(2);
			Ve[j].V12 = x2(0);
			Ve[j].V23 = x2(1);
			Ve[j].V13 = x2(2);
		}
		INL += bbJz;
		for (int i = 0; i < num_pts - node_bdr; i++) {
			unknown_b[i] = INL(node_reorder(i));
		}
		//time[tt++] = SuperLU_timer_();
		//调用superLU三角求解
		dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
		//myTriSolve(1, &L, &U, perm_r, perm_c, &sluB, &info);
		//time[tt++] = SuperLU_timer_();
		//pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);
		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		//读取求解结果
		if (info != 0) {
//			qDebug() << "Error: superlumt.slove";
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
		//误差判断
		double error = norm((A_old - A), 2) / norm(A, 2);
//		qDebug() << "iter: " << count;
//		qDebug() << "error: " << error;

		//绘制计算结果
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
//		qDebug() << i << "\t" << time[i] - time[i - 1];
	}

//	qDebug() << "TLM steps:" << count;
	// 回收空间
	if (rm != NULL) free(rm);
	if (ydot != NULL) free(ydot);
	if (Vs != NULL) free(Vs);
	if (Ve != NULL) free(Ve);
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

bool CFastFEMcore::transientCaseNDDR(){

	return true;
}

bool CFastFEMcore::transientCaseTLM(){

	return true;
}

//2018-12-03
//by Poofee
//solve 2D transient electro-mechnical-circuit 
bool CFastFEMcore::transientCaseNR(){
	//define some variables
	double current_time = 0;//total time from simualation
	double current_position = 0;//the moving distance of mover
	double dtime;//time step length
	double MIN_TIME=0, MAX_TIME;//time constraint
	int total_time_step = 0;//number of time discretization
	int current_time_step = 0;//current time step

	double m;// the mass of the moving body
	double Fr;//the resistant mechanical force acting on the body

	QVector<double> time_list;//store time length per step
	QVector<double> position_list;//store mover's position per step
	QVector<double> Fem_list;//store the electromagnetic force acting on the body per step
	QVector<double> velocity_list;//store mover's velocity per step

	/////--------------start time loop resolution--------------////////
	//genertate discretized time list

	//begin time loop	

	//call mesh generator triangle.exe to remesh the geometry

	//where is the dB/dt???


	////----start to solve magnetostatic field-----//////
	//begin nonlinear loop

	//end nonlinear loop

	//update force, acceleration, velocity and new position
	//1.calcualte electromagnetic force

	//2.calculate acceleration

	//3.calculate new velocity

	//4.calculate new position of mover at netx time step

	//move the mover to new position
	//1.read the initial CAD file

	//2.create new CAD file with mover at new position

	//end time loop

	//output solution

	return true;
}
