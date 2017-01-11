//2016-02-26 by poofee
//add the readfile function to simplify the code.
//add the sourcecode in the git using version control.
//2016-01-20 by poofee
//cheng the Y0 as the diffrent numbers
//as in the papers say
//by poofee
//add the LU fact directly using superLU
//and we don't process the LU fact in every iteration.
#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <stdio.h>
#include <QtAlgorithms>
//#include <stdlib.h>
//#include <omp.h>
#include "datatype.h"
#include "spline.h"
#include <cmath>
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo> 
#include <vector>
#include <ctime>
#include <omp.h>
#include "slu_mt_ddefs.h"
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

MainWindow::MainWindow(QWidget *parent) :
QMainWindow(parent),
ui(new Ui::MainWindow) {
	ui->setupUi(this);
	QAction *start = new QAction("TLM", this);
	addAction(start);
	connect(start, SIGNAL(triggered()), this, SLOT(TLMcalculation()));
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


void MainWindow::TLMcalculation() {
	/*read coarse mesh data from file*/
	double U0 = 1;// PI*4e-7;
	char  coarseFN[] = "E:\\Projects\\axispmmodel\\mesh04.mphtxt";
	CNode * cmeshnode = NULL;	int cnum_pts = 0;
	CElement * cmeshele = NULL;	int cnum_ele = 0;
	readDataFile(coarseFN, cnum_pts, cnum_ele, &cmeshnode, &cmeshele);

	/*using TLM get a coarse miu value*/
	TLM(cnum_pts, cnum_ele, cmeshnode, cmeshele);
	/*read fine mesh data from file*/
	char fineFN[] = "E:\\Projects\\axispmmodel\\mesh04.mphtxt";
	CNode * fmeshnode = NULL;	int fnum_pts = 0;
	CElement * fmeshele = NULL;	int fnum_ele = 0;
	readDataFile(fineFN, fnum_pts, fnum_ele, &fmeshnode, &fmeshele);
	/*set the miu initial value */
	double minx, maxx, miny, maxy;
	QVector <ele> d34;
	ele eletmp;
	for (int i = 0; i < fnum_ele; i++) {
		if (fmeshele[i].domain == 3 || fmeshele[i].domain == 4) {
			eletmp.index = i;
			eletmp.x = fmeshele[i].rc;
			eletmp.y = fmeshele[i].zc;
			d34.push_back(eletmp);
		}
	}
	qSort(d34.begin(), d34.end(), compareele);//ascend

	for (int i = 0; i < cnum_ele; i++) {
		if (cmeshele[i].domain == 3 || cmeshele[i].domain == 4) {
			/*find the minx max miny maxy in the ele*/
			double a, b, c;
			a = cmeshnode[cmeshele[i].n[0]].x;
			b = cmeshnode[cmeshele[i].n[1]].x;
			c = cmeshnode[cmeshele[i].n[2]].x;
			maxx = mymax(mymax(a, b), c);
			minx = mymin(mymin(a, b), c);
			a = cmeshnode[cmeshele[i].n[0]].y;
			b = cmeshnode[cmeshele[i].n[1]].y;
			c = cmeshnode[cmeshele[i].n[2]].y;
			maxy = mymax(mymax(a, b), c);
			miny = mymin(mymin(a, b), c);
			/*binary search to find the index */
			int bottom = 0, top = d34.size() - 1;
			if (top == -1) {
				break;
			}
			int position;
			while (bottom <= top) {
				position = (bottom + top) >> 1;
				if (d34[position].x >= minx && position == 0) {
					break;
				}
				if (d34[position].x >= minx && d34[position - 1].x < minx) {
					break;
				}
				if (d34[position].x < minx) {
					bottom = position + 1;
				} else {
					top = position - 1;
				}
			}
			/*start from the position*/
			for (int j = position; j < d34.size(); j++) {
				/*find the first x between minx~maxx*/
				if (d34[j].x > maxx) {
					break;
				} else if (d34[j].x > minx) {
					if (d34[j].y > miny && d34[j].y < maxy) {
						fmeshele[d34[j].index].miut = cmeshele[i].miu * r;
						d34.remove(j);
						j -= 1;
					}
				}
			}
		}
	}
	/*recall the TLM*/
	TLM(fnum_pts, fnum_ele, fmeshnode, fmeshele);
	CalcForce(fnum_pts, fnum_ele, fmeshnode, fmeshele);
	/*free the space*/
	free(cmeshele); free(cmeshnode);
	free(fmeshele); free(fmeshnode);
}
void MainWindow::TLM(int & num_pts, int & num_ele, CNode *m_n, CElement * m_l) {
	QCustomPlot *customPlot;
	customPlot = ui->widget;
	double U0 = 1;// PI*4e-7;
	std::vector <int> D34;
	D34.empty();
	for (int i = 0; i < num_ele; i++) {
		if (m_l[i].domain == 3 || m_l[i].domain == 4) {
			D34.push_back(i);
		}
	}
	
	long long T1, T2;
	timeval tv1;
	//gettimeofday(&tv1,NULL);
	T1 = tv1.tv_sec * 1000 * 1000 + tv1.tv_usec;
	//------------build C Matrix-----------------------------
	umat locs(2, 9 * num_ele); locs.zeros();
	vec vals = zeros<vec>(9 * num_ele);
	double ce[3][3];
	double cc = 20;
	for (int i = 0; i < num_ele; i++) {
		if (m_l[i].domain == 3 || m_l[i].domain == 4) {
			if (m_l[i].y10 < 0) {
				ce[0][0] = 0;// -abs(m_l[i].y10) / m_l[i].miut;
			} else {
				ce[0][0] = abs(m_l[i].y10) / m_l[i].miut;
			}
			if (m_l[i].y20 < 0) {
				ce[1][1] = 0;//-abs(m_l[i].y20) / m_l[i].miut;
			} else {
				ce[1][1] = abs(m_l[i].y20) / m_l[i].miut;
			}
			if (m_l[i].y30 < 0) {
				ce[2][2] = 0;//-abs(m_l[i].y30) / m_l[i].miut;
			} else {
				ce[2][2] = abs(m_l[i].y30) / m_l[i].miut;
			}
			if (m_l[i].y12 > 0) {
				ce[0][0] += 0;//-abs(m_l[i].y12) / m_l[i].miut;
				ce[1][1] += 0;//-abs(m_l[i].y12) / m_l[i].miut;
				ce[0][1] = 0;//abs(m_l[i].y12) / m_l[i].miut;
			} else {
				ce[0][0] += abs(m_l[i].y12) / m_l[i].miut;
				ce[1][1] += abs(m_l[i].y12) / m_l[i].miut;
				ce[0][1] = -abs(m_l[i].y12) / m_l[i].miut;
			}
			if (m_l[i].y31 > 0) {
				ce[0][0] += 0;//-abs(m_l[i].y31) / m_l[i].miut;
				ce[2][2] += 0;//-abs(m_l[i].y31) / m_l[i].miut;
				ce[0][2] = 0;//abs(m_l[i].y31) / m_l[i].miut;
			} else {
				ce[0][0] += abs(m_l[i].y31) / m_l[i].miut;
				ce[2][2] += abs(m_l[i].y31) / m_l[i].miut;
				ce[0][2] = -abs(m_l[i].y31) / m_l[i].miut;
			}
			if (m_l[i].y23 > 0) {
				ce[1][1] += 0;//-abs(m_l[i].y23) / m_l[i].miut;
				ce[2][2] += 0;//-abs(m_l[i].y23) / m_l[i].miut;
				ce[1][2] = 0;//abs(m_l[i].y23) / m_l[i].miut;
			} else {
				ce[1][1] += abs(m_l[i].y23) / m_l[i].miut;
				ce[2][2] += abs(m_l[i].y23) / m_l[i].miut;
				ce[1][2] = -abs(m_l[i].y23) / m_l[i].miut;
			}
			ce[1][0] = ce[0][1];
			ce[2][0] = ce[0][2];
			ce[2][1] = ce[1][2];
		} else {
			ce[0][0] = m_l[i].y11 / m_l[i].miu;
			ce[0][1] = m_l[i].y12 / m_l[i].miu;
			ce[0][2] = m_l[i].y31 / m_l[i].miu;

			ce[1][0] = ce[0][1];
			ce[1][1] = m_l[i].y22 / m_l[i].miu;
			ce[1][2] = m_l[i].y23 / m_l[i].miu;

			ce[2][0] = ce[0][2];
			ce[2][1] = ce[1][2];
			ce[2][2] = m_l[i].y33 / m_l[i].miu;
		}

		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				locs(0, i * 9 + row * 3 + col) = m_l[i].n[row];
				locs(1, i * 9 + row * 3 + col) = m_l[i].n[col];
				vals(i * 9 + row * 3 + col) = ce[row][col];
			}
		}
	}

	//locs.save("locs1.txt", arma_ascii);
	//vals.save("vals1.txt", arma_ascii);
	//----using armadillo constructor function-----
	sp_mat X(true, locs, vals, num_pts, num_pts, true, true);

	vec bbJz = zeros<vec>(num_pts);	vec b = zeros<vec>(num_pts);//b = bbJz + INL;
	vec A = zeros<vec>(num_pts);	vec A_old = A;
	vec INL = zeros<vec>(num_pts); vec bb = zeros<vec>(num_pts);
	vec rpm = zeros<vec>(num_pts);
	for (int i = 0; i < num_pts; i++) {
		bbJz(i) = m_n[i].I;//set the current
		rpm(i) = m_n[i].pm;
	}
	//---------------------superLU_MT---------------------------------------
	SuperMatrix sluA, L, U;
	SuperMatrix sluB, sluX;
	NCformat    *Astore;
	SCPformat   *Lstore; DNformat	   *Bstore;
	NCPformat   *Ustore;
	int_t         nprocs;
	fact_t      fact;
	trans_t     trans;
	yes_no_t    refact, usepr;
	equed_t     equed;
	double      *a; double t1;
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
	fact = EQUILIBRATE;
	trans = NOTRANS;
	equed = NOEQUIL;
	refact = NO;
	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	u = 1.0;
	usepr = NO;
	drop_tol = 0.0;
	lwork = 0;
	nrhs = 1;

	if (lwork > 0) {
		work = SUPERLU_MALLOC(lwork);
		printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
		if (!work) {
			SUPERLU_ABORT("DLINSOLX: cannot allocate work[]");
		}
	}
	/* create matrix A in Harwell-Boeing format.*/
	m = num_pts; n = num_pts; nnz = X.n_nonzero;
	a = const_cast<double *>(X.values);
	asub = (int*)const_cast<u32 *>(X.row_indices);
	xa = (int*)const_cast<u32 *>(X.col_ptrs);
	dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

	//------create B and X-------------------
	if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&sluX, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
	b = bbJz + INL + rpm;
	rhsb = const_cast<double*>(b.mem);
	dCreate_Dense_Matrix(&sluB, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
	Bstore = (DNformat*)sluB.Store;
	double *Bmat = (double*)Bstore->nzval;
	double *Xmat = (double*)((DNformat*)sluX.Store)->nzval;
	ldx = m;
	if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
	if (!(R = (double *)SUPERLU_MALLOC(sluA.nrow * sizeof(double))))
		SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
	if (!(C = (double *)SUPERLU_MALLOC(sluA.ncol * sizeof(double))))
		SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
	if (!(ferr = (double *)SUPERLU_MALLOC(nrhs * sizeof(double))))
		SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
	if (!(berr = (double *)SUPERLU_MALLOC(nrhs * sizeof(double))))
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
	if (!(superlumt_options.etree = intMalloc(n)))
		SUPERLU_ABORT("Malloc fails for etree[].");
	if (!(superlumt_options.colcnt_h = intMalloc(n)))
		SUPERLU_ABORT("Malloc fails for colcnt_h[].");
	if (!(superlumt_options.part_super_h = intMalloc(n)))
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
	if (info == 0 || info == n + 1) {
		sol = (double*)((DNformat*)sluX.Store)->nzval;
		for (int i = 0; i < num_pts; i++) {
			m_n[i].A = sol[i] * miu0;
			A(i) = sol[i] * miu0;
		}
	} else if (info > 0 && lwork == -1) {

	}
	//---------------------superLU--end----------------------------------
	QVector<double> x1(D34.size()), y1(D34.size());
	QVector<double> x2(num_pts), y2(num_pts);
	for (int i = 0; i < D34.size(); ++i) {
		x1[i] = i;
		x2[i] = i;
	}
	//----------QCustomPlot setting---------------------------
	QCPGraph *graph1 = customPlot->addGraph();
	graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 3));
	graph1->setPen(QPen(QColor(120, 120, 120), 2));
	graph1->setLineStyle(QCPGraph::lsNone);
	QCPGraph *graph2 = customPlot->addGraph();
	graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 1.5), QBrush(Qt::red), 3));
	graph2->setPen(QPen(QColor(120, 120, 120), 2));
	graph2->setLineStyle(QCPGraph::lsNone);
	customPlot->xAxis->setLabel("x");
	customPlot->yAxis->setLabel("y");
	customPlot->rescaleAxes(true);
	customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	//---------openMP setting--------------------------------
	int num_threads = omp_get_num_procs();
	omp_set_num_threads(nprocs);
	//---------the main loop---------------------------------
	double Precision = 1e-6;
	int steps = 250;
	clock_t t2 = clock();

	int rowequ, colequ, notran;
	int ldb = Bstore->lda;
	Gstat_t   Gstat;
	int count;
	//int i,j;
	for (count = 0; count < steps; count++) {
		//#pragma omp parallel for
		//------update miu----------------
		for (int i = 0; i < D34.size(); i++) {
			//meshelement = &m_l[D34[i]];
			int d34 = D34[i];
			m_l[d34].Bx = (m_l[d34].Q[0] * m_n[m_l[d34].n[0]].A
				+ m_l[d34].Q[1] * m_n[m_l[d34].n[1]].A
				+ m_l[d34].Q[2] * m_n[m_l[d34].n[2]].A) / 2. / m_l[d34].AREA;
			m_l[d34].By = (m_l[d34].P[0] * m_n[m_l[d34].n[0]].A
				+ m_l[d34].P[1] * m_n[m_l[d34].n[1]].A
				+ m_l[d34].P[2] * m_n[m_l[d34].n[2]].A) / 2. / m_l[d34].AREA;
			m_l[d34].By += (m_n[m_l[d34].n[0]].A
				+ m_n[m_l[d34].n[1]].A
				+ m_n[m_l[d34].n[2]].A) / 3. / m_l[d34].rc;
			m_l[d34].B = sqrt(m_l[d34].By*m_l[d34].By +
				m_l[d34].Bx*m_l[d34].Bx);
			m_l[d34].miu = HB(m_l[d34].B);// *miu0;
			y1[i] = m_l[d34].B;
			//y2[i] = m_l[D34[i]].miut / r;
		}
		//#pragma omp parallel for
		for (int j = 0; j < D34.size(); j++) {
			int i = D34[j];
			CElement *meshelement = &m_l[i];
			double rtmp;//do this to mark it as private
			rtmp = (meshelement->miu - meshelement->miut) / (meshelement->miu + meshelement->miut);

			meshelement->vr12 = (m_n[meshelement->n[0]].A - m_n[meshelement->n[1]].A) - meshelement->vi12;
			meshelement->vr23 = (m_n[meshelement->n[1]].A - m_n[meshelement->n[2]].A) - meshelement->vi23;
			meshelement->vr31 = (m_n[meshelement->n[2]].A - m_n[meshelement->n[0]].A) - meshelement->vi31;

			meshelement->vr10 = (m_n[meshelement->n[0]].A - 0) - meshelement->vi10;
			meshelement->vr20 = (m_n[meshelement->n[1]].A - 0) - meshelement->vi20;
			meshelement->vr30 = (m_n[meshelement->n[2]].A - 0) - meshelement->vi30;

			if (meshelement->y10 > 0) {//the conductor
				meshelement->vi10 = rtmp*meshelement->vr10;
				INL(m_l[i].n[0]) += 2. / m_l[i].miut*m_l[i].vi10*abs(m_l[i].y10) / miu0;
			} else {//the controlled current source
				meshelement->vi10 = (m_n[meshelement->n[0]].A - 0);//rtmp*meshelement->vr10;//
				INL(m_l[i].n[0]) += 1. / m_l[i].miu*m_l[i].vi10*abs(m_l[i].y10) / miu0;
			}
			if (meshelement->y20 > 0) {
				meshelement->vi20 = rtmp*meshelement->vr20;
				INL(m_l[i].n[1]) += 2. / m_l[i].miut*m_l[i].vi20*abs(m_l[i].y20) / miu0;
			} else {
				meshelement->vi20 = (m_n[meshelement->n[1]].A - 0);//rtmp*meshelement->vr20;
				INL(m_l[i].n[1]) += 1. / m_l[i].miu*m_l[i].vi20*abs(m_l[i].y20) / miu0;
			}
			if (meshelement->y30 > 0) {
				meshelement->vi30 = rtmp*meshelement->vr30;
				INL(m_l[i].n[2]) += 2. / m_l[i].miut*m_l[i].vi30*abs(m_l[i].y30) / miu0;
			} else {
				meshelement->vi30 = (m_n[meshelement->n[2]].A - 0);//rtmp*meshelement->vr30;
				INL(m_l[i].n[2]) += 1. / m_l[i].miu*m_l[i].vi30*abs(m_l[i].y30) / miu0;
			}
			if (meshelement->y12 < 0) {
				meshelement->vi12 = rtmp*meshelement->vr12;
				INL(m_l[i].n[1]) += -2. / m_l[i].miut*m_l[i].vi12*abs(m_l[i].y12) / miu0;
				INL(m_l[i].n[0]) += 2. / m_l[i].miut* m_l[i].vi12 *abs(m_l[i].y12) / miu0;
			} else {
				meshelement->vi12 = (m_n[meshelement->n[0]].A - m_n[meshelement->n[1]].A);//rtmp*meshelement->vr12;
				INL(m_l[i].n[1]) += -1. / m_l[i].miu*m_l[i].vi12*abs(m_l[i].y12) / miu0;
				INL(m_l[i].n[0]) += 1. / m_l[i].miu* m_l[i].vi12 *abs(m_l[i].y12) / miu0;
			}
			if (meshelement->y23 < 0) {
				meshelement->vi23 = rtmp*meshelement->vr23;
				INL(m_l[i].n[1]) += 2. / m_l[i].miut* m_l[i].vi23*abs(m_l[i].y23) / miu0;
				INL(m_l[i].n[2]) += -2. / m_l[i].miut*m_l[i].vi23*abs(m_l[i].y23) / miu0;
			} else {
				meshelement->vi23 = (m_n[meshelement->n[1]].A - m_n[meshelement->n[2]].A);//rtmp*meshelement->vr23;
				INL(m_l[i].n[1]) += 1. / m_l[i].miu* m_l[i].vi23*abs(m_l[i].y23) / miu0;
				INL(m_l[i].n[2]) += -1. / m_l[i].miu*m_l[i].vi23*abs(m_l[i].y23) / miu0;
			}
			if (meshelement->y31 < 0) {
				meshelement->vi31 = rtmp*meshelement->vr31;
				INL(m_l[i].n[2]) += 2. / m_l[i].miut* m_l[i].vi31*abs(m_l[i].y31) / miu0;
				INL(m_l[i].n[0]) += -2.0 / m_l[i].miut*m_l[i].vi31*abs(m_l[i].y31) / miu0;
			} else {
				meshelement->vi31 = (m_n[meshelement->n[2]].A - m_n[meshelement->n[0]].A);//rtmp*meshelement->vr31;
				INL(m_l[i].n[2]) += 1. / m_l[i].miu* m_l[i].vi31*abs(m_l[i].y31) / miu0;
				INL(m_l[i].n[0]) += -1.0 / m_l[i].miu*m_l[i].vi31*abs(m_l[i].y31) / miu0;
			}
			

		}
		b = bbJz + INL + rpm; INL.save("INL.txt", arma_ascii);
		A_old = A;

		//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
		Bstore->nzval = const_cast<double*>(b.mem);
		superlumt_options.fact = FACTORED; /* Indicate the factored form of sluA is supplied. */
		superlumt_options.usepr = YES;
		//get_perm_c(permc_spec, &sluA, perm_c);
		t1 = SuperLU_timer_();
		pdgssvx(nprocs, &superlumt_options, &sluA, perm_c, perm_r,
			&equed, R, C, &L, &U, &sluB, &sluX, &rpg, &rcond,
			ferr, berr, &superlu_memusage, &info);
		t1 = SuperLU_timer_() - t1;

		double error = 0;
		for (int i = 0; i < num_pts; i++) {
			m_n[i].A = sol[i] * miu0;
			A(i) = sol[i] * miu0;
		}

		error = norm((A_old - A), 2) / norm(A, 2);
		if (error < Precision) {
			break;
		}
		INL.zeros();

		//this->setWindowTitle(QString::number(t1));
		graph1->setData(x1, y1);

		customPlot->rescaleAxes(true);
		customPlot->replot();//necessary
		if (QString::number(error) == "nan") {
			return;
		} else {
			this->setWindowTitle(QString::number(count));
		}
	}
	for (int i = 0; i < num_ele; i++) {
		m_l[i].Bx = (m_l[i].Q[0] * m_n[m_l[i].n[0]].A
			+ m_l[i].Q[1] * m_n[m_l[i].n[1]].A
			+ m_l[i].Q[2] * m_n[m_l[i].n[2]].A) / 2. / m_l[i].AREA;
		m_l[i].By = (m_l[i].P[0] * m_n[m_l[i].n[0]].A
			+ m_l[i].P[1] * m_n[m_l[i].n[1]].A
			+ m_l[i].P[2] * m_n[m_l[i].n[2]].A) / 2. / m_l[i].AREA;
		m_l[i].By += (m_n[m_l[i].n[0]].A
			+ m_n[m_l[i].n[1]].A
			+ m_n[m_l[i].n[2]].A) / 3. / m_l[i].rc;
		m_l[i].B = sqrt(m_l[i].By*m_l[i].By +
			m_l[i].Bx*m_l[i].Bx);
	}
	//A.save("A.txt", arma_ascii);
	/*using another loop*/
	//gettimeofday(&tv1,NULL);
	T2 = tv1.tv_sec * 1000 * 1000 + tv1.tv_usec;
	//this->setWindowTitle(QString::number(T2 - T1));
	this->setWindowTitle(QString::number(count) + "," + QString::number(T2 - T1));
	graph1->setData(x1, y1);
	//graph2->setData(x2, y2);
	customPlot->rescaleAxes(true);
	customPlot->replot();//necessary
	//customPlot->removeGraph(graph1);
	//customPlot->removeGraph(graph2);
	/*free the space allocated by superLU*/
	//SUPERLU_FREE(rhsb);
	SUPERLU_FREE(rhsx);
	//SUPERLU_FREE (xact);
	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	SUPERLU_FREE(R);
	SUPERLU_FREE(C);
	SUPERLU_FREE(ferr);
	SUPERLU_FREE(berr);
	//Destroy_CompCol_Matrix(&sluA);
	//Destroy_SuperMatrix_Store(&sluB);
	//Destroy_SuperMatrix_Store(&sluX);
	SUPERLU_FREE(superlumt_options.etree);
	SUPERLU_FREE(superlumt_options.colcnt_h);
	SUPERLU_FREE(superlumt_options.part_super_h);
	if (lwork == 0) {
		Destroy_SuperNode_SCP(&L);
		Destroy_CompCol_NCP(&U);
	} else if (lwork > 0) {
		SUPERLU_FREE(work);
	}
	//customPlot->removeGraph(graph1);
}
//return miur
double MainWindow::HB(double B) {
	double miu;
	double hint;
	double h[10] = { 0, 200, 500, 1000, 2500, 5000, 10000, 15000, 100000, 500000 };
	double b[10] = { 0, 1.2, 1.4, 1.5, 1.62, 1.71, 1.85, 1.851, 1.8511, 5 };
	std::vector<double> BB(b, b + 10);
	std::vector<double> HH(h, h + 10);
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

void MainWindow::readDataFile(char *fileName, int & num_pts, int & num_ele, CNode **ppmeshnode, CElement ** ppmeshele) {
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
	fp = fopen(fileName, "r");
	if (fp == NULL) {
		QMessageBox::warning(NULL, "Error:", "Error: opening file!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	//--------------Read the head-----------------------------
	for (int i = 0; i < 18; i++) {
		fgets(ch, 256, fp);
	}
	//-----------------mesh point-----------------------------
	CNode* meshnode = NULL;
	re = fscanf(fp, "%d # number of mesh points\n", &num_pts);
	if (re == 1) {
		*ppmeshnode = (CNode*)calloc(num_pts, sizeof(CNode));
		meshnode = *ppmeshnode;
		for (int i = 0; i < num_pts; i++) {
			meshnode[i].I = 0;
			meshnode[i].pm = 0;
		}
	} else {
		QMessageBox::warning(NULL, "Error:", "Error: reading num_pts!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	int pts_ind;//the beginning of the points index
	re = fscanf(fp, "%d # lowest mesh point index\n", &pts_ind);
	if (re != 1) {
		QMessageBox::warning(NULL, "Error:", "Error: reading pts_ind!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	fgets(ch, 256, fp);
	for (int i = pts_ind; i < num_pts; i++) {
		re = fscanf(fp, "%lf %lf \n", &(meshnode[i].x), &(meshnode[i].y));
		if (re != 2) {
			QMessageBox::warning(NULL, "Error:", "Error: reading mesh point!",
				QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
	}
	//---------------vertex-------------------------------
	for (int i = 0; i < 7; i++)
		fgets(ch, 256, fp);
	int num_vtx_ns, num_vtx_ele;
	re = fscanf(fp, "%d # number of nodes per element\n", &num_vtx_ns);
	if (re != 1) {
		QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ns!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}

	re = fscanf(fp, "%d # number of elements\n", &num_vtx_ele);
	if (re != 1) {
		QMessageBox::warning(NULL, "Error:", "Error: reading num_vtx_ele!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	fgets(ch, 256, fp);

	int *vtx;
	vtx = (int*)calloc(num_vtx_ele, sizeof(int));
	for (int i = 0; i < num_vtx_ele; i++) {
		re = fscanf(fp, "%d \n", vtx + i);
		if (re != 1) {
			QMessageBox::warning(NULL, "Error:", "Error: reading vertex condition!",
				QMessageBox::Ok, QMessageBox::Ok);
			return;
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
			QMessageBox::warning(NULL, "Error:", "Error: reading vertex condition!",
				QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
	}
	if (vtx2 != NULL) free(vtx2); vtx2 = NULL;
	//--------------boundary--------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int num_bdr_ns, num_bdr_ele;//number of nodes per element;number of elements
	re = fscanf(fp, "%d # number of nodes per element\n", &num_bdr_ns);
	if (re != 1) {
		QMessageBox::warning(NULL, "Error:", "Error: reading num_bdr_ns!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}

	re = fscanf(fp, "%d # number of elements\n", &num_bdr_ele);
	if (re != 1) {
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
			//-----process the A=0 boundary--------------------------
			/*if (abs(meshnode[p1[i]].length() - 0.05) < 5e-3 || abs(meshnode[p1[i]].x) < 5e-5) {
				meshnode[p1[i]].bdr = 1;
				}*/
		} else {
			QMessageBox::warning(NULL, "Error:", "Error: reading boundary condition!",
				QMessageBox::Ok, QMessageBox::Ok);
			return;
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
			QMessageBox::warning(NULL, "Error:", "Error: reading boundary condition!",
				QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
	}
	if (entity != NULL) free(entity); entity = NULL;
	//----------------elements------------------------------
	for (int i = 0; i < 5; i++)
		fgets(ch, 256, fp);
	int ns_per_ele;// num_ele;//number of nodes per element;number of elements
	re = fscanf(fp, "%d # number of nodes per element\n", &ns_per_ele);
	if (re != 1) {
		QMessageBox::warning(NULL, "Error:", "Error: reading ns_per_ele!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}

	CElement* meshele = NULL;
	re = fscanf(fp, "%d # number of elements\n", &num_ele);
	if (re == 1) {
		*ppmeshele = (CElement*)calloc(num_ele, sizeof(CElement));
		meshele = *ppmeshele;
	} else {
		QMessageBox::warning(NULL, "Error:", "Error: reading num_ele!",
			QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	fgets(ch, 256, fp);

	for (int i = 0; i < num_ele; i++) {
		re = fscanf(fp, "%d %d %d \n", &meshele[i].n[0], &meshele[i].n[1], &meshele[i].n[2]);
		if (re == 3) {
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
			meshele[i].zc = (meshnode[meshele[i].n[0]].y +
				meshnode[meshele[i].n[1]].y +
				meshnode[meshele[i].n[2]].y) / 3;

			meshele[i].y12 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[0] * meshele[i].P[1] +
				meshele[i].Q[0] * meshele[i].Q[1]);
			meshele[i].y12 += (meshele[i].P[0] + meshele[i].P[1]) / 6;
			meshele[i].y12 += meshele[i].AREA / 9 / meshele[i].rc;
			meshele[i].y12 *= 2 * PI;

			meshele[i].y23 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[1] * meshele[i].P[2] +
				meshele[i].Q[1] * meshele[i].Q[2]);
			meshele[i].y23 += (meshele[i].P[1] + meshele[i].P[2]) / 6;
			meshele[i].y23 += meshele[i].AREA / 9 / meshele[i].rc;
			meshele[i].y23 *= 2 * PI;

			meshele[i].y31 = meshele[i].rc / 4 / meshele[i].AREA*(meshele[i].P[2] * meshele[i].P[0] +
				meshele[i].Q[2] * meshele[i].Q[0]);
			meshele[i].y31 += (meshele[i].P[2] + meshele[i].P[0]) / 6;
			meshele[i].y31 += meshele[i].AREA / 9 / meshele[i].rc;
			meshele[i].y31 *= 2 * PI;

			meshele[i].y10 = 3 * meshele[i].AREA / 9 / meshele[i].rc;
			meshele[i].y10 += (4 * meshele[i].P[0] + meshele[i].P[1] + meshele[i].P[2]) / 6;
			meshele[i].y10 *= 2 * PI;
			meshele[i].y20 = 3 * meshele[i].AREA / 9 / meshele[i].rc;
			meshele[i].y20 += (meshele[i].P[0] + 4 * meshele[i].P[1] + meshele[i].P[2]) / 6;
			meshele[i].y20 *= 2 * PI;
			meshele[i].y30 = 3 * meshele[i].AREA / 9 / meshele[i].rc;
			meshele[i].y30 += (meshele[i].P[0] + meshele[i].P[1] + 4 * meshele[i].P[2]) / 6;
			meshele[i].y30 *= 2 * PI;


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

			meshele[i].vi10 = 0; meshele[i].vi20 = 0; meshele[i].vi30 = 0;
			//meshele[i].vi40 = 0; meshele[i].vi50 = 0; meshele[i].vi60 = 0;
			meshele[i].vi12 = 0; meshele[i].vi23 = 0; meshele[i].vi31 = 0;
		} else {
			QMessageBox::warning(NULL, "Error:", "Error: reading elements points!",
				QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
	}
	//---------------Domain----------------------------------
	int num_domain;
	re = fscanf(fp, "%d # number of geometric entity indices\n", &num_domain);
	fgets(ch, 256, fp);

	for (int i = 0; i < num_domain; i++) {
		re = fscanf(fp, "%d \n", &meshele[i].domain);
		if (re == 1) {
			switch (meshele[i].domain) {
			case INF:
			case AIR:
				meshele[i].miu = U0;
				meshele[i].miut = U0;
				meshnode[meshele[i].n[0]].I += 0;
				meshnode[meshele[i].n[1]].I += 0;
				meshnode[meshele[i].n[2]].I += 0;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
				break;
			case DOWNCOIL:
				meshele[i].miu = U0;
				meshele[i].miut = U0;
				meshnode[meshele[i].n[0]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*Jdown;
				meshnode[meshele[i].n[1]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*Jdown;
				meshnode[meshele[i].n[2]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*Jdown;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
				break;
			case UPCOIL:
				meshele[i].miu = U0;
				meshele[i].miut = U0;
				meshnode[meshele[i].n[0]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*J;
				meshnode[meshele[i].n[1]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*J;
				meshnode[meshele[i].n[2]].I += 2. / 3 * PI*meshele[i].rc*meshele[i].AREA*J;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
				break;
			case FIX:
			case MOV:
				meshele[i].miu = 1.0 * U0;
				meshele[i].miut = (nonlinear_miu)* U0;
				meshnode[meshele[i].n[0]].I += 0;
				meshnode[meshele[i].n[1]].I += 0;
				meshnode[meshele[i].n[2]].I += 0;
				meshnode[meshele[i].n[0]].pm += 0;
				meshnode[meshele[i].n[1]].pm += 0;
				meshnode[meshele[i].n[2]].pm += 0;
				break;
			case PM:
				int k;
				double r0, K, be[3];
				be[0] = 0; be[1] = 0; be[2] = 0;
				meshele[i].miu = U0; meshele[i].miut = U0;
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
		} else {
			QMessageBox::warning(NULL, "Error:", "Error: reading domain points!",
				QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
	}
	fclose(fp);
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
	customPlot = ui->widget;
	customPlot->xAxis->setLabel("x");
	customPlot->xAxis->setRange(0, 0.09);
	customPlot->xAxis->setAutoTickStep(false);
	customPlot->xAxis->setTicks(false);
	customPlot->yAxis->setLabel("y");
	customPlot->yAxis->setRange(-0.09, 0.09);
	customPlot->xAxis2->setTicks(false);
	customPlot->yAxis->setScaleRatio(ui->widget->xAxis, 1.0);

	//customPlot->yAxis->setTickStep(ui->widget->xAxis->tickStep());
	customPlot->yAxis->setAutoTickStep(false);
	ui->widget->yAxis->setAutoTickLabels(false);
	customPlot->yAxis->setTicks(false);
	customPlot->yAxis->grid()->setVisible(false);
	customPlot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
	customPlot->yAxis2->setTicks(false);
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
	/*cc[0] = QColor(0,0,0);
	cc[1] = QColor(0, 120, 0);
	cc[2] = QColor(0, 0, 120);
	cc[3] = QColor(120, 0, 0);
	cc[4] = QColor(120, 120, 0);
	cc[5] = QColor(0, 120, 120);
	cc[6] = QColor(120, 0, 120);*/
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
			x1[j] = m_n[m_l[i].n[j]].x;
			y1[j] = m_n[m_l[i].n[j]].y;
		}
		x1[3] = x1[0];
		y1[3] = y1[0];
		if (m_l[i].domain == AIR) {
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
			double tmp = PI / 2 * PI * m_l[i].rc * m_l[i].AREA / (4 * PI*1e-7);
			A1 = m_n[m_l[i].n[0]].A;
			A2 = m_n[m_l[i].n[1]].A;
			A3 = m_n[m_l[i].n[2]].A;
			Ac = (A1 + A2 + A3) / 3;
			if (ind != 10) {
				ind = index[n];
				xf1 = m_l[i].B * m_l[i].B / m_l[i].AREA;
				xf1 *= m_l[i].Q[ind] / 2.0;

				beta2 = m_l[i].Q[0] * A1 +
					m_l[i].Q[1] * A2 +
					m_l[i].Q[2] * A3;
				beta2 /= (2. * m_l[i].AREA);
				xf2 = -(2. * beta2 + 2.0 * Ac / m_l[i].rc)*beta2 / (2 * m_l[i].AREA)*(m_l[i].Q[ind]);

				beta3 = m_l[i].P[0] * A1 +
					m_l[i].P[1] * A2 +
					m_l[i].P[2] * A3;
				beta3 /= (2. * m_l[i].AREA);
				if (ind == 0) {
					xf3 = 2. * beta3*((A3 - A2)*(2. * m_l[i].AREA) - beta3*(2. * m_l[i].AREA)*(m_l[i].Q[0]));
					xf3 /= (2. * m_l[i].AREA)*(2. * m_l[i].AREA);
				} else if (ind == 1) {
					xf3 = 2. * beta3*((A1 - A3)*(2. * m_l[i].AREA) - beta3*(2. * m_l[i].AREA)*(m_l[i].Q[1]));
					xf3 /= (2. * m_l[i].AREA)*(2. * m_l[i].AREA);
				} else if (ind == 2) {
					xf3 = 2. * beta3*((A2 - A1)*(2. * m_l[i].AREA) - beta3*(2. * m_l[i].AREA)*(m_l[i].Q[2]));
					xf3 /= (2. * m_l[i].AREA)*(2. * m_l[i].AREA);
				}
				yForce -= (tmp*(xf1 + xf2 + xf3));
				//fprintf(fp, "%d\t%lf\t%lf\n", i, m_l[i].B, tmp*(xf1 + xf2 + xf3));
			}
		}

		//if (m_l[i].domain != AIR && m_l[i].domain != INF) {
		newCurve->setData(x1, y1);
		newCurve->setPen(QPen(cc[m_l[i].domain - 1]));
		//newCurve->setBrush(cc[m_l[i].domain - 1]);
		//}		
		//if (ind != 10) {
		//	fprintf(fp, "%d\n", ind);			
		//} 

	}
	//graph1->setData(gpx1, gpy1);
	//graph2->setData(gpx2, gpy2);
	fclose(fp);
	this->setWindowTitle(QString::number(yForce));
	customPlot->replot();
}


