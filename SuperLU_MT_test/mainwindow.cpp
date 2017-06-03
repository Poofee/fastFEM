#include <QMenu>
#include <QVector>
#include <QAction>
#include <QMenuBar>
#include <QtAlgorithms>
#include <stdio.h>
#include  <math.h>
#include "mainwindow.h"
#include "ui_mainwindow.h"
/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
* -- SuperLU MT routine (version 3.0) --
* Lawrence Berkeley National Lab, Univ. of California Berkeley,
* and Xerox Palo Alto Research Center.
* September 10, 2007
*
*/
#include "slu_mt_ddefs.h"

typedef struct{
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

/* Eat up the rest of the current line */
int_t dDumpLine(FILE *fp) {
    register int_t c;
    while ((c = fgetc(fp)) != '\n');
    return 0;
}

int_t dParseIntFormat(char *buf, int_t *num, int_t *size) {
    char *tmp;

    tmp = buf;
    while (*tmp++ != '(');
    *num = atoi(tmp);
    while (*tmp != 'I' && *tmp != 'i') ++tmp;
    ++tmp;
    *size = atoi(tmp);
    return 0;
}

int_t dParseFloatFormat(char *buf, int_t *num, int_t *size) {
    char *tmp, *period;

    tmp = buf;
    while (*tmp++ != '(');
    *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
    while (*tmp != 'E' && *tmp != 'e' && *tmp != 'D' && *tmp != 'd'
           && *tmp != 'F' && *tmp != 'f') {
        /* May find kP before nE/nD/nF, like (1P6F13.6). In this case the
        num picked up refers to P, which should be skipped. */
        if (*tmp == 'p' || *tmp == 'P') {
            ++tmp;
            *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
        } else {
            ++tmp;
        }
    }
    ++tmp;
    period = tmp;
    while (*period != '.' && *period != ')') ++period;
    *period = '\0';
    *size = atoi(tmp); /*sscanf(tmp, "%2d", size);*/

    return 0;
}

int_t dReadVector(FILE *fp, int_t n, int_t *where, int_t perline, int_t persize) {
    register int_t i, j, item;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
        fgets(buf, 100, fp);    /* read a line at a time */
        for (j = 0; j<perline && i<n; j++) {
            tmp = buf[(j + 1)*persize];     /* save the char at that place */
            buf[(j + 1)*persize] = 0;       /* null terminate */
            item = atoi(&buf[j*persize]);
            buf[(j + 1)*persize] = tmp;     /* recover the char at that place */
            where[i++] = item - 1;
        }
    }

    return 0;
}

int_t dReadValues(FILE *fp, int_t n, double *destination, int_t perline, int_t persize) {
    register int_t i, j, k, s;
    char tmp, buf[100];

    i = 0;
    while (i < n) {
        fgets(buf, 100, fp);    /* read a line at a time */
        for (j = 0; j<perline && i<n; j++) {
            tmp = buf[(j + 1)*persize];     /* save the char at that place */
            buf[(j + 1)*persize] = 0;       /* null terminate */
            s = j*persize;
            for (k = 0; k < persize; ++k) /* No D_ format in C */
                if (buf[s + k] == 'D' || buf[s + k] == 'd') buf[s + k] = 'E';
            destination[i++] = atof(&buf[s]);
            buf[(j + 1)*persize] = tmp;     /* recover the char at that place */
        }
    }

    return 0;
}


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    QAction * calAction = new QAction(tr("Cal"),this);
    connect(calAction,SIGNAL(triggered()),this,SLOT(cal()));

    QMenu * file = menuBar()->addMenu(tr("File"));
    file->addAction(calAction);


    ui->setupUi(this);
    connect(ui->actionCal,SIGNAL(triggered()),this,SLOT(cal()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::cal(){
    SuperMatrix   A;
    NCformat *Astore;
    double   *a;
    int_t      *asub, *xa;
    int_t      *perm_r; /* row permutations from partial pivoting */
    int_t      *perm_c; /* column permutation vector */
    SuperMatrix   L;       /* factor L */
    SCPformat *Lstore;
    SuperMatrix   U;       /* factor U */
    NCPformat *Ustore;
    SuperMatrix   B;
    int_t      nrhs, ldx, info, m, n, nnz, b;
    int_t      nprocs; /* maximum number of processors to use. */
    int_t      panel_size, relax, maxsup;
    int_t      permc_spec;
    trans_t  trans;
    double   *xact, *rhs;
    superlu_memusage_t   superlu_memusage;
    //void   parse_command_line();

    nrhs = 1;
    trans = NOTRANS;
    nprocs = 1;
    n = 1000;
    b = 1;
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    maxsup = sp_ienv(3);

    //parse_command_line(argc, argv, &nprocs, &n, &b, &panel_size,
    //	&relax, &maxsup);
    //dreadhb(int_t *nrow, int_t *ncol, int_t *nonz,
    //	double **nzval, int_t **rowind, int_t **colptr)
    //dreadhb(&m, &n, &nnz, &a, &asub, &xa);
    register int_t i, numer_lines, rhscrd = 0;
    int_t tmp, colnum, colsize, rownum, rowsize, valnum, valsize;
    char buf[100], type[4], key[10];
    FILE *fp;
    //fp = stdin;
    fp = fopen("..\\..\\SuperLU_MT_test\\g5.rua", "r");
    //fp = fopen("..\\..\\SuperLU_MT_test\\big.rua", "r");
    qDebug()<<fp;
    /* Line 1 */
    fscanf(fp, "%72c", buf); buf[72] = 0;
    qDebug()<<"Title: "<< buf;
    fscanf(fp, "%8c", key);  key[8] = 0;
    qDebug()<<"Key: "<<key;
    dDumpLine(fp);
    /* Line 2 */
    for (i = 0; i<5; i++) {
        fscanf(fp, "%14c", buf); buf[14] = 0;
        tmp = atoi(buf); /*sscanf(buf, "%d", &tmp);*/
        if (i == 3) numer_lines = tmp;
        if (i == 4 && tmp) rhscrd = tmp;
    }
    dDumpLine(fp);

    /* Line 3 */
    fscanf(fp, "%3c", type);
    fscanf(fp, "%11c", buf); /* pad */
    type[3] = 0;
#if ( DEBUGlevel>=1 )
    printf("Matrix type %s\n", type);
#endif

    fscanf(fp, "%14c", buf); m = atoi(buf);
    fscanf(fp, "%14c", buf); n = atoi(buf);
    fscanf(fp, "%14c", buf); nnz = atoi(buf);
    fscanf(fp, "%14c", buf); tmp = atoi(buf);

    if (tmp != 0)
        qDebug()<<"This is not an assembled matrix!";
    if (m != n)
        qDebug()<<"Matrix is not square.\n";
    dDumpLine(fp);

    /* Line 4: format statement */
    fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &colnum, &colsize);
    fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &rownum, &rowsize);
    fscanf(fp, "%20c", buf);
    dParseFloatFormat(buf, &valnum, &valsize);
    fscanf(fp, "%20c", buf);
    dDumpLine(fp);

    /* Line 5: right-hand side */
    if (rhscrd) dDumpLine(fp); /* skip RHSFMT */

#if ( DEBUGlevel>=1 )
    printf("%d rows, %d nonzeros\n", *nrow, *nonz);
    printf("colnum %d, colsize %d\n", colnum, colsize);
    printf("rownum %d, rowsize %d\n", rownum, rowsize);
    printf("valnum %d, valsize %d\n", valnum, valsize);
#endif

    /* Allocate storage for the three arrays ( nzval, rowind, colptr ) */
    dallocateA(n, nnz, &a, &asub, &xa);

    dReadVector(fp, m + 1, xa, colnum, colsize);

    dReadVector(fp, nnz, asub, rownum, rowsize);
    //    for(int i = 0;i <nnz;i++){
    //        qDebug()<<asub[i];
    //    }
    if (numer_lines) {
        dReadValues(fp, nnz, a, valnum, valsize);
    }

    fclose(fp);

    //-----------qt plot----------------------
    QCustomPlot * customplot;
    customplot = ui->qplot;
    QCPGraph *graph1 = customplot->addGraph();
    graph1->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 3), QBrush(Qt::black), 3));
    graph1->setPen(QPen(QColor(120, 120, 120), 2));
    graph1->setLineStyle(QCPGraph::lsNone);
    customplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    customplot->xAxis->setLabel("x");
    customplot->xAxis->setRange(0,26);
    customplot->xAxis->setAutoTickStep(false);
    customplot->xAxis->setAutoTickLabels(false);
    customplot->xAxis->setTickStep(1);
    customplot->xAxis->setTickLabels(false);
    customplot->xAxis->grid()->setVisible(false);
    customplot->xAxis->setTicks(false);
    customplot->yAxis->setLabel("y");
    customplot->yAxis->setRange(0,26);
    customplot->yAxis->setAutoTickStep(false);
    customplot->yAxis->setAutoTickLabels(false);
    customplot->yAxis->setTickStep(1);
    customplot->yAxis->setTickLabels(false);
    customplot->yAxis->setTicks(false);
    customplot->yAxis->grid()->setVisible(false);
    customplot->yAxis->setScaleRatio(customplot->xAxis, 1);

    QVector<double> x1(4), y1(4);
    x1[0] = 0;x1[1] = m-1;x1[2]=m-1;x1[3]=0;
    y1[0] = 0;y1[1] = 0;y1[2]=m-1;y1[3]=m-1;
    QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
    newCurve->setBrush(QColor(255, 0, 0,100));
    //newCurve->setData(x1, y1);

    QVector<double> x(nnz), y(nnz);
    for(int i=0;i < m;i++){//col
        for(int j = xa[i];j < xa[i+1];j++){
            x[j] = i;
            y[j] = m-1-asub[j];//row
            //qDebug()<<"y"<<y[j];
            //qDebug()<<"x"<<x[j];
        }
    }
    //graph1->setData(x, y);
    //customplot->rescaleAxes(true);
    customplot->replot();
    //-----------------------------------------------

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = (NCformat*)A.Store;
    qDebug()<<"Dimension " <<A.nrow<< "x"<< A.ncol<< "; # nonzeros "<< Astore->nnz ;

    if (!(rhs = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhs[].");
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);

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
    get_perm_c(permc_spec, &A, perm_c);

    pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);

    if (info == 0) {
        dinf_norm_error(nrhs, &B, xact); /* Inf. norm of the error */

        Lstore = (SCPformat *)L.Store;
        Ustore = (NCPformat *)U.Store;
        //-----------------------------------------
        //绘制U矩阵
        QCPGraph *graphU = customplot->addGraph();
        graphU->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 1), QBrush(Qt::red), 1));
        graphU->setPen(QPen(QColor(120, 120, 120), 2));
        graphU->setLineStyle(QCPGraph::lsNone);

        //绘制边框
        QVector<double> xU(4), yU(4);
        xU[0] = 0;xU[1] = m-1;xU[2]=m-1;xU[3]=0;
        yU[0] = 0;yU[1] = 0;yU[2]=m-1;yU[3]=m-1;
        QCPCurve *newCurveU = new QCPCurve(customplot->xAxis, customplot->yAxis);
        //newCurveU->setBrush(QColor(255, 0, 0,100));


        //绘制L矩阵
        QCPGraph *graphL = customplot->addGraph();
        graphL->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::black), 1));
        graphL->setPen(QPen(QColor(120, 120, 120), 2));
        graphL->setLineStyle(QCPGraph::lsNone);

        QVector<double> xL(2), yL(2);
        xL[0] = 0+m;xL[1] = m-1+m;
        yL[0] = m-1;yL[1] = 0;
        QCPCurve *newCurveL = new QCPCurve(customplot->xAxis, customplot->yAxis);
        //newCurveL->setBrush(QColor(255, 0, 0,100));
        //newCurveL->setData(xL, yL);

        int nnzL = Lstore->nnz;
        //asub = Lstore->rowind;

        qDebug()<<"L"<<nnzL;
        qDebug()<<"nsuper"<<Lstore->nsuper;

        QVector<double> xnL(nnzL), ynL(nnzL);
        QVector<double> valueL(nnzL);

        int nnzU = Ustore->nnz;
        //asubU = Ustore->rowind;

        QVector<double> xnU(nnzU), ynU(nnzU),valueU(nnzU);

        int fsupc,istart,nsupr,nsupc,nrow;//
        int count = 0;int countU = 0;
        double value;
        QVector<int> level(m),sortlevel(m),sortlevelU(m);
        for(int i = 0;i < m;i++){
            level[i] = 0;
        }
        //对Lstore进行读取,nsuper,超级节点个数，从0开始
        for (int ksupno = 0; ksupno <= Lstore->nsuper; ++ksupno) {
            fsupc = Lstore->sup_to_colbeg[ksupno];//第ksupno个超级节点长方形区域的起始列号
            istart = Lstore->rowind_colbeg[fsupc];//第fsupc列的起始行号，这个列必须是超级节点第一列，
            //超级节点是长方形结构，行号一样,但是某些位置是被零填充的
            nsupr = Lstore->rowind_colend[fsupc] - istart;//第fsupc列的结束行号，算出来是超级节点的列宽
            nsupc = Lstore->sup_to_colend[ksupno] - fsupc;//算出来是超级节点的行高度
            nrow = nsupr - nsupc;//算出来是超级节点不在对角区域的行高度

            if ( nsupc == 1 ) {//超级节点只有一列
                for (int j = 0; j < nrhs; j++) {//对B是多列的情况
                    int luptr = Lstore->nzval_colbeg[fsupc];//非零元素的序号
                    for (int iptr=istart; iptr < Lstore->rowind_colend[fsupc]; iptr++){
                        int irow = Lstore->rowind[iptr];//行号
                        value = ((double*)Lstore->nzval)[luptr];//值

                        if(fabs(value)>1e-9){//非零元素,如果不考虑这个，那个level就没法算了
                            if(irow >= fsupc){//下三角
                                xnL[count] = fsupc;
                                ynL[count] = irow;

                                if(irow == fsupc){//对角线上元素
                                    level[irow] += 1;//计算level=max(level)+1;
                                    valueL[count] = 1;
                                    //xnU[countU] = fsupc;
                                    //ynU[countU] = irow;
                                    //valueU[countU] = value;
                                    //countU++;
                                }else if(irow > fsupc){//非对角线元素
                                    //计算level，取该行的最大值
                                    level[irow] = std::max(level[fsupc],level[irow]);
                                    valueL[count] = value;
                                }
                                count++;
                            }else {
                                //xnU[countU] = fsupc;
                                //ynU[countU] = irow;
                                //valueU[countU] = value;
                                //countU++;
                            }
                        }
                        ++luptr;
                    }
                }
            } else {
                for(int j = 0; j < nrhs;j++){
                    for(int i = 0; i < nsupc;i++){
                        int luptr = Lstore->nzval_colbeg[fsupc+i];
                        for(int iptr = istart; iptr < Lstore->rowind_colend[fsupc];iptr++){
                            int irow = Lstore->rowind[iptr];
                            value = ((double*)Lstore->nzval)[luptr];

                            if(fabs(value)>1e-9){
                                if(irow >= fsupc+i){
                                    xnL[count] = fsupc+i;
                                    ynL[count] = irow;

                                    if(irow == fsupc+i){
                                        level[irow] += 1;
                                        valueL[count] = 1;
                                        //xnU[countU] = fsupc+i;
                                        //ynU[countU] = irow;
                                        //valueU[countU] = value;
                                        //countU++;
                                    }else if(irow > fsupc+i){
                                        level[irow] = std::max(level[fsupc+i],level[irow]);
                                        valueL[count] = value;
                                    }
                                    count++;
                                }else{
                                    //xnU[countU] = fsupc+i;
                                    //ynU[countU] = irow;
                                    //valueU[countU] = value;
                                    //countU++;
                                }
                            }
                            ++luptr;
                        }
                    }
                }
            } /* if-else: nsupc == 1 ... */
        } /* for L-solve */

        //对Ustore进行读取
        QVector <int> levelU(m);levelU.fill(0);
        for (int ksupno = Lstore->nsuper; ksupno >= 0; --ksupno) {
            fsupc = Lstore->sup_to_colbeg[ksupno];//第ksupno个超级节点长方形区域的起始列号
            istart = Lstore->rowind_colbeg[fsupc];//第fsupc列的起始行号，这个列必须是超级节点第一列，
            //超级节点是长方形结构，行号一样,但是某些位置是被零填充的
            nsupr = Lstore->rowind_colend[fsupc] - istart;//第fsupc列的结束行号，算出来是超级节点的列宽
            nsupc = Lstore->sup_to_colend[ksupno] - fsupc;//算出来是超级节点的行高度
            nrow = nsupr - nsupc;//算出来是超级节点不在对角区域的行高度
            int iend = Lstore->rowind_colend[fsupc]-1;

            if ( nsupc == 1 ) {//超级节点只有一列
                for (int j = 0; j < nrhs; j++) {//对B是多列的情况
                    int luptr = Lstore->nzval_colend[fsupc]-1;//非零元素的序号
                    for (int iptr=iend; iptr >= Lstore->rowind_colbeg[fsupc]; iptr--){
                        int irow = Lstore->rowind[iptr];//行号
                        value = ((double*)Lstore->nzval)[luptr];//值

                        if(fabs(value)>1e-9){//非零元素,如果不考虑这个，那个level就没法算了
                            if(irow <= fsupc){//上三角
                                xnU[countU] = fsupc;
                                ynU[countU] = irow;
                                valueU[countU] = value;

                                if(irow == fsupc){//对角线上元素
                                    levelU[irow] += 1;//计算level=max(level)+1;
                                }else if(irow < fsupc){//非对角线元素
                                    //计算level，取该行的最大值
                                    levelU[irow] = std::max(levelU[fsupc],levelU[irow]);
                                }
                                countU++;
                            }
                        }
                        --luptr;
                    }
                }
                //读取Ustore中数据
                for(int u = Ustore->colend[fsupc]-1;u >= Ustore->colbeg[fsupc];u--){
                    double value1 = ((double*)Ustore->nzval)[u];
                    int irow = Ustore->rowind[u];//row
                    if(fabs(value1) > 1e-9){
                        xnU[countU] = fsupc;//col
                        ynU[countU] = irow;
                        valueU[countU] = value1;
                        levelU[irow] = std::max(levelU[fsupc],levelU[irow]);

                        countU++;
                    }
                }
            } else {
                for(int j = 0; j < nrhs;j++){
                    for(int i = nsupc-1; i >= 0;i--){
                        int luptr = Lstore->nzval_colend[fsupc+i]-1;
                        for(int iptr = iend; iptr >= Lstore->rowind_colbeg[fsupc];iptr--){
                            int irow = Lstore->rowind[iptr];
                            value = ((double*)Lstore->nzval)[luptr];

                            if(fabs(value)>1e-9){
                                if(irow <= fsupc+i){
                                    xnU[countU] = fsupc+i;
                                    ynU[countU] = irow;
                                    valueU[countU] = value;

                                    if(irow == fsupc+i){
                                        levelU[irow] += 1;
                                    }else if(irow < fsupc+i){
                                        levelU[irow] = std::max(levelU[fsupc+i],levelU[irow]);
                                    }
                                    countU++;
                                }
                            }
                            --luptr;
                        }
                        //读取该列的Ustore数据
                        for(int u = Ustore->colend[fsupc+i]-1;u >= Ustore->colbeg[fsupc+i];u--){
                            double value1 = ((double*)Ustore->nzval)[u];
                            int irow = Ustore->rowind[u];//row
                            if(fabs(value1) > 1e-9){
                                xnU[countU] = fsupc+i;//col
                                ynU[countU] = irow;
                                valueU[countU] = value1;
                                levelU[irow] = std::max(levelU[fsupc+i],levelU[irow]);

                                countU++;
                            }
                        }
                    }
                }
            } /* if-else: nsupc == 1 ... */
        } /* for U-solve */


        qDebug()<<"countL "<<count;
        qDebug()<<"countU "<<countU;
        int maxLevel = 0;
        int maxLevelU = 0;
        QVector <ele> slevel(m);
        QVector <ele> slevelU(m);
        for(int i = 0;i < m;i++){
            //qDebug()<<i<<" "<<level[i];
            //计算最大的level
            maxLevel = maxLevel > level[i] ? maxLevel : level[i];
            maxLevelU = maxLevelU > levelU[i] ? maxLevelU : levelU[i];
            slevel[i].index = i;
            slevel[i].level = level[i];
            slevelU[i].index = i;
            slevelU[i].level = levelU[i];
        }
        qDebug()<<"maxLevel"<<maxLevel;
        qDebug()<<"maxLevelU"<<maxLevelU;


        QVector <double> xxnL(count),yynL(count);
        QVector <double> xxnU(countU),yynU(countU);
        //对level进行排序
        qSort(slevel.begin(), slevel.end(), compareele);//ascend
        qSort(slevelU.begin(), slevelU.end(), compareele);//ascend
        for(int i = 0;i < m;i++){
            //qDebug()<<slevel[i].index<<"  "<<slevel[i].level;
            sortlevel[slevel[i].index] = i;
            sortlevelU[slevelU[i].index] = i;
            //qDebug()<<perm_c[i];
            //qDebug()<<i<<"\t"<<levelU[i];
        }

        for(int i = 0; i < count;i++){
            xxnL[i] = xnL[i];
            yynL[i] = ynL[i];
            //qDebug()<<"L"<<xnL[i]<<"\t"<<ynL[i]<<"\t"<<valueL[i];
        }
        //转换坐标
        for(int i = 0; i < count;i++){
            xxnL[i] = sortlevel[xnL[i]];
            yynL[i] = m-1-sortlevel[ynL[i]];
        }
        for(int i = 0;i < countU;i++){
            xxnU[i] = m-1-sortlevelU[xnU[i]];
            yynU[i] = sortlevelU[ynU[i]];
            //qDebug()<<"U"<<xnU[i]<<"\t"<<ynU[i]<<"\t"<<valueU[i];
        }
        //绘制level分割线等
        double lastline = 0;
        QPen pen;
        pen.setStyle(Qt::DashLine);
        pen.setWidth(2);
        pen.setColor(QColor(0,0,0));
        QVector <int> levelLrow(maxLevel+1);levelLrow[0] = 0;
        for(int i = 0;i < m-1;i++){
            if(slevel[i].level != slevel[i+1].level){
                levelLrow[slevel[i+1].level] = i+1;
                qDebug()<<i+1;
                QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
                QVector <double> xline(2),yline(2);
                xline[0] = 0; xline[1] = i+0.5;
                yline[0] = m-1-(i+0.5); yline[1] = m-1-(i+0.5);
                newCurve->setData(xline, yline);
                newCurve->setPen(pen);
                //newCurve->pen().setStyle(Qt::DotLine);
                //newCurve->setBrush(QColor(0, 120, 0));
                //--绘制并行区域
                QCPCurve *rec1 = new QCPCurve(customplot->xAxis, customplot->yAxis);
                QCPCurve *rec2 = new QCPCurve(customplot->xAxis, customplot->yAxis);
                QVector <double> recdatax(5),recdatay(5);
                recdatax[0] = lastline; recdatax[1] = i;
                recdatax[2] = i; recdatax[3] = lastline; recdatax[4] = lastline;
                recdatay[0] = m-1-i; recdatay[1] = m-1-i;
                recdatay[2] = m-1-lastline; recdatay[3] = m-1-lastline;recdatay[4] = m-1-i;
                rec1->setData(recdatax, recdatay);
                rec1->setPen(Qt::NoPen);
                rec1->setBrush(QColor(255, 0, 0,100));

                recdatax[0] = lastline; recdatax[1] = i;
                recdatax[2] = i; recdatax[3] = lastline; recdatax[4] = lastline;
                recdatay[0] = 0; recdatay[1] = 0;
                recdatay[2] = m-1-i; recdatay[3] = m-1-i;recdatay[4] = 0;
                if(slevel[i].level%2==1){
                    rec2->setData(recdatax, recdatay);
                    rec2->setPen(Qt::NoPen);
                    rec2->setBrush(QColor(0, 255, 0,100));
                }

                lastline = i+1;
            }
        }
        //绘制U矩阵的分割线
        lastline = 0;
        for(int i = 0;i < m-1;i++){
            if(slevelU[i].level != slevelU[i+1].level){
                QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
                QVector <double> xline(2),yline(2);
                xline[0] = m-1; xline[1] = m-1-i-0.5;
                yline[0] = (i+0.5); yline[1] = (i+0.5);
                //newCurve->setData(xline, yline);
                //newCurve->setPen(QPen(QColor(0, 120, 255),1));
                //newCurve->setBrush(QColor(0, 120, 255));
                //--绘制并行区域
                QCPCurve *rec1 = new QCPCurve(customplot->xAxis, customplot->yAxis);
                QCPCurve *rec2 = new QCPCurve(customplot->xAxis, customplot->yAxis);
                QVector <double> recdatax(5),recdatay(5);
                recdatax[0] = m-1-lastline; recdatax[1] = m-1-i;
                recdatax[2] = m-1-i; recdatax[3] = m-1-lastline; recdatax[4] = m-1-lastline;
                recdatay[0] = lastline; recdatay[1] = lastline;
                recdatay[2] = i; recdatay[3] = i;recdatay[4] = lastline;
                //rec1->setData(recdatax, recdatay);
                //rec1->setPen(QPen(QColor(65, 36, 230),2));

                recdatax[0] = m-1-lastline; recdatax[1] = m-1-i;
                recdatax[2] = m-1-i; recdatax[3] = m-1-lastline; recdatax[4] = m-1-lastline;
                recdatay[0] = i; recdatay[1] = i;
                recdatay[2] = m-1; recdatay[3] = m-1;recdatay[4] = i;
                //rec2->setData(recdatax, recdatay);
                //rec2->setPen(QPen(QColor(65, 36, 230),2));
                //rec1->setBrush(QColor(255, 0, 0));
                lastline = i+1;
            }
        }
        //graphU->setData(xxnU, yynU);
        graphL->setData(xxnL, yynL);
        //newCurveU->setData(xU, yU);

        //-----------------------------------------
        qDebug()<<"#NZ in factor L = "<<Lstore->nnz;
        qDebug()<<"#NZ in factor U = "<<Ustore->nnz;
        qDebug()<<"#NZ in L+U = "<<Lstore->nnz + Ustore->nnz - L.ncol;

        superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
        qDebug()<<"L\\U MB "<<superlu_memusage.for_lu / 1024 / 1024<<"\ttotal MB needed "<<superlu_memusage.total_needed / 1024 / 1024<<"\texpansions "<<superlu_memusage.expansions ;

    }
    customplot->replot();

    SUPERLU_FREE(rhs);
    SUPERLU_FREE(xact);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
}

