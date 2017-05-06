#include <QMenu>
#include <QVector>
#include <QAction>
#include <QMenuBar>
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
    //fp = fopen("..\\..\\SuperLU_MT_test\\g5.rua", "r");
    fp = fopen("..\\..\\SuperLU_MT_test\\big.rua", "r");
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
    customplot->yAxis->setLabel("y");
    customplot->yAxis->setRange(0,26);
    customplot->yAxis->setAutoTickStep(false);
    customplot->yAxis->setAutoTickLabels(false);
    customplot->yAxis->setTickStep(1);
    customplot->yAxis->setScaleRatio(customplot->xAxis, 1);

    QVector<double> x1(4), y1(4);
    x1[0] = 0;x1[1] = m-1;x1[2]=m-1;x1[3]=0;
    y1[0] = 0;y1[1] = 0;y1[2]=m-1;y1[3]=m-1;
    QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
    newCurve->setBrush(QColor(255, 0, 0,100));
    newCurve->setData(x1, y1);

    QVector<double> x(nnz), y(nnz);
    for(int i=0;i < m;i++){//col
        for(int j = xa[i];j < xa[i+1];j++){
            x[j] = i;
            y[j] = m-1-asub[j];//row
            //qDebug()<<"y"<<y[j];
            //qDebug()<<"x"<<x[j];
        }
    }
    graph1->setData(x, y);
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
        QCPGraph *graphU = customplot->addGraph();
        graphU->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::red, 3), QBrush(Qt::red), 3));
        graphU->setPen(QPen(QColor(120, 120, 120), 2));
        graphU->setLineStyle(QCPGraph::lsNone);

        QVector<double> xU(4), yU(4);
        xU[0] = 0+m;xU[1] = m-1+m;xU[2]=m-1+m;xU[3]=0+m;
        yU[0] = 0;yU[1] = 0;yU[2]=m-1;yU[3]=m-1;
        QCPCurve *newCurveU = new QCPCurve(customplot->xAxis, customplot->yAxis);
        newCurveU->setBrush(QColor(255, 0, 0,100));
        newCurveU->setData(xU, yU);

        QCPGraph *graphL = customplot->addGraph();
        graphL->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 3), QBrush(Qt::black), 3));
        graphL->setPen(QPen(QColor(120, 120, 120), 2));
        graphL->setLineStyle(QCPGraph::lsNone);

        QVector<double> xL(2), yL(2);
        xL[0] = 0+m;xL[1] = m-1+m;
        yL[0] = m-1;yL[1] = 0;
        QCPCurve *newCurveL = new QCPCurve(customplot->xAxis, customplot->yAxis);
        newCurveL->setBrush(QColor(255, 0, 0,100));
        newCurveL->setData(xL, yL);

        int nnzL = Lstore->nnz;
        //asub = Lstore->rowind;

        qDebug()<<"L"<<nnzL;
        qDebug()<<"nsuper"<<Lstore->nsuper;

        QVector<double> xnL(Lstore->nzval_colend[m-1]), ynL(Lstore->nzval_colend[m-1]);

        int nnzU = Ustore->nnz;
        //asubU = Ustore->rowind;

        QVector<double> xnU(nnzU), ynU(nnzU);

        int fsupc,istart,nsupr,nsupc,nrow;
        int count = 0;
        double value;
        QVector<int> level(m);
        for(int i = 0;i < m;i++){
            level[i] = 0;
        }
        for (int ksupno = 0; ksupno <= Lstore->nsuper; ++ksupno) {
            fsupc = Lstore->sup_to_colbeg[ksupno];
            istart = Lstore->rowind_colbeg[fsupc];
            nsupr = Lstore->rowind_colend[fsupc] - istart;
            nsupc = Lstore->sup_to_colend[ksupno] - fsupc;
            nrow = nsupr - nsupc;

            if ( nsupc == 1 ) {
                for (int j = 0; j < nrhs; j++) {
                    int luptr = Lstore->nzval_colbeg[fsupc];
                    for (int iptr=istart; iptr < Lstore->rowind_colend[fsupc]; iptr++){
                        int irow = Lstore->rowind[iptr];
                        value = ((double*)Lstore->nzval)[luptr];


                        if(fabs(value)>1e-9){
                            xnL[count] = fsupc+m;
                            ynL[count] = m-1 - irow;

                            count++;

                            if(irow >= fsupc){
                                if(irow == fsupc){
                                    level[irow] += 1;
                                }else if(irow > fsupc){
                                    level[irow] = std::max(level[fsupc],level[irow]);
                                }
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
                                xnL[count] = fsupc+i+m;
                                ynL[count] = m-1 - irow;

                                count++;

                                if(irow == fsupc+i){
                                    level[irow] += 1;
                                }else if(irow > fsupc+i){
                                    level[irow] = std::max(level[fsupc+i],level[irow]);
                                }
                            }
                            ++luptr;
                        }
                    }
                }

            } /* if-else: nsupc == 1 ... */
        } /* for L-solve */
        graphU->setData(xnU, ynU);
        graphL->setData(xnL, ynL);
        qDebug()<<"count "<<count;
        int maxLevel = 0;
        for(int i = 0;i < m;i++){
            //qDebug()<<i<<" "<<level[i];
            maxLevel = maxLevel > level[i] ? maxLevel : level[i];
        }
        qDebug()<<maxLevel;


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

