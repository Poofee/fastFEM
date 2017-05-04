#include <QMenu>
#include <QVector>
#include <QAction>
#include <QMenuBar>
#include <stdio.h>
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
    fp = fopen("..\\..\\SuperLU_MT_test\\g5.rua", "r");
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
    if (numer_lines) {
        dReadValues(fp, nnz, a, valnum, valsize);
    }

    fclose(fp);

    //-----------qt plot
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


    QVector<double> x(nnz), y(nnz);
    QVector<double> x1(4), y1(4);
    x1[0] = 0;x1[1] = 24;x1[2]=24;x1[3]=0;
    y1[0] = 0;y1[1] = 0;y1[2]=24;y1[3]=24;
    QCPCurve *newCurve = new QCPCurve(customplot->xAxis, customplot->yAxis);
    newCurve->setBrush(QColor(255, 0, 0,100));
    newCurve->setData(x1, y1);

    for(int i=0;i < m;i++){
        for(int j = xa[i];j < xa[i+1];j++){
            y[j] = 24 - i;
            x[j] = asub[j];
        }
    }
    graph1->setData(x, y);
    //customplot->rescaleAxes(true);
    customplot->replot();


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
        qDebug()<<"#NZ in factor L = "<<Lstore->nnz;
        qDebug()<<"#NZ in factor U = "<<Ustore->nnz;
        qDebug()<<"#NZ in L+U = "<<Lstore->nnz + Ustore->nnz - L.ncol;

        superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
        qDebug()<<"L\\U MB "<<superlu_memusage.for_lu / 1024 / 1024<<"\ttotal MB needed "<<superlu_memusage.total_needed / 1024 / 1024<<"\texpansions "<<superlu_memusage.expansions ;

    }

    SUPERLU_FREE(rhs);
    SUPERLU_FREE(xact);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
}

