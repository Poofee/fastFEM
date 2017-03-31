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
#include "spral.h"
#include "qcustomplot.h"

using namespace std;
using namespace arma;

#define PI 3.14159265358979323846
#define r 2

const double ste = 4;
const double miu0 = PI*4e-7;

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


// 这部分使用TLM进行求解
bool CFastFEMcore::StaticAxisymmetricTLM() {
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
    mat vals(1,9 * num_ele);
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
    double* unknown_b = (double*)calloc(num_pts - node_bdr,sizeof(double));
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
            }
            else{
                ce[0][0] += rm[i].Y12;//在对角项上减去受控源
                ce[1][1] += rm[i].Y12;
                ce[0][1] = 0;//受控源在右侧，所以为0
            }
            if (rm[i].Y13 < 0){
                ce[0][2] = rm[i].Y13;
            }
            else{
                ce[0][0] += rm[i].Y13;
                ce[2][2] += rm[i].Y13;
                ce[0][2] = 0;
            }
            if (rm[i].Y23 < 0){
                ce[1][2] = rm[i].Y23;
            }
            else{
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
                    if(n_row >= n_col){
                        locs(0, pos) = n_row;
                        locs(1, pos) = n_col;
                        vals(0, pos) = ce[row][col];
                        pos++;
                    }
                }
            }
        }
        //计算电流密度//要注意domain会不会越界
        double jr = pmeshele[i].AREA*materialList[pmeshele[i].domain - 1].Jr / 3;
        for (int j = 0; j < 3; j++) {
            bbJz(pmeshele[i].n[j]) += jr;
            // 计算永磁部分
            bbJz(pmeshele[i].n[j]) += materialList[pmeshele[i].domain - 1].H_c / 2.*pmeshele[i].Q[j];
        }
    }
    locs.reshape(2, pos);//重新调整大小
    vals.reshape(1, pos);
    //----using armadillo constructor function-----
    sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

    INL += bbJz;
    for (int i = 0; i < num_pts - node_bdr; i++) {
        unknown_b[i] = INL(node_reorder(i));
    }
    //---------------------spral---------------------------------------
    /* Derived types */
    void *akeep, *fkeep;
    struct spral_ssids_options options;
    struct spral_ssids_inform inform;

    /* Initialize derived types */
    akeep = NULL; fkeep = NULL; /* Important that these are NULL to start with */
    spral_ssids_default_options(&options);
    bool posdef = true;
    int n = num_pts - node_bdr;
    int * ptr   = (int*)const_cast<unsigned int*>(X.col_ptrs);
    long *ptr64 = (long*)malloc((n+1)*sizeof(long));
    for(int iptr=0;iptr<n+1;iptr++){
        ptr64[iptr] = ptr[iptr] + 1;
    }
    int * row    = (int*)malloc(X.n_nonzero*sizeof(int));
    for(int irow=0;irow < X.n_nonzero;irow++){
        row[irow] = X.row_indices[irow] + 1;
    }
    double *val = const_cast<double *>(X.values);

    /* Perform analyse and factorise with data checking */
    bool check = true;
    spral_ssids_analyse(check, n, NULL, ptr64, row, NULL, &akeep, &options,
                        &inform);
    if(inform.flag<0) {
        spral_ssids_free(&akeep, &fkeep);
        exit(1);
    }
    spral_ssids_factor(posdef, NULL, NULL, val, NULL, akeep, &fkeep, &options,
                       &inform);
    if(inform.flag<0) {
        spral_ssids_free(&akeep, &fkeep);
        exit(1);
    }

    /* Solve */
    spral_ssids_solve1(0, unknown_b, akeep, fkeep, &options, &inform);
    if(inform.flag<0) {
        spral_ssids_free(&akeep, &fkeep);
        qDebug() << "Error: spral.slove";
    }else{
        for (int i = 0; i < num_pts - node_bdr; i++) {
            pmeshnode[node_reorder(i)].A = unknown_b[i];// / pmeshnode[i].x;//the A is r*A_real
            A(node_reorder(i)) = unknown_b[i];
        }
    }
    //---------------------spral--end----------------------------------
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
    customplot->yAxis->setRange(0, 0.6);
    //customplot->xAxis2->setTicks(false);
    //customplot->yAxis->setScaleRatio(ui->widget->xAxis, 1.0);
    //---------the main loop---------------------------------
    int steps = 250;
    int count;
    double alpha = 1;
    Voltage *Vr = (Voltage*)calloc(D34.size(), sizeof(Voltage));
    Voltage *Vi = (Voltage*)calloc(D34.size(), sizeof(Voltage));
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

            //y[i] = (A(pmeshele[i].n[0]) + A(pmeshele[i].n[1]) + A(pmeshele[i].n[2]))/3;
            y[i] = pmeshele[i].B;
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

            Vr[j].V12 = (pmeshnode[k].A - pmeshnode[m].A) - Vi[j].V12;
            Vr[j].V23 = (pmeshnode[m].A - pmeshnode[n].A) - Vi[j].V23;
            Vr[j].V13 = (pmeshnode[n].A - pmeshnode[k].A) - Vi[j].V13;

            if (rm[i].Y12 < 0) {
                Vi[j].V12 = Vr[j].V12 * rtmp;
            } else {
                Vi[j].V12 = (1 - alpha)*Vi[j].V12 + alpha*0.5*(pmeshnode[k].A - pmeshnode[m].A);
            }
            INL(k) += -2. *Vi[j].V12*rm[i].Y12;
            INL(m) += 2. * Vi[j].V12 *rm[i].Y12;
            if (rm[i].Y23 < 0) {
                Vi[j].V23 = Vr[j].V23*rtmp;
            } else {
                Vi[j].V23 = (1 - alpha)*Vi[j].V23 + alpha*0.5*(pmeshnode[m].A - pmeshnode[n].A);
            }
            INL(m) += -2. * Vi[j].V23*rm[i].Y23;
            INL(n) += 2. *Vi[j].V23*rm[i].Y23;
            if (rm[i].Y13 < 0) {
                Vi[j].V13 = Vr[j].V13 * rtmp;
            } else {
                Vi[j].V13 = (1 - alpha)*Vi[j].V13 + alpha*0.5*(pmeshnode[n].A - pmeshnode[k].A);
            }
            INL(n) += -2. * Vi[j].V13*rm[i].Y13;
            INL(k) += 2.0 *Vi[j].V13*rm[i].Y13;
        }
        INL += bbJz;
        for (int i = 0; i < num_pts - node_bdr; i++) {
            unknown_b[i] = INL(node_reorder(i));
        }
        /* Solve */
        spral_ssids_solve1(0, unknown_b, akeep, fkeep, &options, &inform);
        if(inform.flag<0) {
            spral_ssids_free(&akeep, &fkeep);
            qDebug() << "Error: spral.slove";
            break;
        }else{
            for (int i = 0; i < num_pts - node_bdr; i++) {
                pmeshnode[node_reorder(i)].A = unknown_b[i];// / pmeshnode[i].x;//the A is r*A_real
                A(node_reorder(i)) = unknown_b[i];
            }
        }
        double error = norm((A_old - A), 2) / norm(A, 2);
        //qDebug() << "iter: " << count;
        //qDebug() << "error: " << error;

        //graph1->setData(x, y);
        //customplot->rescaleAxes(true);
        //customplot->replot();

        if (error < Precision) {
            break;
        }
        INL.zeros();
    }
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
    // 回收空间
    if (rm != NULL) free(rm);
    if (ydot != NULL) free(ydot);
    if (Vi != NULL) free(Vi);
    if (Vr != NULL) free(Vr);
    if (unknown_b != NULL) free(unknown_b);

    int cuda_error = spral_ssids_free(&akeep, &fkeep);
    if(cuda_error!=0) exit(1);
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
            pmeshele[i].miu = 1 * miu0;
            pmeshele[i].miut = 1 * miu0;//must be 1
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
                    //qDebug() << "domainNum = " << numDomain;
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
            //qDebug() << "domainName = " << reader.readElementText();
        } else if (reader.name() == "miu") {
            materialList[i].miu = reader.readElementText().toDouble() * 4 * PI*1e-7;
            //qDebug() << "miu = " << materialList[i].miu;
        } else if (reader.name() == "BH") {
            readBHElement(reader, i);
        } else if (reader.name() == "Jr") {
            materialList[i].Jr = reader.readElementText().toDouble();
            //qDebug() << "Jr = " << materialList[i].Jr;
        } else if (reader.name() == "H_c") {
            materialList[i].H_c = reader.readElementText().toDouble();
            //qDebug() << "H_c = " << materialList[i].H_c;
        }
    }
}


void CFastFEMcore::readBHElement(QXmlStreamReader &reader, int i) {
    reader.readNextStartElement();
    //qDebug()<<reader.name();
    if (reader.name() == "BHpoints") {
        materialList[i].BHpoints = reader.readElementText().toInt();
        //qDebug() << "BHpoints = " << materialList[i].BHpoints;
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
                    bbJz(pmeshele[i].n[j]) += materialList[pmeshele[i].domain - 1].H_c / 2.*pmeshele[i].Q[j];
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
                    tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
                    tmp /= ydot[i] * ydot[i] * ydot[i];
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
                        if(n_row >= n_col){
                            locs(0, pos) = n_row;
                            locs(1, pos) = n_col;
                            vals(0, pos) = ce[row][col];
                            pos++;
                        }
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
        //---------------------spral---------------------------------------
        /* Derived types */
        void *akeep, *fkeep;
        struct spral_ssids_options options;
        struct spral_ssids_inform inform;

        /* Initialize derived types */
        akeep = NULL; fkeep = NULL; /* Important that these are NULL to start with */
        spral_ssids_default_options(&options);
        bool posdef = false;
        int n = num_pts - node_bdr;
        int * ptr   = (int*)const_cast<unsigned int*>(X.col_ptrs);
        long *ptr64 = (long*)malloc((n+1)*sizeof(long));
        for(int iptr=0;iptr<n+1;iptr++){
            ptr64[iptr] = ptr[iptr] + 1;
        }
        int * row    = (int*)malloc(X.n_nonzero*sizeof(int));
        for(int irow=0;irow < X.n_nonzero;irow++){
            row[irow] = X.row_indices[irow] + 1;
        }
        double *val = const_cast<double *>(X.values);

        /* Perform analyse and factorise with data checking */
        bool check = true;
        spral_ssids_analyse(check, n, NULL, ptr64, row, NULL, &akeep, &options,
                            &inform);
        if(inform.flag<0) {
            spral_ssids_free(&akeep, &fkeep);
            exit(1);
        }
        spral_ssids_factor(posdef, NULL, NULL, val, NULL, akeep, &fkeep, &options,
                           &inform);
        if(inform.flag<0) {
            spral_ssids_free(&akeep, &fkeep);
            exit(1);
        }

        /* Solve */
        spral_ssids_solve1(0, unknown_b, akeep, fkeep, &options, &inform);
        if(inform.flag<0) {
            spral_ssids_free(&akeep, &fkeep);
            qDebug() << "Error: spral.slove";
            break;
        }else{
            for (int i = 0; i < num_pts - node_bdr; i++) {
                pmeshnode[node_reorder(i)].A = unknown_b[i];// / pmeshnode[i].x;//the A is r*A_real
                A(node_reorder(i)) = unknown_b[i];
            }
        }

        int cuda_error = spral_ssids_free(&akeep, &fkeep);
        if(cuda_error!=0) exit(1);

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

            y[i] = pmeshele[i].B;
        }
        double error = norm((A_old - A), 2) / norm(A, 2);
        iter++;
        //qDebug() << "iter: " << iter;
        //qDebug() << "error: " << error;
        if (error < Precision || iter > 200) {
            break;
        }
        bn.zeros();
        pos = 0;

        //graph1->setData(x, y);
        //customplot->rescaleAxes(true);
        //customplot->replot();
    }
    if (rm != NULL) free(rm);
    if (ydot != NULL) free(ydot);
    if (unknown_b != NULL) free(unknown_b);
    return 0;
}


