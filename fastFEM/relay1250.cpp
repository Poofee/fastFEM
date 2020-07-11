#include "relay1250.h"

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
#include "qcustomplot.h"
#include "slu_mt_ddefs.h"


using namespace std;
using namespace arma;

#define PI 3.14159265358979323846
#define r 2

const double miu0 = PI*4e-7;

Relay1250::Relay1250()
{

}

Relay1250::~Relay1250()
{

}

/*!
 \brief 参数初始化

*/
void Relay1250::init()
{
    /** 文件名，需不带扩展名 **/
    sprintf(fileName,"%s","JRS1250/JRS1250bgm");
    /** 初始化时间步长 **/
    current_step = 0;
    min_time = 0;
    max_time = 4e-3;
    for(double t = 1e-5; t < 3e-3; t+=1e-5){
        timesteps.emplace_back(1e-5);
    }
    for(double t = 3e-3; t <= 3e-3; t+=5e-5){
        timesteps.emplace_back(5e-5);
    }

    /** 初始化变量 **/
    current_step = 0;
    Ddisplacements.resize(timesteps.size()+1);
    Ddisplacements.at(current_step) = 0;
    displacements.resize(timesteps.size()+1);
    displacements.at(current_step) = 0;
    velocities.resize(timesteps.size()+1);
    velocities.at(current_step) = 0;
    accelerations.resize(timesteps.size()+1);
    accelerations.at(current_step) = 0;
    PhiCoil.resize(timesteps.size()+1);
    PhiCoil.at(current_step) = 0;
    ICoil.resize(timesteps.size()+1);
    ICoil.at(current_step) = 0;
    magforces.resize(timesteps.size()+1);
    magforces.at(current_step) = 0;
    UCoil.resize(timesteps.size()+1);
    UCoil.at(current_step) = 420;

    MAX_NONLINEARSTEPS = 20;
    /** 设置材料参数 **/
    CMaterial* air_material = new CMaterial;
    materialList.emplace_back(air_material);

    CMaterial* pm_material = new CMaterial;
    materialList.emplace_back(pm_material);
    pm_material->H_c = 1e6;
    pm_material->miu = 1.12;

    CMaterial* dt4e_material = new CMaterial;
    materialList.emplace_back(dt4e_material);
    dt4e_material->sigma = 1/(2e-7);
    dt4e_material->BHpoints = 20;
    double* Hdata = (double*)malloc(20*sizeof(double));
    double* Bdata = (double*)malloc(20*sizeof(double));
    /** 1250的BH曲线 **/


    CMaterial* coil_material = new CMaterial;
    materialList.emplace_back(coil_material);
    coil_material->tau = 60/(3e-3 * 17e-3);

    /** 设置映射关系 **/
    materialMap[1] = dt4e_material;/** 衔铁 **/
    materialMap[2] = coil_material;/** 线圈 **/
    materialMap[3] = dt4e_material;/** 铁芯 **/
    materialMap[4] = dt4e_material;/** 铁芯 **/
    materialMap[5] = dt4e_material;/** 铁芯 **/
    materialMap[6] = pm_material;/** 永磁 **/
    materialMap[7] = dt4e_material;/** 铁芯（导磁环） **/
    materialMap[8] = dt4e_material;/** 铁芯 **/
    materialMap[9] = air_material;/** 外部空气 **/
    materialMap[10] = air_material;/** 可压缩空气 **/
    /** 设置形变区域 **/
    tag_xiantie = 1;
    tag_air = 10;

    Precision = 1e-6;
}



/*!
 \brief 计算区域index的电磁力，如果没有指定的话，就按衔铁算。

 \param index
*/
void Relay1250::calcMagForce(int index)
{

}

/*!
 \brief 计算三角形单元的基本项

 \param index
*/
void Relay1250::makeTriangle(int index)
{
    int k,m,n;
    double p0,p1,p2,q0,q1,q2,area;
    k = pmeshele[index].n[0];
    m = pmeshele[index].n[1];
    n = pmeshele[index].n[2];

    p0 = pmeshnode[m].y - pmeshnode[n].y;
    pmeshele[index].P[0] = p0;
    p1 = pmeshnode[n].y - pmeshnode[k].y;
    pmeshele[index].P[1] = p1;
    p2 = pmeshnode[k].y - pmeshnode[m].y;
    pmeshele[index].P[2] = p2;

    q0 = pmeshnode[n].x - pmeshnode[m].x;
    pmeshele[index].Q[0] = q0;
    q1 = pmeshnode[k].x - pmeshnode[n].x;
    pmeshele[index].Q[1] = q1;
    q2 = pmeshnode[m].x - pmeshnode[k].x;
    pmeshele[index].Q[2] = q2;

    area = 0.5*abs(p1 * q2 - q1 * p2);
    pmeshele[index].AREA = area;
    pmeshele[index].rc = (pmeshnode[k].x +
                          pmeshnode[m].x +
                          pmeshnode[n].x) / 3;
    pmeshele[index].zc = (pmeshnode[k].y +
                          pmeshnode[m].y +
                          pmeshnode[n].y) / 3;

    int flag = 0;
    for (int f = 0; f < 3; f++) {
        if (pmeshnode[pmeshele[index].n[f]].x < 1e-7) {
            flag++;
        }
    }
    /** 计算三角形重心半径 **/
    if (flag == 2) {
        pmeshele[index].ydot = pmeshele[index].rc;
    } else {
        pmeshele[index].ydot  = 1 / (pmeshnode[k].x + pmeshnode[m].x);
        pmeshele[index].ydot += 1 / (pmeshnode[k].x + pmeshnode[n].x);
        pmeshele[index].ydot += 1 / (pmeshnode[m].x + pmeshnode[n].x);
        pmeshele[index].ydot = 1.5 / pmeshele[index].ydot;
    }

//    double Y11,Y12,Y13,Y22,Y23,Y33;
    pmeshele[index].Y11 = q0 * q0 + p0 * p0;
    pmeshele[index].Y12 = q0 * q1 + p0 * p1;
    pmeshele[index].Y13 = q0 * q2 + p0 * p2;
    pmeshele[index].Y22 = q1 * q1 + p1 * p1;
    pmeshele[index].Y23 = q1 * q2 + p1 * p2;
    pmeshele[index].Y33 = q2 * q2 + p2 * p2;

    pmeshele[index].Y11 /= 4. * area;
    pmeshele[index].Y12 /= 4. * area;
    pmeshele[index].Y13 /= 4. * area;
    pmeshele[index].Y22 /= 4. * area;
    pmeshele[index].Y23 /= 4. * area;
    pmeshele[index].Y33 /= 4. * area;
}

/*!
 \brief 使用传统的牛顿法进行求解

*/
void Relay1250::AxisNRsolve()
{
    /** 计时 **/
    clock_t time[10];
    int tt = 0;
    time[tt++] = SuperLU_timer_();
    /** 电路参数 **/
    double C1 = 200e-6;
    double C2 = 200e-6;
    double R1 = 1e3;
    double Rc = 28;
    double Y1= 1/R1;
    double Yc= 1/Rc;
    /** 变量 **/
    double beta = 1;
    /** 时间步循环 **/
    current_step = 0;
    double totalTime = 0;
    for(int i = 0; i < timesteps.size();i++){
        totalTime += timesteps.at(i);
        printf("Time step %d, current time is %lf s.\n",i,totalTime);
        /** remesh，要使用增量位移，向下为负，向上为正 **/
        remesh(0,Ddisplacements.at(i));
        /** 读取第i步的分网 **/
        loadMesh();
        /** 查找边界，默认一次边界 **/
        findBoundaryEdges(-1);
        findBoundaryPoints(-1);

        double* unknown_b = (double*)calloc(num_pts - boundaryPoints.size()+1, sizeof(double));
        /** 将边界点排序到末尾 **/
        int node_bdr = boundaryPoints.size();
        int node_all = allPoints.size();
        /** 初始化，提出读取的gmsh多余的节点以及边界点 **/
        for(int i2=0;i2<num_pts;i2++){
            pmeshnode[i2].bdr = 3;
        }
        for(int i2=0;i2<node_all;i2++){
            pmeshnode[allPoints.at(i2)].bdr = 0;
        }
        for(int i2=0;i2<node_bdr;i2++){
            pmeshnode[boundaryPoints.at(i2)].bdr = 3;
        }
        /** 计算固定区域的单元节点数目 **/
        for(int i_tri = 0; i_tri < num_triangle; i_tri++){
            int i1 = num_ele-num_triangle+i_tri;
            if(pmeshele[i1].geometry_tag == tag_air){
                air_position = i_tri;
                break;
            }
        }

        umat locs(2, 9 * num_triangle);
        locs.zeros();
        mat vals(1, 9 * num_triangle);
        vec bbJz = zeros<vec>(num_pts);
        uvec node_reorder = zeros<uvec>(num_pts);
        uvec node_pos = zeros<uvec>(num_pts);
        vec bn = zeros<vec>(num_pts);
        vec A = zeros<vec>(num_pts);
        vec A_old = A;/** 上一步的大小肯定和这一步不一样 **/
        double Y11,Y12,Y13,Y22,Y23,Y33;

        for (int ibdr = 0; ibdr < num_pts; ibdr++) {
            if (pmeshnode[ibdr].bdr == 3) {
                node_bdr++;
                node_reorder(num_pts - node_bdr) = ibdr;
                node_pos(ibdr) = num_pts - node_bdr;
                pmeshnode[ibdr].A = 0;
                A(ibdr) = 0;
            } else {
                node_reorder(ibdr - node_bdr) = ibdr;
                node_pos(ibdr) = ibdr - node_bdr;
            }
        }
        /** 单元初始化，如果是非线性区域，可以考虑优化先初始值设置 **/
        for(int i_tri = 0; i_tri < num_triangle; i_tri++){
            int i1 = num_ele-num_triangle+i_tri;
            pmeshele[i1].B = 0;
            /** 对于涡流区域，没有形变，直接将miu带过来 **/
            if(i_tri < air_position){
                pmeshele[i1].miu = pmeshelelast[i1].miut;
                pmeshele[i1].miut = pmeshelelast[i1].miut;
            }else{
                pmeshele[i1].miu = miu0;
                pmeshele[i1].miut = miu0;
            }
        }
        double miuold = miu0;
        /** 非线性迭代求解 **/
        for(int non_iter = 0;non_iter < MAX_NONLINEARSTEPS;non_iter++){

            /** 有限元装配 **/
            double ce[3][3] = { {0} };/** 雅可比矩阵 **/
            double ce1[3][3] = { {0} };/** 右侧矩阵 **/
            double cn[3][3] = { {0} };/** 用于牛顿迭代 **/
            double Te[3][3] = { {0} };/** 涡流矩阵 **/
            double De[3] = { 0 };/** 电流密度矩阵 **/
            double dt = timesteps.at(non_iter);

            int pos = 0;
            for(int i_tri = 0; i_tri < num_triangle; i_tri++){
                int k,m,n,i1;
                i1 = num_ele-num_triangle+i_tri;
                k = pmeshele[i1].n[0];m = pmeshele[i1].n[1];n = pmeshele[i1].n[2];

                /** 计算除了磁导率部分的系数 **/
                makeTriangle(i1);

                /** 计算雅可比矩阵 **/
                double ydot = pmeshele[i1].ydot;
                double miut = pmeshele[i1].miut;
                /** 相当于C矩阵 **/
                Y11 = pmeshele[i1].Y11;Y22 = pmeshele[i1].Y22;Y33 = pmeshele[i1].Y33;
                Y12 = pmeshele[i1].Y12;Y13 = pmeshele[i1].Y13;Y23 = pmeshele[i1].Y23;

                ce[0][0] = Y11 / miut / ydot;/** 应当是当前迭代的mu **/
                ce[1][1] = Y22 / miut / ydot;
                ce[2][2] = Y33 / miut / ydot;

                ce[0][1] = Y12 / miut / ydot;
                ce[0][2] = Y13 / miut / ydot;
                ce[1][2] = Y23 / miut / ydot;

                double v[3];
                v[0] = Y11*A(k) + Y12*A(m) + Y13*A(n);
                v[1] = Y12*A(k) + Y22*A(m) + Y23*A(n);
                v[2] = Y13*A(k) + Y23*A(m) + Y33*A(n);

                CMaterial* mat = materialMap[pmeshele[i1].geometry_tag];
                if (non_iter != 0) {
                    double tmp;
                    if (mat->BHpoints == 0) {
                        tmp = 0;
                    } else {
                        tmp = mat->getdvdB(pmeshele[i1].B);
                        if (pmeshele[i].B > 1e-9){
                            tmp /= pmeshele[i].B * pmeshele[i].AREA;//B==0?
                            tmp /= ydot * ydot * ydot;
                        }

                    }
                    cn[0][0] = v[0] * v[0] * tmp;cn[1][1] = v[1] * v[1] * tmp;cn[2][2] = v[2] * v[2] * tmp;

                    cn[0][1] = v[0] * v[1] * tmp;cn[0][2] = v[0] * v[2] * tmp;cn[1][2] = v[1] * v[2] * tmp;

                    cn[1][0] = cn[0][1];cn[2][0] = cn[0][2];cn[2][1] = cn[1][2];
                }
                ce[0][0] += cn[0][0];ce[1][1] += cn[1][1];ce[2][2] += cn[2][2];

                ce[0][1] += cn[0][1];ce[0][2] += cn[0][2];ce[1][2] += cn[1][2];

                ce[1][0] = ce[0][1];ce[2][0] = ce[0][2];ce[2][1] = ce[1][2];

                /** 计算右侧向量 **/
                double dtmp = mat->tau * pmeshele[i1].AREA/3;
                De[0] = dtmp;De[1] = dtmp;De[2] = dtmp;
                double dtmp1 = 2*PI*Yc*dtmp*dtmp;
                /** 计算涡流矩阵 **/
                double etmp = mat->sigma/ydot*pmeshele[i1].AREA/12;
                Te[0][0] = etmp * 2+dtmp1;Te[1][1] = etmp+dtmp1;Te[2][2] = etmp+dtmp1;

                Te[0][1] = etmp+dtmp1;Te[0][2] = etmp+dtmp1;Te[1][2] = etmp+dtmp1;

                Te[1][0] = Te[0][1];Te[2][0] = Te[0][2];Te[2][1] = Te[1][2];
                for(int trow=0;trow<3;trow++){
                    for(int tcol=0;tcol<3;tcol++){
                        Te[trow][tcol] /= dt;
                        ce[trow][tcol] += Te[trow][tcol];
                    }
                }

                double jr = pmeshele[i].AREA*mat->Jr / 3;
                for (int j = 0; j < 3; j++) {
                    bbJz(pmeshele[i].n[j]) += jr;
                    /** 计算永磁部分 **/
                    bbJz(pmeshele[i].n[j]) -= mat->H_c / 2.*pmeshele[i].Q[j];
                }

                for (int row = 0; row < 3; row++) {
                    for (int col = 0; col < 3; col++) {
                        /** 判断节点是否在未知节点内 **/
                        /** 得到排序之后的编号 **/
                        int n_row = node_pos(pmeshele[i].n[row]);
                        int n_col = node_pos(pmeshele[i].n[col]);
                        if (n_row < num_pts - node_bdr && n_col < num_pts - node_bdr) {
                            locs(0, pos) = n_row;
                            locs(1, pos) = n_col;
                            vals(0, pos) = ce[row][col];
                            pos++;
                        }
                        /** 与雅可比矩阵相关的项 **/
                        bn(pmeshele[i].n[row]) += cn[row][col] * A(pmeshele[i].n[col]);
                    }
                }
            }/** 单元循环结束 **/
            /** 右下角的点 **/

            if (non_iter == 0) {
                locs.reshape(2, pos);
                vals.reshape(1, pos);
            }
            bn += bbJz;
            /** 调用线性求解器求解 **/
            /** 使用构造函数来生成稀疏矩阵 **/
            sp_mat X(true, locs, vals, num_pts - node_bdr, num_pts - node_bdr, true, true);

            for (int i = 0; i < num_pts - node_bdr; i++) {
                unknown_b[i] = bn(node_reorder(i));
            }
            //---------------------superLU_MT---------------------------------------
            CSuperLU_MT superlumt(num_pts - node_bdr, X, unknown_b);
            if (superlumt.solve() == 1) {
                printf("Error: superlumt.slove. Info:%d\n",superlumt.info);
                break;
            } else {
                double *sol = nullptr;
                A_old = A;
                sol = superlumt.getResult();

                for (int i = 0; i < num_pts - node_bdr; i++) {
                    pmeshnode[node_reorder(i)].A = sol[i];// / pmeshnode[i].x;//the A is r*A_real
                    A(node_reorder(i)) = sol[i];
                }
            }
            /** 更新磁场结果 **/
//            FILE *fp1 = fopen("B_T3_NR.txt", "w");
            for (int i = 0; i < num_triangle; i++) {
                double bx = 0;
                double by = 0;
                int i1 = num_ele-num_triangle+i;
                for (int j = 0; j < 3; j++) {
                    bx += pmeshele[i1].Q[j] * A(pmeshele[i1].n[j]);
                    by += pmeshele[i1].P[j] * A(pmeshele[i1].n[j]);
                }
                CMaterial* mat = materialMap[pmeshele[i1].geometry_tag];
                pmeshele[i].B = sqrt(bx*bx + by*by) / 2. / pmeshele[i].AREA / pmeshele[i1].ydot;
                pmeshele[i].Bx = bx / 2. / pmeshele[i1].AREA / pmeshele[i1].ydot;
                pmeshele[i].By = by / 2. / pmeshele[i1].AREA / pmeshele[i1].ydot;
                pmeshele[i].miut = mat->getMiu(pmeshele[i].B);
//                fprintf(fp1, "%lf \t %lf \t %lf \t %lf \t %lf\n", pmeshele[i].rc, pmeshele[i].zc, pmeshele[i].Bx, pmeshele[i].By, pmeshele[i].B);
//                y[i] = pmeshele[i].miut;
            }
//            fclose(fp1);
            double error = norm((A_old - A), 2) / norm(A, 2);
//            iter++;
            if (error < Precision) {
                break;
            }
            bn.zeros();
            pos = 0;

        }/** 牛顿迭代结束 **/



        /** 更新结果 **/
        current_step += 1;
        Ddisplacements.at(current_step) = 0;
        displacements.at(current_step) = 0;
        velocities.at(current_step) = 0;
        accelerations.at(current_step) = 0;
        PhiCoil.at(current_step) = 0;
        ICoil.at(current_step) = 0;
        magforces.at(current_step) = 0;
        UCoil.at(current_step) = 0;

        /** 回收空间 **/
        if(unknown_b) delete unknown_b;
    }


}

/*!
 \brief 使用传输线法进行求解

*/
void Relay1250::AxisTLMsolve()
{

}

/*!
 \brief 运行测试案例

*/
void Relay1250::run()
{
    init();
    openGeo();

    AxisNRsolve();

    /** 输出信息 **/
    outputResults();
}

/*!
 \brief 将一些结果进行输出。

*/
void Relay1250::outputResults()
{

}
