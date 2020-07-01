#ifndef RELAY1250_H
#define RELAY1250_H

#include "datatype.h"
#include "slu_mt_ddefs.h"

#include "meshGFace.h"
#include "Generator.h"
#include "GModel.h"
#include "gmsh.h"
#include "MLine.h"
#include "Field.h"

#include <vector>
#include <map>
#include<iostream>
#include<algorithm>
#include<cassert>

template <class ForwardIterator>
ForwardIterator Myunique (ForwardIterator first, ForwardIterator last)
{
    if (first==last) return last;

    ForwardIterator result = first;
    while (++first < last)
    {
        if(*result == *first){
            /** 由于在这里已经到达last了，导致下一次加的时候就超过
             * last了，报错 **/
            while (++first != last) {
//                qDebug()<<*first;
                if(*result != *first){
                    *result = *first;
                    break;
                }
            }
        }else{
            *(++result)=*first;
//            qDebug()<<*first;
        }
        /** 最后一个数是不重复的 **/
//        if (first==last-1){
//            break;
//        }
        /** 最后几个数是重复的 **/
        if (first==last){
            return result;
        }
    }
    return ++result;
}

class Relay1250
{
public:
    Relay1250();
    ~Relay1250();

    void init();
    void openGeo();
    void loadMesh();
    void findBoundaryPoints(int index);
    void findBoundaryEdges(int index);
    void deleteFaceMesh(GFace* f);
    void moveFace(GFace *f,double dx,double dy,double dz);
    void remesh(double dx, double dy);
    void AxisNRsolve();
    void AxisTLMsolve();
    void run();
    void outputResults();

protected:


private:
    int num_pts;/** 节点数目 **/
    int num_ele;/** 单元数目 **/
    int numDomain;/** 域数目 **/

    int tag_xiantie; /** 可移动的区域 **/
    int tag_air; /** 重分网区域 **/

    int current_step;
    int MAX_NONLINEARSTEPS;
    int total_steps;

    double  Precision;/** 计算精度 **/
    double  Relax;

    double min_time;
    double max_time;

    /** 结果 **/
    std::vector<double> timesteps;/** 每一步的时间差 **/
    std::vector<double> displacements;/** 总位移 **/
    std::vector<double> Ddisplacements;/** 每一步发生的位移 **/
    std::vector<double> velocities;/** 速度 **/
    std::vector<double> accelerations;/** 加速度 **/
    std::vector<double> PhiCoil;/** 线圈磁通 **/
    std::vector<double> ICoil;/** 线圈电流 **/
    std::vector<double> magforces;/** 电磁力 **/
    std::vector<double> UCoil;/** 线圈电压 **/

    /** 边界 **/
    std::vector<int> boundaryPoints;/** 外部边界点 **/
    std::vector<int> allPoints;/** 所有的分网点编号 **/
    std::vector<FEMedge> boundaryEdges;/** 外部边界边 **/
    std::vector<FEMedge> allEdges;/** 所有的边 **/

    CNode * pmeshnode;
    CElement * pmeshele;
    CNode * pmeshnodelast;
    CElement * pmeshelelast;

    char fileName[256];//

    GModel* model;/** gmsh model **/

//    int		LengthUnits;/** 长度单位变量，整型 **/
    std::vector<CMaterial*> materialList;/** 定义材料类 **/
    std::map<int,CMaterial*> materialMap;/** 建立区域与材料的映射关系 **/
};

#endif // RELAY1250_H
