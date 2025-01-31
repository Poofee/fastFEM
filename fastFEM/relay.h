// Copyright 2020 Poofee (https://github.com/Poofee)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------
/*****************************************************************************
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Poofee                                                          *
 *  Email:   poofee@qq.com                                                   *
 *  Address:                                                                 *
 *  Original Date: 2020-07-15                                                *
 *                                                                           *
 *****************************************************************************/
#ifndef RELAY_H
#define RELAY_H

#include "datatype.h"
#include "slu_mt_ddefs.h"

#include "meshGFace.h"
#include "meshGRegion.h"
#include "Generator.h"
#include "GModel.h"
#include "gmsh.h"
#include "MLine.h"
#include "Field.h"

#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <cassert>

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

class Relay
{
public:
    Relay();
    ~Relay();

    virtual void init();
    /** 几何相关 **/
    void setFileName(char fn[]);
    void openGeo();
    void moveFace(GFace *f,double dx,double dy,double dz);
    void moveRegion(GRegion *r,double dx,double dy,double dz);
    void rotateRegion(GRegion* r,double dangle,double x1, double y1,double z1,
                      double v1, double v2,double v3);
    Point3f RotateByVector(double theta, double old_x, double old_y, double old_z, double vx, double vy, double vz);

    /** 边界相关 **/
    void findBoundaryPoints(int index);
    void findBoundaryEdges(int index);
    void findBoundaryEdges3D(int index);
    void findBoundaryFaces(int index);
    /** 分网相关 **/
    void setXiantieTag(int xiantie);
    void setAirTag(int air);
    void loadMesh();
    void deleteFaceMesh(GFace* f);
    void deleteRegionMesh(GRegion* r);
    void remesh(double dx, double dy);
    void remesh3DParallel(double dx, double dy, double dz);
    void remesh3DRotate(double dangle,double x1, double y1,double z1,
                                      double v1, double v2,double v3);
    void updateField();


    /** 求解器相关 **/
    virtual void run();
    void stepIncrement();
    /** 后处理相关 **/
    virtual void calcMagForce(int index);
    virtual void calcMagTorque(int index);
    virtual void outputResults();
protected:
    int num_pts;/** 读取的节点数目 **/
    int num_ele;/** 读取的单元数目 **/
    int numDomain;/** 域数目 **/

    int num_triangle;/** 三角形单元的数目 **/
    int index_triagle;/** 三角形单元在msh文件保存的起始位置 **/

    int num_tet;/** 四面体单元的数目 **/
    int index_tet;/** 四面体单元在msh文件保存的起始位置 **/

    int tag_xiantie; /** 可移动的区域 **/
    int tag_air; /** 重分网区域，默认geo文件中该区域为最后一个区域。 **/
    int air_position;/** 在分网中，形变区域的起始位置。 **/

    int current_step;/** 当前仿真的步数 **/
    int MAX_NONLINEARSTEPS;/** 牛顿迭代法的最大迭代步数 **/
    int total_steps;/** 总共的迭代步数 **/

    double  Precision;/** 计算精度 **/
    double  Relax;

    double min_time;
    double max_time;

    double min_position;
    double max_position;

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
    std::vector<int> freePoints;/** 所有的自由分网点编号 **/
    std::vector<FEMface> boundaryFaces;/** 外部边界面 **/
    std::vector<FEMface> allFaces;/** 所有的边 **/
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

#endif // RELAY_H
