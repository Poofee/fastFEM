#pragma once

class FEMface{
public:
    int n1,n2,n3;/** 三角形面单元的三个顶点编号 **/

    FEMface():n1(0),n2(0),n3(0){}

    /** 赋值函数 **/
    FEMface& operator=(const FEMface& f);
    /** 比较函数 **/
    bool operator < (const FEMface& f);
    bool operator == (const FEMface& f);
    bool operator != (const FEMface& f);
};

class FEMedge{
public:
    int start;
    int end;

    FEMedge():start(0),end(0){}
    FEMedge& operator=(const FEMedge& e){this->start = e.start;this->end = e.end; return *this;}
    /** 注意两个都要排序 **/
    bool operator < (const FEMedge& e) { return this->start < e.start || (this->start == e.start && this->end < e.end); }
    bool operator == (const FEMedge& e) { return this->start == e.start && this->end == e.end ; }
    bool operator != (const FEMedge& e) { return this->start != e.start || this->end != e.end ; }
};

typedef struct _CNode
{
    double x, y, z;//the co
	double A;//the solution
    double bdr;/** boundary type 0:非边界;1:几何边界;2:物理边界3:既是物理边界，又是几何边界 **/
	double I;//current source
	double pm;
}CNode;

typedef struct _CElement
{
    int n[10];/** 0-3 存储的是节点编号，4-9存储的是棱单元编号 **/
    int ele_type;
    int physic_tag;
    int geometry_tag;/** 与face的编号对应 **/
    int domain;/** 从工程文件会读取一个材料列表，这是当前单元的编号；**/
	double P[3];// Pi, Pj, Pk;
	double Q[3];// Qi, Qj, Qk;
    double AREA;/** 单元的面积变量 **/
    double Bx,By,Bz,B;/** 单元的磁感应强度B及其水平和竖直方向上的分量 **/
    double miu,miut;/** 对于非线性单元，当前步的miu和上一步的miu；**/
    double rc,zc;
    double ydot;
    double Y11,Y12,Y13,Y22,Y23,Y33;
    bool LinearFlag;/** 定义逻辑变量LinearFlag，用来判断具体单元是否处于线性区域 **/
}CElement;

typedef struct _CElement4
{
    int n[4];// ni, nj, nk;//
    double AREA;/*单元的面积变量*/
    double Bx,By,B;/** 单元的磁感应强度B及其水平和竖直方向上的分量 **/
    double miu,miut;/** 对于非线性单元，迭代时的miu和初始时的miu；**/
    int domain;/** 从工程文件会读取一个材料列表，这是当前单元的编号；**/
    double rc,zc;
    bool LinearFlag;/** 定义逻辑变量LinearFlag，用来判断具体单元是否处于线性区域 **/
}CElement4;
class CMaterial{
public:
	double miu;		// permeabilities, relative
	int BHpoints;
	double   *Bdata;
	double *Hdata;
	//double *slope;
	double H_c;				// magnetization, A/m
	double Theta_m;			// orientation of magnetization, degrees
	double Jr;			// applied current density, MA/m^2
    double I;/** 线圈电流 **/
    double tau;/** 线圈匝数比 **/
    double sigma;/** 电导率 **/
	
	CMaterial();
	~CMaterial();
	double getMiu(double B);
    double getdvdB(double B);/** 相对磁阻率对B的偏微分 **/
};
class CBlock{
public:
	int BlockNum;
	
};
class ResistMarix{
public:
	double Y11;
	double Y12;
	double Y13;
	double Y22;
	double Y23;
	double Y33;
};
typedef struct _Resist4Matrix{
	double Y[10];
}Resist4Matrix;

typedef struct _VoltageQ4{
	double V[4];
}VoltageQ4;
class Voltage{
public:
	double V12, V13, V23;
};
