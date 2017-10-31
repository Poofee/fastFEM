#pragma once

typedef struct _CNode
{
	double x, y;//the co
	double A;//the solution
	double bdr;//boundary type 0:非边界;1:几何边界;2:物理边界3:既是物理边界，又是几何边界
	double I;//current source
	double pm;
}CNode;

typedef struct _CElement
{
    int n[3];// ni, nj, nk;//
	double P[3];// Pi, Pj, Pk;
	double Q[3];// Qi, Qj, Qk;
	double AREA;//单元的面积变量
	double Bx,By,B;//单元的磁感应强度B及其水平和竖直方向上的分量
	double miu,miut;//对于非线性单元，迭代时的miu和初始时的miu；
	int domain;//从工程文件会读取一个材料列表，这是当前单元的编号；
    double rc,zc;
	bool LinearFlag;//定义逻辑变量LinearFlag，用来判断具体单元是否处于线性区域
}CElement;

typedef struct _CElement4
{
    int n[4];// ni, nj, nk;//
    double AREA;//单元的面积变量
    double Bx,By,B;//单元的磁感应强度B及其水平和竖直方向上的分量
    double miu,miut;//对于非线性单元，迭代时的miu和初始时的miu；
    int domain;//从工程文件会读取一个材料列表，这是当前单元的编号；
    double rc,zc;
    bool LinearFlag;//定义逻辑变量LinearFlag，用来判断具体单元是否处于线性区域
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
	
	CMaterial();
	~CMaterial();
	double getMiu(double B);
	double getdvdB(double B);//相对磁阻率对B的偏微分
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
    double Y11;
    double Y12;
    double Y13;
    double Y14;
    double Y22;
    double Y23;
    double Y24;
    double Y33;
    double Y34;
    double Y44;
	double Y10;
	double Y20;
	double Y30;
	double Y40;
}Resist4Matrix;

typedef struct _VoltageQ4{
	double V12, V13, V14, V23, V24, V34, V10, V20, V30, V40;
}VoltageQ4;
class Voltage{
public:
	double V12, V13, V23;
};
