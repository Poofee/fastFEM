#pragma once

class CNode
{
public:

	double x, y;//the co
	double A;//the solution
	double bdr;//boundary type 0:constant;
	double I;//current source
	double pm;
	//double length();
private:

};
class CElement
{
public:

	int n[3];// ni, nj, nk;//
	double P[3];// Pi, Pj, Pk;
	double Q[3];// Qi, Qj, Qk;
	double AREA;
	double Bx,By,B;
	double miu,miut;//对于非线性单元，迭代时的miu和初始时的miu；
	int domain;//从工程文件会读取一个材料列表，这是当前单元的编号；
    double rc,zc;
	bool LinearFlag;
	//double y12, y23, y31;
	//double y11, y22, y33;
	//double y10, y20, y30;
	//double y40, y50, y60;
	//double vr12, vr23, vr31;//To save some space...
	//double vi12, vi23, vi31;
	//double vr10, vr20, vr30;
	//double vi10, vi20, vi30;
	//double vr40, vr50, vr60;
	//double vi40, vi50, vi60;

private:
};
class CMaterial{
public:
	double miu;		// permeabilities, relative
	int BHpoints;
	double   *Bdata;
	double *Hdata;
	double *slope;
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

class Voltage{
public:
	double V12, V13, V23;
};
