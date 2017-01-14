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
	double y11,y22,y33;
	double y12, y23, y13;
	double AREA;
	double Bx, By, B;
	double miu,miu_t;
	int domain;
    double rc,zc;
	

private:
};
