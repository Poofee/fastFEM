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
	double miu,miut;
	int domain;
    double rc,zc;
	double y12, y23, y31;
	double y11, y22, y33;
	double y10, y20, y30;
	//double y40, y50, y60;
	double vr12, vr23, vr31;//To save some space...
	double vi12, vi23, vi31;
	double vr10, vr20, vr30;
	double vi10, vi20, vi30;
	//double vr40, vr50, vr60;
	//double vi40, vi50, vi60;

private:
};
