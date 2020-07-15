#pragma once

class FEMface{
public:
    int n1,n2,n3;/** �������浥Ԫ������������ **/

    FEMface():n1(0),n2(0),n3(0){}

    /** ��ֵ���� **/
    FEMface& operator=(const FEMface& f);
    /** �ȽϺ��� **/
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
    /** ע��������Ҫ���� **/
    bool operator < (const FEMedge& e) { return this->start < e.start || (this->start == e.start && this->end < e.end); }
    bool operator == (const FEMedge& e) { return this->start == e.start && this->end == e.end ; }
    bool operator != (const FEMedge& e) { return this->start != e.start || this->end != e.end ; }
};

typedef struct _CNode
{
    double x, y, z;//the co
	double A;//the solution
    double bdr;/** boundary type 0:�Ǳ߽�;1:���α߽�;2:����߽�3:��������߽磬���Ǽ��α߽� **/
	double I;//current source
	double pm;
}CNode;

typedef struct _CElement
{
    int n[10];/** 0-3 �洢���ǽڵ��ţ�4-9�洢�����ⵥԪ��� **/
    int ele_type;
    int physic_tag;
    int geometry_tag;/** ��face�ı�Ŷ�Ӧ **/
    int domain;/** �ӹ����ļ����ȡһ�������б����ǵ�ǰ��Ԫ�ı�ţ�**/
	double P[3];// Pi, Pj, Pk;
	double Q[3];// Qi, Qj, Qk;
    double AREA;/** ��Ԫ��������� **/
    double Bx,By,Bz,B;/** ��Ԫ�ĴŸ�Ӧǿ��B����ˮƽ����ֱ�����ϵķ��� **/
    double miu,miut;/** ���ڷ����Ե�Ԫ����ǰ����miu����һ����miu��**/
    double rc,zc;
    double ydot;
    double Y11,Y12,Y13,Y22,Y23,Y33;
    bool LinearFlag;/** �����߼�����LinearFlag�������жϾ��嵥Ԫ�Ƿ����������� **/
}CElement;

typedef struct _CElement4
{
    int n[4];// ni, nj, nk;//
    double AREA;/*��Ԫ���������*/
    double Bx,By,B;/** ��Ԫ�ĴŸ�Ӧǿ��B����ˮƽ����ֱ�����ϵķ��� **/
    double miu,miut;/** ���ڷ����Ե�Ԫ������ʱ��miu�ͳ�ʼʱ��miu��**/
    int domain;/** �ӹ����ļ����ȡһ�������б����ǵ�ǰ��Ԫ�ı�ţ�**/
    double rc,zc;
    bool LinearFlag;/** �����߼�����LinearFlag�������жϾ��嵥Ԫ�Ƿ����������� **/
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
    double I;/** ��Ȧ���� **/
    double tau;/** ��Ȧ������ **/
    double sigma;/** �絼�� **/
	
	CMaterial();
	~CMaterial();
	double getMiu(double B);
    double getdvdB(double B);/** ��Դ����ʶ�B��ƫ΢�� **/
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
