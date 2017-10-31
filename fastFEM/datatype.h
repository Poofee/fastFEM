#pragma once

typedef struct _CNode
{
	double x, y;//the co
	double A;//the solution
	double bdr;//boundary type 0:�Ǳ߽�;1:���α߽�;2:����߽�3:��������߽磬���Ǽ��α߽�
	double I;//current source
	double pm;
}CNode;

typedef struct _CElement
{
    int n[3];// ni, nj, nk;//
	double P[3];// Pi, Pj, Pk;
	double Q[3];// Qi, Qj, Qk;
	double AREA;//��Ԫ���������
	double Bx,By,B;//��Ԫ�ĴŸ�Ӧǿ��B����ˮƽ����ֱ�����ϵķ���
	double miu,miut;//���ڷ����Ե�Ԫ������ʱ��miu�ͳ�ʼʱ��miu��
	int domain;//�ӹ����ļ����ȡһ�������б����ǵ�ǰ��Ԫ�ı�ţ�
    double rc,zc;
	bool LinearFlag;//�����߼�����LinearFlag�������жϾ��嵥Ԫ�Ƿ�����������
}CElement;

typedef struct _CElement4
{
    int n[4];// ni, nj, nk;//
    double AREA;//��Ԫ���������
    double Bx,By,B;//��Ԫ�ĴŸ�Ӧǿ��B����ˮƽ����ֱ�����ϵķ���
    double miu,miut;//���ڷ����Ե�Ԫ������ʱ��miu�ͳ�ʼʱ��miu��
    int domain;//�ӹ����ļ����ȡһ�������б����ǵ�ǰ��Ԫ�ı�ţ�
    double rc,zc;
    bool LinearFlag;//�����߼�����LinearFlag�������жϾ��嵥Ԫ�Ƿ�����������
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
	double getdvdB(double B);//��Դ����ʶ�B��ƫ΢��
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
