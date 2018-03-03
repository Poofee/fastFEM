#pragma once
#include "datatype.h"
#include "plot.h"
#include "slu_mt_ddefs.h"

class CFastFEMcore
{
public:
	int num_pts;//节点数目
	int num_ele;//单元数目
	int numDomain;//域数目
	CNode * pmeshnode;
    CElement * pmeshele;
    CElement4 * pmeshele4;
	char filename[256];//

	double  Precision;//计算精度
	double  Relax;
	int		LengthUnits;//长度单位变量，整型
	CMaterial* materialList;//定义材料类

	Plot * thePlot;
	CFastFEMcore();
	~CFastFEMcore();	
	// load mesh
	int Load2DMeshCOMSOL(const char fn[]);//载入分网信息函数
    int LoadQ4MeshCOMSOL(const char fn[]);//读入四节点单元
	bool StaticAxisymmetricTLM();//使用TLM静态轴对称磁场的计算函数
	bool StaticAxisTLMNR();//新的迭代方法，NR+TLM
	bool StaticAxisT3TLMgroup();//三角形单元看作为整体
    bool StaticAxisQ4NR();//使用NR求解，四边形分网
    bool StaticAxisQ4TLM();//使用TLM求解，四边形分网
	bool StaticAxisT3NRTLM();//NR迭代的每一步使用TLM求解线性方程组
	bool StaticAxisT3VTM();//采用VTM方法来求解，单元隔离
	double CalcForce();//电磁力计算函数
	int openProject(QString proFile);//打开工程文件函数
	int preCalculation();//预计算函数
	int solve();//求解函数
	void readProjectElement(QXmlStreamReader &reader);//读取project节点
	void readDomainElement(QXmlStreamReader &reader, int i);//读取Domain节点
	void readBHElement(QXmlStreamReader &reader,int i);//读取BH节点
	int StaticAxisymmetricNR();//使用NR静态轴对称磁场的计算函数
    void myTriSolve(int ncore, SuperMatrix *L, SuperMatrix *U,
                    int_t *perm_r, int_t *perm_c, SuperMatrix *B, int_t *info);
    double getLocal4Matrix(int Ki, int Kj, int index);//返回单元矩阵的元素ij
	double getTLMKij(int Ki, int Kj, int index);
    double getP(int Ki, int Kj, double xi, double eta,int index);//被积分的那个函数
	double getPx(int Ki, int Kj, double xi, double eta, int index);
    double getdNidx(int Ki, double xi, double eta, int index);//返回形函数的对x偏导数
    double getdNidy(int Ki,double xi,double eta,int index);//返回形函数的对y的偏导数
    double getdxdxi(double eta, int index);
    double getdxdeta(double xi,int index);
    double getdydxi(double eta, int index);
    double getdydeta(double xi,int index);
    double getdNdxi(int i,double eta);
    double getdNdeta(int i,double xi);
    double getJacobi(double xi, double eta, int index);
	double getJi(int Ki, int index);
	double Ne(double xi, double eta, int index);
	double getx(double xi, double eta, int index);
	double gety(double xi, double eta, int index);
	double getA(double xi, double eta, int index);
	double getrA(double xi, double eta, int index);
	double getDij(int i, int j, int index);
	double getCij(int i, int j, double xi, double eta, int index);
	double getD(int i, int j, double xi, double eta, int index);

};

