#pragma once
#include "datatype.h"
#include "plot.h"

class CFastFEMcore
{
public:
	int num_pts;//节点数目
	int num_ele;//单元数目
	int numDomain;//域数目
	CNode * pmeshnode;
	CElement * pmeshele;
	char filename[256];//

	double  Precision;//计算精度
	double  Relax;
	int		LengthUnits;//长度单位变量，整型
	CMaterial* materialList;//定义材料类

	Plot * thePlot;
	CFastFEMcore();
	~CFastFEMcore();	
	// load mesh
	int LoadMeshCOMSOL(char*fn);//载入分网信息函数
	bool StaticAxisymmetricTLM();//使用TLM静态轴对称磁场的计算函数
	double CalcForce();//电磁力计算函数
	int openProject(QString proFile);//打开工程文件函数
	int preCalculation();//预计算函数
	int solve();//求解函数
	void readProjectElement(QXmlStreamReader &reader);//读取project节点
	void readDomainElement(QXmlStreamReader &reader, int i);//读取Domain节点
	void readBHElement(QXmlStreamReader &reader,int i);//读取BH节点
	int StaticAxisymmetricNR();//使用NR静态轴对称磁场的计算函数
};

