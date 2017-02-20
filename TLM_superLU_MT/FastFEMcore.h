#pragma once
#include "datatype.h"
#include "plot.h"

class CFastFEMcore
{
public:
	int num_pts;//节点数目
	int num_ele;//单元数目
	CNode * pmeshnode;
	CElement * pmeshele;
	char filename[256];//

	double  Precision;
	double  Relax;
	int		LengthUnits;
	CMaterial* materialList;

	Plot * thePlot;
	CFastFEMcore();
	~CFastFEMcore();	
	// load mesh
	int LoadMesh();
	bool StaticAxisymmetricTLM();
	double CalcForce();
	int openProject();
	int preCalculation();
	int solve();
	void readProjectElement(QXmlStreamReader &reader);
	void readDomainElement(QXmlStreamReader &reader);
	void readBHElement(QXmlStreamReader &reader);
};

