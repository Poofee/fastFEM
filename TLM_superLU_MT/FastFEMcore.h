#pragma once
#include "datatype.h"

class CFastFEMcore
{
public:
	int num_pts;
	int num_ele;
	CNode * pmeshnode;
	CElement * pmeshele;
	char filename[256];

	double  Precision;
	double  Relax;
	int		LengthUnits;

	CFastFEMcore();
	~CFastFEMcore();	
	// load mesh
	int LoadMesh();
	bool StaticAxisymmetric();
};

