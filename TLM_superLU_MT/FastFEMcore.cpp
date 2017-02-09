#include "FastFEMcore.h"


CFastFEMcore::CFastFEMcore() {
	Precision = 0.0;
	LengthUnits = 0;
}


CFastFEMcore::~CFastFEMcore() {
}


// load mesh
int CFastFEMcore::LoadMesh() {
	return 0;
}


bool CFastFEMcore::StaticAxisymmetric() {
	return false;
}
