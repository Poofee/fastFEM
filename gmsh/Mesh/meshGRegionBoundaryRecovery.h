// Gmsh - Copyright (C) 1997-2014 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef MESH_GREGION_BOUNDARY_RECOVERY_H
#define MESH_GREGION_BOUNDARY_RECOVERY_H

class GRegion;
class splitQuadRecovery;

bool meshGRegionBoundaryRecovery(GRegion *gr, splitQuadRecovery *sqr = 0);

#endif
