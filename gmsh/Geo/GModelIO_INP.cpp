// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include "GModel.h"
#include "OS.h"
#include "MPoint.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MQuadrangle.h"
#include "MTetrahedron.h"
#include "MHexahedron.h"
#include "MPrism.h"
#include "MPyramid.h"

template <class T>
static void writeElementsINP(FILE *fp, GEntity *ge, std::vector<T *> &elements,
                             bool saveAll)
{
  if(elements.size() && (saveAll || ge->physicals.size())) {
    const char *typ = elements[0]->getStringForINP();
    if(typ) {
      const char *str =
        (ge->dim() == 3) ? "Volume" :
        (ge->dim() == 2) ? "Surface" :
        (ge->dim() == 1) ? "Line" :
        "Point"; // currently unused
      fprintf(fp, "*ELEMENT, type=%s, ELSET=%s%d\n", typ, str, ge->tag());
      for(std::size_t i = 0; i < elements.size(); i++)
        elements[i]->writeINP(fp, elements[i]->getNum());
    }
  }
}

static std::string physicalName(GModel *m, int dim, int num)
{
  std::string name = m->getPhysicalName(dim, num);
  if(name.empty()) {
    char tmp[256];
    sprintf(tmp, "%s%d",
            (dim == 3) ? "PhysicalVolume" :
            (dim == 2) ? "PhysicalSurface" :
            (dim == 1) ? "PhysicalLine" :
            "PhysicalPoint",
            num);
    name = tmp;
  }
  for(std::size_t i = 0; i < name.size(); i++)
    if(name[i] == ' ') name[i] = '_';
  return name;
}

int GModel::writeINP(const std::string &name, bool saveAll,
                     bool saveGroupsOfNodes, double scalingFactor)
{
  FILE *fp = Fopen(name.c_str(), "w");
  if(!fp) {
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  if(noPhysicalGroups()) saveAll = true;

  indexMeshVertices(saveAll);
  std::vector<GEntity *> entities;
  getEntities(entities);

  fprintf(fp, "*Heading\n");
  fprintf(fp, " %s\n", name.c_str());

  fprintf(fp, "*NODE\n");
  for(std::size_t i = 0; i < entities.size(); i++)
    for(std::size_t j = 0; j < entities[i]->mesh_vertices.size(); j++)
      entities[i]->mesh_vertices[j]->writeINP(fp, scalingFactor);

  fprintf(fp, "******* E L E M E N T S *************\n");
  for(viter it = firstVertex(); it != lastVertex(); ++it) {
    writeElementsINP(fp, *it, (*it)->points, saveAll);
  }
  for(eiter it = firstEdge(); it != lastEdge(); ++it) {
    writeElementsINP(fp, *it, (*it)->lines, saveAll);
  }
  for(fiter it = firstFace(); it != lastFace(); ++it) {
    writeElementsINP(fp, *it, (*it)->triangles, saveAll);
    writeElementsINP(fp, *it, (*it)->quadrangles, saveAll);
  }
  for(riter it = firstRegion(); it != lastRegion(); ++it) {
    writeElementsINP(fp, *it, (*it)->tetrahedra, saveAll);
    writeElementsINP(fp, *it, (*it)->hexahedra, saveAll);
    writeElementsINP(fp, *it, (*it)->prisms, saveAll);
    writeElementsINP(fp, *it, (*it)->pyramids, saveAll);
  }

  std::map<int, std::vector<GEntity *> > groups[4];
  getPhysicalGroups(groups);

  // save elements sets for each physical group (currently we don't save point
  // elements: is there this concept in Abaqus?)
  for(int dim = 1; dim <= 3; dim++) {
    for(std::map<int, std::vector<GEntity *> >::iterator it = groups[dim].begin();
        it != groups[dim].end(); it++) {
      std::vector<GEntity *> &entities = it->second;
      fprintf(fp, "*ELSET,ELSET=%s\n",
              physicalName(this, dim, it->first).c_str());
      int n = 0;
      for(std::size_t i = 0; i < entities.size(); i++) {
        for(std::size_t j = 0; j < entities[i]->getNumMeshElements(); j++) {
          MElement *e = entities[i]->getMeshElement(j);
          if(n && !(n % 10)) fprintf(fp, "\n");
          fprintf(fp, "%lu, ", e->getNum());
          n++;
        }
      }
      fprintf(fp, "\n");
    }
  }

  // save node sets for each physical group (here we include node sets on
  // physical points)
  if(saveGroupsOfNodes) {
    for(int dim = 0; dim <= 3; dim++) {
      for(std::map<int, std::vector<GEntity *> >::iterator it = groups[dim].begin();
          it != groups[dim].end(); it++) {
        std::set<MVertex *> nodes;
        std::vector<GEntity *> &entities = it->second;
        for(std::size_t i = 0; i < entities.size(); i++) {
          for(std::size_t j = 0; j < entities[i]->getNumMeshElements(); j++) {
            MElement *e = entities[i]->getMeshElement(j);
            for(std::size_t k = 0; k < e->getNumVertices(); k++)
              nodes.insert(e->getVertex(k));
          }
        }
        fprintf(fp, "*NSET,NSET=%s\n",
                physicalName(this, dim, it->first).c_str());
        int n = 0;
        for(std::set<MVertex *>::iterator it2 = nodes.begin();
            it2 != nodes.end(); it2++) {
          if(n && !(n % 10)) fprintf(fp, "\n");
          fprintf(fp, "%ld, ", (*it2)->getIndex());
          n++;
        }
        fprintf(fp, "\n");
      }
    }
  }

  fclose(fp);
  return 1;
}
