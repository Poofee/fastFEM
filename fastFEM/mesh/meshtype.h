// Copyright 2020 Poofee (https://github.com/Poofee)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------
/*****************************************************************************
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Poofee                                                          *
 *  Email:   poofee@qq.com                                                   *
 *  Address:                                                                 *
 *  Original Date: 2020-07-15                                                *
 *                                                                           *
 *****************************************************************************/

#ifndef MESHTYPE_H
#define MESHTYPE_H
enum GenTypes {
	GEN_UNKNOWN,
	GEN_TETLIB,
	GEN_NGLIB,
	GEN_ELMERGRID
};

enum PdeTypes {
	PDE_UNKNOWN,
	PDE_BOUNDARY,
	PDE_BULK
};
// node class
class node_t {
public:
    node_t();
    ~node_t();

    void setX(int, double);
    double getX(int) const;
    void setXvec(double*);
    double* getXvec();
    void setIndex(int);
    int getIndex() const;

private:
    double x[3];                     // 3d-coordinates
    int index;                       // optional tag
};

// base element class
class element_t {
public:
    element_t();
    ~element_t();

    void setNature(int);
    int getNature() const;
    void setCode(int);
    int getCode() const;
    void setNodes(int);
    int getNodes() const;
    void setIndex(int);
    int getIndex() const;
    void setSelected(int);
    int getSelected() const;
    int getNodeIndex(int) const;
    void setNodeIndex(int, int);
    int* getNodeIndexes() const;
    void newNodeIndexes(int);
    void deleteNodeIndexes();

private:
    int nature;                      // PDE_BULK, ...
    int code;                        // element code for Elmer (504, 808, ...)
    int nodes;                       // number of nodes
    int index;                       // bc/mat index as defined in input file
    int selected;                    // element is selected or not
    int* node;                       // list of nodes
};

// zero dimensional elements
class point_t: public element_t {
public:
    point_t();
    ~point_t();

    void setSharp(bool);
    bool isSharp() const;
    void setEdges(int);
    int getEdges() const;
    void setEdgeIndex(int, int);
    int getEdgeIndex(int) const;
    void newEdgeIndexes(int);
    void deleteEdgeIndexes();

private:
    bool sharp_point;                // marker
    int edges;                       // number of parent edges
    int* edge;                       // list of parent edges
};

// one dimensional elements
class edge_t: public element_t {
public:
    edge_t();
    ~edge_t();

    void setSharp(bool);
    bool isSharp() const;
    void setPoints(int);
    int getPoints() const;
    void setPointIndex(int, int);
    int getPointIndex(int) const;
    void newPointIndexes(int);
    void deletePointIndexes();
    void setSurfaces(int);
    int getSurfaces() const;
    void setSurfaceIndex(int, int);
    int getSurfaceIndex(int) const;
    void newSurfaceIndexes(int);
    void deleteSurfaceIndexes();

private:
    bool sharp_edge;                 // marker
    int points;                      // number of child points
    int* point;                      // list of points
    int surfaces;                    // number of parent surfaces
    int* surface;                    // list of parent surfaces
};

// two dimensional elements
class surface_t: public element_t {
public:
    surface_t();
    ~surface_t();

    void setEdges(int);
    int getEdges() const;
    void setEdgeIndex(int, int);
    int getEdgeIndex(int) const;
    void newEdgeIndexes(int);
    void deleteEdgeIndexes();
    void setElements(int);
    int getElements() const;
    void setElementIndex(int, int);
    int getElementIndex(int) const;
    void newElementIndexes(int);
    void deleteElementIndexes();
    void setNormalVec(double*);
    double* getNormalVec();
    double getNormal(int) const;
    void setNormal(int, double);
    void setVertexNormalVec(int, double*);
    void addVertexNormalVec(int, double*);
    void subVertexNormalVec(int, double*);
    double* getVertexNormalVec(int);

private:
    int edges;                       // number of child edges
    int* edge;                       // list of child edges
    int elements;                    // number of parent elements
    int* element;                    // list of parent elements
    double normal[3];                // unit (outward) normal
    double vertex_normals[4][3];     // unit (outward) normal on corner points
};

// mesh class
class mesh_t {
public:
    mesh_t();
    ~mesh_t();

    bool isUndefined() const;
    void clear();
    bool load(char*);
    bool save(char*);
    double* boundingBox();
    void setCdim(int);
    int getCdim() const;
    void setDim(int);
    int getDim() const;
    void setNodes(int);
    int getNodes() const;
    void setPoints(int);
    int getPoints() const;
    void setEdges(int);
    int getEdges() const;
    void setSurfaces(int);
    int getSurfaces() const;
    void setElements(int);
    int getElements() const;
    node_t* getNode(int);
    void setNodeArray(node_t*);
    void newNodeArray(int);
    void deleteNodeArray();
    point_t* getPoint(int);
    void setPointArray(point_t*);
    void newPointArray(int);
    void deletePointArray();
    edge_t* getEdge(int);
    void setEdgeArray(edge_t*);
    void newEdgeArray(int);
    void deleteEdgeArray();
    surface_t* getSurface(int);
    void setSurfaceArray(surface_t*);
    void newSurfaceArray(int);
    void deleteSurfaceArray();
    element_t* getElement(int);
    void setElementArray(element_t*);
    void newElementArray(int);
    void deleteElementArray();

private:
    void setDefaults();

    int cdim;                        // model coordinate dimension
    int dim;                         // model max element dimension
    int nodes;                       // number of nodes
    int points;                      // number of point elements
    int edges;                       // number of edge elements
    int surfaces;                    // number of surface elements
    int elements;                    // number of volume elements

    node_t* node;                    // array of nodes
    point_t* point;                  // array of point elements
    edge_t* edge;                    // array of edge elements
    surface_t* surface;              // array of surface elements
    element_t* element;              // array of volume elements

};

#endif // MESHTYPE_H
