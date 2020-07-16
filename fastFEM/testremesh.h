#ifndef TESTREMESH_H
#define TESTREMESH_H

#include "relay.h"


/*!
 \brief 对gmsh的分网的remesh进行测试。

*/
class testRemesh
{
public:
    testRemesh();

    void testRemesh2D();
    void testRemesh3DParallel();
    void testRemesh3DRotate();
    void testMove3DMeshParallel();
    void testRotate3DMesh();
};

#endif // TESTREMESH_H
