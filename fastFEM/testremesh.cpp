#include "testremesh.h"


testRemesh::testRemesh()
{

}

void testRemesh::testRemesh2D()
{
    Relay relay;
    char fileName[] = "JRS1250/JRS1250bgmband";
    relay.setFileName(fileName);
    relay.openGeo();
    relay.setAirTag(13);
    relay.setXiantieTag(1);

    double xdisp = 0;
    double ydisp = 0.1;
    double xingcheng = 2.2;

    for(int i = 1; i <= 20;i++){
        relay.remesh(xdisp,ydisp);
        relay.stepIncrement();
    }
}

void testRemesh::testRemesh3DParallel()
{
    Relay relay;
    char fileName[] = "valve/valve3dnewnoairnoocc";
    relay.setFileName(fileName);
    relay.openGeo();
}

/**
 * @brief 测试to5的remesh
 *
 */
void testRemesh::testRemesh3DRotate()
{
    Relay relay;
    char fileName[] = "to5/to5";
    relay.setFileName(fileName);
    relay.openGeo();
}
