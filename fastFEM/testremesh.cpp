#include "testremesh.h"


testRemesh::testRemesh()
{

}

void testRemesh::testRemesh2D()
{
    Relay relay;
    char fileName[] = "JRS1250/JRS1250bgm";
    relay.setFileName(fileName);
    relay.openGeo();
    relay.setAirTag(10);
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

}

void testRemesh::testRemesh3DRotate()
{

}
