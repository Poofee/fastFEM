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
    double ydisp = -0.1;
    double xingcheng = 2.2;

    for(int i = 1; i <= 20;i++){
        relay.remesh(xdisp,ydisp);
        relay.stepIncrement();
    }
}

void testRemesh::testRemesh3DParallel()
{
    Relay relay;
//    char fileName[] = "valve/valve3dnewnoairnoocc";
    char fileName[] = "valve/bb";
    relay.setFileName(fileName);
    relay.openGeo();
    relay.setAirTag(1);
    relay.setXiantieTag(2);

    double xdisp = 0;
    double ydisp = -1;
    double zdisp = 0;
    double xingcheng = 2.2;

    relay.stepIncrement();
    relay.remesh3DParallel(xdisp,ydisp,zdisp);
}

/**
 * @brief 测试to5的remesh
 *
 */
void testRemesh::testRemesh3DRotate()
{
    Relay relay;
    char fileName[] = "to5/to5new";
    relay.setFileName(fileName);
    relay.openGeo();
    relay.setAirTag(13);
    relay.setXiantieTag(1);

    double angle = 0.5;
    /** TO5的轴其实与x轴平行 **/
    for(double a = 0; a <= 3.4;a += angle){
        relay.stepIncrement();
        relay.remesh3DRotate(angle,1.74579962392938, -4.10996380979215, 0.00232963976730521,1,0,0);
    }
}

void testRemesh::testMove3DMeshParallel()
{
    Relay relay;
    char fileName[] = "valve/bb";
    relay.setFileName(fileName);
    relay.openGeo();
    relay.setAirTag(1);
    relay.setXiantieTag(2);

    double xdisp = -1;
    double ydisp = -5;
    double zdisp = -1;
    double xingcheng = 2.2;

    relay.stepIncrement();
    relay.remesh3DParallel(xdisp,ydisp,zdisp);
}

void testRemesh::testRotate3DMesh()
{
    Relay relay;
    char fileName[] = "valve/bb";
    relay.setFileName(fileName);
    relay.openGeo();
//    relay.setAirTag(1);
    relay.setXiantieTag(1);

    double angle = 90;
    double xingcheng = 2.2;

    for(int i = 0;i < 3;i++){
        relay.stepIncrement();
        relay.remesh3DRotate(angle,0,-1.8-1.25,-2.1+1.05+0.4/2,0,1,0);
    }
}
