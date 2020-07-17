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
#ifndef MAG2DTIME_H
#define MAG2DTIME_H

#include "libstructure.h"

#include <vector>

typedef struct {
  int   *nop ;
  int   *flag ;
  double *crd ;
  int   *apid2phiid ;
  double *ini ;
  COIL  *coil ;
  COIL  *pmagnet ;
  double (*make_nu)( int nl, int flg, double *xe, double D,
                     double *X, double *Y, double *Z, double mem_nu,
                     MTRL mtrl ) ;
  double mem_nu ;
  void (*make_sigma)( int el, int flg, MTRL mtrl, double *vec,
                      int *_rs, double *_sigma ) ;
  double *hetero_conductivity ;
  MTRL  *mtrl ;
  double starttime;
  double endtime;
  std::vector<double> timesteps;
//#ifdef __ADVMAG_NEW_SOLVER__
//  void **solvmat ;
//  SOLVER_FUNC *sfunc ;
//#else
//  SOLVER_MAT *smat ;
//  void (*SetMat)( SOLVER_MAT *smat, int szRow, int *indRow,
//                  int szCol, int *indCol, double *ae, double *fe ) ;
//#endif
} DATA_NS_EDDY_2D ;

class Mag2DTime
{
public:
    Mag2DTime();

    void init();
    void doMesh();
    void reMeshParallelShift(double dx,double dy);
    void solve();
    void run();

private:
    int MAX_NONLINEARSTEPS;
    DATA_NS_EDDY_2D data;
};

#endif // MAG2DTIME_H
