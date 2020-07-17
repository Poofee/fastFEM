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
#ifndef MAG3DSTATIC_H
#define MAG3DSTATIC_H

#include "libstructure.h"

typedef struct {
  int   *nop ;/**  **/
  int   *flag ;/**  **/
  double *crd ;/**  **/
  int   *apid2phiid ;/**  **/
  COIL  *coil ;/**  **/
  COIL  *pmagnet ;/**  **/
  MTRL  *mtrl ;/**  **/
  double weight ;/**  **/
  double *x_d ;/**  **/
  double *nl_nu ;/**  **/
} DATA_STATIC ;

class Mag3DStatic
{
public:
    Mag3DStatic();

    void domesh();
    void solve();
    void assemble();
    void run();
    void tet_MakeElement_NL_Mag_Static_Newton_A_p(int el);
    void make_nu( int el, int nl, double *nu, double *dnu, double B2, double *nl_nu,MTRL mtrl, double weight );
    void make_SR( double *ae, double *fe, double *xi, double *yi, double *zi,
                  int *si, double *l, double D, double *xe, double B2, double dnu );
    void MakeElement_IP_Get_nuB2( int nl, int *n, double **a, double **b, double **v,
                                  MTRL mtrl );
    void tet_MakeElement_NL_Mag_Static_IP_1_sub( double *nu, double *dnu, double B2,
                                int n, double *a, double *b, double *v, MTRL mtrl ) ;
private:
    DATA_STATIC data;
};

#endif // MAG3DSTATIC_H
