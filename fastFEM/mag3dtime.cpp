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
#include "tetralheral.h"

#include "mag3dtime.h"

Mag3DTime::Mag3DTime()
{

}

void Mag3DTime::doMesh()
{

}

void Mag3DTime::reMeshRotate(double theta)
{

}

void Mag3DTime::reMeshParallelShift(double dx, double dy, double dz)
{

}

void Mag3DTime::assemble()
{

}

void Mag3DTime::solve()
{

}

void Mag3DTime::run()
{

}
/*+ Function: make element
                 conventional piecewise linear tetrahedral
                 Non-Steady Eddy Current Analysis
                 A method +*/
double Mag3DTime::tet_MakeElement_NS_Eddy_A(int el)
{

    int    i, j, k ;

    int*   nop = data.nop ;
    double* crd = data.crd ;
    int*   flag = data.flag ;
    double* ini = data.ini ;
    COIL*  coil = data.coil ;
    COIL*  pmagnet = data.pmagnet ;
    MTRL*  mtrl = data.mtrl ;
    double mem_nu = data.mem_nu ;
    double* hetero_conductivity = data.hetero_conductivity ;

    int    *ne ;
    int    flg = flag[el] ;
    double delta_t = 1e-3;//OPTIONS_delta_t_get() ;
    double nu = 0.0 ;
    double sigma = 0.0 ;
    int    rs = 1 ;
    int    nl = -1 ;
    int    nc = -1 ;
    int    np = -1 ;
    double x[ap_elem], y[ap_elem], z[ap_elem] ;
    double inie[mp_elem] ;
    double D ;
    double l[mp_elem] ;
    double tetg[3] ;
    int    si[mp_elem] ;
    double xi[mp_elem], yi[mp_elem], zi[mp_elem] ;
    double X[mp_elem], Y[mp_elem], Z[mp_elem] ;
    double xy[mp_elem], yz[mp_elem], zx[mp_elem] ;
    double rr[mp_elem][mp_elem] ;
#ifdef szRow_def
#undef szRow_def
#endif
#define szRow_def nd_elem
    int    szRow = szRow_def ;
    int    indRow[szRow_def] ;

    double ae[szRow_def*szRow_def] ;
    double be[szRow_def] ;


    data.make_sigma( el, flg, *mtrl, hetero_conductivity, &rs, &sigma ) ;
    for( i=0 ; i<mtrl->ncoil ; i++ ) {
        if( flg == mtrl->coil_def[i].flag ) {
            nc = i ;
            break ;
        }
    }
    for( i=0 ; i<mtrl->npmagnet ; i++ ) {
        if( flg == mtrl->pmagnet_def[i].flag ) {
            np = i ;
            break ;
        }
    }


    ne = &(nop[el*nd_elem]) ;
    for( i=0,j=ap_elem ; i<mp_elem ; i++,j++ ) {
        indRow[i] = ne[j] ;
    }
    szRow = mp_elem ;

    for( i=0 ; i<mp_elem ; i++ ) inie[i] = ini[indRow[i]] ;

    for( i=0,k=0 ; i<szRow ; i++ ) {
        for( j=0 ; j<szRow ; j++,k++ ) {
            ae[k] = 0.0 ;
        }
        be[i] = 0.0 ;
    }


    tet_pickup_coordinate_4vertex( ne, crd, x, y, z ) ;
    D = tet_Volume6( x, y, z ) ;
    tet_SideLength( x, y, z, l ) ;
    tet_Center( x, y, z, tetg ) ;
    tetNedelec_Direction( x, y, z, si ) ;
    tetNedelec_mk( x, y, z, xi, yi, zi ) ;
    tetNedelec_XYZ( l, si, xi, yi, zi, X, Y, Z ) ;
    tetNedelec_km_mk( x, y, z, xy, yz, zx ) ;


    nu = data.make_nu( nl, flg, inie, D, X, Y, Z, mem_nu, *mtrl ) ;


    tetNedelec_rotrot( D, X, Y, Z, rr ) ;
    for( i=0 ; i<mp_elem ; i++ ) {
        int    ii = i * szRow ;
        for( j=0 ; j<mp_elem ; j++,ii++ ) {
            ae[ii] += nu * rr[i][j] ;
        }
    }


    if( rs == 2 ) {
        double ss[mp_elem][mp_elem] ;
        double be_ss[mp_elem] ;

        double ee = sigma / delta_t ;

        tetNedelec_sideside( D, l, tetg, si, x, y, z, xi, yi, zi,
                             xy, yz, zx, ss ) ;
        for( i=0 ; i<mp_elem ; i++ ) {
            int    ii = i * szRow ;
            for( j=0 ; j<mp_elem ; j++,ii++ ) {
                ae[ii] += ee * ss[i][j] ;
            }
        }


        tetNedelec_be_sideside( inie, D, l, tetg, si, x, y, z,
                                xi, yi, zi, xy, yz, zx, be_ss ) ;
        for( i=0 ; i<mp_elem ; i++ ) {
            be[i] += ee * be_ss[i] ;
        }
    }


    if( nc != -1 ) {
        double b[ap_elem], c[ap_elem], d[ap_elem] ;
        double Fx[mp_elem], Fy[mp_elem], Fz[mp_elem] ;

        double Joe[ap_elem][dimension] ;
        double Ihe[ap_elem] ;

        double be_JoA[mp_elem] ;
        double be_gs[mp_elem] ;


        for( i=0 ; i<ap_elem ; i++ ) {
            int    ii = coil[nc].enf[ne[i]] * dimension ;
            for( j=0 ; j<dimension ; j++,ii++ ) {
                Joe[i][j] = coil[nc].Jo[ii] ;
            }
            Ihe[i] = coil[nc].Ihr[coil[nc].enf[ne[i]]] ;
        }


        tet_simple_b( x, y, z, b ) ;
        tet_simple_c( x, y, z, c ) ;
        tet_simple_d( x, y, z, d ) ;
        tetNedelec_Fxyz( l, si, tetg, xi, yi, zi, xy, yz, zx, Fx, Fy, Fz ) ;
        tetNedelec_JoA( Joe, si, l, b, c, d, be_JoA ) ;
        tetNedelec_be_gradside( Ihe, D, b, c, d, Fx, Fy, Fz, be_gs ) ;
        for( i=0 ; i<mp_elem ; i++ ) {
            be[i] += be_JoA[i] ;
            be[i] -= be_gs[i] ;
        }
    }


    if( np != -1 ) {
        double a[ap_elem], b[ap_elem], c[ap_elem], d[ap_elem] ;
        double Me[ap_elem][dimension] ;
        double be_MrotA[mp_elem] ;

        tet_simple_a( x, y, z, a ) ;
        tet_simple_b( x, y, z, b ) ;
        tet_simple_c( x, y, z, c ) ;
        tet_simple_d( x, y, z, d ) ;

        for( i=0 ; i<ap_elem ; i++ ) {
            int    ii = pmagnet[np].enf[ne[i]] * dimension ;
            for( j=0 ; j<dimension ; j++,ii++ ) Me[i][j] = pmagnet[np].Jo[ii] ;
        }
        tetNedelec_MrotA( Me, D, si, l, tetg, a, b, c, d, xi, yi, zi, be_MrotA ) ;
        for( i=0 ; i<mp_elem ; i++ ) be[i] += nu * be_MrotA[i] ;
    }


//#if 0
//    {
//        int iii,jjj;
//        fprintf(stdout,"\n%d\n",el);
//        for(iii=0;iii<szRow;iii++){
//            for(jjj=0;jjj<szRow;jjj++) fprintf(stdout," % le",ae[iii*szRow+jjj]);
//            fprintf(stdout,"\n");
//        }
//        for(iii=0;iii<szRow;iii++) fprintf(stdout," % le\n",be[iii]);
//    }
//#endif


//#ifdef __ADVMAG_NEW_SOLVER__
//    *data.solvmat = data.sfunc->_Set( *data.solvmat,
//                                      szRow, indRow, szRow, indRow, ae, be ) ;
//#else
//    data.SetMat( data.smat, szRow, indRow, szRow, indRow, ae, be ) ;
//#endif


    return nu ;
}
