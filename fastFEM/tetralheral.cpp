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

#include <math.h>

/*+ Function: pickup coordinates of 4 vertexes of tetrahedral +*/
void tet_pickup_coordinate_4vertex( int *ne, double *crd,
                                    double *x, double *y, double *z )
{
    int    i ;

    for( i=0 ; i<4 ; i++ ) {
        size_t ii = (size_t)ne[i] * (size_t)3 ;
        x[i] = crd[ii] ;
        y[i] = crd[++ii] ;
        z[i] = crd[++ii] ;
    }
}

/*+ Function: calculate lengths of 6 sides of tetrahedral +*/
void tet_SideLength( double *x, double *y, double *z, double *l )
{
    int    i ;

    for( i=0 ; i<6 ; i++ ){
        int   ii = tet_n[i][0] ;
        int   jj = tet_n[i][1] ;
        l[i] = sqrt( pow( x[ii] - x[jj], 2)
                     +pow( y[ii] - y[jj], 2)
                     +pow( z[ii] - z[jj], 2) ) ;
    }
}

/*+ Function: calculate 6 times volume of tetrahedral +*/
double tet_Volume6( double *x, double *y, double *z )
{
    double D ;

    D = - x[0]*y[1]*z[2] + x[0]*y[1]*z[3] + x[0]*y[2]*z[1]
            - x[0]*y[2]*z[3] - x[0]*y[3]*z[1] + x[0]*y[3]*z[2]
            + x[1]*y[0]*z[2] - x[1]*y[0]*z[3] - x[1]*y[2]*z[0]
            + x[1]*y[2]*z[3] + x[1]*y[3]*z[0] - x[1]*y[3]*z[2]
            - x[2]*y[0]*z[1] + x[2]*y[0]*z[3] + x[2]*y[1]*z[0]
            - x[2]*y[1]*z[3] - x[2]*y[3]*z[0] + x[2]*y[3]*z[1]
            + x[3]*y[0]*z[1] - x[3]*y[0]*z[2] - x[3]*y[1]*z[0]
            + x[3]*y[1]*z[2] + x[3]*y[2]*z[0] - x[3]*y[2]*z[1] ;


    return D ;
}

/*+ Function: calculate volume of tetrahedral +*/
double tet_Volume( double *x, double *y, double *z )
{
    double V ;


    V = tet_Volume6( x, y, z ) / 6.0 ;


    return V  ;
}

/*+ Function: calculate center of gravity of tetrahedral +*/
void tet_Center( double *x, double *y, double *z, double *g )
{
    g[0] = ( x[0] + x[1] + x[2] + x[3] ) * 0.25 ;
    g[1] = ( y[0] + y[1] + y[2] + y[3] ) * 0.25 ;
    g[2] = ( z[0] + z[1] + z[2] + z[3] ) * 0.25 ;
}

/*+ Function: calculate determinant `a'
                  of conventional piecewise linear tetrahedral +*/
/*+     | xj yj zj |
  a = | xk yk zk |
      | xm ym zm | +*/
void tet_simple_a( double *x, double *y, double *z, double *a )
{
    int i;

    for( i=0 ; i<4 ; i++ ) {
        int    ii = tet_m[i][1] ;
        int    jj = tet_m[i][2] ;
        int    kk = tet_m[i][3] ;
        a[i] = x[ii]*(y[jj]*z[kk] - y[kk]*z[jj])
                +x[jj]*(y[kk]*z[ii] - y[ii]*z[kk])
                +x[kk]*(y[ii]*z[jj] - y[jj]*z[ii]) ;
    }
}


/*+ Function: calculate determinant `b'
                  of conventional piecewise linear tetrahedral +*/
/*+       | 1 yj zj |
  b = - | 1 yk zk |
        | 1 ym zm | +*/
void tet_simple_b( double *x, double *y, double *z, double *b )
{

    int i;


    for( i=0 ; i<4 ; i++ ) {
        int    ii = tet_m[i][1] ;
        int    jj = tet_m[i][2] ;
        int    kk = tet_m[i][3] ;
        b[i] = -1.0*( y[ii]*(z[jj] - z[kk])
                      +y[jj]*(z[kk] - z[ii])
                      +y[kk]*(z[ii] - z[jj]) ) ;
    }
}


/*+ Function: calculate determinant `c'
                  of conventional piecewise linear tetrahedral +*/
/*+       | xj 1 zj |
  c = - | xk 1 zk |
        | xm 1 zm | +*/
void tet_simple_c( double *x, double *y, double *z, double *c )
{

    int i;


    for( i=0 ; i<4 ; i++ ) {
        int    ii = tet_m[i][1] ;
        int    jj = tet_m[i][2] ;
        int    kk = tet_m[i][3] ;
        c[i] = -1.0*( z[ii]*(x[jj] - x[kk])
                      +z[jj]*(x[kk] - x[ii])
                      +z[kk]*(x[ii] - x[jj]) ) ;
    }
}


/*+ Function: calculate determinant `d'
                  of conventional piecewise linear tetrahedral +*/
/*+       | xj yj 1 |
  d = - | xk yk 1 |
        | xm ym 1 | +*/
void tet_simple_d( double *x, double *y, double *z, double *d )
{

    int i;


    for( i=0 ; i<4 ; i++ ) {
        int    ii = tet_m[i][1] ;
        int    jj = tet_m[i][2] ;
        int    kk = tet_m[i][3] ;
        d[i] = -1.0*( x[ii]*(y[jj] - y[kk])
                      +x[jj]*(y[kk] - y[ii])
                      +x[kk]*(y[ii] - y[jj]) ) ;
    }
}


/*+ Function: calculate determinant `abcd'
                  of conventional piecewise linear tetrahedral +*/
/*+ abcd[0][0] = a[0]ux[0]+a[1]ux[1]+a[2]ux[2]+a[3]ux[3]
  abcd[0][1] = a[0]uy[0]+a[1]uy[1]+a[2]uy[2]+a[3]uy[3]
  abcd[0][2] = a[0]uz[0]+a[1]uz[1]+a[2]uz[2]+a[3]uz[3] +*/
/*+ abcd[1][0] = b[0]ux[0]+b[1]ux[1]+b[2]ux[2]+b[3]ux[3]
  abcd[1][1] = b[0]uy[0]+b[1]uy[1]+b[2]uy[2]+b[3]uy[3]
  abcd[1][2] = b[0]uz[0]+b[1]uz[1]+b[2]uz[2]+b[3]uz[3] +*/
/*+ abcd[2][0] = c[0]ux[0]+c[1]ux[1]+c[2]ux[2]+c[3]ux[3]
  abcd[2][1] = c[0]uy[0]+c[1]uy[1]+c[2]uy[2]+c[3]uy[3]
  abcd[2][2] = c[0]uz[0]+c[1]uz[1]+c[2]uz[2]+c[3]uz[3] +*/
/*+ abcd[3][0] = d[0]ux[0]+d[1]ux[1]+d[2]ux[2]+d[3]ux[3]
  abcd[3][1] = d[0]uy[0]+d[1]uy[1]+d[2]uy[2]+d[3]uy[3]
  abcd[3][2] = d[0]uz[0]+d[1]uz[1]+d[2]uz[2]+d[3]uz[3] +*/
void tet_simple_abcd( double *a, double *b, double *c, double *d,
                      double ue[][3], double abcd[][3] )
{

    int i, j ;


    for( j=0 ; j<3 ; j++ ) {
        abcd[0][j] = a[0] * ue[0][j] ;
        abcd[1][j] = b[0] * ue[0][j] ;
        abcd[2][j] = c[0] * ue[0][j] ;
        abcd[3][j] = d[0] * ue[0][j] ;
    }
    for( i=1 ; i<4 ; i++ ) {
        for( j=0 ; j<3 ; j++ ) {
            abcd[0][j] += a[i] * ue[i][j] ;
            abcd[1][j] += b[i] * ue[i][j] ;
            abcd[2][j] += c[i] * ue[i][j] ;
            abcd[3][j] += d[i] * ue[i][j] ;
        }
    }
}


    /*+ Function: calculate (grad u, grad v) of tetrahedral +*/
void tet_simple_gradgrad( double D, double *b, double *c, double *d,
                          double gg[][4] )
{

    int    i, j ;

    double dd = 1.0 / (6.0 * D) ;


    for( i=0 ; i<4 ; i++ ) {
        for( j=0 ; j<4 ; j++ ) {
            gg[i][j] = (b[i]*b[j] + c[i]*c[j] + d[i]*d[j]) * dd ;
        }
    }
}


/*+ Function: calculate (u(ap), grad v) of tetrahedral for be +*/
void tet_simple_be_apgrad( double D, double *b, double *c, double *d,
                           double u[][3], double *ag )
{

    int    i, j ;

    double dd = 1.0 / 24.0 ;


    for( i=0 ; i<4 ; i++ ) {
        ag[i] = 0.0 ;
        for( j=0 ; j<4 ; j++ ) {
            ag[i] += (b[i]*u[j][0] + c[i]*u[j][1] + d[i]*u[j][2]) * dd ;
        }
    }
}


/*+ Function: determine direction of Nedelec by coordinates +*/
void tetNedelec_Direction( double *x, double *y, double *z, int *si )
{

    int    i ;

    double pp ;


    for( i=0 ; i<6 ; i++ ){
        /* try x-define */
        int    ii = tet_n[i][1] ;
        int    jj = tet_n[i][0] ;
        pp = x[ii] - x[jj] ;
        if( pp > 0.0 ) {
            si[i] =  1 ;
        } else if( pp < 0.0 ) {
            si[i] = -1 ;
        } else {
            /* try y-define */
            pp = y[ii] - y[jj] ;
            if( pp > 0.0 ) {
                si[i] =  1 ;
            } else if( pp < 0.0 ) {
                si[i] = -1 ;
            } else {
                /* try z-define */
                pp = z[ii] - z[jj] ;
                if( pp > 0.0 ) {
                    si[i] =  1 ;
                } else if( pp < 0.0 ) {
                    si[i] = -1 ;
                } else {
                    /* error */
//                    emessage( 2050, __ABORT_LIBFEM__, NULL ) ;
                }
            }
        }
    }
}


/*+ Function: calculate xm-xk, ym-yk, zm-zk of Nedelec +*/
void tetNedelec_mk( double *x, double *y, double *z,
                    double *xi, double *yi, double *zi )
{

    int    i ;


    for( i=0 ; i<6 ; i++ ) {
        int    ii = tet_n[i][2] ;
        int    jj = tet_n[i][3] ;
        xi[i] = x[jj] - x[ii] ;
        yi[i] = y[jj] - y[ii] ;
        zi[i] = z[jj] - z[ii] ;
    }
}


/*+ Function: calculate X = l*si*xi,
                      Y = l*si*yi,
                      Z = l*si*zi of Nedelec +*/
void tetNedelec_XYZ( double *l, int *si,
                     double *xi, double *yi, double *zi,
                     double *X, double *Y, double *Z )
{

    int    i ;


    for( i=0 ; i<6 ; i++ ) {
        double d = l[i] * (double)si[i] ;
        X[i] = d * xi[i] ;
        Y[i] = d * yi[i] ;
        Z[i] = d * zi[i] ;
    }
}


/*+ Function: calculate xkym-xmyk, ykzm-ymzk, zkxm-zmxk of Nedelec +*/
void tetNedelec_km_mk( double *x, double *y, double *z,
                       double *xy, double *yz, double *zx )
{

    int    i ;


    for( i=0 ; i<6 ; i++ ) {
        int    ii = tet_n[i][2] ;
        int    jj = tet_n[i][3] ;
        xy[i] = x[ii]*y[jj] - x[jj]*y[ii] ;
        yz[i] = y[ii]*z[jj] - y[jj]*z[ii] ;
        zx[i] = z[ii]*x[jj] - z[jj]*x[ii] ;
    }
}


/*+ Function: calculate Fx = l*si*(yz-zi*gy+yi*gz),
                      Fy = l*si*(zx-xi*gz+zi*gx),
                      Fz = l*si*(xy-yi*gx+xi*gy) of Nedelec +*/
void tetNedelec_Fxyz( double *l, int *si, double *g,
                      double *xi, double *yi, double *zi,
                      double *xy, double *yz, double *zx,
                      double *Fx, double *Fy, double *Fz )
{

    int    i ;


    for( i=0 ; i<6 ; i++ ) {
        double d = l[i] * (double)si[i] ;
        Fx[i] = d * (yz[i] - zi[i]*g[1] + yi[i]*g[2]) ;
        Fy[i] = d * (zx[i] - xi[i]*g[2] + zi[i]*g[0]) ;
        Fz[i] = d * (xy[i] - yi[i]*g[0] + xi[i]*g[1]) ;
    }
}


/*+ Function: calculate C [n] = abcd[0][0]X[n]
                             +abcd[0][1]Y[n]+abcd[0][2]Z[n],
                      Cx[n] = abcd[1][0]X[n]
                             +abcd[1][1]Y[n]+abcd[1][2]Z[n],
                      Cy[n] = abcd[2][0]X[n]
                             +abcd[2][1]Y[n]+abcd[2][2]Z[n],
                      Cz[n] = abcd[3][0]X[n]
                             +abcd[3][1]Y[n]+abcd[3][2]Z[n],
                                                         of Nedelec +*/
void tetNedelec_C( double abcd[][3], double *X, double *Y, double *Z,
double C[][6] )
{

    int    i, j ;


    for( i=0 ; i<4 ; i++ ) {
        for( j=0 ; j<6 ; j++ ) {
            C[i][j] = abcd[i][0]*X[j] + abcd[i][1]*Y[j] + abcd[i][2]*Z[j] ;
        }
    }
}


/*+ Function: calculate rot u of Nedelec +*/
void tetNedelec_rot( double D, double *X, double *Y, double *Z,
                     double *u, double *rot )
{

    int    i ;

    double d = 2.0 / D ;


    for( i=0 ; i<3 ; i++ ) {
        rot[i] = 0.0 ;
    }
    for( i=0 ; i<6 ; i++ ) {
        rot[0] += u[i] * X[i] ;
        rot[1] += u[i] * Y[i] ;
        rot[2] += u[i] * Z[i] ;
    }
    for( i=0 ; i<3 ; i++ ) {
        rot[i] *= d ;
    }
}


/*+ Function: calculate (rot u, rot v) of Nedelec +*/
void tetNedelec_rotrot( double D, double *X, double *Y, double *Z,
                        double rr[][6] )
{

    int    i, j ;

    double d = 2.0 / (3.0 * D) ;


    for( i=0 ; i<6 ; i++ ) {
        for( j=0 ; j<6 ; j++ ) {
            rr[i][j] = d * (X[i]*X[j] + Y[i]*Y[j] + Z[i]*Z[j]) ;
        }
    }
}


/*+ Function: calculate (u(side), v(side)) of Nedelec +*/
void tetNedelec_sideside( double D, double *l, double *tetg, int *si,
                          double *x, double *y, double *z,
                          double *xi, double *yi, double *zi,
                          double *xy, double *yz, double *zx,
                          double ss[][6] )
{

    int    i, j ;

    double d = 1.0 / (6.0 * D) ;
    double gxx[6] ;
    double dx = 0.0, dy = 0.0, dz = 0.0 ;
    double lsi[6] ;


    for( i=0 ; i<4 ; i++ ) {
        dx += x[i] ;
        dy += y[i] ;
        dz += z[i] ;
    }
    gxx[0] = dx * dx ;
    gxx[1] = dy * dy ;
    gxx[2] = dz * dz ;
    gxx[3] = dx * dy ;
    gxx[4] = dy * dz ;
    gxx[5] = dz * dx ;
    for( i=0 ; i<4 ; i++ ) {
        gxx[0] += x[i] * x[i] ;
        gxx[1] += y[i] * y[i] ;
        gxx[2] += z[i] * z[i] ;
        gxx[3] += x[i] * y[i] ;
        gxx[4] += y[i] * z[i] ;
        gxx[5] += z[i] * x[i] ;
    }
    for( i=0 ; i<6 ; i++ ) {
        gxx[i] *= 0.05 ;
    }


    for( i=0 ; i<6 ; i++ ) {
        lsi[i] = l[i] * (double)si[i] ;
    }


    for( i=0 ; i<6 ; i++ ) {
        for( j=0 ; j<6 ; j++ ) {
            ss[i][j] = lsi[i] * lsi[j] * d
                    *(-zi[i]*(-zi[j]*gxx[1]  + yi[j]*gxx[4]  + yz[j]*tetg[1])
                    +yi[i]*(-zi[j]*gxx[4]  + yi[j]*gxx[2]  + yz[j]*tetg[2])
                    +yz[i]*(-zi[j]*tetg[1] + yi[j]*tetg[2] + yz[j])
                    -xi[i]*(-xi[j]*gxx[2]  + zi[j]*gxx[5]  + zx[j]*tetg[2])
                    +zi[i]*(-xi[j]*gxx[5]  + zi[j]*gxx[0]  + zx[j]*tetg[0])
                    +zx[i]*(-xi[j]*tetg[2] + zi[j]*tetg[0] + zx[j])
                    -yi[i]*(-yi[j]*gxx[0]  + xi[j]*gxx[3]  + xy[j]*tetg[0])
                    +xi[i]*(-yi[j]*gxx[3]  + xi[j]*gxx[1]  + xy[j]*tetg[1])
                    +xy[i]*(-yi[j]*tetg[0] + xi[j]*tetg[1] + xy[j])) ;
        }
    }
}


/*+ Function: calculate (u(side), v(side)) of Nedelec +*/
void tetNedelec_be_sideside( double *u,
                             double D, double *l, double *tetg, int *si,
                             double *x, double *y, double *z,
                             double *xi, double *yi, double *zi,
                             double *xy, double *yz, double *zx,
                             double *be_ss )
{

    int    i, j ;

    double ss[6][6] ;


    tetNedelec_sideside( D, l, tetg, si, x, y, z, xi, yi, zi,
                         xy, yz, zx, ss ) ;
    for( i=0 ; i<6 ; i++ ) {
        be_ss[i] = 0.0 ;
        for( j=0 ; j<6 ; j++ ) {
            be_ss[i] += u[j] * ss[i][j] ;
        }
    }
}


/*+ Function: calculate (u(element), rot v) of Nedelec for be +*/
void tetNedelec_be_elrot( double *u, double *X, double *Y, double *Z,
                          double *er )
{

    int    i ;

    double d = 1.0 / 3.0 ;


    for( i=0 ; i<6 ; i++ ) {
        er[i] = d * (u[0]*X[i] + u[1]*Y[i] + u[2]*Z[i] ) ;
    }
}


/*+ Function: calculate (grad u, v(side)) of Nedelec +*/
void tetNedelec_gradside( double D,
                          double *b, double *c, double *d,
                          double *Fx, double *Fy, double *Fz,
                          double gs[][4] )
{

    int    i, j ;

    double dd = 1.0 / (6.0 * D) ;


    for( i=0 ; i<6 ; i++ ) {
        for( j=0 ; j<4 ; j++ ) {
            gs[i][j] = dd * (Fx[i]*b[j] + Fy[i]*c[j] + Fz[i]*d[j]) ;
        }
    }
}


/*+ Function: calculate (grad u, v(side)) of Nedelec +*/
void tetNedelec_be_gradside( double *u, double D,
                             double *b, double *c, double *d,
                             double *Fx, double *Fy, double *Fz,
                             double *be_gs )
{

    int    i, j ;

    double gs[6][4] ;


    tetNedelec_gradside( D, b, c, d, Fx, Fy, Fz, gs ) ;
    for( i=0 ; i<6 ; i++ ) {
        be_gs[i] = 0.0 ;
        for( j=0 ; j<4 ; j++ ) {
            be_gs[i] += u[j] * gs[i][j] ;
        }
    }
}


/*+ Function: calculate (u(side), grad v) of Nedelec +*/
void tetNedelec_sidegrad( double D,
                          double *b, double *c, double *d,
                          double *Fx, double *Fy, double *Fz,
                          double sg[][6] )
{

    int    i, j ;

    double dd = 1.0 / (6.0 * D) ;


    for( i=0 ; i<4 ; i++ ) {
        for( j=0 ; j<6 ; j++ ) {
            sg[i][j] = dd * (b[i]*Fx[j] + c[i]*Fy[j] + d[i]*Fz[j]) ;
        }
    }
}


/*+ Function: calculate (u(side), grad v) of Nedelec +*/
void tetNedelec_be_sidegrad( double *u, double D,
                             double *b, double *c, double *d,
                             double *Fx, double *Fy, double *Fz,
                             double *be_sg )
{

    int    i, j ;

    double sg[4][6] ;


    tetNedelec_sidegrad( D, b, c, d, Fx, Fy, Fz, sg ) ;
    for( i=0 ; i<4 ; i++ ) {
        be_sg[i] = 0.0 ;
        for( j=0 ; j<6 ; j++ ) {
            be_sg[i] += u[j] * sg[i][j] ;
        }
    }
}


/*+ Function: calculate (Jo(apex*3), A*) of Nedelec +*/
void tetNedelec_JoA( double Joe[][3], int *si, double *l,
                    double *b, double *c, double *d,
                    double *be_JoA )
{

    int    i, j, k ;

    double PN[6][12] ;
    double d0[][4] = { {  2.0,  1.0,  1.0,  1.0},
                       {  2.0,  1.0,  1.0,  1.0},
                       {  2.0,  1.0,  1.0,  1.0},
                       {  1.0,  2.0,  1.0,  1.0},
                       {  1.0,  2.0,  1.0,  1.0},
                       {  1.0,  1.0,  2.0,  1.0} } ;
    double d1[][4] = { { -1.0, -2.0, -1.0, -1.0},
                       { -1.0, -1.0, -2.0, -1.0},
                       { -1.0, -1.0, -1.0, -2.0},
                       { -1.0, -1.0, -2.0, -1.0},
                       { -1.0, -1.0, -1.0, -2.0},
                       { -1.0, -1.0, -1.0, -2.0} } ;


    for( i=0 ; i<6 ; i++ ) {
        double dd = (double)si[i] * l[i] / 120.0 ;
        for( j=0 ; j<12 ; j++ ) {
            PN[i][j] = dd ;
        }
    }

    for( i=0 ; i<6 ; i++ ) {
        for( j=0 ; j<4 ; j++ ) {
            PN[i][j] *= d0[i][j]*b[tet_n[i][1]] + d1[i][j]*b[tet_n[i][0]] ;
        }
    }
    for( i=0 ; i<6 ; i++ ) {
        for( j=0,k=4 ; j<4 ; j++,k++ ) {
            PN[i][k] *= d0[i][j]*c[tet_n[i][1]] + d1[i][j]*c[tet_n[i][0]] ;
        }
    }
    for( i=0 ; i<6 ; i++ ) {
        for( j=0,k=8 ; j<4 ; j++,k++ ) {
            PN[i][k] *= d0[i][j]*d[tet_n[i][1]] + d1[i][j]*d[tet_n[i][0]] ;
        }
    }

    for( i=0 ; i<6 ; i++ ) {
        be_JoA[i] = 0.0 ;
        for( j=0 ; j<4 ; j++ ) {
            be_JoA[i] += PN[i][j] * Joe[j][0] ;
        }
    }
    for( i=0 ; i<6 ; i++ ) {
        for( j=0,k=4 ; j<4 ; j++,k++ ) {
            be_JoA[i] += PN[i][k] * Joe[j][1] ;
        }
    }
    for( i=0 ; i<6 ; i++ ) {
        for( j=0,k=8 ; j<4 ; j++,k++ ) {
            be_JoA[i] += PN[i][k] * Joe[j][2] ;
        }
    }
}



/*+ Function: calculate (M(apex*3), rot A*) of Nedelec +*/
void tetNedelec_MrotA( double Me[][3], double D, int *si, double *l,
                        double *tetg,
                        double *a, double *b, double *c, double *d,
                        double *xi, double *yi, double *zi,
                        double *be_MrotA )
{

    int    i ;

    double da = 0.0 ;
    double db = 0.0 ;
    double dc = 0.0 ;
    double dd = 0.0 ;

    double dx = 0.0 ;
    double dy = 0.0 ;
    double dz = 0.0 ;

    double ee = 1.0 / (D * 3.0) ;


    for( i=0 ; i<4 ; i++ ) {
        da += a[i] * Me[i][0] ;
        db += b[i] * Me[i][0] ;
        dc += c[i] * Me[i][0] ;
        dd += d[i] * Me[i][0] ;
    }
    db *= tetg[0] ;
    dc *= tetg[1] ;
    dd *= tetg[2] ;
    dx = da + db + dc + dd ;

    da = db = dc = dd = 0 ;
    for( i=0 ; i<4 ; i++ ) {
        da += a[i] * Me[i][1] ;
        db += b[i] * Me[i][1] ;
        dc += c[i] * Me[i][1] ;
        dd += d[i] * Me[i][1] ;
    }
    db *= tetg[0] ;
    dc *= tetg[1] ;
    dd *= tetg[2] ;
    dy = da + db + dc + dd ;

    da = db = dc = dd = 0 ;
    for( i=0 ; i<4 ; i++ ) {
        da += a[i] * Me[i][2] ;
        db += b[i] * Me[i][2] ;
        dc += c[i] * Me[i][2] ;
        dd += d[i] * Me[i][2] ;
    }
    db *= tetg[0] ;
    dc *= tetg[1] ;
    dd *= tetg[2] ;
    dz = da + db + dc + dd ;


    for( i=0 ; i<6 ; i++ ) {
        be_MrotA[i] = l[i] * (double)si[i] * ee
                * (dx*xi[i] + dy*yi[i] + dz*zi[i]) ;
    }
}


    /*+ Function: calculate ((rotA.u)u, rot A*) of Nedelec +*/
void tetNedelec_rotAu_u_rotA( double *x, double *y, double *z,
                              double *tetg, double D, double C[][6],
double rotAu_u_rotA[][6] )
{

    int    i, j ;

    double xx = 0.0 ;
    double yy = 0.0 ;
    double zz = 0.0 ;
    double xy = 0.0 ;
    double yz = 0.0 ;
    double zx = 0.0 ;

    double dd = 1.0 / D ;


    dd = 2.0 / 15 * dd * dd * dd ;

    xx = 16.0 * tetg[0] * tetg[0] ;
    yy = 16.0 * tetg[1] * tetg[1] ;
    zz = 16.0 * tetg[2] * tetg[2] ;
    xy = 16.0 * tetg[0] * tetg[1] ;
    yz = 16.0 * tetg[1] * tetg[2] ;
    zx = 16.0 * tetg[2] * tetg[0] ;
    for( i=0 ; i<4 ; i++ ) {
        xx += x[i] * x[i] ;
        yy += y[i] * y[i] ;
        zz += z[i] * z[i] ;
        xy += x[i] * y[i] ;
        yz += y[i] * z[i] ;
        zx += z[i] * x[i] ;
    }
    xx *= 0.25 ;
    yy *= 0.25 ;
    zz *= 0.25 ;
    xy *= 0.25 ;
    yz *= 0.25 ;
    zx *= 0.25 ;


    for( i=0 ; i<6 ; i++ ) {
        for( j=0 ; j<6 ; j++ ) {
            rotAu_u_rotA[i][j] = dd * (5.0*(C[0][i]*C[0][j]
                    + C[0][i]*(C[1][j]+C[2][j]+C[3][j])
                    + (C[1][i]+C[2][i]+C[3][i])*C[0][j])
                    + (C[1][i]*C[2][j]+C[2][i]*C[1][j])*xy
                    + (C[2][i]*C[3][j]+C[3][i]*C[2][j])*yz
                    + (C[3][i]*C[1][j]+C[1][i]*C[3][j])*zx
                    + C[1][i]*C[1][j]*xx
                    + C[2][i]*C[2][j]*yy
                    + C[3][i]*C[3][j]*zz) ;
        }
    }
}
