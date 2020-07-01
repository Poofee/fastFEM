#include "tetralheral.h"

#include "mag3dstatic.h"

#include <math.h>

Mag3DStatic::Mag3DStatic()
{

}

/*!
 \brief Function: make element
        conventional piecewise linear tetrahedral
        Non-linear Magnetostatic Analysis
        Newton Method
        A-p method
        计算得到一个10 乘 10 的单元矩阵

 \param el
*/
void Mag3DStatic::tet_MakeElement_NL_Mag_Static_Newton_A_p(int el)
{
    int    i, j, k ;

    int*   nop = data.nop ;
    double* crd = data.crd ;
    int*   flag = data.flag ;
    int*   apid2phiid = data.apid2phiid ;
    double* x_d = data.x_d ;
    double* nl_nu = data.nl_nu ;
    COIL*  pmagnet = data.pmagnet ;
    MTRL*  mtrl = data.mtrl ;
    double weight = data.weight ;

    int    *ne ;
    int    flg = flag[el] ;
    double nu = 0.0 ;
    int    nl = -1 ;
    int    np = -1 ;
    double x[ap_elem], y[ap_elem], z[ap_elem] ;
    double D ;
    double tetg[ap_elem] ;
    double l[mp_elem] ;
    int    si[mp_elem] ;
    double xi[mp_elem], yi[mp_elem], zi[mp_elem] ;
    double xy[mp_elem], yz[mp_elem], zx[mp_elem] ;
    double X[mp_elem], Y[mp_elem], Z[mp_elem] ;
    double Fx[mp_elem], Fy[mp_elem], Fz[mp_elem] ;
    double b[ap_elem], c[ap_elem], d[ap_elem] ;
    double rr[mp_elem][mp_elem], gs[mp_elem][ap_elem], sg[ap_elem][mp_elem] ;

#ifdef szRow_def
#undef szRow_def
#endif
#define szRow_def nd_elem
    int    szRow = szRow_def ;
    int    indRow[szRow_def] ;

    double rotA[dimension] ;
    double B2 = 0.0 ;

    double ae[szRow_def*szRow_def] ;
    double fe[szRow_def] ;


    ne = &(nop[el*nd_elem]) ;
    /** 6个棱的编号 **/
    for( i=0,j=ap_elem ; i<mp_elem ; i++,j++ )
        indRow[i] = ne[j] ;
    /** 4个节点的编号 **/
    for( i=mp_elem,j=0 ; i<szRow ; i++,j++ )
        indRow[i] = apid2phiid[ne[j]] ;

    /** 初始化 **/
    for( i=0,k=0 ; i<szRow ; i++ ) {
        for( j=0 ; j<szRow ; j++,k++ ) {
            ae[k] = 0.0 ;
        }
        fe[i] = 0.0 ;
    }


    tet_pickup_coordinate_4vertex( ne, crd, x, y, z ) ;/** 单元坐标 **/
    D = tet_Volume6( x, y, z ) ;/** 计算体积 **/
    tet_Center( x, y, z, tetg ) ;/** 计算中心 **/
    tet_SideLength( x, y, z, l ) ;/** 计算棱长 **/
    tetNedelec_Direction( x, y, z, si ) ;/** 计算棱的方向 **/
    tetNedelec_mk( x, y, z, xi, yi, zi ) ;/**  **/
    tetNedelec_XYZ( l, si, xi, yi, zi, X, Y, Z ) ;/**  **/
    tetNedelec_km_mk( x, y, z, xy, yz, zx ) ;/**  **/
    tetNedelec_Fxyz( l, si, tetg, xi, yi, zi, xy, yz, zx, Fx, Fy, Fz ) ;/**  **/
    tet_simple_b( x, y, z, b ) ;/**  **/
    tet_simple_c( x, y, z, c ) ;/**  **/
    tet_simple_d( x, y, z, d ) ;/**  **/

    /** 材料编号？ **/
    for( i=0 ; i<mtrl->nnl ; i++ ) {
        if( flg == mtrl->nl_i[i] ) {
            nl = i ;
            break ;
        }
    }
    /**  **/
    for( i=0 ; i<mtrl->npmagnet ; i++ ) {
        if( flg == mtrl->pmagnet_def[i].flag ) {
            np = i ;
//            if( (mtrl->pmagnet_def[i].mode == MTRL_DEF_NL_RF_NUMBER)
//                    || (mtrl->pmagnet_def[i].mode == MTRL_DEF_NL_MD_NUMBER) ) {
//                /* Nothing to do */
//            } else {
//                np = -1 ;
//            }
            break ;
        }
    }


    /** 计算非线性部分 **/
    if( nl > -1 ) {
        double dnu = 0.0 ;
        double xe[mp_elem] ;
        for( i=0,j=ap_elem ; i<mp_elem ; i++,j++ ) {
            xe[i] = x_d[ne[j]] ;
        }
        /* |B[T]| */
        for( i=0 ; i<dimension ; i++ ) {
            rotA[i] = 0.0 ;
        }
        for( i=0 ; i<mp_elem ; i++ ) {
            rotA[0] += xe[i] * X[i] ;
            rotA[1] += xe[i] * Y[i] ;
            rotA[2] += xe[i] * Z[i] ;
        }
        {
            double dd =2.0 / D ;
            for( i=0 ; i<dimension ; i++ ) {
                rotA[i] *= dd ;
            }
        }
        B2 = rotA[0]*rotA[0] + rotA[1]*rotA[1] + rotA[2]*rotA[2] ;
        make_nu( el, nl, &nu, &dnu, B2, nl_nu, *mtrl, weight ) ;

        {
            double sr[mp_elem] ;
            for( i=0 ; i<mp_elem ; i++ ) {
                sr[i] = 0.0 ;
            }
            make_SR( ae, sr, xi, yi, zi, si, l, D, xe, B2, dnu ) ;
            for( i=0 ; i<mp_elem ; i++ ) {
                fe[i] += sr[i] ;
            }
        }


    } else {/** 线性部分 **/
        for( i=0 ; i<mtrl->nflag ; i++ ) {
            if( flg == mtrl->flag_i[i] ) {
                nu = mtrl->flag_nu[i] ;
                break ;
            }
        }

        for( i=0 ; i<mtrl->npmagnet ; i++ ) {
            if( flg == mtrl->pmagnet_def[i].flag ) {
//                if( mtrl->pmagnet_def[i].mode == MTRL_DEF_AL_NUMBER ) {
//                    nu = mtrl->pmagnet_def[i].ialpha ;
//                }
                break ;
            }
        }
    }


    /** 计算电流和永磁 **/
    if( np != -1 ) {
        double a[ap_elem] ;
        double Me[ap_elem][dimension] ;
        double be_MrotA[mp_elem] ;

        tet_simple_a( x, y, z, a ) ;

        for( i=0 ; i<ap_elem ; i++ ) {
            int    ii = pmagnet[np].enf[ne[i]] * dimension ;
            for( j=0 ; j<dimension ; j++,ii++ )
                Me[i][j] = pmagnet[np].Jo[ii] ;
        }
        tetNedelec_MrotA( Me, D, si, l, tetg, a, b, c, d, xi, yi, zi, be_MrotA ) ;
        for( i=0 ; i<mp_elem ; i++ )
            fe[i] += nu * be_MrotA[i] ;
    }


    /** 计算双旋度矩阵 **/
    tetNedelec_rotrot( D, X, Y, Z, rr ) ;
    for( i=0 ; i<mp_elem ; i++ ) {
        int    ii = i * szRow ;
        for( j=0 ; j<mp_elem ; j++,ii++ ) {
            ae[ii] += nu * rr[i][j] ;
        }
    }


    /** 计算拉格朗日算子矩阵 **/
    tetNedelec_gradside( D, b, c, d, Fx, Fy, Fz, gs ) ;
    for( i=0 ; i<mp_elem ; i++ ) {
        int    ii = i*szRow + mp_elem ;
        for( j=mp_elem,k=0 ; j<szRow ; j++,k++,ii++)
            ae[ii] += gs[i][k] ;
    }
    for( i=0 ; i<mp_elem ; i++ ) {
        for( j=0 ; j<ap_elem ; j++ ) {
            sg[j][i] = gs[i][j] ;
        }
    }
    for( i=mp_elem,k=0 ; i<szRow ; i++,k++ ) {
        int    ii = i * szRow ;
        for( j=0 ; j<mp_elem ; j++,ii++ )
            ae[ii] += sg[k][j] ;
    }
    /*for( i=mp_elem,k=0 ; i<szRow ; i++,k++) {
          int    ii = i * szRow ;
          for( j=0 ; j<mp_elem ; j++,ii++ ) ae[ii] += gs[j][k] ;
        }*/


    //#ifdef __ADVMAG_NEW_SOLVER__
    //    *data.solvmat = data.sfunc->_Set( *data.solvmat,
    //                                      szRow, indRow, szRow, indRow, ae, fe ) ;
    //#else
    //    data.SetMat( data.smat, szRow, indRow, szRow, indRow, ae, fe ) ;
    //#endif
}

/*!
 \brief tet_MakeElement_NL_Mag_Static_IP_1_nuB2
        Function: make element
        make nu
        conventional piecewise linear tetrahedral
        Non-linear Magnetostatic Analysis
        first order interpolation

 \param el
 \param nl
 \param nu
 \param dnu
 \param B2
 \param nl_nu
 \param mtrl
 \param weight
*/
void Mag3DStatic::make_nu(int el, int nl, double *nu, double *dnu, double B2, double *nl_nu, MTRL mtrl, double weight)
{
    /*+  +*/
    int    n = 0 ;
    double* a = nullptr ;
    double* b = nullptr ;
    double* v = nullptr ;

    MakeElement_IP_Get_nuB2( nl, &n, &a, &b, &v, mtrl ) ;
    tet_MakeElement_NL_Mag_Static_IP_1_sub( nu, dnu, B2,
                                            n, a, b, v, mtrl ) ;
}

/*!
 \brief tet_MakeElement_NL_Mag_Static_Newton_A_SR_IP_1
        Function: make element
        conventional piecewise linear tetrahedral
        Non-linear Magnetostatic Analysis
        Newton Method
        A method
        make side-rot
        first order interpolation

 \param ae
 \param fe
 \param xi
 \param yi
 \param zi
 \param si
 \param l
 \param D
 \param xe
 \param B2
 \param dnu
*/
void Mag3DStatic::make_SR(double *ae, double *fe, double *xi, double *yi, double *zi, int *si, double *l, double D, double *xe, double B2, double dnu)
{
    /*+  +*/
    int    i, j, k ;

    double dtmp[mp_elem] ;
    double sr[mp_elem][mp_elem] ;
    double rD = 1.0 / D ;
    double rB = 1.0 / sqrt( B2 ) ;


    for( i=0 ; i<mp_elem ; i++ ) {
        dtmp[i] = 0.0 ;
        for( j=0 ; j<mp_elem ; j++ ) {
            dtmp[i] += (xi[i]*xi[j]  + yi[i]*yi[j] + zi[i]*zi[j])
                    * si[i] * si[j] * l[i] * l[j] * 2.0 / 3.0 * rD * xe[j] ;
            sr[i][j] = 0.0 ;
        }
    }
    for( i=0 ; i<mp_elem ; i++ ) {
        for( j=0 ; j<mp_elem ; j++ ) {
            sr[i][j] = 6.0 * dnu * rD * dtmp[i] * dtmp[j] * rB ;
        }
    }
    for( i=0,k=0 ; i<mp_elem ; i++ ) {
        for( j=0 ; j<mp_elem ; j++,k++ ) {
            ae[k] += sr[i][j] ;
            fe[i] += sr[i][j] * xe[j] ;
        }
    }
}


/*!
 \brief Function: make element
        make nu
        get nu-B2 curve
        conventional piecewise linear tetrahedral
        Non-linear Magnetostatic Analysis

 \param nl
 \param n
 \param a
 \param b
 \param v
 \param mtrl
*/
void Mag3DStatic::MakeElement_IP_Get_nuB2(int nl, int *n, double **a, double **b, double **v, MTRL mtrl)
{
    (*n) = mtrl.nl_mag[nl].n ;
    (*a) = mtrl.nl_mag[nl].B2 ;
    (*b) = mtrl.nl_mag[nl].nu ;
    (*v) = mtrl.nl_mag[nl].v_newton ;
}


/*!
 \brief Function: make element
        make nu
        conventional piecewise linear tetrahedral
        Non-linear Magnetostatic Analysis
        first order interpolation

 \param nu
 \param dnu
 \param B2
 \param n
 \param a
 \param b
 \param v
 \param mtrl
*/
void Mag3DStatic::tet_MakeElement_NL_Mag_Static_IP_1_sub(double *nu, double *dnu, double B2, int n, double *a, double *b, double *v, MTRL mtrl)
{
    int    i ;


    {
        int    ii = n - 1 ;
        if( a[ii] <= B2 ) {
            int    jj = ii - 1 ;
            (*dnu) = (b[ii] - b[jj]) / (a[ii] - a[jj]) ;
            (*nu) = b[ii] + (*dnu)*(B2 - a[ii]) ;
        } else if( B2 < a[0] ) {
            (*dnu) = (b[1] - b[0]) / (a[1] - a[0]) ;
            (*nu) = b[0] + (*dnu)*(B2 - a[0]) ;
            /*} else if( B2 < a[1] ) {
                  (*dnu) = (b[2] - b[1]) / (a[2] - a[1]) ;
                  (*nu) = b[1] + (*dnu)*(B2 - a[1]) ;
                } else {
                  for( i=2 ; i<n ; i++ ) {*/
        } else {
            for( i=1 ; i<n ; i++ ) {
                if( B2 < a[i] ) {
                    int    jj = i - 1 ;
                    (*dnu) = (b[i] - b[jj]) / (a[i] - a[jj]) ;
                    (*nu) = b[jj] + (*dnu)*(B2 - a[jj]) ;
                    break ;
                }
            }
        }
    }
}
