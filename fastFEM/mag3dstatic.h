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
