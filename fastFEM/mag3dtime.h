#ifndef MAG3DTIME_H
#define MAG3DTIME_H

#include "libstructure.h"

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
//#ifdef __ADVMAG_NEW_SOLVER__
//  void **solvmat ;
//  SOLVER_FUNC *sfunc ;
//#else
//  SOLVER_MAT *smat ;
//  void (*SetMat)( SOLVER_MAT *smat, int szRow, int *indRow,
//                  int szCol, int *indCol, double *ae, double *fe ) ;
//#endif
} DATA_NS_EDDY ;

class Mag3DTime
{
public:
    Mag3DTime();
    ~Mag3DTime();

    void doMesh();
    void reMeshRotate(double theta);
    void reMeshParallelShift(double dx,double dy,double dz);
    void assemble();
    void solve();
    void run();
    double tet_MakeElement_NS_Eddy_A(int el);

private:
    DATA_NS_EDDY data;
};

#endif // MAG3DTIME_H
