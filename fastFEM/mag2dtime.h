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
