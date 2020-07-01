#ifndef LIBSTRUCTURE_H
#define LIBSTRUCTURE_H


#define N_INI 2


/* dirichlet B.C. */
typedef struct {
    int    need ;
    int    nbc ;
    int    *ibc ;
    int    *direction ;
    double *vbc ;
    //		Complex *vbc_c ;
} DirichletBC ;


typedef struct {
    int    need ;
    int    type ;
    int    nface ;
    int    *el ;
    int    *face ;
    int    *direction ;
    double *vbc ;
//    Complex *vbc_c ;
} DirichletBC_EF ;


typedef struct {
    int    n ;
    double ap[5][3] ;
    int    analysis ;
    double Jor[3] ;
    double Joi[3] ;

    double abc[3][3] ;
    double d[2][3] ;
    double tolerance ;
} DEF_FORM_PP ;

typedef struct {
    int    n ;
    int    xyz ; /* 0: parallel to x  1: parallel to y  2: parallel to z */
    double ap[4] ;
    double theta[2] ;
    double r[2] ;
    int    analysis ;
    double Jor ;
    double Joi ;
} DEF_FORM_DSC ;

typedef struct {
    int    n ;
    double base0[3] ;
    double base1[3] ;
    double ref[3] ;
    double theta[2] ;
    double r[2] ;
    int    analysis ;
    double Jor ;
    double Joi ;
} DEF_FORM_DSC2 ;

typedef struct {
    int    npp ;
    DEF_FORM_PP  *pp ;

    int    ndsc ;
    DEF_FORM_DSC *dsc ;

    int    ndsc2 ;
    DEF_FORM_DSC2 *dsc2 ;
} DEF_FORM ;

typedef struct {
    double af ; /* angular frequency */
    double phase ;
    double magnification ;
    double constant ;
} DEF_TIME_SINUSOIDAL ;

typedef struct {
    double s_magnification ;
    double e_magnification ;
    double constant ;
} DEF_TIME_LINEAR ;

typedef struct {
    double t ;

    int    nsin ;
    DEF_TIME_SINUSOIDAL* sin ;

    int    nlin ;
    DEF_TIME_LINEAR* lin ;
} DEF_TIME ;

typedef struct {
//    char    fname[FILENAME_MAX] ;
    int     n ;
    double  *H ;
    double  *B ;
    double  *B2 ;
    double  *nu ;
    double  *v ;
    double  *v_newton ;
    double  *M ;
} NONLINEAR_MAG ;

typedef struct {
    int    flag ;
    int    mode ;
//    char   fname[FILENAME_MAX] ;

    DEF_FORM *form ;

    int    ntime ;
    DEF_TIME *time ;

    double remanent_magnetic_flux_density ;
    double coercivity ;
    double ialpha ;

    NONLINEAR_MAG *nl ;
} DEFINITION ;

typedef struct {
    char   label[1024] ;
    int    id ;
    int    num ;
    int    size ;
    double *index ;
    double **data ;
} DATA_ARRAY ;

typedef struct {
    int    nvol ;
    int    *vol ;

    int    nbc ;
    int    *bci ;


    int    nflag ;
    int    *flag_i ;
    double *flag_nu ;

    int    npermittivity ;
    int    *permittivity_i ;
    double *permittivity ;

    int    ncoil ;
    DEFINITION *coil_def ;
    int    ncoil_omega ;
    double coil_omega ;

    int    npmagnet ;
    int    npmagnet_nl ;
    DEFINITION *pmagnet_def ;

    /* Non-linear Analysis */
    int    nnl ;
    int    *nl_i ;
    NONLINEAR_MAG *nl_mag ;

    int    ntime_series ;
    DATA_ARRAY *time_series ;

    /* Conductor domain (Eddy Current) */
    int    nR ;
    int    *R_i ;
    double *R_sigma ;
} MTRL ;


typedef struct {
    int     snid ;
    int     nfid ;
    int     coordinate ;
    int     charge_part ;
    int     ln_chpar ;  /* local interface B.C. number in charge part */
    int     ln_part ; /* local interface B.C. number in part */
} Inbc ;




/* Coil或者永磁 */
typedef struct {
    int    coil ;
    int    nf ;
    int    *enf ;/**  **/
    int    *apid2gnid ;
    double *Jo ;/**  **/
    double *Jor ;
    double *Ihr ;
    double *Ihi ;
    double *Ih_mem ;
    DirichletBC dbc ;
    //#ifdef __ADVMAG_NEW_SOLVER__
    //	void*  solvmat ;
    //#else
    //	//TMP_AIJ tmp_aij ;
    //	SOLVER_MAT smat ;
    //#endif
    int     ninbc ;
    Inbc    *inbc ;
    DirichletBC dbc_inbc ;

    /* NS-Eddy */
    double omega ;
    double phase ;
} COIL ;

#endif // LIBSTRUCTURE_H
