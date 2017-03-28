#pragma once
#include "slu_mt_ddefs.h"
//define 32 bit word, must be consistant with other solver
//must be placed before include file
#if !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD
#endif
#include <armadillo>
#include <vector>



class CSuperLU_MT
{
public:
	SuperMatrix sluA, L, U;
	SuperMatrix sluB, sluX;
	NCformat    *Astore;
	SCPformat   *Lstore; DNformat	   *Bstore;
	NCPformat   *Ustore;
	int_t         nprocs;
	fact_t      fact;
	trans_t     trans;
	yes_no_t    refact, usepr;
	equed_t     equed;
	double      *a; double t1;
	int_t         *asub, *xa;
	int_t         *perm_c; /* column permutation vector */
	int_t         *perm_r; /* row permutations from partial pivoting */
	void        *work;
	superlumt_options_t superlumt_options;
	int_t         info, lwork, nrhs, ldx, panel_size, relax;
	int_t         m, n, nnz, permc_spec;
	double      *rhsb, *rhsx, *xact;
	double      *R, *C;
	double      *ferr, *berr;
	double      u, drop_tol, rpg, rcond;
	superlu_memusage_t superlu_memusage;

	CSuperLU_MT();
	CSuperLU_MT(int mm, arma::sp_mat& X,double *b);
	int solve();
	int LUsolve();
	double * getResult();
	~CSuperLU_MT();
};

