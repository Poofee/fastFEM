#include "SuperLU_MT.h"


CSuperLU_MT::CSuperLU_MT() {
	
}


CSuperLU_MT::CSuperLU_MT(int mm, arma::sp_mat &X, double * b) {
	/* Default parameters to control factorization. */
	nprocs = 1;
	fact = EQUILIBRATE;
	trans = NOTRANS;
	equed = NOEQUIL;
	refact = NO;
	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	u = 1.0;
	usepr = NO;
	drop_tol = 0.0;
	lwork = 0;
	nrhs = 1;

	if (lwork > 0) {
		work = SUPERLU_MALLOC(lwork);
		//printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
		if (!work) {
			SUPERLU_ABORT("DLINSOLX: cannot allocate work[]");
		}
	}

	/* create matrix A in Harwell-Boeing format.*/
	m = mm; n = mm; nnz = X.n_nonzero;
    a = const_cast<double *>(X.values);
    //bug:arma,uint;superlu,int
    //this place, use <> deconst the const variable,
    //and use () to force change the type,
    //note there exist data loss
    asub = (int*)const_cast<unsigned int*>(X.row_indices);
    xa = (int*)const_cast<unsigned int*>(X.col_ptrs);
	dCreate_CompCol_Matrix(&sluA, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

	StatAlloc(n, nprocs, panel_size, relax, &Gstat1);
	StatInit(n, nprocs, &Gstat1);
	//------create B and X-------------------
	if (!(rhsx = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
	dCreate_Dense_Matrix(&sluX, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
	//

	rhsb = b;
	dCreate_Dense_Matrix(&sluB, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
	Bstore = (DNformat*)sluB.Store;
	//double *Bmat = (double*)Bstore->nzval;
	//double *Xmat = (double*)((DNformat*)sluX.Store)->nzval;
	//ldx = m;
	if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
	if (!(R = (double *)SUPERLU_MALLOC(sluA.nrow * sizeof(double))))
		SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
	if (!(C = (double *)SUPERLU_MALLOC(sluA.ncol * sizeof(double))))
		SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
	if (!(ferr = (double *)SUPERLU_MALLOC(nrhs * sizeof(double))))
		SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
	if (!(berr = (double *)SUPERLU_MALLOC(nrhs * sizeof(double))))
		SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");
	/*   * Get column permutation vector perm_c[], according to permc_spec:
	*   permc_spec = 0: natural ordering
	*   permc_spec = 1: minimum degree ordering on structure of A'*A
	*   permc_spec = 2: minimum degree ordering on structure of A'+A
	*   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	*/
	permc_spec = 1;
	get_perm_c(permc_spec, &sluA, perm_c);

	superlumt_options.nprocs = nprocs;
	superlumt_options.fact = EQUILIBRATE;
	superlumt_options.trans = NOTRANS;
	superlumt_options.refact = NO;
	superlumt_options.panel_size = panel_size;
	superlumt_options.relax = relax;
	superlumt_options.usepr = NO;
	superlumt_options.drop_tol = drop_tol;
	superlumt_options.diag_pivot_thresh = u;
	superlumt_options.SymmetricMode = NO;
	superlumt_options.PrintStat = NO;
	superlumt_options.perm_c = perm_c;
	superlumt_options.perm_r = perm_r;
	superlumt_options.work = work;
	superlumt_options.lwork = 0;
	if (!(superlumt_options.etree = intMalloc(n)))
		SUPERLU_ABORT("Malloc fails for etree[].");
	if (!(superlumt_options.colcnt_h = intMalloc(n)))
		SUPERLU_ABORT("Malloc fails for colcnt_h[].");
	if (!(superlumt_options.part_super_h = intMalloc(n)))
		SUPERLU_ABORT("Malloc fails for colcnt_h[].");
}


int CSuperLU_MT::solve() {
	/*   * Solve the system and compute the condition number
	* and error bounds using pdgssvx.
	*/
	//t1 = SuperLU_timer_();
	pdgssvx(nprocs, &superlumt_options, &sluA, perm_c, perm_r,
		&equed, R, C, &L, &U, &sluB, &sluX, &rpg, &rcond,
		ferr, berr, &superlu_memusage, &info);
	//t1 = SuperLU_timer_() - t1;
	//printf("%lf\n",t1);
	/* This is how you could access the solution matrix. */
	//double *sol = NULL;
	if (info == 0 || info == n + 1) {
		//sol = (double*)((DNformat*)sluX.Store)->nzval;
		//for (int i = 0; i < num_pts; i++) {
		//	pmeshnode[i].A = sol[i] * miu0;
		//	A(i) = sol[i] * miu0;
		//}
		return 0;
	} else if (info > 0 && lwork == -1) {
		return 1;
	}
	return 1;
}
int CSuperLU_MT::solve1() {
	pdgssv(nprocs, &sluA, perm_c, perm_r, &L, &U, &sluB, &info);	
	if (info == 0 || info == n + 1) {
		//sol = (double*)((DNformat*)sluX.Store)->nzval;
		//for (int i = 0; i < num_pts; i++) {
		//	pmeshnode[i].A = sol[i] * miu0;
		//	A(i) = sol[i] * miu0;
		//}
		return 0;
	} else if (info > 0 && lwork == -1) {
		return 1;
	}
	return 1;
}
int CSuperLU_MT::triangleSolve(){
	dgstrs(trans, &L, &U, perm_r, perm_c, &sluB, &Gstat1, &info);
	if (info == 0 || info == n + 1) {
		//sol = (double*)((DNformat*)sluX.Store)->nzval;
		//for (int i = 0; i < num_pts; i++) {
		//	pmeshnode[i].A = sol[i] * miu0;
		//	A(i) = sol[i] * miu0;
		//}
		return 0;
	} else if (info > 0 && lwork == -1) {
		return 1;
	}
	return 1;
}
int CSuperLU_MT::LUsolve() {
	//NOW WE SOLVE THE LINEAR SYSTEM USING THE FACTORED FORM OF sluA.
	//Bstore->nzval = const_cast<double*>(b.mem);
	superlumt_options.fact = FACTORED; /* Indicate the factored form of sluA is supplied. */
	superlumt_options.usepr = YES;
	//get_perm_c(permc_spec, &sluA, perm_c);
	//t1 = SuperLU_timer_();
	pdgssvx(nprocs, &superlumt_options, &sluA, perm_c, perm_r,
		&equed, R, C, &L, &U, &sluB, &sluX, &rpg, &rcond,
		ferr, berr, &superlu_memusage, &info);
	//t1 = SuperLU_timer_() - t1;
	if (info == 0 || info == n + 1) {
		//sol = (double*)((DNformat*)sluX.Store)->nzval;
		//for (int i = 0; i < num_pts; i++) {
		//	pmeshnode[i].A = sol[i] * miu0;
		//	A(i) = sol[i] * miu0;
		//}
		return 0;
	} else if (info > 0 && lwork == -1) {
		return 1;
	}
	return 1;
}

double * CSuperLU_MT::getResult() {
	if (info == 0 || info == n + 1) {
		//sol = (double*)((DNformat*)sluX.Store)->nzval;
		//for (int i = 0; i < num_pts; i++) {
		//	pmeshnode[i].A = sol[i] * miu0;
		//	A(i) = sol[i] * miu0;
		//}
		return (double*)((DNformat*)sluX.Store)->nzval;
	} else if (info > 0 && lwork == -1) {
		return NULL;
	}	
	return NULL;
}

double * CSuperLU_MT::get1Result() {
	if (info == 0 || info == n + 1) {
		//sol = (double*)((DNformat*)sluX.Store)->nzval;
		//for (int i = 0; i < num_pts; i++) {
		//	pmeshnode[i].A = sol[i] * miu0;
		//	A(i) = sol[i] * miu0;
		//}
		return (double*)((DNformat*)sluB.Store)->nzval;
	} else if (info > 0 && lwork == -1) {
		return NULL;
	}
	return NULL;
}
CSuperLU_MT::~CSuperLU_MT() {
	/*free the space allocated by superLU*/
	//SUPERLU_FREE(rhsb);//not allocate memory, don't free space!
	SUPERLU_FREE(rhsx);
	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	SUPERLU_FREE(R);
	SUPERLU_FREE(C);
	SUPERLU_FREE(ferr);
	SUPERLU_FREE(berr);
	Destroy_SuperMatrix_Store(&sluA);
	Destroy_SuperMatrix_Store(&sluB);
	Destroy_SuperMatrix_Store(&sluX);
	SUPERLU_FREE(superlumt_options.etree);
	SUPERLU_FREE(superlumt_options.colcnt_h);
	SUPERLU_FREE(superlumt_options.part_super_h);
	if (lwork == 0) {
		Destroy_SuperNode_SCP(&L);
		Destroy_CompCol_NCP(&U);
	} else if (lwork > 0) {
		SUPERLU_FREE(work);
	}
}
