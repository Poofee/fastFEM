/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
* -- SuperLU MT routine (version 3.0) --
* Lawrence Berkeley National Lab, Univ. of California Berkeley,
* and Xerox Palo Alto Research Center.
* September 10, 2007
*
*/
#include "slu_mt_ddefs.h"

/* Eat up the rest of the current line */
int_t dDumpLine(FILE *fp) {
	register int_t c;
	while ((c = fgetc(fp)) != '\n');
	return 0;
}

int_t dParseIntFormat(char *buf, int_t *num, int_t *size) {
	char *tmp;

	tmp = buf;
	while (*tmp++ != '(');
	*num = atoi(tmp);
	while (*tmp != 'I' && *tmp != 'i') ++tmp;
	++tmp;
	*size = atoi(tmp);
	return 0;
}

int_t dParseFloatFormat(char *buf, int_t *num, int_t *size) {
	char *tmp, *period;

	tmp = buf;
	while (*tmp++ != '(');
	*num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
	while (*tmp != 'E' && *tmp != 'e' && *tmp != 'D' && *tmp != 'd'
		&& *tmp != 'F' && *tmp != 'f') {
		/* May find kP before nE/nD/nF, like (1P6F13.6). In this case the
		num picked up refers to P, which should be skipped. */
		if (*tmp == 'p' || *tmp == 'P') {
			++tmp;
			*num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
		} else {
			++tmp;
		}
	}
	++tmp;
	period = tmp;
	while (*period != '.' && *period != ')') ++period;
	*period = '\0';
	*size = atoi(tmp); /*sscanf(tmp, "%2d", size);*/

	return 0;
}

int_t dReadVector(FILE *fp, int_t n, int_t *where, int_t perline, int_t persize) {
	register int_t i, j, item;
	char tmp, buf[100];

	i = 0;
	while (i < n) {
		fgets(buf, 100, fp);    /* read a line at a time */
		for (j = 0; j<perline && i<n; j++) {
			tmp = buf[(j + 1)*persize];     /* save the char at that place */
			buf[(j + 1)*persize] = 0;       /* null terminate */
			item = atoi(&buf[j*persize]);
			buf[(j + 1)*persize] = tmp;     /* recover the char at that place */
			where[i++] = item - 1;
		}
	}

	return 0;
}

int_t dReadValues(FILE *fp, int_t n, double *destination, int_t perline, int_t persize) {
	register int_t i, j, k, s;
	char tmp, buf[100];

	i = 0;
	while (i < n) {
		fgets(buf, 100, fp);    /* read a line at a time */
		for (j = 0; j<perline && i<n; j++) {
			tmp = buf[(j + 1)*persize];     /* save the char at that place */
			buf[(j + 1)*persize] = 0;       /* null terminate */
			s = j*persize;
			for (k = 0; k < persize; ++k) /* No D_ format in C */
				if (buf[s + k] == 'D' || buf[s + k] == 'd') buf[s + k] = 'E';
			destination[i++] = atof(&buf[s]);
			buf[(j + 1)*persize] = tmp;     /* recover the char at that place */
		}
	}

	return 0;
}

int main(int argc, char *argv[]) {
	SuperMatrix   A;
	NCformat *Astore;
	double   *a;
	int_t      *asub, *xa;
	int_t      *perm_r; /* row permutations from partial pivoting */
	int_t      *perm_c; /* column permutation vector */
	SuperMatrix   L;       /* factor L */
	SCPformat *Lstore;
	SuperMatrix   U;       /* factor U */
	NCPformat *Ustore;
	SuperMatrix   B;
	int_t      nrhs, ldx, info, m, n, nnz, b;
	int_t      nprocs; /* maximum number of processors to use. */
	int_t      panel_size, relax, maxsup;
	int_t      permc_spec;
	trans_t  trans;
	double   *xact, *rhs;
	superlu_memusage_t   superlu_memusage;
	//void   parse_command_line();

	nrhs = 1;
	trans = NOTRANS;
	nprocs = 1;
	n = 1000;
	b = 1;
	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
	maxsup = sp_ienv(3);

	//parse_command_line(argc, argv, &nprocs, &n, &b, &panel_size,
	//	&relax, &maxsup);
	//dreadhb(int_t *nrow, int_t *ncol, int_t *nonz,
	//	double **nzval, int_t **rowind, int_t **colptr)
	//dreadhb(&m, &n, &nnz, &a, &asub, &xa);
	register int_t i, numer_lines, rhscrd = 0;
	int_t tmp, colnum, colsize, rownum, rowsize, valnum, valsize;
	char buf[100], type[4], key[10];
	FILE *fp;
	//fp = stdin;
	fp = fopen("E:\\Projects\\cplusplus\\SuperLU_MT_3.1\\EXAMPLE\\big.rua", "r");

	/* Line 1 */
	fscanf(fp, "%72c", buf); buf[72] = 0;
	printf("Title: %s", buf);
	fscanf(fp, "%8c", key);  key[8] = 0;
	printf("Key: %s\n", key);
	dDumpLine(fp);

	/* Line 2 */
	for (i = 0; i<5; i++) {
		fscanf(fp, "%14c", buf); buf[14] = 0;
		tmp = atoi(buf); /*sscanf(buf, "%d", &tmp);*/
		if (i == 3) numer_lines = tmp;
		if (i == 4 && tmp) rhscrd = tmp;
	}
	dDumpLine(fp);

	/* Line 3 */
	fscanf(fp, "%3c", type);
	fscanf(fp, "%11c", buf); /* pad */
	type[3] = 0;
#if ( DEBUGlevel>=1 )
	printf("Matrix type %s\n", type);
#endif

	fscanf(fp, "%14c", buf); m = atoi(buf);
	fscanf(fp, "%14c", buf); n = atoi(buf);
	fscanf(fp, "%14c", buf); nnz = atoi(buf);
	fscanf(fp, "%14c", buf); tmp = atoi(buf);

	if (tmp != 0)
		printf("This is not an assembled matrix!\n");
	if (m != n)
		printf("Matrix is not square.\n");
	dDumpLine(fp);

	/* Line 4: format statement */
	fscanf(fp, "%16c", buf);
	dParseIntFormat(buf, &colnum, &colsize);
	fscanf(fp, "%16c", buf);
	dParseIntFormat(buf, &rownum, &rowsize);
	fscanf(fp, "%20c", buf);
	dParseFloatFormat(buf, &valnum, &valsize);
	fscanf(fp, "%20c", buf);
	dDumpLine(fp);

	/* Line 5: right-hand side */
	if (rhscrd) dDumpLine(fp); /* skip RHSFMT */

#if ( DEBUGlevel>=1 )
	printf("%d rows, %d nonzeros\n", *nrow, *nonz);
	printf("colnum %d, colsize %d\n", colnum, colsize);
	printf("rownum %d, rowsize %d\n", rownum, rowsize);
	printf("valnum %d, valsize %d\n", valnum, valsize);
#endif

	/* Allocate storage for the three arrays ( nzval, rowind, colptr ) */
	dallocateA(n, nnz, &a, &asub, &xa);

	dReadVector(fp, m + 1, xa, colnum, colsize);
	dReadVector(fp, nnz, asub, rownum, rowsize);
	if (numer_lines) {
		dReadValues(fp, nnz, a, valnum, valsize);
	}

	fclose(fp);

	dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
	Astore = (NCformat*)A.Store;
	printf("Dimension " IFMT "x" IFMT "; # nonzeros " IFMT "\n", A.nrow, A.ncol, Astore->nnz);

	if (!(rhs = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhs[].");
	dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
	xact = doubleMalloc(n * nrhs);
	ldx = n;
	dGenXtrue(n, nrhs, xact, ldx);
	dFillRHS(trans, nrhs, xact, ldx, &A, &B);

	if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
	if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");

	/*
	* Get column permutation vector perm_c[], according to permc_spec:
	*   permc_spec = 0: natural ordering
	*   permc_spec = 1: minimum degree ordering on structure of A'*A
	*   permc_spec = 2: minimum degree ordering on structure of A'+A
	*   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	*/
	permc_spec = 1;
	get_perm_c(permc_spec, &A, perm_c);

	pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);

	if (info == 0) {
		dinf_norm_error(nrhs, &B, xact); /* Inf. norm of the error */

		Lstore = (SCPformat *)L.Store;
		Ustore = (NCPformat *)U.Store;
		printf("#NZ in factor L = " IFMT "\n", Lstore->nnz);
		printf("#NZ in factor U = " IFMT "\n", Ustore->nnz);
		printf("#NZ in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - L.ncol);

		superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
			superlu_memusage.for_lu / 1024 / 1024,
			superlu_memusage.total_needed / 1024 / 1024,
			superlu_memusage.expansions);

	}

	SUPERLU_FREE(rhs);
	SUPERLU_FREE(xact);
	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_SCP(&L);
	Destroy_CompCol_NCP(&U);
	getchar();
}

/*
* Parse command line to get relaxed snode size, panel size, etc.
*/
//void
//parse_command_line(int argc, char *argv[], int_t *procs, int_t *n,
//int_t *b, int_t *w, int_t *r, int_t *maxsup) {
//	register int c;
//	extern char *optarg;
//
//	while ((c = getopt(argc, argv, "ht:p:n:b:w:x:s:")) != EOF) {
//		switch (c) {
//		case 'h':
//			printf("Options: (default values are in parenthesis)\n");
//			printf("\t-p <int> - number of processes     ( " IFMT " )\n", *procs);
//			printf("\t-n <int> - dimension               ( " IFMT " )\n", *n);
//			printf("\t-b <int> - semi-bandwidth          ( " IFMT " )\n", *b);
//			printf("\t-w <int> - panel size              ( " IFMT " )\n", *w);
//			printf("\t-x <int> - relax                   ( " IFMT " )\n", *r);
//			printf("\t-s <int> - maximum supernode size  ( " IFMT " )\n", *maxsup);
//			exit(1);
//			break;
//		case 'p': *procs = atoi(optarg);
//			break;
//		case 'n': *n = atoi(optarg);
//			break;
//		case 'b': *b = atoi(optarg);
//			break;
//		case 'w': *w = atoi(optarg);
//			break;
//		case 'x': *r = atoi(optarg);
//			break;
//		case 's': *maxsup = atoi(optarg);
//			break;
//		}
//	}
//}

