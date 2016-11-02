/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef __SURFPACK_LAPACK_WRAPPERS_H__
#define __SURFPACK_LAPACK_WRAPPERS_H__

#ifdef HAVE_CONFIG_H
#include "surfpack_config.h"
#else
#include "surf77_config.h"
#endif

/***************************************************************************/
/***************************************************************************/
/**** The BLAS and LAPACK wrappers should be the same whichever version ****/
/**** of SurfMat we use, and since they are _ONLY_ wrappers they should ****/
/**** the should be inline functions so they should be in a header file ****/
/**** so it makes sense to put them here.                               ****/
/***************************************************************************/
/***************************************************************************/

#ifdef HAVE_CONFIG_H
#define PIVOTCHOL_F77 F77_FUNC(pivotchol,PIVOTCHOL)
#define BLOCKPIVOTCHOL_F77 F77_FUNC(blockpivotchol,BLOCKPIVOTCHOL)


#define DGETRF_F77 F77_FUNC(dgetrf,DGETRF)
#define DGETRI_F77 F77_FUNC(dgetri,DGETRI)
#define DGEMV_F77  F77_FUNC(dgemv,DGEMV)
#define DGEMM_F77  F77_FUNC(dgemm,DGEMM)
#define DDOT_F77   F77_FUNC(ddot, DDOT)
#define DGELS_F77  F77_FUNC(dgels,DGELS)
#define DGESVD_F77 F77_FUNC(dgesvd,DGESVD)

#define DPOTRF_F77 F77_FUNC(dpotrf,DPOTRF)
#define DPOTRI_F77 F77_FUNC(dpotri,DPOTRI)
#define DPOTRS_F77 F77_FUNC(dpotrs,DPOTRS)
#define DPOCON_F77 F77_FUNC(dpocon,DPOCON)
#define DGETRS_F77 F77_FUNC(dgetrs,DGETRS)
#define DLANGE_F77 F77_FUNC(dlange,DLANGE)
#define DGECON_F77 F77_FUNC(dgecon,DGECON)
#define DGGLSE_F77 F77_FUNC(dgglse,DGGLSE)
#define DSYEV_F77  F77_FUNC(dsyev,DSYEV)
#define DSYTRF_F77 F77_FUNC(dsytrf,DSYTRF)
#define DSYTRI_F77 F77_FUNC(dsytri,DSYTRI)
#define DSYTRS_F77 F77_FUNC(dsytrs,DSYTRS)
#define DSYCON_F77 F77_FUNC(dsycon,DSYCON)
//
#else
// Use the CMake generated fortran name mangling macros (eliminate warnings)
#define PIVOTCHOL_F77 SURF77_GLOBAL(pivotchol,PIVOTCHOL)
#define BLOCKPIVOTCHOL_F77 SURF77_GLOBAL(blockpivotchol,BLOCKPIVOTCHOL)


#define DGETRF_F77 SURF77_GLOBAL(dgetrf,DGETRF) 
#define DGETRI_F77 SURF77_GLOBAL(dgetri,DGETRI) 
#define DGEMV_F77  SURF77_GLOBAL(dgemv,DGEMV) 
#define DGEMM_F77  SURF77_GLOBAL(dgemm,DGEMM) 
#define DDOT_F77   SURF77_GLOBAL(ddot, DDOT) 
#define DGELS_F77  SURF77_GLOBAL(dgels,DGELS)
#define DGESVD_F77 SURF77_GLOBAL(dgesvd,DGESVD)

#define DPOTRF_F77 SURF77_GLOBAL(dpotrf,DPOTRF)
#define DPOTRI_F77 SURF77_GLOBAL(dpotri,DPOTRI)
#define DPOTRS_F77 SURF77_GLOBAL(dpotrs,DPOTRS)
#define DPOCON_F77 SURF77_GLOBAL(dpocon,DPOCON)
#define DGETRS_F77 SURF77_GLOBAL(dgetrs,DGETRS)
#define DLANGE_F77 SURF77_GLOBAL(dlange,DLANGE)
#define DGECON_F77 SURF77_GLOBAL(dgecon,DGECON)
#define DGGLSE_F77 SURF77_GLOBAL(dgglse,DGGLSE)
#define DSYEV_F77  SURF77_GLOBAL(dsyev,DSYEV)
#define DSYTRF_F77 SURF77_GLOBAL(dsytrf,DSYTRF)
#define DSYTRI_F77 SURF77_GLOBAL(dsytri,DSYTRI)
#define DSYTRS_F77 SURF77_GLOBAL(dsytrs,DSYTRS)
#define DSYCON_F77 SURF77_GLOBAL(dsycon,DSYCON)
//
#endif

/***************************************************************************/
/**** Fortran to C name mangling                                        ****/
/***************************************************************************/

extern "C" { // prevent C++ name mangling

//wrapper for KRD's implementation of the pivoting cholesky algorithm of: C. Lucas, "LAPACK-style Codes for LEvel2 and 3 Pivoted Cholesky Factorizations", Numerical Analysis Report No. 442, February 2004, from the Manchester Center for Computational Mathematics, I downloaded it from http://www.maths.manchester.ac.uk/~nareports/narep442.pdf, I think that Lucas's lev3pchol.f was ALWAYS defaulting to his lev2pchol.f, but this implementation was between 5% and 10% faster than Lucas's for the test problem (a 5500x5500 R matrix paviani 10D, 500 pts)
void PIVOTCHOL_F77(const char* uplo, const int* n, double* A, const int* lda,
		   int* piv, int* rank, const double* tol, int* info);
void BLOCKPIVOTCHOL_F77(const char* uplo, const int* n, double* A, 
			const int* lda, const int *blocksize, int* piv, 
			int* rank, double *dwork, const double* tol, 
			int* info);


/***************************************************************************/
/**** BLAS Fortran to C name mangling                                   ****/
/***************************************************************************/

// Vector-vector inner product
double DDOT_F77(const int* n, const double* x, const int* incx,
		const double* y, const int* incy);


// Matrix-vector multiplication
void DGEMV_F77(char* trans, const int* m, const int* n, const double* alpha, 
	       const double* A, const int* lda, const double* x,
	       const int* incx, const double* beta, double* y, const int* incy);

// Matrix-matrix multiplication
void DGEMM_F77(char* transa, char* transb, const int* m, const int* n,
	       const int* k, const double* alpha, const double* A,
	       const int* lda, const double* B, const int* ldb, 
	       const double* beta, double* C, const int* ldc);

/***************************************************************************/
/**** LAPACK Fortran to C name mangling                                 ****/
/***************************************************************************/

// Perform Cholesky factorization
void DPOTRF_F77(const char* uplo, const int* n, double* AChol, const int* lda, int* info);

// Compute the inverse of a matrix expressed as an cholesky decomposition (i.e., call dpotrf on the matrix first)
void DPOTRI_F77(const char* uplo, const int* n, double* ACholInv, const int* lda, int* info);

// solve A*X=B for X, after A has been Cholesky factorized (i.e., call dptorf on the matrix first)
void DPOTRS_F77(const char* uplo, const int* n, const int* nRHS, 
		const double* AChol,
		const int* ldAChol , double* RHS, 
		const int* ldRHS, int* info);

// function to compute the condition number of a matrix from the Cholesky factorization
void DPOCON_F77(const char* uplo, const int* n, const double* AChol, const int* lda, const double* anorm, double* rconda, double* work, int* iwork, int* info);

  
// Performs an L*D*L^T  (or U*D*U^T, either can be used) factorization 
void DSYTRF_F77(const char* uplo, const int* n, double* ALDLT, const int* lda, int* ipiv, double* work, const int* lwork, const int* info);

// Compute the inverse of a matrix expressed as an L*D*L^T (or U*D*U^T, either can be used) factorization (i.e. call dsytrf on the matrix first)
void DSYTRI_F77(const char* uplo, const int* n, double* ALDLTINV, const int* lda, const int* ipiv, double* work, int* info);

// solve A*X=B for X, after A has been L*D*L^T (or U*D*U^T, either can be used) factorized (i.e. call dsytrf on the matrix first)
  void DSYTRS_F77(const char* uplo, const int* n, const int* nRHS, const double* ALDLT, const int* lda, const int* ipiv, double* RHS, const int* ldRHS, int* info);

// function to compute the condition number of a matrix from the L*D*L^T (or U*D*U^T, either can be used) factorization
void DSYCON_F77(const char* uplo, const int* n, const double* ALDLT, const int* lda, const int* ipiv, const double* anorm, double* rconda, double* work, int* iwork, int* info);


// Perform SVD factorization
void DGESVD_F77(const char* jobu, const char* jobvt, const int* m, const int* n, double* A, const int* lda, double* S, double* U, const int* ldu, double* VT, const int* ldvt, double* work, const int* lwork, int* info);


// Perform LU factorization
void DGETRF_F77(const int* m, const int* n, double* a, const int* lda,
		int* ipiv, int* info);

// Compute the inverse of a matrix expressed as an LU decomposition (i.e., call dgetrf on the matrix first)
void DGETRI_F77(const int* n, double* a, const int* lda, const int* ipiv,
		double* work, const int* lwork, int* info);

// solve A*X=B for X, where A={A || A^T} after A has been LU factorized (i.e., call dgetrf on the matrix first)
void DGETRS_F77(const char* transLU, const int* n, const int* nRHS,
		const double* LU, const int* ldLU , const int* ipiv, 
		double* RHS, const int* ldRHS, int* info);

//function to compute the norm of a matrix A, choices are
//M max(abs(A(i,j))) this is not a consistent matrix norm, 
//1 one norm of a matrix, maximum column sum, 
//I infinity norm of matrix, maximum row sum,or 
//F frobenius norm of a matrix, square root of sum of squares
double DLANGE_F77(char *whichnorm, int *M, int *N, const double *A, int *LDA,
		  double *work);

//function to compute the condition number of a matrix
void DGECON_F77(const char *whichnorm, const int *N, const double *ALU,
        const int *LDA, const double *anorm,
        double *rcond, double *work, int *iwork, int *info);

// Least-squares solution to linear system of equations
void DGELS_F77(const char* trans, const int* nrows, const int* ncols,
	       const int* nrhs, double* A, const int* lda, double* b,
	       const int* ldb, double* work, const int* lwork, int* info);

// Performs least-squares solve subject to equality constraints
void DGGLSE_F77(const int* m, const int* n, const int* p, double* A,
		const int* lda, double* B, const int* ldb, double* c,
		double* d, double* x, double* work, const int* lwork,
		int* info);

// determines eigenvalues and (optionally) eigenvectors for a real symmetric matrix
void DSYEV_F77(const char* jobz, const char* uplo, const int* N, 
	       double *A_EIGVECT, const int* lda, double* eigval, 
	       double* work, const int* lwork, int* info);

} // extern "C" (prevent C++ name mangling)


#endif // __SURFPACK_LAPACK_WRAPPERS_H__

