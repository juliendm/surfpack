#ifndef __SURFMAT_HPP__
//do not #define __SURFMAT_HPP__ here, that should/must only be done in either CustomSurfMat.hpp OR TeuchosSurfMat.hpp, to keep them from being included directly when the other has already been included

//#define __SURFMAT_ERR_CHECK__

#include "surfpack_LAPACK_wrappers.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iostream>

#include "NKM_SurfMat_Native.hpp" //a native implementation
//#include "SurfMat_Teuchos.hpp" //a wrapper for the Teuchos Serial Dense Matrix class


namespace nkm {

typedef SurfMat<double> MtxDbl;
typedef SurfMat<int> MtxInt;


/***************************************************************************/
/**** The BLAS wrappers start here                                      ****/ 
/***************************************************************************/

/// sum of element by element products of a and b (works for matrices as well as vectors and that functionality is required for the Kriging/GP implementation), wraps DDOT
double dot_product(const MtxDbl& a, const MtxDbl& b);

/// matrix matrix Mult OR Matrix Vector Mult: C=scaleAB*A*B+scaleC*C, where A={A || A^T} and B={B || B^T}, wraps and chooses between DGEMV and DGEMM
MtxDbl& matrix_mult(MtxDbl& C, const MtxDbl& A, const MtxDbl& B, 
		    double scaleC=0.0, double scaleAB=1.0, 
		    char transA='N', char transB='N');

/***************************************************************************/
/**** The LAPACK wrappers start here                                    ****/ 
/***************************************************************************/


/// computes the L*D*L^T factorization of "matrix" with partial pivoting ("matrix" must be real and symmetric but does not need to be positive definite) wraps LAPACK subroutine DSYTRF (and also DLANGE and DSYCON to return the 1-norm rcond of "matrix", doing them together makes the rcond cheap)
MtxDbl& LDLT_fact(MtxDbl& matrix, MtxInt& ipvt, MtxDbl& scale, int& info, double& rcond);

/// inverts a real symmetric positive definite matrix in place after the L*D*L^T factorization has already been done, requires the lower triangular portion as input returns the full symmetric inverse (as opposed to just the lower triangular part), wraps LAPACK subroutine DSYTRI
MtxDbl& inverse_after_LDLT_fact(MtxDbl& matrix, const MtxInt& ipvt, const MtxDbl& scale);

//computes the (1 norm) reciprocal of the condition number of A from the L*D*L^T factorization of A (A must be real and symmetric but does not need to be positive definite) wraps LAPACK subroutines DLANGE and DSYCON
double rcond_after_LDLT_fact(const MtxDbl& A, const MtxDbl& ALDLT, const MtxInt& ipvt);

/// solves A*X=B for X, where A is real symmetric and B={B || B^T}, without changing the contents of B, after A has been L*D*L^T factorized, ALDLT must contain the lower triangular portion of the factorization of A, wraps LAPACK subroutine DSYTRS 
MtxDbl& solve_after_LDLT_fact(MtxDbl& result, const MtxDbl& ALDLT, const MtxInt& ipvt, const MtxDbl& scale, const MtxDbl& BRHS,char transB='N');

///compute the pseudo inverse of A, this function wraps the LAPACK subroutine DGESVD in this context rcond is considered to be the ratio of the smallest to largest singular values, A becomes A_pinv
MtxDbl& pseudo_inverse(MtxDbl& A, double min_allowed_rcond, double& rcondA, double& log_abs_det, int& if_det_eq_zero);

/// perform a numerically optimally preconditioned Cholesky factorization of a real symmetric positive semi-definite matrix, wraps DPOTRF, returns the lower triangular portion, the preconditioning is "undone" so the L matrix is the same (minus some rounding error due to poor conditioning) as it would be without precondtioning, if info>0 then the preconditined matrix is singular, rcondprecond is the standard fast estimate of the reciprocal of the condition number of the numerically optimally preconditioned real symmetric positive definite matrix (that's what affects the round off error)
MtxDbl& Chol_fact(MtxDbl& matrix, int& info, double& rcondprecond);

MtxDbl& Chol_fact_workspace(MtxDbl& matrix, MtxDbl& scalefactor, 
			    MtxDbl& work, MtxInt& iwork, int& info, 
			    double& rcondprecond);

/// inverts a real symmetric positive definite matrix in place after the Cholesky factorization has already been done, requires the lower triangular portion as input returns the full symmetric inverse (as opposed to just the lower triangular part), wraps DPOTRF
MtxDbl& inverse_after_Chol_fact(MtxDbl& matrix);

//computes the (1 norm) reciprocal of the condition number of A from the Cholesky factorization of A (A must be real symmetric and positive semi-definite)
double rcond_after_Chol_fact(const MtxDbl& A, const MtxDbl& AChol);

/// solves A*X=B for X, where A is symmetric positive definite and B={B || B^T}, without changing the contents of B, after A has been Cholesky factorized, AChol must contain the lower triangular portion of the factorization of A, wraps DPOTRS 
MtxDbl& solve_after_Chol_fact(MtxDbl& result, const MtxDbl& AChol, const MtxDbl& BRHS,char transB='N');

/// perform an LU factorization with partial pivoting, wraps LAPACK subroutine DGETRF
MtxDbl& LU_fact(MtxDbl& matrix, MtxInt& ipvt);

/// inverts a matrix, in place, after an LU factorization has already been done, wraps DGETRI
MtxDbl& inverse_after_LU_fact(MtxDbl& matrix, MtxInt& ipvt);

//computes the (1 norm) reciprocal of the condition number of A from the LU factorization of A
double rcond_after_LU_fact(const MtxDbl& A, const MtxDbl& ALU);

/// solves A*X=B for X, where A={A || A^T} and B={B || B^T}, without changing the contents of B, after A has been LU factorized, wraps DGETRS 
MtxDbl& solve_after_LU_fact(MtxDbl& result, const MtxDbl& ALU, const MtxInt& ipvt, 
			    const MtxDbl& BRHS, char transA='N', 
			    char transB='N');

/// wraps DGELS
void least_squares(MtxDbl& A, MtxDbl& x, MtxDbl& b);

/// wrapse DGGLSE
void least_squares_with_equality_constraints(MtxDbl& A, 
     MtxDbl& x, MtxDbl& c, MtxDbl& B, MtxDbl& d);

///finds the eigenvalues and optionally (by default) eigenvectors of a real symmetric matrix, returns a reference to the vector of eigenvalues, this function wraps the LAPACK subroutine DSYEV
MtxDbl& eig_sym(MtxDbl& eigvect, MtxDbl& eigval, const MtxDbl& A, char jobz='V');

///finds the eigenvalues (but not eigenvectors) of a real symmetric matrix, returns a reference to the vector of eigenvalues, this function wraps the function inline MtxDbl& eig_sym(MtxDbl& eigvect, MtxDbl& eigval, const MtxDbl& A, char jobz='V') which inturn wraps LAPACK subroutine DSYEV
inline MtxDbl& eig_sym(MtxDbl& eigval, const MtxDbl& A) {
  MtxDbl eigvect;
  return (eig_sym(eigvect,eigval,A,'N'));
}


/*****************************************************************************/
/* extra functions added for convenience                                     */
/*****************************************************************************/


MtxDbl& pseudo_inverse_sym(MtxDbl& A, double min_allowed_rcond, double& rcond_A, double& log_abs_det_A, double& sign_of_det_A);

/// computes log(fabs(det(A))) and sign(det(A)) from the L*D*L^T factorization of A by adding logs of absolute value of determinants of component diagonal blocks to prevent underflow/overflow errors that would by taking product of diagonal block determinants before taking the log.  You must have called LDLT_fact() before calling this function
double log_det_after_LDLT_fact(const MtxDbl& ALDLT, const MtxInt& ipvt, const MtxDbl& scale, double& sign_of_det);

/// computes the determinant of A from its LDLT factorization as sign_of_det_A*exp(log_fabs_det_A) to prevent intermediate underflow/overflow errors (i.e. give large and small components the best chance of balancing), you must have called LDLT_fact() first.
inline double det_after_LDLT_fact(const MtxDbl& ALDLT, const MtxInt& ipvt, const MtxDbl& scale) {
  double sign_of_det;
  double log_det_A=log_det_after_LDLT_fact(ALDLT, ipvt, scale, sign_of_det);
  return (sign_of_det*std::exp(log_det_A));
}

} // end namespace nkm

#endif
