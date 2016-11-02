#include "NKM_SurfMat.hpp" 


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
// export for each template type for now (currently only double and int)
BOOST_CLASS_EXPORT(nkm::SurfMat<double>)
BOOST_CLASS_EXPORT(nkm::SurfMat<int>)
#endif

namespace nkm {

/***************************************************************************/
/**** The BLAS wrappers start here                                      ****/ 
/***************************************************************************/

/// sum of element by element products of a and b (works for matrices as well as vectors and that functionality is required for the Kriging/GP implementation), wraps DDOT
double dot_product(const MtxDbl& a, const MtxDbl& b)
{
  int ncolsa=a.getNCols();
  int nrowsa=a.getNRows();
  int nrowsb=b.getNRows();
  int ncolsb=b.getNCols();
  int nelem=nrowsa*ncolsa;
#ifdef __SURFMAT_ERR_CHECK__
  assert(nelem==nrowsb*b.getNCols());
#endif
  int inc=1;
  if(((nrowsa==a.getNRowsAct())&&(nrowsb==b.getNRowsAct()))||
     ((ncolsa==1)&&(ncolsb==1))
     ) {
    // ddot will not violate the constness
    return DDOT_F77(&nelem, a.ptr(0,0), &inc, b.ptr(0,0), &inc);
  }
  else if(nrowsa==nrowsb) {
    double dotprod=DDOT_F77(&nrowsa, a.ptr(0,0), &inc, b.ptr(0,0), &inc);
    for(int j=1; j<ncolsa; ++j)
      dotprod+=DDOT_F77(&nrowsa, a.ptr(0,j), &inc, b.ptr(0,j), &inc);

    return dotprod;
  }
  else{
    double dotprod=a(0,0)*b(0,0);
    if((nrowsa==ncolsb)&&(ncolsa==1)&&(nrowsb==1)) 
      for(int i=1; i<nrowsa; ++i)
	dotprod+=a(i,0)*b(0,i);
    else if((nrowsb==ncolsa)&&(ncolsb==1)&&(nrowsa==1)) 
      for(int i=1; i<ncolsa; ++i)
	dotprod+=a(0,i)*b(i,0);
    else
      assert(false);
    return dotprod;
  }
}

/// matrix matrix Mult OR Matrix Vector Mult: C=scaleAB*A*B+scaleC*C, where A={A || A^T} and B={B || B^T}, wraps and chooses between DGEMV and DGEMM
MtxDbl& matrix_mult(MtxDbl& C, const MtxDbl& A, const MtxDbl& B, 
		    double scaleC, double scaleAB, char transA, char transB)
//#ifdef __SURFMAT_ERR_CHECK__
//{
//  return matrix_mult_debug(C, A, B, scaleC, scaleAB, transA, transB);
//} 
//#else
{
#ifdef __SURFMAT_ERR_CHECK__
  assert(((transA=='N')||(transA=='T'))&&((transB=='N')||(transB=='T')));  
#endif
  int nrowsC, ncolsC, ninnerA, ninnerB;
  if(transA=='N') {
    nrowsC =A.getNRows();
    ninnerA=A.getNCols();
  }
  else{
    nrowsC =A.getNCols();
    ninnerA=A.getNRows();
  }
  if(transB=='N') {
    ncolsC =B.getNCols();
    ninnerB=B.getNRows();
  }
  else{
    ncolsC =B.getNRows();
    ninnerB=B.getNCols();
  }
#ifdef __SURFMAT_ERR_CHECK__
  if(!(ninnerA==ninnerB)){
    printf("ninnerA=%d ninnerB=%d",ninnerA,ninnerB);
    printf("\n");
    assert(ninnerA==ninnerB);
  }
#endif
  C.newSize(nrowsC,ncolsC);
  C.putTol(A.getTol()); //revise this once nonzero-tol starts to be used
  
  int nrowsA=static_cast<int>(A.getNRows()); //transpose does not affect
  int ncolsA=static_cast<int>(A.getNCols()); //transpose does not affect
  //int nrowsB=static_cast<int>(B.getNRows()); //transpose does not affect
  int lda=static_cast<int>(A.getNRowsAct()); //transpose does not affect
  int ldb=static_cast<int>(B.getNRowsAct()); //transpose does not affect
  int ldc=static_cast<int>(C.getNRowsAct()); //transpose does not affect
  int inc=1;
  if(ncolsC==1)//BLAS2 matrix vector multiply
    { 
      
      //printf("transA=%c transB=%c nrowsC=%d ninnnerA=%d ninnerB=%d A.getNRows()=%d\n",transA,transB,nrowsC,ninnerA,ninnerB,A.getNRows());
      
      //DGEMV_F77(&transA,&nrowsC,&ninnerA,&scaleAB,A.ptr(0,0),&nrowsC,B.ptr(0),&inc,&scaleC,C.ptr(0),&inc); 
      DGEMV_F77(&transA,&nrowsA,&ncolsA,&scaleAB,A.ptr(0,0),&lda,B.ptr(0,0),&inc,&scaleC,C.ptr(0,0),&inc);
      
    }
  else //BLAS3 matrix matrix multiply
    DGEMM_F77(&transA,&transB,&nrowsC,&ncolsC,&ninnerA,&scaleAB,A.ptr(0,0),
	      &lda,B.ptr(0,0),&ldb,&scaleC,C.ptr(0,0),&ldc);
  return C;
}


/***************************************************************************/
/**** The LAPACK wrappers start here                                    ****/ 
/***************************************************************************/


/// computes the L*D*L^T factorization of "matrix" with partial pivoting ("matrix" must be real and symmetric but does not need to be positive definite) wraps LAPACK subroutine DSYTRF (and also DLANGE and DSYCON to return the 1-norm rcond of "matrix", doing them together makes the rcond cheap)
MtxDbl& LDLT_fact(MtxDbl& matrix, MtxInt& ipvt, MtxDbl& scale, int& info, double& rcond)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  //printf("LDLT_fact() nrows=%d ncols=%d\n",nrows,ncols);
  char uplo='L';
  int lda=static_cast<int>(matrix.getNRowsAct());
  int info_local=0;
  ipvt.newSize(nrows,1);
  scale.newSize(nrows,1);

  int abspower;
  double log_of_2=std::log(2.0);
  for(int i=0; i<nrows; ++i) {
    abspower=static_cast<int>(std::floor(0.5+std::log(std::sqrt(matrix(i,i)))/log_of_2));
    scale(i,0)=std::pow(2.0,static_cast<double>(-abspower)); //this is the "numerically optimal" equilibration of a real symmetric positive definite matrix by "numerically optimal" I meant the analytically optimal (for reducing condition number) scaling has been rounded to the nearest power of 2 so that we don't lose any bits of accuracy due to rounding error due to equilibration (and yes I know that the matrix may not be positive definite, but we want it to work at least as good as equilibrated Cholesky for rspd matrix)
  }

  for(int j=0; j<nrows; ++j)
    for(int i=0; i<nrows; ++i)
      matrix(i,j)*=scale(i,0)*scale(j,0);


  //determine the needed amount of work space, and allocate it
  int lwork=-1; //lwork=-1 => tell me how much work space I need
  double workneeded;
  DSYTRF_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,ipvt.ptr(0,0),&workneeded,&lwork,&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local>=0);
#endif
  lwork = static_cast<int>(workneeded);
  if(lwork<2*nrows) lwork=2*nrows;
  MtxDbl work(lwork,1);
  MtxInt iwork(nrows,1);



  //calculate orignorm now so we don't have to store a copy of the unfactored matrix
  char whichnorm='1';  //DSYCON can only return the 1-norm rcond
  double orignorm=DLANGE_F77(&whichnorm,&nrows,&ncols,matrix.ptr(0,0),&lda,work.ptr(0,0));


  //perform the L*D*L^T factorization
  info_local=0;
  DSYTRF_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,ipvt.ptr(0,0),work.ptr(0,0),&lwork,&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local>=0);
#endif
  info=info_local;


  //the rcond of the real symmetric matrix, feed it orignorm, note that 
  info_local=0;
  DSYCON_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,ipvt.ptr(0,0),&orignorm,&rcond,work.ptr(0,0),iwork.ptr(0,0),&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local==0);
#endif

  return matrix;
}

/// inverts a real symmetric positive definite matrix in place after the L*D*L^T factorization has already been done, requires the lower triangular portion as input returns the full symmetric inverse (as opposed to just the lower triangular part), wraps LAPACK subroutine DSYTRI
MtxDbl& inverse_after_LDLT_fact(MtxDbl& matrix, const MtxInt& ipvt, const MtxDbl& scale)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  char uplo='L';
  int lda=static_cast<int>(matrix.getNRowsAct());
  int info=0;
  MtxDbl work(nrows,1);

  DSYTRI_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,ipvt.ptr(0,0),work.ptr(0,0),&info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  //fill in the top half of the inverse
  for(int j=0; j<ncols-1; ++j)
    for(int i=j+1; i<nrows; ++i) 
      matrix(j,i)=(matrix(i,j)*=scale(i,0)*scale(j,0));

  return matrix;
}

//computes the (1 norm) reciprocal of the condition number of A from the L*D*L^T factorization of A (A must be real and symmetric but does not need to be positive definite) wraps LAPACK subroutines DLANGE and DSYCON
double rcond_after_LDLT_fact(const MtxDbl& A, const MtxDbl& ALDLT, const MtxInt& ipvt)
{
  std::cerr << "rcond_after_LDLT_fact doesn't work because ALDLT was scaled during LDLT_fact" << std::endl;
  assert(false);
  double rconda;
  char whichnorm='1'; //DSYCON can only return the 1-norm rcond
  char uplo='L';
  int nrows=A.getNRows();
  int ncols=A.getNCols();
  int lda=ALDLT.getNRowsAct();
  MtxDbl work(2*nrows,1);
  MtxInt iwork(nrows,1);
  int info=0;
#ifdef __SURFMAT_ERR_CHECK__
  int nrows_ALDLT=ALDLT.getNRows();
  assert((nrows_ALDLT==nrows)&&(nrows_ALDLT==ALDLT.getNCols())&&(nrows==ncols));
#endif
  double anorm=DLANGE_F77(&whichnorm,&nrows,&ncols,A.ptr(0,0),&lda,work.ptr(0,0));

  DSYCON_F77(&uplo,&nrows,ALDLT.ptr(0,0),&lda,ipvt.ptr(0,0),&anorm,&rconda,work.ptr(0,0),iwork.ptr(0,0),&info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  return rconda;
}


/// solves A*X=B for X, where A is real symmetric and B={B || B^T}, without changing the contents of B, after A has been L*D*L^T factorized, ALDLT must contain the lower triangular portion of the factorization of A, wraps LAPACK subroutine DSYTRS 
MtxDbl& solve_after_LDLT_fact(MtxDbl& result, const MtxDbl& ALDLT, const MtxInt& ipvt, const MtxDbl& scale, const MtxDbl& BRHS,char transB)
{
  int n   = static_cast<int>(ALDLT.getNRows());
#ifdef __SURFMAT_ERR_CHECK__
  assert((((transB=='N') && (BRHS.getNRows()==n)) ||
	  ((transB=='T') && (BRHS.getNCols()==n))) &&
	 (scale.getNRows()==n) && (scale.getNCols()==1));
#endif
  int lda = static_cast<int>(ALDLT.getNRowsAct());
  int ldb = static_cast<int>(BRHS.getNRows());;
  char uplo='L';
  if(transB=='N'){ //copy B into the work matrix "result"
    /*
    result.newSize(BRHS.getNRows(),BRHS.getNCols());
    result.putTol(BRHS.getTol());
    for(int j=0; j<BRHS.getNCols(); j++)
      for(int i=0; i<BRHS.getNRows(); i++)
	result(i,j)=BRHS(i,j);
    */
    result = BRHS;
  }
  else{ //copy the transpose of B into work matrix "result"
    result.newSize(BRHS.getNCols(),BRHS.getNRows());
    result.putTol(BRHS.getTol());
    for(int i=0; i<BRHS.getNRows(); i++)
      for(int j=0; j<BRHS.getNCols(); j++)
	result(j,i)=BRHS(i,j);    
  }
  int nrhs= static_cast<int>(result.getNCols());

  for(int j=0; j<nrhs; ++j)
    for(int i=0; i<n; ++i)
      result(i,j)*=scale(i,0);

  //printf("solve_after_LDLT_fact: n=%d nrhs=%d\n",n,nrhs);
  int info=0;
  DSYTRS_F77(&uplo, &n, &nrhs, ALDLT.ptr(0,0), &lda, ipvt.ptr(0,0), result.ptr(0,0), &ldb, &info);

  for(int j=0; j<nrhs; ++j)
    for(int i=0; i<n; ++i)
      result(i,j)*=scale(i,0);

  return result;
}


///compute the pseudo inverse of A, this function wraps the LAPACK subroutine DGESVD in this context rcond is considered to be the ratio of the smallest to largest singular values, A becomes A_pinv
MtxDbl& pseudo_inverse(MtxDbl& A, double min_allowed_rcond, double& rcondA, double& log_abs_det, int& if_det_eq_zero) {
  int nrowsA = static_cast<int>(A.getNRows());
  int ncolsA = static_cast<int>(A.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrowsA>0)&&(ncolsA>0));
#endif  
  int lda = static_cast<int>(A.getNRowsAct());
  char jobu='S';
  char jobvt='S';
  int min_nrows_ncols_A=(nrowsA<ncolsA)?nrowsA:ncolsA;
  MtxDbl U(nrowsA,min_nrows_ncols_A);
  MtxDbl S(min_nrows_ncols_A,1);
  MtxDbl VT(min_nrows_ncols_A,ncolsA);
  int ldu =static_cast<int>(U.getNRowsAct());
  int ldvt=static_cast<int>(VT.getNRowsAct());
  int local_info;

  double size_of_work;
  int lwork=-1; //tell me how much work space I need;
  DGESVD_F77(&jobu,&jobvt,&nrowsA,&ncolsA,A.ptr(0,0),&lda,S.ptr(0,0),U.ptr(0,0),&ldu,VT.ptr(0,0),&ldvt,&size_of_work,&lwork,&local_info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(local_info==0);
#endif
  lwork=static_cast<int>(size_of_work);
  MtxDbl work(lwork,1);
  
  DGESVD_F77(&jobu,&jobvt,&nrowsA,&ncolsA,A.ptr(0,0),&lda,S.ptr(0,0),U.ptr(0,0),&ldu,VT.ptr(0,0),&ldvt,work.ptr(0,0),&lwork,&local_info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(local_info==0);
#endif

  rcondA=S(min_nrows_ncols_A-1,0)/S(0,0); //this is actually rcond if A is symmetric

  log_abs_det=0.0;
  if_det_eq_zero=0;
    
  if(S(0,0)==0.0) {
    A.zero();
    if_det_eq_zero=1;
  }
  else{
    double min_allowed_sing=S(0,0)*min_allowed_rcond;
    double inv_sing;
    for(int j=0; j<min_nrows_ncols_A; ++j) {

      if(S(j,0)>0.0)
	log_abs_det+=std::log(S(j,0));
      else
	if_det_eq_zero=1;

      if(min_allowed_sing<=S(j,0)) {
	inv_sing=1.0/S(j,0);
	for(int i=0; i<nrowsA; ++i)
	  U(i,j)*=inv_sing;
      }
      else{
	for(int i=0; i<nrowsA; ++i)
	  U(i,j)=0.0;
      }
    }
    matrix_mult(A,VT,U,0.0,1.0,'T','T');
  }

  return A;
}



/// perform a numerically optimally equilibrated Cholesky factorization of a real symmetric positive semi-definite matrix, wraps DPOTRF, returns the lower triangular portion, the pequilibration is "undone" so the L matrix is the same (minus some rounding error due to poor conditioning) as it would be without equilibration, if info>0 then the equilibrated matrix is singular, rcond is the standard fast estimate of the reciprocal of the condition number of the numerically optimally equilibrated real symmetric positive definite matrix (that's what affects the round off error)
MtxDbl& Chol_fact(MtxDbl& matrix, int& info, double& rcond)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  //printf("Chol_fact() nrows=%d ncols=%d\n",nrows,ncols);
  char uplo='L';
  int lda=static_cast<int>(matrix.getNRowsAct());;
  int info_local=0;

  MtxDbl work(3*nrows,1);
  MtxInt iwork(nrows,1);

  int abspower;
  int minabspower;
  int maxabspower;
  MtxDbl scalefactor(nrows,1);
  double log_of_2=std::log(2.0);
  abspower=static_cast<int>(std::floor(0.5+std::log(std::sqrt(matrix(0,0)))/log_of_2));
  scalefactor(0,0)=std::pow(2.0,static_cast<double>(-abspower));
  //abspower=log2(sqrt(matrix(0,0)));
  //scalefactor(0,0)=1.0/sqrt(matrix(0,0));
  minabspower=maxabspower=abspower;
  for(int i=1; i<nrows; ++i) {
    abspower=static_cast<int>(std::floor(0.5+std::log(std::sqrt(matrix(i,i)))/log_of_2));
    scalefactor(i,0)=std::pow(2.0,static_cast<double>(-abspower)); //this is the "numerically optimal" equilibration of a real symmetric positive definite matrix by "numerically optimal" I meant the analytically optimal (for reducing condition number) scaling has been rounded to the nearest power of 2 so that we don't lose any bits of accuracy due to rounding error due to equilibration
    //abspower=log2(sqrt(matrix(i,i)));
    //scalefactor(i,0)=1.0/sqrt(matrix(i,i));
    minabspower=(abspower<minabspower)?abspower:minabspower;
    maxabspower=(abspower>maxabspower)?abspower:maxabspower;   
  }
  
  if(maxabspower!=minabspower) {
    //only do the equilibration if the maximum and minimum numerically optimal scaling factors are different (by a factor of 2 or more) because otherwise the real symmetric positive definite matrix is already numerically optimally equilibrated.
    for(int j=0; j<nrows; ++j)
      for(int i=0; i<nrows; ++i)
	matrix(i,j)*=scalefactor(i,0)*scalefactor(j,0);
  }

  //calculate orignorm now so we don't have to store a copy of the unfactored matrix
  char whichnorm='1';
  double orignorm=DLANGE_F77(&whichnorm,&nrows,&ncols,matrix.ptr(0,0),&lda,work.ptr(0,0));

  DPOTRF_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local>=0);
#endif
  info=info_local;

  //the rcond of the "numerically optimally" equilibrated real symmetric positive definte matrix, feed it orignorm
  DPOCON_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,&orignorm,&rcond,work.ptr(0,0),iwork.ptr(0,0),&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local==0);
#endif

  if(maxabspower!=minabspower) {
    //undo the "numerically optimal" equilibration of the real symmetric positive definite matrix, that is other than possibly avoiding rounding error due to poorly condition matrix this function produces the same cholesky "L" matrix that you would get without equilibration.
    for(int i=0; i<nrows; ++i)
      scalefactor(i,0)=1.0/scalefactor(i,0); //multiplication can be faster than division

    for(int j=0; j<nrows; ++j)
      for(int i=j; i<nrows; ++i)  //it's lower triangular
	matrix(i,j)*=scalefactor(i,0);
  }

  return matrix;
}

MtxDbl& Chol_fact_workspace(MtxDbl& matrix, MtxDbl& scalefactor, 
			    MtxDbl& work, MtxInt& iwork, int& info, 
			    double& rcond)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  //printf("Chol_fact() nrows=%d ncols=%d\n",nrows,ncols);
  char uplo='L';
  int lda=static_cast<int>(matrix.getNRowsAct());;
  int info_local=0;

  work.newSize(3*nrows,1);
  iwork.newSize(nrows,1);

  int abspower;
  int minabspower;
  int maxabspower;
  scalefactor.newSize(nrows,1);
  double log_of_2=std::log(2.0);
  abspower=static_cast<int>(std::floor(0.5+std::log(std::sqrt(matrix(0,0)))/log_of_2));
  scalefactor(0,0)=std::pow(2.0,static_cast<double>(-abspower));
  //abspower=log2(sqrt(matrix(0,0)));
  //scalefactor(0,0)=1.0/sqrt(matrix(0,0));
  minabspower=maxabspower=abspower;
  for(int i=1; i<nrows; ++i) {
    abspower=static_cast<int>(std::floor(0.5+std::log(std::sqrt(matrix(i,i)))/log_of_2));
    scalefactor(i,0)=std::pow(2.0,static_cast<double>(-abspower)); //this is the "numerically optimal" equilibration of a real symmetric positive definite matrix by "numerically optimal" I meant the analytically optimal (for reducing condition number) scaling has been rounded to the nearest power of 2 so that we don't lose any bits of accuracy due to rounding error due to equilibration
    //abspower=log2(sqrt(matrix(i,i)));
    //scalefactor(i,0)=1.0/sqrt(matrix(i,i));
    minabspower=(abspower<minabspower)?abspower:minabspower;
    maxabspower=(abspower>maxabspower)?abspower:maxabspower;   
  }
  
  if(maxabspower!=minabspower) {
    //only do the equilibration if the maximum and minimum numerically optimal scaling factors are different (by a factor of 2 or more) because otherwise the real symmetric positive definite matrix is already numerically optimally equilibrated.
    for(int j=0; j<nrows; ++j)
      for(int i=0; i<nrows; ++i)
	matrix(i,j)*=scalefactor(i,0)*scalefactor(j,0);
  }

  //calculate orignorm now so we don't have to store a copy of the unfactored matrix
  char whichnorm='1';
  double orignorm=DLANGE_F77(&whichnorm,&nrows,&ncols,matrix.ptr(0,0),&lda,work.ptr(0,0));

  DPOTRF_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local>=0);
#endif
  info=info_local;

  //the rcond of the "numerically optimally" equilibrated real symmetric positive definte matrix, feed it orignorm
  DPOCON_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,&orignorm,&rcond,work.ptr(0,0),iwork.ptr(0,0),&info_local);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info_local==0);
#endif

  if(maxabspower!=minabspower) {
    //undo the "numerically optimal" equilibration of the real symmetric positive definite matrix, that is other than possibly avoiding rounding error due to poorly condition matrix this function produces the same cholesky "L" matrix that you would get without equilibration.
    for(int i=0; i<nrows; ++i)
      scalefactor(i,0)=1.0/scalefactor(i,0); //multiplication can be faster than division

    for(int j=0; j<nrows; ++j)
      for(int i=j; i<nrows; ++i)  //it's lower triangular
	matrix(i,j)*=scalefactor(i,0);
  }

  return matrix;
}

/// inverts a real symmetric positive definite matrix in place after the Cholesky factorization has already been done, requires the lower triangular portion as input returns the full symmetric inverse (as opposed to just the lower triangular part), wraps DPOTRF
MtxDbl& inverse_after_Chol_fact(MtxDbl& matrix)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  char uplo='L';
  int lda=static_cast<int>(matrix.getNRowsAct());
  int info=0;

  DPOTRI_F77(&uplo,&nrows,matrix.ptr(0,0),&lda,&info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  //fill in the top half of the inverse
  for(int j=0; j<ncols-1; ++j)
    for(int i=j+1; i<nrows; ++i)
      matrix(j,i)=matrix(i,j);

  return matrix;
}

//computes the (1 norm) reciprocal of the condition number of A from the Cholesky factorization of A (A must be real symmetric and positive semi-definite)
 double rcond_after_Chol_fact(const MtxDbl& A, const MtxDbl& AChol)
{
  double rconda;
  char whichnorm='1';
  char uplo='L';
  int nrows=A.getNRows();
  int ncols=A.getNCols();
  int lda=static_cast<int>(A.getNRowsAct());
  int ld_AChol=static_cast<int>(AChol.getNRowsAct());
  MtxDbl work(3*nrows,1);
  MtxInt iwork(nrows,1);
  int info;
#ifdef __SURFMAT_ERR_CHECK__
  int nrows_AChol=AChol.getNRows();
  assert((nrows_AChol==nrows)&&(nrows_AChol==AChol.getNCols())&&(nrows==ncols));
#endif
  double anorm=DLANGE_F77(&whichnorm,&nrows,&ncols,A.ptr(0,0),&lda,work.ptr(0,0));
  DPOCON_F77(&uplo,&nrows,AChol.ptr(0,0),&ld_AChol,&anorm,&rconda,work.ptr(0,0),iwork.ptr(0,0),&info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  return rconda;
}

  //MtxDbl& debug_solve_after_Chol_fact(MtxDbl& result, const MtxDbl& AChol, const MtxDbl& BRHS,char transB='N');

/// solves A*X=B for X, where A is symmetric positive definite and B={B || B^T}, without changing the contents of B, after A has been Cholesky factorized, AChol must contain the lower triangular portion of the factorization of A, wraps DPOTRS 
MtxDbl& solve_after_Chol_fact(MtxDbl& result, const MtxDbl& AChol, const MtxDbl& BRHS,char transB)
//#ifdef __SURFMAT_ERR_CHECK__
//{
//return debug_solve_after_Chol_fact(result, AChol, BRHS,transB);
//}
//#else
{
  int n   = static_cast<int>(AChol.getNRows());
#ifdef __SURFMAT_ERR_CHECK__
  if(!(((transB=='N') && (BRHS.getNRows()==n)) ||
       ((transB=='T') && (BRHS.getNCols()==n)))) {
    printf("solve_after_Chol_fact transB='%c' size(Achol)=[%d x %d] size(BRHS)=[%d x %d]",
	   transB,n,n,BRHS.getNRows(),BRHS.getNCols());
    fflush(stdout);
    printf("\n");
    assert(((transB=='N') && (BRHS.getNRows()==n)) ||
	   ((transB=='T') && (BRHS.getNCols()==n)));
  }
#endif
  int lda = static_cast<int>(AChol.getNRowsAct());
  char uplo='L';
  if(transB=='N'){ //copy B into the work matrix "result"
    /*
    result.newSize(BRHS.getNRows(),BRHS.getNCols());
    result.putTol(BRHS.getTol());
    for(int j=0; j<BRHS.getNCols(); j++)
      for(int i=0; i<BRHS.getNRows(); i++)
	result(i,j)=BRHS(i,j);
    */
    result = BRHS;
  }
  else{ //copy the transpose of B into work matrix "result"
    result.newSize(BRHS.getNCols(),BRHS.getNRows());
    result.putTol(BRHS.getTol());
    for(int i=0; i<BRHS.getNRows(); i++)
      for(int j=0; j<BRHS.getNCols(); j++)
	result(j,i)=BRHS(i,j);    
  }
  int nrhs= static_cast<int>(result.getNCols());
  int ldb = static_cast<int>(result.getNRowsAct());

  
  //printf("solve_after_Chol_fact: n=%d nrhs=%d\n",n,nrhs);
  int info=0;
  DPOTRS_F77(&uplo, &n, &nrhs, AChol.ptr(0,0), &lda, result.ptr(0,0), &ldb, &info);
  return result;
}
  //#endif

/// perform an LU factorization with partial pivoting, wraps LAPACK subroutine DGETRF
MtxDbl& LU_fact(MtxDbl& matrix, MtxInt& ipvt)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  // dgetrf may reorder the rows, the mapping between old and new rows
  // is returned in ipvt.
  if((ipvt.getNRows()!=nrows)||(ipvt.getNCols()!=1))
    ipvt.newSize(nrows,1); 
  
  int lda = static_cast<int>(matrix.getNRowsAct());
  int info = 0;
  //printf("nrows=%d ncols=%d\n",nrows,ncols);
  //std::cout << "Matrix size: " << n_rows << " " << n_cols << std::endl;
  DGETRF_F77(&nrows,&ncols,matrix.ptr(0,0),&lda,ipvt.ptr(0,0),&info);
  //std::cout << "Done with dgetrf" << std::endl;
  return matrix;
}

/// inverts a matrix, in place, after an LU factorization has already been done, wraps DGETRI
MtxDbl& inverse_after_LU_fact(MtxDbl& matrix, MtxInt& ipvt)
{
  int nrows = static_cast<int>(matrix.getNRows());
  int ncols = static_cast<int>(matrix.getNCols());
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrows==ncols);
#endif
  int lwork = ncols;  // should be optimal blocksize
  MtxDbl work(lwork,1);
  int lda = static_cast<int>(matrix.getNRowsAct());
  int info = 0;
  //std::cout << "Matrix size: " << n_rows << " " << n_cols << std::endl;
  DGETRI_F77(&nrows,matrix.ptr(0,0),&lda,ipvt.ptr(0,0),work.ptr(0,0),&lwork,&info);
  //std::cout << "Done with getri" << std::endl;
  return matrix;
}



//computes the (1 norm) reciprocal of the condition number of A from the LU factorization of A
double rcond_after_LU_fact(const MtxDbl& A, const MtxDbl& ALU) {
  double rcond;
  char whichnorm='1';
  int nrows_A=A.getNRows();
  int ncols_A=A.getNCols();
  int LDA=ALU.getNRowsAct();
  int nrows_ALU=ALU.getNRows();
  int LD_ALU=ALU.getNRowsAct();

  MtxDbl work(4*ncols_A,1); 
  MtxInt iwork(ncols_A,1);
  int info;
#ifdef __SURFMAT_ERR_CHECK__
  int ncols_ALU=ALU.getNCols();
  assert((nrows_ALU==ncols_ALU)&&(nrows_ALU==nrows_A)&&(nrows_A==ncols_A));
#endif
  double anorm=DLANGE_F77(&whichnorm,&nrows_A,&ncols_A,A.ptr(0,0),&LDA,
			  work.ptr(0,0));
  DGECON_F77(&whichnorm,&nrows_ALU,ALU.ptr(0,0),&LD_ALU,&anorm,&rcond,
	     work.ptr(0,0), iwork.ptr(0,0),&info);
  return rcond;
}

/// solves A*X=B for X, where A={A || A^T} and B={B || B^T}, without changing the contents of B, after A has been LU factorized, wraps DGETRS 
MtxDbl& solve_after_LU_fact(MtxDbl& result, const MtxDbl& ALU, const MtxInt& ipvt, 
			    const MtxDbl& BRHS, char transA, char transB)
{
  int nrows_ALU   = static_cast<int>(ALU.getNRows());
#ifdef __SURFMAT_ERR_CHECK__
  assert((transA=='N' || transA=='T')&&(nrows_ALU == ALU.getNCols())&&
	 (nrows_ALU=ipvt.getNRows()));
  assert(((transB=='N') && (BRHS.getNRows()==nrows_ALU)) ||
	 ((transB=='T') && (BRHS.getNCols()==nrows_ALU)));
#endif
  int lda = static_cast<int>(ALU.getNRowsAct());
  int ldb = static_cast<int>(BRHS.getNRowsAct());
  if(transB=='N'){ //copy B into the work matrix "result"
    /*
    result.newSize(BRHS.getNRows(),BRHS.getNCols());
    result.putTol(BRHS.getTol());
    for(int j=0; j<BRHS.getNCols(); j++)
      for(int i=0; i<BRHS.getNRows(); i++)
	result(i,j)=BRHS(i,j);
    */
    result = BRHS;
  }
  else{ //copy the transpose of B into work matrix "result"
    result.newSize(BRHS.getNCols(),BRHS.getNRows());
    result.putTol(BRHS.getTol());
    for(int i=0; i<BRHS.getNRows(); i++)
      for(int j=0; j<BRHS.getNCols(); j++)
	result(j,i)=BRHS(i,j);    
  }
  int nrhs= static_cast<int>(result.getNCols());
  
  //printf("solve_after_LU_fact: n=%d nrhs=%d\n",n,nrhs);
  int info=0;
  DGETRS_F77(&transA, &nrows_ALU, &nrhs, ALU.ptr(0,0), &lda, ipvt.ptr(0,0), result.ptr(0,0), &ldb, &info);
  return result;
}


void least_squares(MtxDbl& A, MtxDbl& x, MtxDbl& b)
{
  // # Rows in A must == # Rows in b
#ifdef __SURFMAT_ERR_CHECK__
  assert(A.getNRows() == b.getNElems()); 
  // System must be square or over-constrained
  assert(A.getNRows() >= A.getNCols());
#endif
  int nrowsA = static_cast<int>(A.getNRows());
  int ncolsA = static_cast<int>(A.getNCols());
  int lda = static_cast<int>(A.getNRowsAct());
  // Client may supply a "blank" initialized vector for x
  int lwork = nrowsA*ncolsA * 2;
  MtxDbl work(lwork,1);
  // values must be passed by reference to Fortran, so variables must be 
  // declared for info, nrhs, trans
  int info;
  int nrhs=1;
  char trans = 'N';
  x = b; //preserve the original b
  int ldb=x.getNRowsAct();
  DGELS_F77(&trans,&nrowsA,&ncolsA,&nrhs,A.ptr(0,0),&lda,x.ptr(0,0),
	    &ldb,work.ptr(0,0),&lwork,&info);
  x.reshape(ncolsA,1);
}


void least_squares_with_equality_constraints(MtxDbl& A, 
     MtxDbl& x, MtxDbl& c, MtxDbl& B, MtxDbl& d)
{
  int nrowsA = static_cast<int>(A.getNRows());
  int ncolsA = static_cast<int>(A.getNCols());
  int nrowsB = static_cast<int>(B.getNRows());
  int lda = static_cast<int>(A.getNRowsAct());
  int ldb = static_cast<int>(B.getNRowsAct());
#ifdef __SURFMAT_ERR_CHECK__
  int ncolsB = static_cast<int>(B.getNCols());
  assert(ncolsB == ncolsA);
  assert(nrowsB <= ncolsA);
  assert(ncolsA <= nrowsB + nrowsA);
  assert(x.getNElems() == ncolsA);
#endif
  int lwork = (nrowsA + ncolsA + nrowsB);
  lwork *= lwork;
  ///\todo Compute optimal blocksize before running dgglse
    MtxDbl work(lwork,1);
  int info = 0;
  DGGLSE_F77(&nrowsA,&ncolsA,&nrowsB,A.ptr(0,0),&lda,B.ptr(0,0),&ldb,c.ptr(0,0),d.ptr(0,0),x.ptr(0,0),work.ptr(0,0),&lwork,
	     &info);
#ifdef __SURFMAT_ERR_CHECK__
  if (info != 0) {
    printf("least_squares_with_equality_constraints: Error encountered in DGGLSE\n");
    assert(info == 0);
  }
#endif

  //if (info != 0) throw string("Error in dgglse\n");
}

///finds the eigenvalues and optionally (by default) eigenvectors of a real symmetric matrix, returns a reference to the vector of eigenvalues, this function wraps the LAPACK subroutine DSYEV
MtxDbl& eig_sym(MtxDbl& eigvect, MtxDbl& eigval, const MtxDbl& A, char jobz) {
  char uplo='L';
  eigvect.copy(A);
  int nrowsA=static_cast<int>(eigvect.getNRows());
  int lda=static_cast<int>(eigvect.getNRowsAct()); //because a matrix's apparent and actual (memory footprint) size can be different and "copy" only ensures that the apparent size is the same, we need to get lda from eigvect instead of A, getting lda from was causing a segfault when Pivoting Cholesky is used to discard points
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<nrowsA)&&(nrowsA==eigvect.getNCols()));
#endif
  eigval.newSize(nrowsA,1); eigval.zero(); //zero is for sanity check
  int info;
  int lwork=-1;
  double work_size;
  DSYEV_F77(&jobz, &uplo, &nrowsA, eigvect.ptr(0,0), &lda, eigval.ptr(0,0), 
	    &work_size, &lwork, &info);
#ifdef __SURFMAT_ERR_CHECK__
  assert(info==0);
#endif
  lwork=static_cast<int>(work_size);
  MtxDbl work(lwork,1);
  DSYEV_F77(&jobz, &uplo, &nrowsA, eigvect.ptr(0,0), &lda, eigval.ptr(0,0), 
	    work.ptr(0,0), &lwork, &info);

  return eigval;
}

/*****************************************************************************/
/* extra functions added for convenience                                     */
/*****************************************************************************/


MtxDbl& pseudo_inverse_sym(MtxDbl& A, double min_allowed_rcond, double& rcond_A, double& log_abs_det_A, double& sign_of_det_A) {
  int nrowsA=A.getNRows();
#ifdef __SURFMAT_ERR_CHECK__
  assert(nrowsA==A.getNCols());
#endif
  MtxDbl scale(nrowsA,1);
  int abspower;
  double log_of_2=std::log(2.0);
  for(int i=0; i<nrowsA; ++i) {
    abspower=static_cast<int>(std::floor(0.5+std::log(std::sqrt(std::fabs(A(i,i))))/log_of_2));
    scale(i,0)=std::pow(2.0,static_cast<double>(-abspower));
    log_abs_det_A-=std::log(scale(i,0));
  }
  log_abs_det_A*=2.0;
  for(int j=0; j<nrowsA; ++j)
    for(int i=0; i<nrowsA; ++i)
      A(i,j)*=scale(i,0)*scale(j,0);

  MtxDbl eigvect(nrowsA,nrowsA), eigval(nrowsA,1);  
  eig_sym(eigvect, eigval, A);

  sign_of_det_A=1.0;
  MtxDbl eigvect_eigval_inv(nrowsA,nrowsA);
  double smallest_eigval=std::fabs(eigval(0,0));
  for(int j=1; j<nrowsA; ++j) {
    sign_of_det_A*=(static_cast<double> ((eigval(j,0)>0.0) - (eigval(j,0)<0.0)));
    if(std::fabs(eigval(j,0))>smallest_eigval) {
      log_abs_det_A+=std::log(std::fabs(eigval(j,0)));
      smallest_eigval=std::fabs(eigval(j,0));  
    }
  }
  smallest_eigval*=min_allowed_rcond;
  
  double eigval_inv;
  for(int j=0; j<nrowsA; ++j) 
    if(smallest_eigval<=std::fabs(eigval(j,0))) {
      eigval_inv=1.0/eigval(j,0);
      for(int i=0; i<nrowsA; ++i) 
	eigvect_eigval_inv(i,j)=eigvect(i,j)*eigval_inv;
    }
    else
      for(int i=0; i<nrowsA; ++i) 	
	eigvect_eigval_inv(i,j)=eigvect(i,j)=0.0;
  
  matrix_mult(A,eigvect_eigval_inv,eigvect,0.0,1.0,'N','T');

  for(int j=0; j<nrowsA; ++j)
    for(int i=0; i<nrowsA; ++i)
      A(i,j)*=scale(i,0)*scale(j,0);

  return A;
}

/// computes log(fabs(det(A))) and sign(det(A)) from the L*D*L^T factorization of A by adding logs of absolute value of determinants of component diagonal blocks to prevent underflow/overflow errors that would by taking product of diagonal block determinants before taking the log.  You must have called LDLT_fact() before calling this function
double log_det_after_LDLT_fact(const MtxDbl& ALDLT, const MtxInt& ipvt, const MtxDbl& scale, double& sign_of_det) {
  //if A=L*D*L^T then det(A)=det(L)*det(D)*det(L^T)=det(L)^2*det(D), the LDLT factorization assures that det(L)=1 (because it is lower triangular and has ones on the diagonal) so det(A) = det(D).  D is block diagonal so its determinant is the product of the determinant of its diagonal blocks.  The blocks are all 1x1 or 2x2, and the block sizes are indicated by the contents of ipvt (see the documentation for the LAPACK subroutine DSYTRF for more details), note that A and (therfore) D are symmetric, and only the lower triangular parts of them are stored in ALDLT.

  int nrowsA   = static_cast<int>(ALDLT.getNRows());
#ifdef __SURFMAT_ERR_CHECK__
  assert((ALDLT.getNCols()==nrowsA)&&(scale.getNRows()==nrowsA)&&
	 (scale.getNCols()==1));
#endif
  double log_det_A=0.0;
  for(int i=0; i<nrowsA; ++i)
    log_det_A-=std::log(scale(i,0));
  log_det_A*=2.0;

  sign_of_det=1.0;
  double block2_det;
  for(int i=0; i<nrowsA; ) {
#ifdef __SURFMAT_ERR_CHECK__
    assert(ipvt(i,0)!=0);
#endif
    if(ipvt(i,0)>0) {
      //sign_of_det*=sign(ALDLT(i,i))
      sign_of_det*=(static_cast<double> ((ALDLT(i,i)>0.0) - (ALDLT(i,i)<0.0)));       
      log_det_A+=std::log(std::fabs(ALDLT(i,i)));
      i++;
    }
    else if(ipvt(i,0)==ipvt(i+1,0)){
#ifdef __SURFMAT_ERR_CHECK__
      assert(i<nrowsA-1);
#endif
      block2_det=ALDLT(i,i)*ALDLT(i+1,i+1)-ALDLT(i+1,i)*ALDLT(i+1,i);
      sign_of_det*=(static_cast<double> ((block2_det>0.0) - (block2_det<0.0)));       
      log_det_A+=std::log(std::fabs(block2_det));      
      i+=2;
    }
    else{
      std::cerr << "in log_det_after_LDLT_fact must have ipvt(i,0)>0 or ipvt(i,0)==ipvt(i+1,0) (in latter case we do i+=2, so we never see the 'second negative' i.e. don't see ipvt(i+1)<0 when ipvt(i)==ipvt(i+1)<0.  See LAPACK DSYTRF for more details about diagonal block size of 2 when UPLO='L'" << std::endl;
      assert((ipvt(i,0)>0)||(ipvt(i,0)==ipvt(i+1,0)));
    }
  }
  return log_det_A;
}


} // end namespace nkm
