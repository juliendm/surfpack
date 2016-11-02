      SUBROUTINE PIVOTCHOL( UPLO, N, A, LDA, PIV, RANK, TOL, INFO )
      IMPLICIT NONE
C     
C     Keith Dalbey, Sandia National Laboratories, June, 2011 some 
C     code taken from LEV3PCHOL by Craig Lucas, University of 
C     Manchester. January, 2004 who took some code from LAPACK 
C     routine DPOTF2, This is Keith Dalbey's implementation of 
c     Craig Lucas's algorithm (intended to speed up the 
C     decomposition, for my test problem which was a 5500 by 5500 
C     matrix, my code took 53.6001 seconds, Craig's took 57.1924,
C     LAPACK's non-pivoting cholesky took 8.74516, I think that
C     Craig's level 3 code may have defaulted to his level 2 code 
C     for the times I quoted)
C     
C     .. Scalar Arguments ..
      DOUBLE PRECISION   TOL
      INTEGER            INFO, LDA, N, RANK
      CHARACTER          UPLO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
      INTEGER            PIV( N )
C     ..
C
C Purpose
C =======
C
C LEV3PCHOL computes the Cholesky factorization with complete
C pivoting of a real symmetric positive semi-definite matrix A.
C
C The factorization has the form
C     P' * A * P = U' * U , if UPLO = 'U',
C     P' * A * P = L * L', if UPLO = 'L',
C where U is an upper triangular matrix and L is lower triangular, and
C P is stored as vector PIV.
C
C This algorithm does not attempt to check that A is positive
C semi-definite. This version of the algorithm calls level 3 BLAS.
C
C Arguments
C =========
C
C UPLO     (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored.
C          = 'U': Upper triangular
C          = 'L': Lower triangular
C
C N        (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C A        (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the symmetric matrix A. If UPLO = 'U', the leading
C          n by n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced. If UPLO = 'L', the
C          leading n by n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C
C          On exit, if INFO = 0, the factor U or L from the Cholesky
C          factorization as above.
C
C PIV      (output) INTEGER array, dimension (N)
C          PIV is such that the nonzero entries are P( PIV( K ), K ) = 1.
C
C RANK     (INPUT/output) INTEGER
C          On output, this is the rank of A given by the number of 
c          steps the algorithm completed.
C          RANK BEING NEGATIVE ON ENTRY IS A FLAG THAT MEANS THE USER
C          ONLY WANTS THE FIRST -RANK VALUES (FOR SUBSET SELECTION)
C
C TOL      (input) DOUBLE PRECISION
C          User defined tolerance. If TOL < 0, then N*U*MAX( A(k,k) )
C          will be used. The algorithm terminates at the (k-1)st step
C          if the pivot <= TOL.
C
C LDA      (input) INTEGER
C          The leading dimension of the array A.   LDA >= max(1,N).
C
C INFO     (output) INTEGER
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          = 0 algorithm completed successfully.
C
C =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION    ONE, ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION    PIVOT, TMPDBL, MXPVIN, TMPDBL2, ALPHA,
     $     BETA 
      INTEGER             IPIV, ITOADD, I, J, TMPINT, INC1, M,
     $     NSTOP
      LOGICAL             UPNLOW
      CHARACTER           TRANS

C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DDOT
      LOGICAL            LSAME
      EXTERNAL           DDOT, LSAME

C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, DSWAP, DGEMV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX0, DSQRT
C     ..
C     
C     Test the input parameters.
C     

c      PRINT*, 'RANK=', RANK 

C     RANK BEING NEGATIVE ON ENTRY IS A FLAG THAT MEANS THE USER
C     ONLY WANTS THE FIRST -RANK VALUES (FOR SUBSET SELECTION)
      NSTOP=-RANK
      IF((RANK.GE.0).OR.(NSTOP.GT.N)) THEN
c         PRINT*, 'RANK SET TO N'
         NSTOP=N
      ENDIF
c      NSTOP=N
c      call flush(8)
      INFO = 0
C
C     MAKE SURE WE HAVE BOTH HALVES OF THE SYMMETRIC MATRIX SO WE CAN
C     SPEED UP MEMORY ACCESS LATER
      UPNLOW = LSAME( UPLO,'B')
      IF(LSAME(UPLO,'L')) THEN
C     IF WE ONLY HAVE THE LOWER HALF COPY IT TO THE UPPER HALF
         DO J=2,N
            DO I=1,J-1
               A(I,J)=A(J,I)
            ENDDO
         ENDDO
         UPNLOW=.TRUE.
      ELSEIF(LSAME(UPLO,'U')) THEN
C     IF WE ONLY HAVE THE UPPER HALF COPY IT TO THE LOWER HALF
         DO J=1,N-1
            DO I=J+1,N
               A(I,J)=A(J,I)
            ENDDO
         ENDDO
         UPNLOW=.TRUE.
      ENDIF
C
      IF( .NOT.UPNLOW) THEN
C     THEY ENTERED SOMETHING OTHER THAN 'B', 'L', OR 'U' FOR UPLO
         INFO = -1
      ELSEIF( N.LT.0 ) THEN
         INFO = -2
      ELSEIF( LDA.LT.MAX0( 1, N ) ) THEN
         INFO = -4
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'LEV3PCHOL', -INFO )
         RETURN
      ENDIF
C     
C     Quick return if possible
C     
      INC1=1
      TRANS='T'

      IF( N.EQ.0 )
     $     RETURN
C      
C     THE FIRST COLUMN OF PIVOTING CHOLESKY IS DIFFERENT THAN ALL 
C     OTHERS SO DO IT SEPARATELY (DOING IT SEPARATELY IS FASTER
C     THAN USING IF'S INSIDE THE LOOP)    
C
C     FILL UP THE PIV (PIVOT) INTEGER ARRAY AND FIND THE FIRST
C     PIVOT
      PIVOT=A(1,1)
      IPIV=1
      PIV(1)=1
      DO J=2,N
         PIV(J)=J
         IF(PIVOT.LT.A(J,J)) THEN
            IPIV=J
            PIVOT=A(J,J)
         ENDIF
      ENDDO
C
      IF(PIVOT.LE.ZERO) THEN
C     ALL OF THE DIAGONALS ARE LESS THAN OR EQUAL ZERO SO THIS MATRIX
C     IS NOT POSITIVE DEFINITE SO EXIT EARLY
         RANK=0
         RETURN
      ENDIF
      RANK=1
      ALPHA=-ONE
C
      IF(IPIV.NE.1) THEN
C     WE NEED TO SWAP ROWS/COLUMNS 1 AND IPIV
         PIV(IPIV)=1
         PIV(1)=IPIV
C
C     SWAP COLUMNS 1 AND IPIV 
         CALL DSWAP(N,A(1,1),INC1,A(1,IPIV),INC1)
C
C     SWAP ROWS 1 AND IPIV 
         CALL DSWAP(N,A(1,1),LDA,A(IPIV,1),LDA)
      ENDIF
C
C     WE NEED 1.0/MAXPIV TO HELP US ESTIMATE RCOND LATER
      MXPVIN=ONE/PIVOT
      BETA=DSQRT(PIVOT)
      A(1,1)=BETA
      BETA=ONE/BETA
      PIVOT=ZERO
      IPIV=2
      DO I=2,N
         TMPDBL2=A(I,1)*BETA
         A(I,1)=TMPDBL2
         A(1,I)=TMPDBL2
         TMPDBL=A(I,I)-TMPDBL2**2
         A(I,I)=TMPDBL
         IF(PIVOT.LT.TMPDBL) THEN
C     FIND THE NEXT LARGEST PIVOT
            PIVOT=TMPDBL
            IPIV=I
         ENDIF
      ENDDO
C
C     NOW DO THE PIVOTING CHOLESKY FOR THE REST OF THE ROWS/
C     COLUMNS
      DO ITOADD=2,NSTOP
C     IF OUR PIVOT ESTIMATE OF "RCOND" IS LESS THAN THE TOLERANCE
C     THEN EXIT THE SUBROUTINE
         TMPDBL=PIVOT*MXPVIN
         IF(TMPDBL.LE.TOL) 
     $        RETURN

         IF(IPIV.NE.ITOADD) THEN
C     WE NEED TO PIVOT: SWAP ROWS/COLUMNS IPIV AND ITOADD
            TMPINT     =PIV(ITOADD)
            PIV(ITOADD)=PIV(IPIV  )
            PIV(IPIV  )=TMPINT
C
C     SWAP COLUMNS IPIV AND ITOADD
            CALL DSWAP(N,A(1,ITOADD),INC1,A(1,IPIV),INC1)
C
C     SWAP ROWS IPIV AND ITOADD
            CALL DSWAP(N,A(ITOADD,1),LDA,A(IPIV,1),LDA)
         ENDIF
C         
C     UPDATE THE CHOLESKY DECOMPOSITION
         BETA=DSQRT(PIVOT)
         A(ITOADD,ITOADD)=BETA
         BETA=ONE/BETA
         PIVOT=ZERO
         IPIV=ITOADD+1
         M=N-ITOADD

C     TECHNICALLY THE DGEMV SHOULD "LOOP" ACROSS ROWS NOT DOWN 
C     COLUMNS BUT WE CAN DO THIS BECAUSE WE ALSO STORE THE UPPER 
C     TRIANGULAR REFLECTION OF THE LOWER TRIANGULAR CHOLESKY 
C     DECOMPOSITION AND THE REASON WE DO THIS IS SO WE CAN "LOOP" 
C     DOWN COLUMNS WHICH HAS FASTER MEMORY ACCESS
C     REMEMBER THAT RANK=ITOADD-1
         CALL DGEMV(TRANS,RANK,M,ALPHA,A(1,ITOADD+1),LDA,A(1,ITOADD),
     $        INC1,ONE,A(ITOADD,ITOADD+1),LDA)

         DO I=ITOADD+1,N
            TMPDBL2=A(ITOADD,I)*BETA
            A(ITOADD,I)=TMPDBL2
            TMPDBL=A(I,I)-TMPDBL2**2
            A(I,I)=TMPDBL
            A(I,ITOADD)=TMPDBL2
            IF(PIVOT.LT.TMPDBL) THEN
C     FIND THE NEXT LARGEST PIVOT
               PIVOT=TMPDBL
               IPIV=I
            ENDIF
         ENDDO
C
         RANK=ITOADD
      ENDDO
C
      RETURN
C
      END
