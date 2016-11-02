      SUBROUTINE BLOCKPIVOTCHOL( UPLO, N, A, LDA, BLCKSZ, PIV, RANK, 
     $     DWORK, TOL, INFO )
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
      INTEGER            BLCKSZ,INFO, LDA, N, RANK
      CHARACTER          UPLO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION   DWORK(BLCKSZ, BLCKSZ, 3 )
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
      DOUBLE PRECISION    ONE, ZERO, NEGONE
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0, 
     $     NEGONE=-1.0d+0 )
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION    PIVOT, TMPDBL, MAXPIV, MINPIV
      INTEGER             IPIV, I, J, TMPINT, NCOL, INFOF, INFOI, 
     $     NBSTOP,IBLOCK, JBLOCK, NBLOCK, IBLOFF, JBLOFF
      LOGICAL             UPNLOW, ANYCON
C     LAST LETTER IN NAME OF CHARACTER VARIABLE IS THE VALUE THAT THE
C     VARIABLE HOLDS, E.G. UPLOL='L', DIAGN='N'
      CHARACTER           TRANSN, TRANST, UPLOL, DIAGN, SIDEL      
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, DGEMM, DPOTRF, DPOTRI, DTRSM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX0, DABS
C     ..
C     
C     Test the input parameters.
C     
C
C     THE NUMBER OF BLOCKS IS THE NUMBER OF ROWS DIVIDED BY THE
C     BLOCK SIZE
      NBLOCK = N/BLCKSZ
C
C     RANK BEING NEGATIVE ON ENTRY IS A FLAG THAT MEANS THE USER
C     ONLY WANTS THE FIRST -RANK VALUES (FOR SUBSET SELECTION)
      NBSTOP=-RANK
      IF((NBSTOP.LT.0).OR.(NBSTOP.GT.NBLOCK)) THEN
c         PRINT*, 'RANK SET TO N'
         NBSTOP=NBLOCK
      ENDIF
C      PRINT*, 'YADA0'
C      call flush(8)
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
      ELSEIF(BLCKSZ.LE.0) THEN
         INFO = -8
      ELSEIF(NBLOCK*BLCKSZ.NE.N) THEN
         INFO =-16
      ENDIF
C      PRINT*, 'YADA0.5'
C      call flush(8)
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'LEV3PCHOL', -INFO )
         RETURN
      ENDIF
C     
C     Quick return if possible
C     
      TRANST='T'
      TRANSN='N'
      UPLOL='L'
      DIAGN='N'
      SIDEL='L'
C
C      PRINT*, 'YADA0.6 N=', N
C      call flush(8)
      IF( N.EQ.0 )
     $     RETURN
C      PRINT*, 'YADA0.7'
C      call flush(8)
C      
C     THE FIRST COLUMN OF PIVOTING CHOLESKY IS DIFFERENT THAN ALL 
C     OTHERS SO DO IT SEPARATELY (DOING IT SEPARATELY IS FASTER
C     THAN USING IF'S INSIDE THE LOOP)    
C
C     FILL UP THE PIV (PIVOT) INTEGER ARRAY AND FIND THE FIRST
C     PIVOT
C
C      PRINT*, 'YADA1'
C      call flush(8)

      DO JBLOCK=1,NBLOCK
         PIV(JBLOCK)=JBLOCK
      ENDDO
C         
C      PRINT*, 'YADA2'
C      call flush(8)
      DO JBLOCK=1,NBLOCK
C         PRINT*, 'JBLOCK=', JBLOCK
C         call flush(8)
         ANYCON=.FALSE.
         JBLOFF=(JBLOCK-1)*BLCKSZ
C     THE NEXT DIAGONAL BLOCK WE WANT IS THE ONE WITH THE SMALLEST
C     "PIVOT" (THE PIVOT IS THE ONE NORM OF THE INVERSE OF THE 
C     DIAGONAL BLOCK)
C
C     WE START BY SAYING THE JBLOCK-TH BLOCK IS THE BEST BLOCK AND
C     THEN WE'LL LOOK FOR A BETTER BLOCK
         IPIV=JBLOCK
C
C     FIGURE OUT THE JBLOCK-TH BLOCK'S PIVOT
C     COPY THE DIAGONAL BLOCK TO WORKSPACE (BEST LOCATION)
         DO J=1,BLCKSZ
            DO I=1,BLCKSZ
               DWORK(I,J,1)=A(I+JBLOFF,J+JBLOFF)
            ENDDO
         ENDDO

C      PRINT*, 'YADA3'
C      call flush(8)

         
C
C     DO A CHOLESKY FACTORIZATION OF THE JBLOCK-TH DIAGONAL BLOCK
C     IN THE BEST BLOCK'S WORKSPACE
         CALL DPOTRF(UPLOL,BLCKSZ,DWORK(1,1,1),BLCKSZ,INFOF)

C      PRINT*, 'YADA4'
C      call flush(8)

C
C     MAKE A COPY OF THE CHOLESKY FACTORIZATION
         DO J=1,BLCKSZ
            DO I=J,BLCKSZ
               DWORK(I,J,3)=DWORK(I,J,1)
            ENDDO
         ENDDO
C
C     INVERT THE DIAGONAL BLOCK (STARTING FROM THE COPY OF THE
C     CHOLESKY FACTORIZATION OF THE BEST (SO-FAR) DIAGONAL BLOCK)
         CALL DPOTRI(UPLOL,BLCKSZ,DWORK(1,1,3),BLCKSZ,INFOI)

C      PRINT*, 'YADA5'
C      call flush(8)

C
C     IF WE WERE ABLE TO INVERT ANY DIAGONAL BLOCK THEN WE WILL BE 
C     ABLE TO CONTINUE ON TO THE NEXT JBLOCK LOOP ITERATION
         IF((INFOF.EQ.0).AND.(INFOI.EQ.0))
     $        ANYCON=.TRUE.
C
C     COMPUTE THE ONE NORM OF THE INVERSE OF THE DIAGONAL BLOCK
C     AND STORE IT IN "PIVOT"
         DO J=2,BLCKSZ
            DWORK(1,J,3)=DABS(DWORK(J,J,3));
         ENDDO
         DO J=1,BLCKSZ-1
            DO I=J+1,BLCKSZ
               TMPDBL=DABS(DWORK(I,J,3))
               DWORK(1,J,3)=DWORK(1,J,3)+TMPDBL
               DWORK(1,I,3)=DWORK(1,I,3)+TMPDBL
            ENDDO
         ENDDO
         PIVOT=DWORK(1,1,3)
         DO J=2,BLCKSZ
            IF(PIVOT.LT.DWORK(1,J,3)) THEN
               PIVOT=DWORK(1,J,3)
            ENDIF
         ENDDO

C      PRINT*, 'YADA6 PIVOT=',PIVOT
C      call flush(8)

C
C     IF THE JBLOCK-TH BLOCK IS NOT THE LAST BLOCK
         IF(JBLOCK.LT.NBLOCK) THEN      
C     THEN CHECK IF ANY LATER "CANDIDATE" DIAGONAL BLOCKS HAVE A 
C     SMALLER/BETTER PIVOT VALUE... LOOP OVER THE REST OF THE 
C     DIAGONAL BLOCKS TO FIND THE ONE WITH THE SMALLEST PIVOT
            DO IBLOCK=JBLOCK+1,NBLOCK
               IBLOFF=(IBLOCK-1)*BLCKSZ
C     COPY THE DIAGONAL BLOCK TO "CANDIDATE" BLOCK WORKSPACE
               DO J=1,BLCKSZ
                  DO I=1,BLCKSZ
                     DWORK(I,J,2)=A(I+IBLOFF,J+IBLOFF)
                  ENDDO
               ENDDO
C
C     DO A CHOLESKY FACTORIZATION OF THE "CANDIDATE" DIAGONAL BLOCK IN
C     ITS WORKSPACE
               CALL DPOTRF(UPLOL,BLCKSZ,DWORK(1,1,2),BLCKSZ,INFOF)
C
C     MAKE A COPY OF THE CHOLESKY FACTORIZATION OF THE "CANDIDATE" 
C     DIAGONAL BLOCK
               DO J=1,BLCKSZ
                  DO I=J,BLCKSZ
                     DWORK(I,J,3)=DWORK(I,J,2)
                  ENDDO
               ENDDO
C
C     INVERT THE "CANDIDATE" DIAGONAL BLOCK (STARTING FROM THE COPY OF
C     THE CANDIDATE DIAGONAL BLOCK'S CHOLESKY FACTORIZATION)
               CALL DPOTRI(UPLOL,BLCKSZ,DWORK(1,1,3),BLCKSZ,INFOI)
C
C     COMPUTE THE ONE NORM OF THE INVERSE OF THE "CANDIDATE" DIAGONAL 
C     BLOCK AND STORE THE FINAL ANSWER IN TMPDBL
               DO J=2,BLCKSZ
                  DWORK(1,J,3)=DABS(DWORK(J,J,3));
               ENDDO
               DO J=1,BLCKSZ-1
                  DO I=J+1,BLCKSZ
                     TMPDBL=DABS(DWORK(I,J,3))
                     DWORK(1,J,3)=DWORK(1,J,3)+TMPDBL
                     DWORK(1,I,3)=DWORK(1,I,3)+TMPDBL
                  ENDDO
               ENDDO
               TMPDBL=DWORK(1,1,3)
               DO J=2,BLCKSZ
                  IF(TMPDBL.LT.DWORK(1,J,3)) THEN
                     TMPDBL=DWORK(1,J,3)
                  ENDIF
               ENDDO
C
C     IF THE "CANDIDATE" DIAGONAL BLOCK HAS A SMALLER PIVOT (THE ONE
C     NORM OF IT'S INVERSE) THAN THE PREVIOUS BEST (SMALLEST) THEN 
C     SAY THE "CANDIDATE" IS THE BEST AND COPY ITS CHOLESKY 
C     FACTORIZATION TO THE "BEST" LOCATION SO IT WON'T GET OVERWRITTEN
               IF(((INFOF.EQ.0).AND.(INFOI.EQ.0)).AND.
     $              ((.NOT.ANYCON).OR.(TMPDBL.LT.PIVOT))) THEN
C     IF WE WERE ABLE TO INVERT ANY DIAGONAL BLOCK THEN WE WILL BE 
C     ABLE TO CONTINUE ON TO THE NEXT JBLOCK LOOP ITERATION
                  ANYCON=.TRUE.
                  IPIV=IBLOCK
                  PIVOT=TMPDBL
                  DO J=1,BLCKSZ
                     DWORK(J,J,1)=DWORK(J,J,2)               
                     DO I=J+1,BLCKSZ
                        DWORK(I,J,1)=DWORK(I,J,2)
                        DWORK(J,I,1)=DWORK(I,J,1)
                     ENDDO
                  ENDDO
               ENDIF
C
C     MOVE ON TO THE NEXT CANDIDATE DIAGONAL BLOCK
C     DO IBLOCK=JBLOCK+1,NBLOCK
            ENDDO
C
C     IF(JBLOCK.LT.NBLOCK) THEN    
         ENDIF
C      PRINT*, 'YADA7 PIVOT=',PIVOT
C      call flush(8)


C
C     CAN USE THIS TO ESTIMATE RCOND, BUT NEED TO FIGURE OUT HOW TO
C     USE IT TO ESTIMATE RCOND... TODO FOR NOW IT'S A PLACEHOLDER
         IF(JBLOCK.EQ.1) THEN
            MAXPIV=PIVOT
            MINPIV=PIVOT
         ELSE
            IF(MAXPIV.LT.PIVOT) 
     $           MAXPIV=PIVOT
            IF(PIVOT.LT.MINPIV)
     $           MINPIV=PIVOT
         ENDIF
C         
C     IF THE CURRENT BLOCK IS NOT THE LAST BLOCK THEN WE MAY NEED TO 
C     PIVOT

C      PRINT*, 'YADA9'
C      call flush(8)

         IF((JBLOCK.LT.NBLOCK).AND.(IPIV.NE.JBLOCK)) THEN
C            PRINT*, 'PIVOTING'
C            call flush(8)
C     WE NEED TO PIVOT
            TMPINT=PIV(JBLOCK)
            PIV(JBLOCK)=PIV(IPIV)
            PIV(IPIV)=TMPINT
            IBLOFF=(IPIV-1)*BLCKSZ
C     SWAP BLOCK COLUMNS
            DO J=1,BLCKSZ
               DO I=1,N
                  TMPDBL=A(I,J+JBLOFF)
                  A(I,J+JBLOFF)=A(I,J+IBLOFF)
                  A(I,J+IBLOFF)=TMPDBL
               ENDDO
            ENDDO
C     SWAP BLOCK ROWS
            DO J=1,N
               DO I=1,BLCKSZ
                  TMPDBL=A(I+JBLOFF,J)
                  A(I+JBLOFF,J)=A(I+IBLOFF,J)
                  A(I+IBLOFF,J)=TMPDBL
               ENDDO
            ENDDO
C     IF((JBLOCK.LT.NBLOCK).AND.(IPIV.NE.JBLOCK)) THEN
         ENDIF

C      PRINT*, 'YADA10'
C      call flush(8)

C     
C     COPY THE CHOLESKY FACTORIZATION OF THE BEST BLOCK TO THE
C     JBLOCK-TH DIAGONAL BLOCK
         DO J=1,BLCKSZ
            DO I=1,BLCKSZ
               A(I+JBLOFF,J+JBLOFF)=DWORK(I,J,1)
            ENDDO
         ENDDO

C      PRINT*, 'YADA11'
C      call flush(8)

C     IF WE WERE NOT ABLE TO INVERT ANY DIAGONAL BLOCK THEN WE 
C     CAN'T CONTINUE (WE NEED TO STOP NOW)
         IF(.NOT.ANYCON) THEN
C            PRINT*, 'ANYCON=.FALSE.'
C            call flush(8)            
            IF(INFOF.NE.0) THEN
               RANK=(JBLOCK-1)*BLCKSZ+ABS(INFOF)
            ELSE
               RANK=JBLOCK*BLCKSZ
            ENDIF
            GOTO 10
         ENDIF

C         PRINT*, 'YADA12'
C         call flush(8)
C
C     IF THE CURRENT BLOCK IS NOT THE LAST BLOCK
         IF(JBLOCK.LT.NBLOCK) THEN
C     THEN MODIFY THE ELEMENTS TO THE RIGHT AND BELOW THE DIAGONAL 
C     BLOCK
            IBLOFF=JBLOFF+BLCKSZ            
            NCOL=N-IBLOFF
C     A(1:BLCKSZ,JBLOFF:END)=L\A(1:BLCKSZ,JBLOFF:END);
            CALL DTRSM(SIDEL,UPLOL,TRANSN,DIAGN,BLCKSZ,NCOL,
     $           ONE,DWORK(1,1,1),BLCKSZ,A(JBLOFF+1,IBLOFF+1),LDA)

C            PRINT*, 'YADA13'
C            call flush(8)

C            
C     A(JBLOFF:END,1:BLCKSZ)=A(1:BLCKSZ,JBLOFF:END)';
            DO I=1,BLCKSZ
               DO J=IBLOFF+1,N
                  A(J,I+JBLOFF)=A(I+JBLOFF,J)
               ENDDO
            ENDDO

C      PRINT*, 'YADA14'
C      call flush(8)


C            
C     A(JBLOFF:END,JBLOFF:END)=1.0*A(JBLOFF:END,JBLOFF:END)...
C     -1.0*A(1:BLCKSZ,JBLOFF:END)'*A(1:BLCKSZ,JBLOFF:END);
            CALL DGEMM(TRANST,TRANSN,NCOL,NCOL,BLCKSZ,NEGONE,
     $           A(JBLOFF+1,IBLOFF+1),LDA,A(JBLOFF+1,IBLOFF+1),LDA,
     $           ONE,A(IBLOFF+1,IBLOFF+1),LDA)
C      PRINT*, 'YADA15'
C      call flush(8)


C     IF(JBLOCK.LT.NBLOCK) THEN
         ENDIF
C     DO JBLOCK=1,NBLOCK
      ENDDO
C      PRINT*, 'YADA16'
C      call flush(8)

C
      RANK=N
 10   CONTINUE
C     CONVERT PIV FROM BLOCK INDICES TO EQN/ROW/COLUMN INDICES
      DO JBLOCK=NBLOCK,1,-1
         JBLOFF=(JBLOCK-1)*BLCKSZ
         IBLOFF=(PIV(JBLOCK)-1)*BLCKSZ
         DO J=1,BLCKSZ
            PIV(J+JBLOFF)=J+IBLOFF
         ENDDO
      ENDDO
C      PRINT*, 'YADA17'
C      call flush(8)
C
C      DO J=1,N
C         PRINT*, 'PIV(', J , ')=', PIV(J)
C      ENDDO
C
C      PRINT*, 'YADA18'
C      call flush(8)


C
      RETURN
C
      END
