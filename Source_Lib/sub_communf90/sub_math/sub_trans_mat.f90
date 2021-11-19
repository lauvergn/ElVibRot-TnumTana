!================================================================
!    inversion de la matrice m1 : m2=m1^-1
!================================================================
      SUBROUTINE inv_m1_TO_m2(m1,m2,n,inv_type,epsi)
      USE mod_system
      IMPLICIT NONE

       integer          :: n
       real(kind=Rkind) :: m1(n,n)
       real(kind=Rkind) :: m2(n,n)
       integer          :: inv_type
       real(kind=Rkind) :: epsi

       integer          :: indx(n)
       real(kind=Rkind) :: trav(n),m1w(n,n)
       real(kind=Rkind) :: vv(n,n)
       real(kind=Rkind) :: b(n)

       real(kind=Rkind) :: wmax,wmin
       real(kind=Rkind) :: d
       integer          :: j



       CALL mat_id(m2,n,n)
       m1w = m1

       SELECT CASE (inv_type)
       CASE (0) ! ludcmp ...
         CALL ludcmp(m1w,n,trav,indx,d)
         DO j=1,n
           CALL lubksb(m1w,n,indx,m2(:,j))
         END DO
       CASE (1) ! svd
         CALL SVDCMP(m1w,n,n,trav,vv,n)
         ! Find maximum singular value
         !write(out_unitp,*) 'SVD : epsi',epsi
         !write(out_unitp,*) 'SVD : trav',trav

         wmax = maxval(trav(:))
         wmin = wmax * epsi
         !write(out_unitp,*) 'SVD : count non zero',count(trav >= wmin)
         ! Zero the "small" singular values
         WHERE (trav < WMIN) trav = ZERO

         DO j=1,n
           b(:) = m2(:,j)
           CALL SVBKSB(m1w,trav,vv,n,n,b,m2(:,j),n)
         END DO
       CASE Default ! ludcmp ...
          CALL ludcmp(m1w,n,trav,indx,d)
          DO j=1,n
            CALL lubksb(m1w,n,indx,m2(:,j))
          END DO
       END SELECT

       END SUBROUTINE inv_m1_TO_m2
!================================================================
!    inversion de la matrice m1 : m2=m1^-1
!================================================================
      SUBROUTINE inv_m1_TO_m2_cplx(m1,m2,n,inv_type,epsi)
      USE mod_system
      IMPLICIT NONE

       integer          :: n
       complex(kind=Rkind) :: m1(n,n)
       complex(kind=Rkind) :: m2(n,n)
       integer          :: inv_type
       real(kind=Rkind) :: epsi

       integer          :: indx(n)
       complex(kind=Rkind) :: trav(n),m1w(n,n)
       complex(kind=Rkind) :: vv(n,n)
       complex(kind=Rkind) :: b(n)

       complex(kind=Rkind) :: wmax,wmin
       complex(kind=Rkind) :: d
       integer          :: j



       CALL Cplx_mat_id(m2,n,n)
       m1w = m1

       SELECT CASE (inv_type)
       CASE (0) ! ludcmp ...
         CALL ludcmp_cplx(m1w,n,trav,indx,d)
         DO j=1,n
           CALL lubksb_cplx(m1w,n,indx,m2(:,j))
         END DO

       CASE (1) ! svd

          STOP 'SVD not yet in complex'

       CASE Default ! ludcmp ...
         CALL ludcmp_cplx(m1w,n,trav,indx,d)
         DO j=1,n
           CALL lubksb_cplx(m1w,n,indx,m2(:,j))
         END DO
       END SELECT

       END SUBROUTINE inv_m1_TO_m2_cplx
!================================================================
!    Dertermniant of m1
!================================================================
      SUBROUTINE Det_OF_m1(m1,det,n)
      USE mod_system
      IMPLICIT NONE

       integer          :: n
       real(kind=Rkind) :: m1(n,n)
       real(kind=Rkind) :: det

       integer          :: index(n)
       real(kind=Rkind) :: trav(n),m1w(n,n)

       real(kind=Rkind) :: d
       integer          :: j

       m1w = m1

       CALL ludcmp(m1w,n,trav,index,d)

       det = d
       DO j=1,n
         det = det * m1w(j,j)
       END DO

       END SUBROUTINE Det_OF_m1
!================================================================
!    inversion de la matrice a : c=1/a
!================================================================
      SUBROUTINE inversion(c,a,trav,index,n)
      USE mod_system
      IMPLICIT NONE

       integer n
       real(kind=Rkind) a(n,n),d
       real(kind=Rkind) c(n,n)

       integer index(n)
       real(kind=Rkind) trav(n)

       integer i,j

       DO i=1,n
         DO j=1,n
           c(i,j)=ZERO
         END DO
         c(i,i)=ONE
       END DO


       CALL ludcmp(a,n,trav,index,d)

       DO j=1,n
         CALL lubksb(a,n,index,c(1,j))
       END DO

       RETURN
       end subroutine inversion
!================================================================
!    Solve,x: a.x=b
!================================================================
      SUBROUTINE Linear_Sys(a,b,x,n)
      USE mod_system
      IMPLICIT NONE

       integer          :: n
       real(kind=Rkind) :: a(n,n),vv(n,n)
       real(kind=Rkind) :: b(n),x(n)

       integer          :: indx(n)
       real(kind=Rkind) :: trav(n),d
       real(kind=Rkind) :: aa(n,n)
       real(kind=Rkind) :: wmax,wmin,epsi=ONETENTH**10
       integer          :: k

       logical          :: svd = .TRUE.
       !logical          :: svd = .FALSE.

       x  = b
       aa = a

       IF (svd) THEN
           ! une facon.... SVD
           CALL SVDCMP(aa,n,n,trav,vv,n)
           ! Find maximum singular value
           wmax = maxval(trav(:))
           wmin = wmax * epsi
           ! Zero the "small" singular values
           DO k=1,n
              IF (trav(k) < WMIN) trav(k) = ZERO
           END DO

           CALL SVBKSB(aa,trav,vv,n,n,b,x,n)


           !write(out_unitp,*) 'solve?',sum(abs(matmul(a,x)-b))
           !STOP
        ELSE
          ! une autre ...
          CALL ludcmp(aa,n,trav,indx,d)
          CALL lubksb(aa,n,indx,x)
          !IF (mpro) CALL CALL mprove(a,aa,n,indx,b,x)

        END IF

       END SUBROUTINE Linear_Sys
!================================================================
!    ameliore la solution d un systeme d equations
!    par une iteration
!================================================================
      SUBROUTINE mprove(A,ALUD,N,INDX,B,X)
      USE mod_system
      IMPLICIT NONE

      integer          :: n
      real(kind=Rkind) :: a(n,n)
      real(kind=Rkind) :: alud(n,n)
      real(kind=Rkind) :: b(n),x(n),r(n)
      integer          :: indx(n)
      integer          :: i,j

      DO I=1,N
        R(I) = -B(I) + dot_product(A(I,:),X(:))
      END DO
      CALL LUBKSB(ALUD,N,INDX,R)

      X(:) = X(:) - R(:)

      END SUBROUTINE mprove

!
!================================================================
!    resolution de a*x=b apres la procedure ludcmp
!
!================================================================

      SUBROUTINE lubksb(a,n,index,b)
      USE mod_system
      IMPLICIT NONE

       integer n
       real(kind=Rkind) a(n,n),b(n)
       integer index(n)
       real(kind=Rkind)  sum

       integer i,j,ii,ll

       ii=0
       DO 12 i=1,n
         ll=index(i)
         sum=b(ll)
         b(ll)=b(i)
         IF (II .NE. 0) THEN
            DO 11 j=ii,i-1
              sum=sum-a(i,j)*b(j)
 11         CONTINUE
         ELSE IF (sum .NE. ZERO) THEN
                ii=i
              ENDIF
         b(i)=sum
 12    CONTINUE
       DO 14 i=n,1,-1
         sum=b(i)
         DO 13 j=i+1,n
           sum=sum-a(i,j)*b(j)
 13      CONTINUE
         b(i)=sum/a(i,i)
 14    CONTINUE

       RETURN
       end subroutine lubksb
!================================================================
!    decomposition de a=l*u (pour la resolution d un systeme d equations
!     l matrice triangulaire inferieur
!     u matrice triangulaire superieur
!
!    a l u matrices n*n
!
!================================================================

      SUBROUTINE ludcmp(a,n,vv,index,d)
      USE mod_system
      IMPLICIT NONE

       integer n
       real(kind=Rkind)   tiny
       parameter (tiny=ONETENTH**20)
       real(kind=Rkind) a(n,n),vv(n)
       real(kind=Rkind) aamax,sum,dum,d
       integer index(n)

       integer i,j,k,imax

       d=ONE
       DO 12 i=1,n
        aamax=ZERO
        DO 11 j=1,n
          IF (abs(a(i,j)) .GT. aamax) aamax=abs(a(i,j))
 11     CONTINUE
        IF (aamax < tiny) STOP "matrice singuliere"
        vv(i)=ONE/aamax
 12    CONTINUE


       DO 19 j=1,n

        DO 14 i=1,j-1
         sum=a(i,j)
         DO 13 k=1,i-1
          sum=sum-a(i,k)*a(k,j)
 13      CONTINUE
         a(i,j)=sum
 14     CONTINUE

        aamax=ZERO
        imax=0
        DO 16 i=j,n
         sum=a(i,j)
         DO 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
 15      CONTINUE
         a(i,j)=sum
         dum=vv(i)*abs(sum)
         IF (dum .GE. aamax) THEN
           imax=i
           aamax=dum
         ENDIF
 16     CONTINUE
        IF (imax ==0) THEN
          write(out_unitp,*) ' ERROR in ludcmp'
          write(out_unitp,*) ' imax = 0 !!!'
          write(out_unitp,*) ' matrix a:'
          CALL ecriture(a,n,n,4,.TRUE.,n)
          STOP
        END IF

        IF (j .NE. imax) THEN
          DO 17 k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
 17       CONTINUE
          d=-d
          vv(imax)=vv(j)
        ENDIF

        index(j)=imax
        IF (a(j,j) .EQ. ZERO) a(j,j)=tiny
        IF (j .NE. n) THEN
          dum=ONE/a(j,j)
          DO 18 i=j+1,n
            a(i,j)=a(i,j)*dum
 18       CONTINUE
        ENDIF

 19    CONTINUE


       RETURN
       END SUBROUTINE ludcmp

      SUBROUTINE SVDCMP(A,M,N,W,V,max_n)
      USE mod_system
      IMPLICIT NONE

      integer max_n,N,M
      real (kind=Rkind) :: A(max_n,max_n),V(max_n,max_n)
      real (kind=Rkind) :: W(max_n),RV1(max_n)
      real (kind=Rkind) :: G,SCALE,ANORM,S,F,H,C,Y,Z,X

      integer I,K,J,NM,JJ,L,ITS


      G=ZERO
      SCALE=ZERO
      ANORM=ZERO
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.ZERO) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(sqrt(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=ZERO
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.ZERO) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(sqrt(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=ZERO
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.ZERO) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=ZERO
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=ZERO
            V(J,I)=ZERO
31        CONTINUE
        ENDIF
        V(I,I)=ONE
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=ZERO
33        CONTINUE
        ENDIF
        IF (G.NE.ZERO) THEN
          G=ONE/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=ZERO
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=ZERO
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+ONE
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=ZERO
          S=ONE
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=sqrt(F*F+G*G)
              W(I)=H
              H=ONE/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.ZERO) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.50) STOP 'No convergence in 50 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
          G=sqrt(F*F+ONE)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=ONE
          S=ONE
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=sqrt(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 JJ=1,N
              X=V(JJ,J)
              Z=V(JJ,I)
              V(JJ,J)= (X*C)+(Z*S)
              V(JJ,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=sqrt(F*F+H*H)
            W(J)=Z
            IF (Z.NE.ZERO) THEN
              Z=ONE/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 JJ=1,M
              Y=A(JJ,J)
              Z=A(JJ,I)
              A(JJ,J)= (Y*C)+(Z*S)
              A(JJ,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=ZERO
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END SUBROUTINE SVDCMP


      SUBROUTINE SVBKSB(U,W,V,M,N,B,X,max_n)
      USE mod_system
      IMPLICIT NONE

      integer max_n,M,N
      real (kind=Rkind) :: U(max_n,max_n),V(max_n,max_n)
      real (kind=Rkind) :: W(max_n),B(max_n),X(max_n),TMP(max_n)
      real (kind=Rkind) :: s
      integer I,J,JJ

      DO 12 J=1,N
        S=ZERO
        IF(W(J).NE.ZERO)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=ZERO
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN

      END SUBROUTINE SVBKSB


!================================================================
!    inversion de la matrice a : c=1/a
!================================================================

      SUBROUTINE inversion_cplx(c,a,trav,index,n)
      USE mod_system
      IMPLICIT NONE

       integer n
       complex(kind=Rkind) a(n,n),d
       complex(kind=Rkind) c(n,n)

       integer index(n)
       complex(kind=Rkind) trav(n)

       integer i,j

       DO i=1,n
         DO j=1,n
           c(i,j)=CZERO
         END DO
         c(i,i)=CONE
       END DO
       CALL ludcmp_cplx(a,n,trav,index,d)

       DO j=1,n
         CALL lubksb_cplx(a,n,index,c(1,j))
       END DO

       RETURN
       end subroutine inversion_cplx
!================================================================
!    resolution de a*x=b apres la procedure ludcmp
!================================================================
  SUBROUTINE Driver_LU_solve_cplx(a,n,LU_index,b,type_lu)
  USE mod_system
  IMPLICIT NONE

  integer,             intent(in)    :: n,type_lu
  complex(kind=Rkind), intent(inout) :: a(n,n),b(n)
  integer,             intent(in)    :: LU_index(n)

  integer               :: err,type_lu_loc
  integer, parameter    :: type_lu_default = 1
  integer(kind=I4kind)  :: n4,ierr4



    type_lu_loc = type_lu

    !when lapack is used and Rkind /= real64 (not a double)
    IF (Rkind /= R8kind .AND. type_lu_loc == 3) type_lu_loc = type_lu_default

#if __LAPACK != 1
    IF ( type_lu_loc == 3) type_lu_loc = type_lu_default
#endif

    SELECT CASE (type_lu)
    CASE(1) ! ori
      CALL lubksb_cplx(a,n,LU_index,b)
    CASE(3) ! lapack
#if __LAPACK == 1
      n4     = int(n,kind=I4kind)
      CALL ZGETRS('No transpose',n4,1,a,n4,LU_index,b,n4,ierr4)
      err = int(ierr4)
      IF (err /= 0) STOP 'LU Driver_LU_solve_cplx'
#else
      write(out_unitp,*) ' ERROR in Driver_LU_solve_cplx'
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makefile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in Driver_LU_solve_cplx: LAPACK case impossible'
#endif
    CASE Default
      CALL lubksb_cplx(a,n,LU_index,b)
    END SELECT

  END SUBROUTINE Driver_LU_solve_cplx
  SUBROUTINE Driver_LU_decomp_cplx(a,n,LU_index,d,type_lu)
  USE mod_system
  IMPLICIT NONE

  integer,             intent(in)    :: n,type_lu
  complex(kind=Rkind), intent(inout) :: d,a(n,n)
  integer,             intent(in)    :: LU_index(n)

  integer               :: err,type_lu_loc
  integer, parameter    :: type_lu_default = 1
  integer(kind=I4kind)  :: n4,ierr4
  complex(kind=Rkind), allocatable :: work(:)



    type_lu_loc = type_lu

    !when lapack is used and Rkind /= real64 (not a double)
    IF (Rkind /= R8kind .AND. type_lu_loc == 3) type_lu_loc = type_lu_default

#if __LAPACK != 1
    IF ( type_lu_loc == 3) type_lu_loc = type_lu_default
#endif

    SELECT CASE (type_lu)
    CASE(1) ! ori
      allocate(work(n))
      CALL ludcmp_cplx(a,n,work,LU_index,d)
      deallocate(work)
    CASE(3) ! lapack
#if __LAPACK == 1
      n4     = int(n,kind=I4kind)
      CALL ZGETRF(n4,n4,a,n4,LU_index,ierr4)
      err = int(ierr4)
      IF (err /= 0) STOP 'Driver_LU_decomp_cplx'
#else
      write(out_unitp,*) ' ERROR in Driver_LU_decomp_cplx'
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makefile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in Driver_LU_decomp_cplx: LAPACK case impossible'
#endif
    CASE Default
      allocate(work(n))
      CALL ludcmp_cplx(a,n,work,LU_index,d)
      deallocate(work)
    END SELECT

  END SUBROUTINE Driver_LU_decomp_cplx
      SUBROUTINE lubksb_cplx(a,n,index,b)
      USE mod_system
      IMPLICIT NONE

       integer n
       complex(kind=Rkind) a(n,n),b(n)
       integer index(n)
       complex(kind=Rkind)  sum

       integer i,j,ii,ll

       ii=0
       DO 12 i=1,n
         ll=index(i)
         sum=b(ll)
         b(ll)=b(i)
         IF (II .NE. 0) THEN
            DO 11 j=ii,i-1
              sum=sum-a(i,j)*b(j)
 11         CONTINUE
         ELSE IF (abs(sum) .NE. ZERO) THEN
                ii=i
              ENDIF
         b(i)=sum
 12    CONTINUE
       DO 14 i=n,1,-1
         sum=b(i)
         DO 13 j=i+1,n
           sum=sum-a(i,j)*b(j)
 13      CONTINUE
         b(i)=sum/a(i,i)
 14    CONTINUE

       RETURN
       end subroutine lubksb_cplx
!================================================================
!    decomposition de a=l*u (pour la resolution d un systeme d equations
!     l matrice triangulaire inferieur
!     u matrice triangulaire superieur
!
!    a l u matrices n*n
!
!================================================================

      SUBROUTINE ludcmp_cplx(a,n,vv,index,d)
      USE mod_system
      IMPLICIT NONE

       integer n
       real(kind=Rkind)   tiny
       parameter (tiny=ONETENTH**20)
       complex(kind=Rkind) a(n,n),vv(n)
       complex(kind=Rkind) aamax,sum,dum,d
       integer index(n)

       integer i,j,k,imax

       d=CONE
       DO 12 i=1,n
        aamax=CZERO
        DO 11 j=1,n
          IF (abs(a(i,j)) .GT. abs(aamax)) aamax=cmplx(abs(a(i,j)),kind=Rkind)
 11     CONTINUE
        IF (abs(aamax) < tiny) STOP "matrice singuliere"
        vv(i)=CONE/aamax
 12    CONTINUE


       DO 19 j=1,n

        DO 14 i=1,j-1
         sum=a(i,j)
         DO 13 k=1,i-1
          sum=sum-a(i,k)*a(k,j)
 13      CONTINUE
         a(i,j)=sum
 14     CONTINUE

        aamax=CZERO
        imax = 0
        DO 16 i=j,n
         sum=a(i,j)
         DO 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
 15      CONTINUE
         a(i,j)=sum
         dum=vv(i)*cmplx(abs(sum),kind=Rkind)
         IF (abs(dum) .GE. abs(aamax)) THEN
           imax=i
           aamax=dum
         ENDIF
 16     CONTINUE

        IF (j .NE. imax) THEN
          DO 17 k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
 17       CONTINUE
          d=-d
          vv(imax)=vv(j)
        ENDIF

        index(j)=imax
        IF (abs(a(j,j)) .EQ. ZERO) a(j,j)=cmplx(tiny,kind=Rkind)
        IF (j .NE. n) THEN
          dum=CONE/a(j,j)
          DO 18 i=j+1,n
            a(i,j)=a(i,j)*dum
 18       CONTINUE
        ENDIF

 19    CONTINUE


       RETURN
       end subroutine ludcmp_cplx
!
!=====================================================================
!
! ++   A = coef*A
!      A : square matrix
!
!=====================================================================
!
      SUBROUTINE mult_mat(a,coef,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer i,j,n,max_niv
       real(kind=Rkind) a(max_niv,max_niv)
       real(kind=Rkind) coef

       DO i=1,n
         DO j=1,n
           a(i,j) = coef*a(i,j)
         END DO
       END DO

       END SUBROUTINE mult_mat
!
!=====================================================================
!
! ++   A = 0
!      A : square matrix
!
!=====================================================================
!
      SUBROUTINE mat0(a,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer i,j,n,max_niv
       real(kind=Rkind) a(max_niv,max_niv)

       a(1:n,1:n) = ZERO

       END SUBROUTINE mat0
!
!=====================================================================
!
! ++   V(n+n) (filled with n real) => V(n) (filled with n cplx)
!
!=====================================================================
!
      SUBROUTINE vectRealTOCplx(v,n)
      USE mod_system
      IMPLICIT NONE

       integer i,j,n
       real(kind=Rkind) v(n+n)

!      DO i=1,2*n
!      write(out_unitp,*) i,v(i)
!      END DO
       DO i=0,n-1
         v(n+n-i-i) = ZERO
         v(n+n-i-i-1) = v(n-i)
       END DO
!      DO i=1,n
!      write(out_unitp,*) i,v(i+i-1),v(i+i)
!      END DO

       END SUBROUTINE vectRealTOCplx
!
!=====================================================================
!
! ++   A = 0
!      A : complex square matrix
!
!=====================================================================
!
      SUBROUTINE cplx_mat0(a,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer             :: n,max_niv
       complex(kind=Rkind) :: a(max_niv,max_niv)

       a(1:n,1:n) = CZERO

       END SUBROUTINE cplx_mat0
!
!=====================================================================
!
! ++   A(i,j) = 0 if |A(i,j)| <epsi
!      A : square matrix
!
!=====================================================================
!
      SUBROUTINE mat_epsiTOzero(a,n,epsi,max_niv)
      USE mod_system
      IMPLICIT NONE

      integer          :: max_niv,n
      real(kind=Rkind) :: a(max_niv,max_niv)
      real(kind=Rkind) :: epsi

      integer   i,j

!---------------------------------------------------------------------
      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING mat_epsiTOzero'
      write(out_unitp,*) 'epsi',epsi
      write(out_unitp,*) 'a',n
      CALL ecriture(a,n,n,5,.TRUE.,n)
      END IF
!---------------------------------------------------------------------


      DO i=1,n
      DO j=1,n

        IF ( abs(a(i,j)) .LT. epsi) a(i,j) = ZERO

      END DO
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'new a'
      CALL ecriture(a,n,n,5,.TRUE.,n)
      write(out_unitp,*) 'END mat_epsiTOzero'
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE mat_epsiTOzero
!
!=====================================================================
!
! ++   A = Id =>  A(i,i)=ONE
!      A : square matrix
!
!=====================================================================
!
      SUBROUTINE mat_id(a,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: i,n,max_niv
       real(kind=Rkind) :: a(max_niv,max_niv)

       a(1:n,1:n) = ZERO

       DO i=1,n
         a(i,i) = ONE
       END DO

       END SUBROUTINE mat_id
      SUBROUTINE Cplx_mat_id(a,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: i,n,max_niv
       complex(kind=Rkind) :: a(max_niv,max_niv)

       a(1:n,1:n) = CZERO

       DO i=1,n
         a(i,i) = CONE
       END DO

       END SUBROUTINE Cplx_mat_id
!
!=====================================================================
!
! ++   copy of 2 matrix A=B
!      square matrix
!
!=====================================================================
!
      SUBROUTINE copy_mat(a,b,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: n,max_niv
       real(kind=Rkind) :: a(max_niv,max_niv)
       real(kind=Rkind) :: b(max_niv,max_niv)

       a(1:n,1:n) = b(1:n,1:n)

       END SUBROUTINE copy_mat
!
!=====================================================================
!
! ++   add  of 2 matrix A = A + B
!      square matrix
!
!=====================================================================
!
      SUBROUTINE mat_p(a,b,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: n,max_niv
       real(kind=Rkind) :: a(max_niv,max_niv)
       real(kind=Rkind) :: b(max_niv,max_niv)

       a(1:n,1:n) = a(1:n,1:n) + b(1:n,1:n)

       END SUBROUTINE mat_p
!
!=====================================================================
!
! ++   add  of 2 matrix A = A - B
!      square matrix
!
!=====================================================================
!
      SUBROUTINE mat_m(a,b,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: n,max_niv
       real(kind=Rkind) :: a(max_niv,max_niv)
       real(kind=Rkind) :: b(max_niv,max_niv)

       a(1:n,1:n) = a(1:n,1:n) - b(1:n,1:n)

       END SUBROUTINE mat_m
!=====================================================================
!
! ++   add  a cte : A = A + cte
!      square matrix
!
!=====================================================================
!
      SUBROUTINE mat_p_cte(a,cte,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: n,max_niv
       real(kind=Rkind) :: a(max_niv,max_niv)
       real(kind=Rkind) :: cte

       a(1:n,1:n) = a(1:n,1:n) + cte

       END SUBROUTINE mat_p_cte
!=====================================================================
!
! ++   add  a cte : A = A + I*cte
!      square matrix
!
!=====================================================================
!
      SUBROUTINE mat_p_Icte(a,cte,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: i,n,max_niv
       real(kind=Rkind) :: a(max_niv,max_niv)
       real(kind=Rkind) :: cte

       DO i=1,n
         a(i,i) = a(i,i) + cte
       END DO

       END SUBROUTINE mat_p_Icte
!
!=====================================================================
!
! ++   add  of 2 matrix A = A + coef*B
!      square matrix
!
!=====================================================================
!
      SUBROUTINE add_mat(a,coef,b,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer          :: n,max_niv
       real(kind=Rkind) :: coef
       real(kind=Rkind) :: a(max_niv,max_niv)
       real(kind=Rkind) :: b(max_niv,max_niv)

       a(1:n,1:n) = a(1:n,1:n) + coef*b(1:n,1:n)

       END SUBROUTINE add_mat
!
!=====================================================================
!
! ++   change v such that v(i) = x
!
!=====================================================================
!
      SUBROUTINE const_d0v(v,x,i,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       real(kind=Rkind)  x,v(n)

       v(i) = x

!      write(out_unitp,*) 'const_d0v : n i v(i)',n,i,v(i)

       end subroutine const_d0v
!
!=====================================================================
!
! ++   change v such that d1v(i,.) = d1x(.)
!
!=====================================================================
!
      SUBROUTINE const_d1v(d1v,d1x,i,n,nd)
      USE mod_system
      IMPLICIT NONE

       integer i,k,n,nd
       real(kind=Rkind)  d1x(nd),d1v(n,nd)

       DO k=1,nd
        d1v(i,k) = d1x(k)
       END DO

!      write(out_unitp,*) 'const_d1v : i d1v(i.)',i,(d1v(i,k),k=1,nd)

       end subroutine const_d1v
!
!=====================================================================
!
! ++   change v such that d2v(i,..) = d1x(..)
!
!=====================================================================
!
      SUBROUTINE const_d2v(d2v,d2x,i,n,nd)
      USE mod_system
      IMPLICIT NONE

       integer i,k,l,n,nd
       real(kind=Rkind)  d2x(nd,nd),d2v(n,nd,nd)

       DO k=1,nd
       DO l=1,nd
        d2v(i,k,l) = d2x(k,k)
       END DO
       END DO

!      write(out_unitp,*) 'const_d1v : i d1v(i.)',i,(d1v(i,k),k=1,nd)

       end subroutine const_d2v
!
!=====================================================================
!
! ++   vector norm
!
!=====================================================================
!
      FUNCTION norme_v(v,n)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: norme_v

       integer i,n
       real(kind=Rkind) v(n)
       real(kind=Rkind) val

       val = ZERO
       DO i=1,n
         val = val + v(i)*v(i)
       END DO

       norme_v = sqrt(val)

       RETURN
       end function norme_v
!=====================================================================
!
! ++   copy of 2 vectors v1=v2
!
!=====================================================================
!
      SUBROUTINE copy_vect(v1,v2,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       real(kind=Rkind) v1(n)
       real(kind=Rkind) v2(n)

       DO i=1,n
         v1(i)= v2(i)
       END DO
!      write(out_unitp,*) v1

       RETURN
       end subroutine copy_vect
!
!=====================================================================
!
! ++   copy of 2 vectors v1=v2 (complex)
!
!=====================================================================
!
      SUBROUTINE copy_vect_cplx(v1,v2,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       complex(kind=Rkind) v1(n),v2(n)

       DO i=1,n
         v1(i)= v2(i)
       END DO
!      write(out_unitp,*) v1

       RETURN
       end subroutine copy_vect_cplx
!
!=====================================================================
!
! ++   s = v1 . v2  (dot product)
!
!=====================================================================
!
      FUNCTION dot_pro(v1,v2,n)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: dot_pro

       integer i,n
       real(kind=Rkind) v1(n)
       real(kind=Rkind) v2(n)
       real(kind=Rkind) s

       s = ZERO
       DO i=1,n
         S = S + v1(i) * v2(i)
       END DO

       dot_pro = S

       RETURN
       end function dot_pro
!
!=====================================================================
!
! ++   (i_s,r_s) = <r_v1|(i_v2,r_v2)> (complex dot product)
!
!=====================================================================
!
      FUNCTION cplx_dot_pro(v1,i_v2,r_v2,n)
      USE mod_system
      IMPLICIT NONE
      complex(kind=Rkind) :: cplx_dot_pro

       integer i,n
       real(kind=Rkind) v1(n)
       real(kind=Rkind) i_v2(n)
       real(kind=Rkind) r_v2(n)
       complex(kind=Rkind) s

       s = CZERO
       DO i=1,n
         S = S + cmplx( v1(i)*i_v2(i), v1(i)*r_v2(i),kind=Rkind )
       END DO

       cplx_dot_pro = S

       RETURN
       end function cplx_dot_pro
!
!=====================================================================
!
! ++   v1 = v1 + v2
!
!=====================================================================
!
      SUBROUTINE vect_p(v1,v2,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       real(kind=Rkind) v1(n)
       real(kind=Rkind) v2(n)

       DO i=1,n
         v1(i)= v1(i) + v2(i)
       END DO

       RETURN
       end subroutine vect_p
!=====================================================================
!
! ++   v(i) = v(i)*v(i)
!
!=====================================================================
!
      SUBROUTINE vect_sqr(v,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       real(kind=Rkind) v(n)

       DO i=1,n
         v(i)= v(i) * v(i)
       END DO

       RETURN
       end subroutine vect_sqr
!
!=====================================================================
!
! ++   v1 = v1 - v2
!
!=====================================================================
!
      SUBROUTINE vect_m(v1,v2,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       real(kind=Rkind) v1(n)
       real(kind=Rkind) v2(n)

       DO i=1,n
         v1(i)= v1(i) - v2(i)
       END DO

       RETURN
       end subroutine vect_m
!
!=====================================================================
!
! ++   copy of 2 vectors v1 = v1 + coef*v2
!
!=====================================================================
!
      SUBROUTINE add_vect(v1,coef,v2,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       real(kind=Rkind) coef
       real(kind=Rkind) v1(n)
       real(kind=Rkind) v2(n)

       DO i=1,n
         v1(i)= v1(i) + coef*v2(i)
       END DO

       RETURN
       end subroutine add_vect
!
!=====================================================================
!
! ++   v1 = coef*v1
!      vector : v1
!
!=====================================================================
!
      SUBROUTINE mult_vect(v1,coef,n)
      USE mod_system
      IMPLICIT NONE

       integer i,n
       real(kind=Rkind) coef
       real(kind=Rkind) v1(n)

       DO i=1,n
         v1(i)= coef*v1(i)
       END DO

       RETURN
       end subroutine mult_vect
!
!============================================================
!
!   multiplication de matrice C = A*B
!   sur des matrice rectangulaire
!
!============================================================
!
      SUBROUTINE matmult_r(c,nlc,ncc,a,nla,nca,b,nlb,ncb,ndim)
      USE mod_system
      IMPLICIT NONE

        integer nlc,ncc,nla,nca,nlb,ncb,ndim
        real(kind=Rkind) a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)

        integer i,j,k

        IF (nca .NE. nlb) THEN
          write(out_unitp,*) ' ERREUR : stop dans matmult'
          write(out_unitp,*) 'le nombre de colonnes de A (',nca,') doit'
          write(out_unitp,*) 'etre egal au nombre de lignes de B (',nlb,')'
          STOP
        END IF

        nlc = nla
        ncc = ncb
        DO i=1,nlc
          DO j=1,ncc
            c(i,j) = ZERO
            DO k=1,nca
              c(i,j) = c(i,j) + a(i,k)*b(k,j)
            END DO
          END DO
        END DO

        RETURN
        end subroutine matmult_r
!
!============================================================
!
!   multiplication de matrice C = A*B
!   sur des matrice rectangulaire
!
!   matrices complexes
!
!============================================================
!
      SUBROUTINE matmult_cr(c,nlc,ncc,a,nla,nca,b,nlb,ncb)
      USE mod_system
      IMPLICIT NONE

        integer      nla,nca,   nlb,ncb,   nlc,ncc
        complex(kind=Rkind) a(nla,nca),b(nlb,ncb),c(nlc,ncc)

        integer i,j,k

        IF (nca .NE. nlb .OR. ncc .NE. ncb .OR. nla .NE. nlc) THEN
          write(out_unitp,*) ' ERROR : stop in matmult_cr'
          write(out_unitp,*) 'nla,nca',nla,nca
          write(out_unitp,*) 'nlb,ncb',nlb,ncb
          write(out_unitp,*) 'nlc,ncc',nlc,ncc
          write(out_unitp,*) 'You MUST have: nca=nlb ncb=ncc nla=nlc'
          STOP
        END IF

        DO i=1,nlc
          DO j=1,ncc
            c(i,j) = cmplx(ZERO,ZERO,kind=Rkind)
            DO k=1,nca
              c(i,j) = c(i,j) + a(i,k)*b(k,j)
            END DO
          END DO
        END DO

        RETURN
        end subroutine matmult_cr
!
!============================================================
!
!   multiplication de matrice C = A*B
!   sur des matrice rectangulaire
!
!   matrices complexes
!
!============================================================
!
      SUBROUTINE fmatmult_cr(c,nlc,ncc,a,nla,nca,b,nlb,ncb,             &
                              ind_a,dim_a)
      USE mod_system
      IMPLICIT NONE

      integer      nla,nca,   nlb,ncb,   nlc,ncc
      complex(kind=Rkind) a(nla,nca),b(nlb,ncb),c(nlc,ncc)
      integer      ind_a(nla,nca),dim_a(nla)

      integer i,j,k,ki

      IF (nca .NE. nlb .OR. ncc .NE. ncb .OR. nla .NE. nlc) THEN
          write(out_unitp,*) ' ERROR : stop in matmult_cr'
          write(out_unitp,*) 'nla,nca',nla,nca
          write(out_unitp,*) 'nlb,ncb',nlb,ncb
          write(out_unitp,*) 'nlc,ncc',nlc,ncc
          write(out_unitp,*) 'You MUST have: nca=nlb ncb=ncc nla=nlc'
          STOP
      END IF


      DO i=1,nlc
        DO j=1,ncc
          c(i,j) = cmplx(ZERO,ZERO,kind=Rkind)
          DO k=1,dim_a(i)
            ki = ind_a(k,i)
            c(i,j) = c(i,j) + a(i,ki)*b(ki,j)
          END DO
        END DO
      END DO

      RETURN
      end subroutine fmatmult_cr
!
!============================================================
!
!   multiplication de matrice C = A*B
!   pour des matrices carrees
!
!============================================================
!
      SUBROUTINE matmult(c,a,b,n,ndim)
      USE mod_system
      IMPLICIT NONE

        integer n,ndim
        real(kind=Rkind) a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)
        integer i,j,k

        DO i=1,n
          DO j=1,n
            c(i,j) = ZERO
            DO k=1,n
              c(i,j) = c(i,j) + a(i,k)*b(k,j)
            END DO
          END DO
        END DO

        RETURN
        end subroutine matmult
!
!============================================================
!
!   multiplication de matrice C = tA*B
!   pour des matrices carrees
!
!============================================================
!
      SUBROUTINE matmult2(c,a,b,n,ndim)
      USE mod_system
      IMPLICIT NONE

        integer n,ndim
        real(kind=Rkind) a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)
        integer i,j,k

        DO i=1,n
          DO j=1,n
            c(i,j) = ZERO
            DO k=1,n
              c(i,j) = c(i,j) + a(k,i)*b(k,j)
            END DO
          END DO
        END DO

        RETURN
        end subroutine matmult2
!
!=====================================================================
!
! ++   calcul de mm=tc*m*c pour i,j
!      vecteur i :   c(.,i)
!      vecteur j :   c(.,j)
!
!=====================================================================
!
      FUNCTION tcmc_ij(i,j,m,c,n,max_niv)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) tcmc_ij

        integer i,j,k,l,n,max_niv
        real(kind=Rkind) m(max_niv,max_niv)
        real(kind=Rkind) val
        real(kind=Rkind) c(max_niv,max_niv)

        val = ZERO
        DO k=1,n
          DO l=1,n
            val = val + c(k,i)*c(l,j)*m(k,l)
          END DO
        END DO
        tcmc_ij = val

       end function tcmc_ij
!
!=====================================================================
!
!  ++  calcul de mm=tc*m*c
!      vecteur i :   c(.,i)
!      vecteur j :   c(.,j)
!
!=====================================================================
!
      SUBROUTINE tcmc(mm,m,c,mtrav,n,max_niv)
      USE mod_system
      IMPLICIT NONE

       integer n,max_niv
       real(kind=Rkind) m(max_niv,max_niv)
       real(kind=Rkind) mm(max_niv,max_niv)
       real(kind=Rkind) c(max_niv,max_niv)
       real(kind=Rkind) mtrav(max_niv,max_niv)


       CALL copy_mat(mtrav,m,n,max_niv)

       CALL matmult(mm,mtrav,c,n,max_niv)
       CALL matmult2(mtrav,c,mm,n,max_niv)

       CALL copy_mat(mm,mtrav,n,max_niv)

       RETURN
       end subroutine tcmc

      SUBROUTINE tcmc_old(mm,m,c,n,max_niv)
      USE mod_system
      IMPLICIT NONE

        integer i,j,k,l,n,max_niv
        real(kind=Rkind) m(max_niv,max_niv)
        real(kind=Rkind) mm(max_niv,max_niv)
        real(kind=Rkind) c(max_niv,max_niv)

!      write(out_unitp,*) 'tcmc_old m c'
!      CALL ecriture(m,n,n,5,.TRUE.,max_niv)
!      CALL ecriture(c,n,n,5,.TRUE.,max_niv)
       DO i=1,n
         DO j=1,n
           mm(i,j) = ZERO
           DO k=1,n
             DO l=1,n
               mm(i,j) = mm(i,j) + c(k,i)*c(l,j)*m(k,l)
             END DO
           END DO
         END DO
       END DO
!      write(out_unitp,*) 'tcmc_old mm'
!      CALL ecriture(mm,n,n,5,.TRUE.,max_niv)

       RETURN
       end subroutine tcmc_old
!
!
!=====================================================================
!
! ++   recouvrement de deux vecteurs i j
!      vecteur i :   c(.,i)
!
!=====================================================================
!
      FUNCTION recouvrement_ij(i,j,s,c,n,max_niv)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) recouvrement_ij

        integer i,j,k,l,n,max_niv
        real(kind=Rkind) val
        real(kind=Rkind) c(max_niv,max_niv)
        real(kind=Rkind) s(max_niv,max_niv)

        val = ZERO
        DO k=1,n
          DO l=1,n
!           val = val + c(k,i)*c(l,j)*s(k,l)
            val = val + c(k,i)*c(l,j)
          END DO
        END DO
        recouvrement_ij = val

       end function recouvrement_ij
!
!=====================================================================
!
! ++   calcul de mm=c1*m*c2 pour i,j
!      vecteur i :   c1(.,i)
!      vecteur j :   c2(.,j)
!
!=====================================================================
!
      FUNCTION tcmc2_ij(i,j,m,c1,c2,n,max_niv)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) tcmc2_ij

        integer i,j,k,l,n,max_niv
        real(kind=Rkind) m(max_niv,max_niv)
        real(kind=Rkind) val
        real(kind=Rkind) c1(max_niv,max_niv)
        real(kind=Rkind) c2(max_niv,max_niv)

        val = ZERO
        DO k=1,n
          DO l=1,n
            val = val + c1(k,i)*c2(l,j)*m(k,l)
          END DO
        END DO
        tcmc2_ij = val

       end function tcmc2_ij
