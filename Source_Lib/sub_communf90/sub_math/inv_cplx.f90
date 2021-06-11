      USE mod_system
      IMPLICIT NONE

       integer n
       parameter (n=3)
       complex(kind=Rkind) a(n,n),c(n,n),a_save(n,n)
       complex(kind=Rkind) trav(n)
       integer    index(n)


       a = 0.d0

       a(1,1) = dcmplx(1.d0,2.d0)
       a(2,2) = dcmplx(2.d0,2.d0)
       a(3,3) = dcmplx(3.d0,2.d0)

       a(1,2) = dcmplx(1.d0,1.d0)
       a(2,1) = dcmplx(1.d0,1.d0)
       a(2,3) = dcmplx(2.d0,1.d0)
       a(1,3) = dcmplx(0.d0,1.d0)
       a(3,1) = dcmplx(1.d0,0.d0)


       a_save = a

       CALL inversion_cplx(c,a_save,trav,index,n)

       write(out_unitp,*) c
       write(out_unitp,*)

       write(out_unitp,*) matmul(c,a)





       END
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
           c(i,j)=0.D0
         END DO
         c(i,i)=1.D0
       END DO
       CALL ludcmp_cplx(a,n,trav,index,d)

       DO j=1,n
         CALL lubksb_cplx(a,n,index,c(1,j))
       END DO

       RETURN
       end subroutine inversion_cplx
!
!================================================================
!    resolution de a*x=b apres la procedure ludcmp
!
!================================================================

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
         ELSE IF (sum .NE. 0.) THEN
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
       parameter (tiny=1.0d-20)
       complex(kind=Rkind) a(n,n),vv(n)
       complex(kind=Rkind) aamax,sum,dum,d
       integer index(n)

       integer i,j,k,imax

       d=1.0
       DO 12 i=1,n
        aamax=0.0
        DO 11 j=1,n
          IF (abs(a(i,j)) .GT. abs(aamax)) aamax=abs(a(i,j))
 11     CONTINUE
        IF (abs(aamax) .EQ. tiny) STOP "matrice singuliere"
        vv(i)=1/aamax
 12    CONTINUE


       DO 19 j=1,n

        DO 14 i=1,j-1
         sum=a(i,j)
         DO 13 k=1,i-1
          sum=sum-a(i,k)*a(k,j)
 13      CONTINUE
         a(i,j)=sum
 14     CONTINUE

        aamax=0.0
        DO 16 i=j,n
         sum=a(i,j)
         DO 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
 15      CONTINUE
         a(i,j)=sum
         dum=vv(i)*abs(sum)
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
        IF (a(j,j) .EQ. 0.0) a(j,j)=tiny
        IF (j .NE. n) THEN
          dum=1./a(j,j)
          DO 18 i=j+1,n
            a(i,j)=a(i,j)*dum
 18       CONTINUE
        ENDIF

 19    CONTINUE


       RETURN
       end subroutine ludcmp_cplx
