
!================================================================
! ++    ecriture
!================================================================
      SUBROUTINE ecriture(f,nbom,nboa,nbcol,tnb,ndim)
      USE mod_system
      IMPLICIT NONE

         integer nbom,nboa,nbcol,ndim
         real(kind=Rkind) f(ndim,ndim)
         logical tnb

         integer i,j,nb,nbblocs,nfin
         character (len=Name_longlen) :: wformat


         nbblocs=int(nbom/nbcol)

         IF (nbblocs*nbcol .EQ. nbom) nbblocs=nbblocs-1

         IF (tnb) THEN

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,nboa
               write(6,51) j,(f(j,i+nb*nbcol),i=1,nbcol)
 51            format(i5,1x,10(f18.10,1x))
!51            format(i5,1x,10(e15.8,1x))
 50          CONTINUE
             write(6,*)
 40        CONTINUE
           DO 60 j=1,nboa
             nfin=nbom-nbcol*nbblocs
             write(6,51) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE
         ELSE

           DO 10 nb=0,nbblocs-1
             DO 20 j=1,nboa
               write(6,21) (f(j,i+nb*nbcol),i=1,nbcol)
 21            format(5x,10(f18.10,1x))
!21            format(5x,10(e15.8,1x))
 20          CONTINUE
             write(6,*)
 10        CONTINUE
           DO 30 j=1,nboa
             nfin=nbom-nbcol*nbblocs
             write(6,21) (f(j,i+nbcol*nbblocs),i=1,nfin)
 30         CONTINUE

         END IF

         RETURN
         end subroutine ecriture
!================================================================
! ++    ecriture complex
!================================================================
      SUBROUTINE ecriture_cplx(f,nbom,nboa,nbcol,tnb,ndim)
      USE mod_system
      IMPLICIT NONE

         integer nbom,nboa,nbcol,ndim
         complex(kind=Rkind) f(ndim,ndim)
         logical tnb

         integer i,j,nb,nbblocs,nfin
         character (len=Name_longlen) :: wformat


         nbblocs=int(nbom/nbcol)

         IF (nbblocs*nbcol .EQ. nbom) nbblocs=nbblocs-1

         IF (tnb) THEN

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,nboa
               write(6,51) j,(f(j,i+nb*nbcol),i=1,nbcol)
 51            format(i3,2x,10('(',f15.7,' +i',f15.7,')'))
!51            format(i3,2x,10(e15.8,3x))
 50          CONTINUE
             write(6,*)
 40        CONTINUE
           DO 60 j=1,nboa
             nfin=nbom-nbcol*nbblocs
             write(6,51) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE
         ELSE

           DO 10 nb=0,nbblocs-1
             DO 20 j=1,nboa
               write(6,21) (f(j,i+nb*nbcol),i=1,nbcol)
 21            format(5x,10('(',f15.7,' +i',f15.7,')'))
!21            format(5x,10(e15.8,3x))
 20          CONTINUE
             write(6,*)
 10        CONTINUE
           DO 30 j=1,nboa
             nfin=nbom-nbcol*nbblocs
             write(6,21) (f(j,i+nbcol*nbblocs),i=1,nfin)
 30         CONTINUE

         END IF

         RETURN
         end subroutine ecriture_cplx
!================================================================
! ++    ecriture lecture
!================================================================
      SUBROUTINE ecriture_f2(f,n,nio,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nio,ndim
         real(kind=Rkind) f(ndim,ndim)
         integer i,j

         DO i=1,n
          write(nio,*) (f(i,j),j=1,n)
         END DO
         end subroutine ecriture_f2
!        --------------------------------------------------------
      SUBROUTINE lecture_f2(f,n,nio,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nio,ndim
         real(kind=Rkind) f(ndim,ndim)
         integer i,j

         DO i=1,n
          read(nio,*) (f(i,j),j=1,n)
         END DO
         end subroutine lecture_f2
!================================================================
! ++    ecriture
!================================================================
      SUBROUTINE ecriture_f(x,f,n,nio,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nio,ndim
         real(kind=Rkind) x,f(ndim,ndim)

         integer i,j

         write(nio,*) x
         DO i=1,n
          write(nio,*) (f(i,j),j=1,n)
         END DO
         end subroutine ecriture_f
!================================================================
! ++    lecture
!================================================================
      SUBROUTINE lecture_f(x,f,n,nio,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nio,ndim
         real(kind=Rkind) x,f(ndim,ndim)

         integer i,j

         read(nio,*) x
         DO i=1,n
          read(nio,*) (f(i,j),j=1,n)
         END DO
         end subroutine lecture_f
!================================================================
! ++    lecture turtle
!================================================================

      SUBROUTINE lecture(v,nbom,nboa,nbcol,tnb,ndim)
      USE mod_system
      IMPLICIT NONE

         integer nbom,nboa,nbcol,ndim
         logical tnb
         real(kind=Rkind) v(ndim,ndim)

         integer i,j,nb,nbblocs,nfin,ifa

         character*132 ligne


!        write(6,*) ' lecture',nbom,nboa,tnb

         nbblocs=int(nbom/nbcol)
         IF (nbblocs*nbcol .EQ. nbom) nbblocs=nbblocs-1


         IF (tnb) THEN
          DO 10 nb=0,nbblocs-1
           DO 20 j=1,nboa
!              read(5,'(132a)') ligne
!              write(6,*) ligne
!              read(ligne,*) ifa,(v(j,i+nb*nbcol),i=1,nbcol)
               read(5,*) ifa,(v(j,i+nb*nbcol),i=1,nbcol)
!              write(6,*) 'a',(v(j,i+nb*nbcol),i=1,nbcol)
 20        CONTINUE
           read(5,*)
 10       CONTINUE
          nfin=nbom-nbcol*nbblocs
            DO 30 j=1,nboa
               read(5,*) ifa,(v(j,i+nbcol*nbblocs),i=1,nfin)
!              write(6,*) 'a',(v(j,i+nbcol*nbblocs),i=1,nfin)
 30          CONTINUE
         ELSE
          DO 40 nb=0,nbblocs-1
           DO 50 j=1,nboa
               read(5,*) (v(j,i+nb*nbcol),i=1,nbcol)
 50        CONTINUE
           read(5,*)
 40       CONTINUE
          nfin=nbom-nbcol*nbblocs
            DO 60 j=1,nboa
               read(5,*) (v(j,i+nbcol*nbblocs),i=1,nfin)
 60          CONTINUE
         ENDIF

         RETURN
         end subroutine lecture

!================================================================
! ++    lecture d'une matrice sous forme triangulaire (hessian)
!================================================================

      SUBROUTINE lecture_tri(v,nbligne,nbcol,nio)
      USE mod_system
      IMPLICIT NONE

         integer nbligne,nbcol
         real(kind=Rkind) v(nbligne,nbligne)

         integer i,j,nb,nbblocs,nfin,ifa
         integer nbcol2,nio


!       write(6,*) ' lecture_tri',nbligne,nbcol,nio

        nbblocs=int(nbligne/nbcol)+1
        IF (nbblocs*nbcol .EQ. nbligne) nbblocs=nbblocs-1

        DO nb=0,nbblocs-1
         DO j=nb*nbcol+1,nbligne
             nbcol2 = min(j-nb*nbcol,nbcol)
             read(nio,11) (v(i+nb*nbcol,j),i=1,nbcol2)
 11          format(21x,5f10.5)
!            write(6,11) (v(i+nb*nbcol,j),i=1,nbcol2)
         END DO
         IF (nb .NE. nbblocs-1) read(nio,*)
!        write(6,*) 'a',nb
        END DO

        DO i=1,nbligne
        DO j=i+1,nbligne
            v(j,i) = v(i,j)
        END DO
        END DO

!        CALL ecriture(v,nbligne,nbligne,5,.true.,nbligne)

         RETURN
         end subroutine lecture_tri
      SUBROUTINE lecture_tri2(v,nbligne,nbcol,nio)
      USE mod_system
      IMPLICIT NONE

         integer nbligne,nbcol
         real(kind=Rkind) v(nbligne,nbligne)

         integer i,j,nb,nbblocs,nfin,ifa
         integer nbcol2,nio


!       write(6,*) ' lecture_tri',nbligne,nbcol,nio

        nbblocs=int(nbligne/nbcol)+1
        IF (nbblocs*nbcol .EQ. nbligne) nbblocs=nbblocs-1

        DO nb=0,nbblocs-1
         DO j=nb*nbcol+1,nbligne
             nbcol2 = min(j-nb*nbcol,nbcol)
             read(nio,11) (v(i+nb*nbcol,j),i=1,nbcol2)
 11          format(4x,5D14.6)
!            write(6,11) (v(i+nb*nbcol,j),i=1,nbcol2)
         END DO
         IF (nb .NE. nbblocs-1) read(nio,*)
!        write(6,*) 'a',nb
        END DO

        DO i=1,nbligne
        DO j=i+1,nbligne
            v(j,i) = v(i,j)
        END DO
        END DO

!        CALL ecriture(v,nbligne,nbligne,5,.true.,nbligne)

         RETURN
         end subroutine lecture_tri2
!================================================================
! ++    ecriture d'une matrice rectangulaire
!================================================================
      SUBROUTINE ecriture_r4(f,nio,nl,nc,nbcol1)
      USE mod_system
      IMPLICIT NONE

         integer nl,nc,nbcol,nbcol1
         real(kind=Rkind) f(nl,nc)

         integer i,j,nb,nbblocs,nfin,nio

         character(len=:), allocatable     :: wformat


         !write(6,*) 'nl,nc,nbcol',nl,nc,nbcol
         nbcol = nbcol1
         IF (nbcol .GT. 10) nbcol=10
         nbblocs=int(nc/nbcol)

         !CALL sub_LineOFmatFormat(wformat,nl,nbcol,.FALSE.)

         wformat = String_TO_String( '(i' //                            &
                      int_TO_char(int(log10(real(nl,kind=Rkind)))+1) // &
                      ',2x,10(f16.10,1x))' )

         IF (nbblocs*nbcol .EQ. nc) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,nl
               write(nio,wformat) j,(f(j,i+nb*nbcol),i=1,nbcol)
 50          CONTINUE
             IF (nl .GT. 1 ) write(nio,*)
 40        CONTINUE
           DO 60 j=1,nl
             nfin=nc-nbcol*nbblocs
             write(nio,wformat) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE

           if (allocated(wformat)) deallocate(wformat)

         RETURN
         end subroutine ecriture_r4
!================================================================
! ++    ecriture d'une matrice rectangulaire
!================================================================
      SUBROUTINE ecriture_c4(f,nio,nl,nc,nbcol1)
      USE mod_system
      IMPLICIT NONE

         integer nl,nc,nbcol1,nbcol
         complex(kind=Rkind) f(nl,nc)

         integer i,j,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


!        write(6,*) 'nl,nc,nbcol',nl,nc,nbcol
         nbcol = nbcol1
         IF (nbcol .GT. 10) nbcol=10
         nbblocs=int(nc/nbcol)


         IF (nbblocs*nbcol .EQ. nc) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,nl
               write(nio,51) j,(f(j,i+nb*nbcol),i=1,nbcol)
 51   format(i3,2x,'( ',f10.4,' +i',f10.4,9(' )( ',f10.4,' i ',f10.4))
 50          CONTINUE
             IF (nl .GT. 1 ) write(nio,*)
 40        CONTINUE
           DO 60 j=1,nl
             nfin=nc-nbcol*nbblocs
             write(nio,51) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE

         RETURN
         end subroutine ecriture_c4
!================================================================
! ++    ecriture d'une matrice rectangulaire
!================================================================
      SUBROUTINE ecriture_intr4(f,nio,nl,nc,nbcol1)
      USE mod_system
      IMPLICIT NONE

         integer nl,nc,nbcol,nbcol1
         integer f(nl,nc)

         integer i,j,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


!        write(6,*) 'nl,nc,nbcol',nl,nc,nbcol
         nbcol = nbcol1
         IF (nbcol .GT. 20) nbcol=20
         nbblocs=int(nc/nbcol)


         IF (nbblocs*nbcol .EQ. nc) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,nl
               write(nio,51) j,(f(j,i+nb*nbcol),i=1,nbcol)
 51            format(i3,2x,20(i5,1x))
 50          CONTINUE
             IF (nl .GT. 1 ) write(nio,*)
 40        CONTINUE
           DO 60 j=1,nl
             nfin=nc-nbcol*nbblocs
             write(nio,51) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE

         RETURN
         end subroutine ecriture_intr4
!================================================================
! ++    lecture d'une matrice rectangulaire
!================================================================
      SUBROUTINE lecture_r4(f,nio,nl,nc,nbcol)
      USE mod_system
      IMPLICIT NONE

         integer nl,nc,nbcol
         real(kind=Rkind) f(nl,nc)

         integer i,j,jj,nb,nbblocs,nfin,nio

!        write(6,*) 'nl,nc,nbcol',nl,nc,nbcol


         nbblocs=int(nc/nbcol)

         IF (nbblocs*nbcol .EQ. nc) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,nl
               read(nio,*) jj,(f(j,i+nb*nbcol),i=1,nbcol)
 50          CONTINUE
             IF (nl .GT. 1 ) read(nio,*)
 40        CONTINUE
           DO 60 j=1,nl
             nfin=nc-nbcol*nbblocs
             read(nio,*) jj,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE

      end subroutine lecture_r4
      SUBROUTINE lecture_c4(f,nio,nl,nc,nbcol)
      USE mod_system
      IMPLICIT NONE

         integer :: nl,nc,nbcol
         complex(kind=Rkind) :: f(nl,nc)

         integer :: i,j,jj,nb,nbblocs,nfin,nio

!        write(6,*) 'nl,nc,nbcol',nl,nc,nbcol


         nbblocs=int(nc/nbcol)

         IF (nbblocs*nbcol .EQ. nc) nbblocs=nbblocs-1

           DO nb=0,nbblocs-1
             DO j=1,nl
               read(nio,*) jj,(f(j,i+nb*nbcol),i=1,nbcol)
             END DO
             IF (nl .GT. 1 ) read(nio,*)
           END DO

           nfin=nc-nbcol*nbblocs
           DO j=1,nl
             read(nio,*) jj,(f(j,i+nbcol*nbblocs),i=1,nfin)
           END DO

      end subroutine lecture_c4
!================================================================
! ++    write a vector in line
!================================================================
      SUBROUTINE ecriture_l3(l,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim,nbcol1
         real(kind=Rkind) l(ndim)

         integer i,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


         nbcol = nbcol1
         IF (nbcol .GT. 5) nbcol=5
         nbblocs=int(n/nbcol)


         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             write(nio,51) (l(i+nb*nbcol),i=1,nbcol)
 51          format(5(e30.23,1x))
 40        CONTINUE
           nfin=n-nbcol*nbblocs
           write(nio,51) (l(i+nbcol*nbblocs),i=1,nfin)

         RETURN
         end subroutine ecriture_l3
!================================================================
! ++    write a vector in line
!================================================================
      SUBROUTINE ecriture_cplxl3(l,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim,nbcol1
         complex(kind=Rkind) l(ndim)

         integer i,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


         nbcol = nbcol1
         IF (nbcol .GT. 5) nbcol=5
         nbblocs=int(n/nbcol)

         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             write(nio,51) (l(i+nb*nbcol),i=1,nbcol)
 51          format(5('(',f15.7,' +i',f15.7,')'))
 40        CONTINUE
           nfin=n-nbcol*nbblocs
           write(nio,51) (l(i+nbcol*nbblocs),i=1,nfin)

         RETURN
         end subroutine ecriture_cplxl3
!================================================================
! ++    write a vector in line
!================================================================
      SUBROUTINE ecriture_intl3(l,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim,nbcol1
         integer l(ndim)

         integer i,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat



         nbcol = nbcol1
         IF (nbcol .GT. 20) nbcol=20
         nbblocs=int(n/nbcol)

         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             write(nio,51) (l(i+nb*nbcol),i=1,nbcol)
 51          format(20(i5))
 40        CONTINUE
           nfin=n-nbcol*nbblocs
           write(nio,51) (l(i+nbcol*nbblocs),i=1,nfin)

         RETURN
         end subroutine ecriture_intl3
!================================================================
! ++    write a vector in line
!================================================================
      SUBROUTINE ecriture_l31(name,l,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim,nbcol1
         real(kind=Rkind) l(ndim)
         character (len=*) :: name

         integer i,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


         nbcol = nbcol1
         IF (nbcol .GT. 5) nbcol=5
         nbblocs=int(n/nbcol)


         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             write(nio,51) name,(l(i+nb*nbcol),i=1,nbcol)
 51          format(a,5(e30.23,1x))
 40        CONTINUE
           nfin=n-nbcol*nbblocs
           write(nio,51) name,(l(i+nbcol*nbblocs),i=1,nfin)

      end subroutine ecriture_l31
      SUBROUTINE ecriture_c31(name,l,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer              :: n,nbcol,ndim,nbcol1
         complex (kind=Rkind) :: l(ndim)
         character (len=*)    :: name

         integer              :: i,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


         nbcol = nbcol1
         IF (nbcol .GT. 5) nbcol=5
         nbblocs=int(n/nbcol)


         IF (nbblocs*nbcol == n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             write(nio,51) trim(name),(l(i+nb*nbcol),i=1,nbcol)
 51          format(a,5(' (',f10.4,1x,f10.4,')'))
 40        CONTINUE
           nfin=n-nbcol*nbblocs
           write(nio,51) trim(name),(l(i+nbcol*nbblocs),i=1,nfin)

         RETURN
         end subroutine ecriture_c31


!================================================================
! ++    read a vector in line
!================================================================
      SUBROUTINE lecture_l3(l,nio,n,nbcol,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim
         real(kind=Rkind) l(ndim)

         integer i,nb,nbblocs,nfin,nio


         nbblocs=int(n/nbcol)


         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             read(nio,51) (l(i+nb*nbcol),i=1,nbcol)
 51          format(5(e30.23,1x))
 40        CONTINUE
           nfin=n-nbcol*nbblocs
           read(nio,51) (l(i+nbcol*nbblocs),i=1,nfin)

         RETURN
         end subroutine lecture_l3
!================================================================
! ++    read a vector in line
!================================================================
      SUBROUTINE lecture_l31(name,l,nio,n,nbcol,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim
         real(kind=Rkind) l(ndim)
         character*10 name

         integer i,nb,nbblocs,nfin,nio


         nbblocs=int(n/nbcol)

         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             read(nio,51) name,(l(i+nb*nbcol),i=1,nbcol)
 51          format(a10,5(e30.23,1x))
 40        CONTINUE
           nfin=n-nbcol*nbblocs
           read(nio,51) name,(l(i+nbcol*nbblocs),i=1,nfin)

         RETURN
         end subroutine lecture_l31
!================================================================
! ++    ecriture
!================================================================
      SUBROUTINE ecriture_f3(f,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim,nbcol1
         real(kind=Rkind) f(ndim,ndim)

         integer i,j,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


         nbcol = nbcol1
         IF (nbcol .GT. 5) nbcol=5
         nbblocs=int(n/nbcol)

         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,n
               write(nio,51) j,(f(j,i+nb*nbcol),i=1,nbcol)
 51            format(i5,2x,5(e30.23,1x))
 50          CONTINUE
             write(nio,*)
 40        CONTINUE
           DO 60 j=1,n
             nfin=n-nbcol*nbblocs
             write(nio,51) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE

         RETURN
         end subroutine ecriture_f3
!================================================================
! ++    ecriture
!================================================================
      SUBROUTINE lecture_f3(f,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim,nbcol1
         real(kind=Rkind) f(ndim,ndim)

         integer i,j,nb,nbblocs,nfin,nio
         integer jdum
         character(len=Name_len) cdum

         nbcol = nbcol1
         IF (nbcol .GT. 5) nbcol=5
         nbblocs=int(n/nbcol)


         nbblocs=int(n/nbcol)

         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,n
               read(nio,*) jdum,(f(j,i+nb*nbcol),i=1,nbcol)
!              write(6,*) jdum,(f(j,i+nb*nbcol),i=1,nbcol)
 50          CONTINUE
             read(nio,*)
 40        CONTINUE
           DO 60 j=1,n
             nfin=n-nbcol*nbblocs
             read(nio,*) jdum,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE

         RETURN
         end subroutine lecture_f3
!================================================================
! ++    ecriture
!================================================================
      SUBROUTINE ecriture_f4(f,nio,n,nbcol1,ndim)
      USE mod_system
      IMPLICIT NONE

         integer n,nbcol,ndim,nbcol1
         real(kind=Rkind) f(ndim,ndim)

         integer i,j,nb,nbblocs,nfin,nio
         character (len=Name_longlen) :: wformat


         nbcol = nbcol1
         IF (nbcol .GT. 5) nbcol=5
         nbblocs=int(n/nbcol)

         IF (nbblocs*nbcol .EQ. n) nbblocs=nbblocs-1

           DO 40 nb=0,nbblocs-1
             DO 50 j=1,n
               write(nio,51) j,(f(j,i+nb*nbcol),i=1,nbcol)
 51            format(i3,2x,10(e15.8,3x))
 50          CONTINUE
             write(nio,*)
 40        CONTINUE
           DO 60 j=1,n
             nfin=n-nbcol*nbblocs
             write(nio,51) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
 60        CONTINUE

         RETURN
         end subroutine ecriture_f4
