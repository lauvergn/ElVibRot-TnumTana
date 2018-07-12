!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong, Josep Maria Luis
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!===========================================================================
!===========================================================================
      MODULE mod_poly
      USE mod_system
      IMPLICIT NONE

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        TYPE para_poly
          logical       :: cplx = .FALSE.
          integer       :: npoly = -1            ! degree of the polynomial
          integer       :: ndim  = 0            ! dimension of the polynomial
          integer       :: nb_coef = 0          ! number of coef.
          real (kind=Rkind),    pointer :: Rcoef(:) => null() ! Rcoef(nb_coef) : coef. of the polynomial
          complex (kind=Rkind), pointer :: Ccoef(:) => null() ! Ccoef(nb_coef) : coef. of the polynomial
          integer      , pointer :: ind(:,:) => null() ! ind(ndim,nb_coef) : index of ndim-poly

          ! in 1D (ndim=1), nb_coef = npoly+1 and ind(1,i)=i-1 (exponent of q^i)
          ! in ndimD => coef(i) * Prod_k=(1...ndim) [q_k^ind(k,i)]
          logical :: alloc_poly = .FALSE.
        END TYPE para_poly
        CONTAINS

!       ==============================================================
!         alloc_poly   : allocation of poly
!         dealloc_poly : deallocation of poly
!         build_ind_poly : build the list of exponent (RECURCIVE)
!         FUNCTION locate(poly,ind)
!
!         write_poly
!
!         p1TOp2 :  p2=p1 or p2=conjugate(p1) if conjug=t
!         p1PLUSp2TOp3
!         p1TIMEp2TOp3
!         d1p1TOp2
!         p1_linearTransfoTOp2
!       ==============================================================

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE check_alloc_poly(A,name_A,name_sub)
        TYPE (para_poly), intent(in)  :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF (.NOT. A%alloc_poly) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been allocated with "alloc_poly"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
        END SUBROUTINE check_alloc_poly
!----- ---------------------------------------------------------------
!----- ---------------------------------------------------------------
!----- ---------------------------------------------------------------
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE alloc_poly(poly,npoly,ndim,cplx)
           TYPE (para_poly), intent(inout)     :: poly
           integer, intent(in)                 :: npoly
           integer, intent(in), optional       :: ndim
           logical, intent(in), optional       :: cplx

           integer :: i,ip

           integer, pointer :: ind_exp(:)

!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'alloc_poly'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             CALL write_poly(poly)
           END IF


           IF (present(ndim) ) THEN
             poly%ndim = ndim
           ELSE
             poly%ndim = 1
           END IF

           IF (present(cplx) ) THEN
             poly%cplx = cplx
           ELSE
             poly%cplx = .FALSE.
           END IF

           IF (npoly < 0) THEN
             write(out_unitp,*) ' ERROR in alloc_poly'
             write(out_unitp,*) ' the degree of the polynomial is <0',npoly
             STOP
           END IF

           poly%npoly   = npoly

           allocate(ind_exp(0:poly%ndim))

           IF (poly%ndim == 1 ) THEN
             poly%nb_coef = npoly + 1
             allocate(poly%ind(poly%ndim,poly%nb_coef))
             DO i=0,npoly
               poly%ind(1,i+1)   = i
             END DO
           ELSE
             poly%nb_coef = 0
             DO i=0,npoly
               ind_exp(:) = 0
               ind_exp(0) = i
               CALL build_ind_poly(poly,ind_exp,1,poly%nb_coef,.TRUE.)
             END DO
!            write(out_unitp,*) 'nb_coef',poly%nb_coef
             allocate(poly%ind(poly%ndim,poly%nb_coef))
             ip = 0
             DO i=0,npoly
               ind_exp(:) = 0
               ind_exp(0) = i
               CALL build_ind_poly(poly,ind_exp,1,ip,.FALSE.)
             END DO
             deallocate(ind_exp)
           END IF

           IF (poly%cplx) THEN
             allocate(poly%Ccoef(poly%nb_coef))
             poly%Ccoef(:)    = ZERO
           ELSE
             allocate(poly%Rcoef(poly%nb_coef))
             poly%Rcoef(:)    = ZERO
           END IF


           poly%alloc_poly = .TRUE.
           IF (debug) THEN
             CALL write_poly(poly)
             write(out_unitp,*) ' END ',name_sub
           END IF

        END SUBROUTINE alloc_poly

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE dealloc_poly(poly)
           TYPE (para_poly), intent(inout)     :: poly

!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'dealloc_poly'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             CALL write_poly(poly)
           END IF


           IF (poly%alloc_poly) THEN
             IF (poly%cplx) THEN
                deallocate(poly%Ccoef)
             ELSE
                deallocate(poly%Rcoef)
             END IF
             deallocate(poly%ind)
             poly%alloc_poly = .FALSE.
           END IF

           IF (debug) THEN
             CALL write_poly(poly)
             write(out_unitp,*) ' END ',name_sub
           END IF

        END SUBROUTINE dealloc_poly
!      ==============================================================
!         build_ind_poly : make the table ind
!         recursive SUBROUTINE
!
!         ind_exp(0) = should be initiated to the degree of the polynomial
!
!         calc_nb_coef_poly: calculation of the number of coef of an nD-polynomial
!         recursive SUBROUTINE
!
!         ind_exp(0) = should be initiated to the degree of the polynomial
!      ==============================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
       RECURSIVE SUBROUTINE build_ind_poly(poly,ind_exp,iq,ip,check)

           TYPE (para_poly), intent(inout)  :: poly
           integer, intent(in)              :: iq
           logical, intent(in)              :: check
           integer, intent(inout)           :: ip
           integer, intent(inout)           :: ind_exp(0:poly%ndim)


!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'build_ind_poly'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             write(out_unitp,*) ' iq,ip,check',iq,ip,check
             write(out_unitp,*) ' ind_exp',ind_exp
           END IF


           IF (.NOT. check .AND. .NOT. associated(poly%ind)) THEN
             write(out_unitp,*) ' ERROR in build_ind_poly'
             write(out_unitp,*) ' poly%ind has NOT been allocated'
             write(out_unitp,*) 'check,associated(poly%ind)',                   &
                         check,associated(poly%ind)
             STOP
           END IF
           IF (iq <0 .OR. iq > poly%ndim) THEN
             write(out_unitp,*) ' ERROR in build_ind_poly'
             write(out_unitp,*) ' iq MUST >= 0 or < ndim',iq,poly%ndim
             write(out_unitp,*) ' Probably, iq MUST be iniatiated to 0'
             STOP
           END IF


!          write(out_unitp,*) 'build_ind_poly: ind_exp',iq,ind_exp
           IF (iq == poly%ndim) THEN
             ip = ip + 1
             ind_exp(iq) = ind_exp(0) - SUM(ind_exp(1:iq-1))
!            write(out_unitp,*) 'ip',ip,'ind_exp',ind_exp(1:poly%ndim)
             IF (.NOT. check) poly%ind(:,ip) = ind_exp(1:poly%ndim)
           ELSE
             ind_exp(iq)=-1
             DO WHILE (ind_exp(iq) < (ind_exp(0)-SUM(ind_exp(1:iq-1))) )
               ind_exp(iq) = ind_exp(iq) + 1
               CALL build_ind_poly(poly,ind_exp,iq+1,ip,check)
             END DO
           END IF

           IF (debug) THEN
             write(out_unitp,*) ' END ',name_sub
           END IF

       END SUBROUTINE build_ind_poly
!      ==============================================================
!         locate function
!      ==============================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
       FUNCTION locate(poly,ind)
           integer                          :: locate
           TYPE (para_poly), intent(in)     :: poly
           integer,intent(in)               :: ind(poly%ndim)

           integer :: i

!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'locate'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             write(out_unitp,*) ' ind',ind
             CALL write_poly(poly)
           END IF

           CALL check_alloc_poly(poly,'poly',name_sub)

           locate = -1


           DO i=1,poly%nb_coef

             IF ( sum(abs(ind(:)-poly%ind(:,i))) == 0 ) THEN
                locate = i
                EXIT
             END IF

           END DO


           IF (debug) THEN
             write(out_unitp,*) 'ind,i',ind,i
             write(out_unitp,*) ' END ',name_sub
           END IF


        END FUNCTION locate
!       ==============================================================
!         write_poly : write a polynomial
!       ==============================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE write_poly(poly)
           TYPE (para_poly), intent(in) :: poly

           integer :: i

!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'write_poly'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
           END IF

           CALL check_alloc_poly(poly,'poly',name_sub)

           write(out_unitp,*) '--------------------------------------'
           write(out_unitp,*) 'ndim,npoly,nb_coef',                     &
                                       poly%ndim,poly%npoly,poly%nb_coef
           write(out_unitp,*) 'alloc_poly',poly%alloc_poly
           write(out_unitp,*) 'cplx',poly%cplx
           write(out_unitp,*)

           IF (poly%cplx) THEN
             DO i=1,poly%nb_coef
               write(out_unitp,*) 'i,ind,Ccoef',i,poly%ind(:,i),poly%Ccoef(i)
             END DO
           ELSE
             DO i=1,poly%nb_coef
               write(out_unitp,*) 'i,ind,Rcoef',i,poly%ind(:,i),poly%Rcoef(i)
             END DO
           END IF

           write(out_unitp,*) '--------------------------------------'

           IF (debug) THEN
             write(out_unitp,*) ' END ',name_sub
           END IF


        END SUBROUTINE write_poly
!       ==============================================================
!         p1TOp2 : p2 = p1
!       ==============================================================

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE p1TOp2(p1,p2,conjug)
           TYPE (para_poly), intent(in)     :: p1
           TYPE (para_poly), intent(inout)    :: p2
           logical, optional, intent(in)    :: conjug

           logical   :: cconjug
!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'p1TOp2'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------

           CALL check_alloc_poly(p1,'p1',name_sub)

           cconjug = .FALSE.
           IF (present(conjug)) cconjug = conjug

           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             write(out_unitp,*) ' conjug',cconjug
             CALL write_poly(p1)
           END IF


           IF (p2%alloc_poly) THEN
             IF (p2%ndim /= p1%ndim) THEN
               write(out_unitp,*) ' ERROR in p1TOp2'
               write(out_unitp,*) ' p2 is allocated so'
               write(out_unitp,*) ' ndim of p2 and p1 MUST be identical'
               write(out_unitp,*) ' ndim of p1,p2:',p1%ndim,p2%ndim
               STOP
             END IF

             IF (p2%npoly < p1%npoly) THEN
               write(out_unitp,*) ' ERROR in p1TOp2'
               write(out_unitp,*) ' p2 is allocated so'
               write(out_unitp,*) ' p2%npoly MUST be >= p1%npoly'
               write(out_unitp,*) ' npoly p1,p2:',p1%npoly,p2%npoly
               STOP
             END IF
           ELSE
             CALL dealloc_poly(p2)
             CALL alloc_poly(p2,ndim=p1%ndim,npoly=p1%npoly,            &
                             cplx=p1%cplx)
           END IF

           IF (p2%cplx) THEN
             p2%Ccoef(:)  = ZERO
             IF (p1%cplx) THEN
               IF (cconjug) THEN
                 p2%Ccoef(1:p1%nb_coef)  = conjg(p1%Ccoef(:))
               ELSE
                 p2%Ccoef(1:p1%nb_coef)  = p1%Ccoef(:)
               END IF
             ELSE
               p2%Ccoef(1:p1%nb_coef)  = p1%Rcoef(:)
             END IF
           ELSE
             p2%Rcoef(:)  = ZERO
             IF (p1%cplx) THEN
               p2%Rcoef(1:p1%nb_coef)  = real(p1%Ccoef(:),kind=Rkind)
             ELSE
               p2%Rcoef(1:p1%nb_coef)  = p1%Rcoef(:)
             END IF
           END IF

           IF (debug) THEN
             CALL write_poly(p1)
             write(out_unitp,*) ' END ',name_sub
           END IF

        END SUBROUTINE p1TOp2

!       ==============================================================
!         p1PLUSp2TOp3 : p3 = p1 + p2
!       ==============================================================

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE p1PLUSp2TOp3(p1,p2,p3)
           TYPE (para_poly), intent(in)     :: p1,p2
           TYPE (para_poly), intent(inout)    :: p3

!          working variables
           integer :: npoly3 ! parameters of p3
           logical :: cplx

!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'p1PLUSp2TOp3'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           CALL check_alloc_poly(p1,'p1',name_sub)

           CALL check_alloc_poly(p2,'p2',name_sub)


           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             CALL write_poly(p1)
             CALL write_poly(p2)
           END IF

           IF (p1%ndim /= p2%ndim) THEN
             write(out_unitp,*) ' ERROR in p1PLUSp2TOp3'
             write(out_unitp,*) ' ndim of p1 and p2 MUST be identical'
             write(out_unitp,*) ' ndim of p1 and p2:',p1%ndim,p2%ndim
             STOP
           END IF


           IF (p3%alloc_poly) THEN
             IF (p3%ndim /= p1%ndim) THEN
               write(out_unitp,*) ' ERROR in p1PLUSp2TOp3'
               write(out_unitp,*) ' p3 is allocated so'
               write(out_unitp,*) ' ndim of p3 and p1 (or p2) MUST be identical'
               write(out_unitp,*) ' ndim of p1,p2 p3:',p1%ndim,p2%ndim,p3%ndim
               STOP
             END IF

             npoly3 = max(p1%npoly,p2%npoly)
             IF (p3%npoly < npoly3) THEN
               write(out_unitp,*) ' ERROR in p1PLUSp2TOp3'
               write(out_unitp,*) ' p3 is allocated so'
               write(out_unitp,*) ' p3%npoly MUST be >= p1%npoly and p2%npoly'
               write(out_unitp,*) ' npoly p1,p2 p3:',p1%npoly,p2%npoly,p3%npoly
               STOP
             END IF
           ELSE
             cplx = p1%cplx .OR. p2%cplx
             CALL dealloc_poly(p3)
             npoly3 = max(p1%npoly,p2%npoly)
             CALL alloc_poly(p3,ndim=p1%ndim,npoly=npoly3,cplx=cplx)
           END IF


           IF (p3%cplx) THEN
             p3%Ccoef(:)  = CZERO
             IF (p1%cplx) THEN
               p3%Ccoef(1:p1%nb_coef)  = p1%Ccoef(:)
             ELSE
               p3%Ccoef(1:p1%nb_coef)  = cmplx(p1%Rcoef(:),kind=Rkind)
             END IF
             IF (p2%cplx) THEN
               p3%Ccoef(1:p2%nb_coef)  = p3%Ccoef(1:p2%nb_coef) +       &
                                         p2%Ccoef(:)
             ELSE
               p3%Ccoef(1:p2%nb_coef)  = p3%Ccoef(1:p2%nb_coef) +       &
                                         cmplx(p2%Rcoef(:),kind=Rkind)
             END IF
           ELSE
             p3%Rcoef(:)  = ZERO
             IF (p1%cplx) THEN
               p3%Rcoef(1:p1%nb_coef)  = real(p1%Ccoef(:),kind=Rkind)
             ELSE
               p3%Rcoef(1:p1%nb_coef)  = p1%Rcoef(:)
             END IF
             IF (p2%cplx) THEN
               p3%Rcoef(1:p2%nb_coef)  = p3%Rcoef(1:p2%nb_coef) +       &
                                         real(p2%Ccoef(:),kind=Rkind)
             ELSE
               p3%Rcoef(1:p2%nb_coef)  = p3%Rcoef(1:p2%nb_coef) +       &
                                         p2%Rcoef(:)
             END IF
           END IF

           IF (debug) THEN
             CALL write_poly(p3)
             write(out_unitp,*) ' END ',name_sub
           END IF

        END SUBROUTINE p1PLUSp2TOp3
!       ==============================================================
!         p1TIMEp2TOp3 : p3 = p1 * p2
!       ==============================================================

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE p1TIMEp2TOp3(p1,p2,p3)
           TYPE (para_poly), intent(in)     :: p1,p2
           TYPE (para_poly), intent(inout)    :: p3

!          working variables
           integer :: npoly3 ! parameters of p3
           integer :: ind3(p1%ndim)
           integer :: i1,i2,i3
           logical :: cplx

!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'p1TIMEp2TOp3'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           CALL check_alloc_poly(p1,'p1',name_sub)

           CALL check_alloc_poly(p2,'p2',name_sub)


           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             CALL write_poly(p1)
             CALL write_poly(p2)
           END IF

           IF (p1%ndim /= p2%ndim) THEN
             write(out_unitp,*) ' ERROR in p1TIMEp2TOp3'
             write(out_unitp,*) ' ndim of p1 and p2 MUST be identical'
             write(out_unitp,*) ' ndim of p1 and p2:',p1%ndim,p2%ndim
             STOP
           END IF


           IF (p3%alloc_poly) THEN
             IF (p3%ndim /= p1%ndim) THEN
               write(out_unitp,*) ' ERROR in p1TIMEp2TOp3'
               write(out_unitp,*) ' p3 is allocated so'
               write(out_unitp,*) ' ndim of p3 and p1 (or p2) MUST be identical'
               write(out_unitp,*) ' ndim of p1,p2 p3:',p1%ndim,p2%ndim,p3%ndim
               STOP
             END IF

             npoly3 = p1%npoly + p2%npoly
             IF (p3%npoly < npoly3) THEN
               write(out_unitp,*) ' ERROR in p1TIMEp2TOp3'
               write(out_unitp,*) ' p3 is allocated so'
               write(out_unitp,*) ' p3%npoly MUST be >= p1%npoly + p2%npoly'
               write(out_unitp,*) ' npoly p1,p2 p3:',p1%npoly,p2%npoly,p3%npoly
               STOP
             END IF
           ELSE
             cplx = p1%cplx .OR. p2%cplx
             CALL dealloc_poly(p3)
             npoly3 = p1%npoly+p2%npoly
             CALL alloc_poly(p3,ndim=p1%ndim,npoly=npoly3,cplx=cplx)
           END IF



           IF (p3%cplx) THEN
             p3%Ccoef(:)  = CZERO
             IF (p1%cplx .AND. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)
!                  write(out_unitp,*) 'i3', i3
                 p3%Ccoef(i3) = p3%Ccoef(i3) + p1%Ccoef(i1)*p2%Ccoef(i2)

               END DO
               END DO


             ELSE IF (.NOT. p1%cplx .AND. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Ccoef(i3) = p3%Ccoef(i3) + cmplx(p1%Rcoef(i1),kind=Rkind)*p2%Ccoef(i2)
               END DO
               END DO
             ELSE IF (p1%cplx .AND. .NOT. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Ccoef(i3) = p3%Ccoef(i3) + p1%Ccoef(i1)*cmplx(p2%Rcoef(i2),kind=Rkind)
               END DO
               END DO
             ELSE
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Ccoef(i3) = p3%Ccoef(i3) + p1%Rcoef(i1)*p2%Rcoef(i2)
               END DO
               END DO
             END IF
           ELSE
             p3%Rcoef(:)  = ZERO

             IF (p1%cplx .AND. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Rcoef(i3) = p3%Rcoef(i3) + real(p1%Ccoef(i1)*p2%Ccoef(i2),kind=Rkind)
               END DO
               END DO
             ELSE IF (.NOT. p1%cplx .AND. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Rcoef(i3) = p3%Rcoef(i3) + p1%Rcoef(i1)*real(p2%Ccoef(i2),kind=Rkind)
               END DO
               END DO
             ELSE IF (p1%cplx .AND. .NOT. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Rcoef(i3) = p3%Rcoef(i3) + real(p1%Ccoef(i1),kind=Rkind)*p2%Rcoef(i2)
               END DO
               END DO
             ELSE
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Rcoef(i3) = p3%Rcoef(i3) + p1%Rcoef(i1)*p2%Rcoef(i2)
               END DO
               END DO
             END IF
           END IF

           IF (debug) THEN
             CALL write_poly(p3)
             write(out_unitp,*) ' END ',name_sub
           END IF


        END SUBROUTINE p1TIMEp2TOp3

!       ==============================================================
!         p2 = d1p1 (first derivative with respect to the variable id)
!       ==============================================================
      !!@description: p2 = d1p1 (first derivative with respect to the variable id)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE d1p1TOp2(p1,p2,id)

           TYPE (para_poly), intent(in)     :: p1
           TYPE (para_poly), intent(inout)    :: p2
           integer,intent(in)               :: id

!          - working variables ----------------
           integer :: npoly2 ! parameters of p3
           logical :: cplx
           integer :: i1,i2,exp_q(p1%ndim)
!          - working variables ----------------

!          - for debuging ----------------------------------
           character (len=*), parameter :: name_sub = 'd1p1TOp2'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------

           CALL check_alloc_poly(p1,'p1',name_sub)


           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub
             CALL write_poly(p1)
           END IF

           IF (p2%alloc_poly) THEN
             IF (p2%ndim /= p1%ndim) THEN
               write(out_unitp,*) ' ERROR in d1p1TOp2'
               write(out_unitp,*) ' p2 is allocated so'
               write(out_unitp,*) ' ndim of p2 and p1 MUST be identical'
               write(out_unitp,*) ' ndim of p1,p2:',p1%ndim,p2%ndim
               STOP
             END IF

             npoly2 = max(0,p1%npoly-1)
             IF (p2%npoly < npoly2) THEN
               write(out_unitp,*) ' ERROR in d1p1TOp2'
               write(out_unitp,*) ' p2 is allocated so'
               write(out_unitp,*) ' p2%npoly MUST be >= p1%npoly-1'
               write(out_unitp,*) ' npoly p1,p2:',p1%npoly,p2%npoly
               STOP
             END IF
             IF (p1%cplx .AND. .NOT. p2%cplx) THEN
               write(out_unitp,*) ' ERROR in d1p1TOp2'
               write(out_unitp,*) ' p1 is complex and not p2 !'
               write(out_unitp,*) ' p1%cplx, p2%cplx',p1%cplx,p2%cplx
               STOP
             END IF
           ELSE
             cplx = p1%cplx
             CALL dealloc_poly(p2)
             npoly2 = max(0,p1%npoly-1)
             CALL alloc_poly(p2,ndim=p1%ndim,npoly=npoly2,cplx=cplx)
           END IF



           IF (p2%cplx .AND. p1%cplx) THEN

             p2%Ccoef(:) = ZERO

             DO i1=2,p1%nb_coef
               IF (p1%ind(id,i1) == 0) CYCLE
               exp_q(:) = p1%ind(:,i1)
               exp_q(id) = exp_q(id) - 1

               i2 = locate(p2,exp_q)
               p2%Ccoef(i2) = p1%Ccoef(i1) *                            &
                   cmplx( p1%ind(id,i1),kind=Rkind)

             END DO

           ELSE IF (p2%cplx .AND. .NOT. p1%cplx) THEN

             p2%Ccoef(:) = ZERO

             DO i1=2,p1%nb_coef
               IF (p1%ind(id,i1) == 0) CYCLE
               exp_q(:) = p1%ind(:,i1)
               exp_q(id) = exp_q(id) - 1

               i2 = locate(p2,exp_q)
               p2%Ccoef(i2) = cmplx(p1%Rcoef(i1),kind=Rkind) *          &
                     cmplx( p1%ind(id,i1),kind=Rkind)

             END DO

           ELSE ! p1 and p2 are real

             p2%Rcoef(:) = ZERO

             DO i1=2,p1%nb_coef
               IF (p1%ind(id,i1) == 0) CYCLE
               exp_q(:) = p1%ind(:,i1)
               exp_q(id) = exp_q(id) - 1

               i2 = locate(p2,exp_q)
               p2%Rcoef(i2) = p1%Rcoef(i1) *                            &
                     real( p1%ind(id,i1),kind=Rkind)

             END DO

           END IF

           IF (debug) THEN
             CALL write_poly(p2)
             write(out_unitp,*) ' END ',name_sub
           END IF

        END SUBROUTINE d1p1TOp2
!       ==============================================================
!       p1_O_p2TOp3   : p3 = p1(+linear transfo) Rq linearTransfo : aX+b
!       ivar : X associated with the index ivar
!       ==============================================================

      !!@description: p1_O_p2TOp3   : p3 = p1(+linear transfo) Rq linearTransfo:
      !!              aX+b ivar : X associated with the index ivar
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE p1_linearTransfoTOp2(p1,ivar,Ra,Ca,Rb,Cb,p2,cplx)
           TYPE (para_poly), intent(in)     :: p1
           real (kind=Rkind), intent(in)        :: Ra,Rb
           complex (kind=Rkind), intent(in)     :: Ca,Cb
           TYPE (para_poly), intent(inout)    :: p2
           logical, intent(in) :: cplx
           integer, intent(in) :: ivar


!          working variables
           integer :: npoly3 ! parameters of p2
           integer, pointer :: ind_exp(:)
           complex (kind=Rkind)    :: CCa,CCb
           integer :: exp_ivar,i1,i3,k

           real (kind=Rkind) :: binomial ! function

!          - for debuging ----------------------------------
           character (len=*), parameter ::                              &
                     name_sub = 'p1_linearTransfoTOp2'
           logical,parameter :: debug = .FALSE.
!          logical,parameter :: debug = .TRUE.
!          -------------------------------------------------
           CALL check_alloc_poly(p1,'p1',name_sub)

           IF (debug) THEN
             write(out_unitp,*) ' BEGINNING ',name_sub

             write(out_unitp,*) ' ivar,cplx',ivar,cplx
             IF (cplx) THEN
               write(out_unitp,*) ' Ca,Cb',Ca,Cb
             ELSE
               write(out_unitp,*) ' Ra,Rb',Ra,Rb
             END IF
             CALL write_poly(p1)
           END IF

           IF (p2%alloc_poly) THEN
             IF (p2%ndim /= p1%ndim) THEN
               write(out_unitp,*) ' ERROR in p1_linearTransfoTOp2'
               write(out_unitp,*) ' p2 is allocated so'
               write(out_unitp,*) ' ndim of p2 and p1 MUST be identical'
               write(out_unitp,*) ' ndim of p1 p2:',p1%ndim,p2%ndim
               STOP
             END IF

             npoly3 = p1%npoly
             IF (p2%npoly < npoly3) THEN
               write(out_unitp,*) ' ERROR in p1_linearTransfoTOp2'
               write(out_unitp,*) ' p2 is allocated so'
               write(out_unitp,*) ' p2%npoly MUST be >= p1%npoly'
               write(out_unitp,*) ' npoly p1 p2:',p1%npoly,p2%npoly
               STOP
             END IF
           ELSE
             CALL dealloc_poly(p2)
             npoly3 = p1%npoly
             CALL alloc_poly(p2,ndim=p1%ndim,npoly=npoly3,              &
                                 cplx=(p1%cplx .OR. cplx))
           END IF

           IF (p1%cplx .OR. cplx) THEN
             IF (.NOT. p2%cplx) THEN
               write(out_unitp,*) ' ERROR in p1_linearTransfoTOp2'
               write(out_unitp,*) ' p2 MUST be complex'
               write(out_unitp,*) ' p2%cplx',p2%cplx
               STOP
             END IF
           END IF


           IF (ivar < 0 .OR. ivar > p1%ndim) THEN
             write(out_unitp,*) ' ERROR in p1_linearTransfoTOp2'
             write(out_unitp,*) ' ivar MUST be > 0 or < ndim'
             write(out_unitp,*) 'ivar,ndim',ivar,p1%ndim
             STOP
           END IF

           allocate(ind_exp(0:p1%ndim))

           IF (cplx .AND. p2%cplx) THEN
            CCa = Ca
            CCb = Cb
           ELSE IF (.NOT. cplx .AND. p2%cplx) THEN
            CCa = cmplx(Ra,kind=Rkind)
            CCb = cmplx(Rb,kind=Rkind)
           END IF

           IF (p2%cplx) THEN
             p2%Ccoef(:)  = CZERO
             IF (p1%cplx) THEN
               DO i1=1,p1%nb_coef
                 ind_exp(0) = 0
                 ind_exp(1:p1%ndim) = p1%ind(:,i1)
                 exp_ivar = p1%ind(ivar,i1)

!                write(out_unitp,*) 'i1,exp_ivar,p1%ind',
!    *                       i1,exp_ivar,p1%ind(:,i1)
!                devlopment of coef*(aX+b)^exp_ivar
                 DO k=0,exp_ivar
                   ind_exp(ivar) = k
                   ind_exp(0) = sum(ind_exp(1:p1%ndim))
                   i3 = locate(p2,ind_exp(1:p1%ndim))
!                  write(out_unitp,*) 'k,ind_exp,i3,p2%ind',
!    *                k,ind_exp(1:p1%ndim),i3,p2%ind(:,i3)
                   p2%Ccoef(i3) = p2%Ccoef(i3) +                        &
                     p1%Ccoef(i1)*cmplx(binomial(exp_ivar,k),kind=Rkind) *&
                             CCa**k * CCb**(exp_ivar-k)
                 END DO
               END DO
             ELSE
               DO i1=1,p1%nb_coef
                 ind_exp(0) = 0
                 ind_exp(1:p1%ndim) = p1%ind(:,i1)
                 exp_ivar = p1%ind(ivar,i1)
!                devlopment of coef*(aX+b)^exp_ivar
                 DO k=0,exp_ivar
                   ind_exp(ivar) = k
                   ind_exp(0) = sum(ind_exp(1:p1%ndim))
                   i3 = locate(p2,ind_exp(1:p1%ndim))
                   p2%Ccoef(i3) = p2%Ccoef(i3) +                        &
                       p1%Rcoef(i1)*cmplx(binomial(exp_ivar,k),kind=Rkind)*&
                                CCa**k * CCb**(exp_ivar-k)
                 END DO
               END DO
             END IF
           ELSE
             p2%Rcoef(:)  = ZERO
             DO i1=1,p1%nb_coef
               ind_exp(0) = 0
               ind_exp(1:p1%ndim) = p1%ind(:,i1)
               exp_ivar = p1%ind(ivar,i1)
!              devlopment of coef*(aX+b)^exp_ivar
               DO k=0,exp_ivar
                 ind_exp(ivar) = k
                 ind_exp(0) = sum(ind_exp(1:p1%ndim))
                 i3 = locate(p2,ind_exp(1:p1%ndim))
                 p2%Rcoef(i3) = p2%Rcoef(i3) +                          &
              p1%Rcoef(i1)*binomial(exp_ivar,k)*Ra**k*Rb**(exp_ivar-k)
                 END DO
               END DO
           END IF

           deallocate(ind_exp)

           IF (debug) THEN
             CALL write_poly(p2)
             write(out_unitp,*) ' END ',name_sub
           END IF
        END SUBROUTINE p1_linearTransfoTOp2

      END MODULE mod_poly

