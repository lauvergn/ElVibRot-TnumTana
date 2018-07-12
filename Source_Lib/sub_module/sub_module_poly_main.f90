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

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      IMPLICIT NONE
        TYPE para_poly
          logical       :: cplx
          integer       :: npoly            ! degree of the polynomial
          integer       :: ndim             ! dimension of the polynomial
          integer       :: nb_coef          ! number of coef.
          real (kind=Rkind),    pointer :: Rcoef(:) ! Rcoef(nb_coef) : coef. of the polynomial
          complex (kind=Rkind), pointer :: Ccoef(:) ! Ccoef(nb_coef) : coef. of the polynomial
          integer      , pointer :: ind(:,:)! ind(ndim,nb_coef) : index of ndim-poly

          ! in 1D (ndim=1), nb_coef = npoly+1 and ind(1,i)=i-1 (exponent of q^i)
          ! in ndimD => coef(i) * Prod_k=(1...ndim) [q_k^ind(k,i)]
          logical :: init0
          logical :: alloc_poly
        END TYPE para_poly
        CONTAINS

!       ==============================================================
!         init0_poly   : initalization of poly
!         alloc_poly   : allocation of poly
!         dealloc_poly : deallocation of poly
!       ==============================================================

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE init0_poly(poly)
          TYPE (para_poly), intent(inout) :: poly

          poly%cplx       = .FALSE.
          poly%npoly      = 0
          poly%ndim       = 0
          poly%nb_coef    = 0

          poly%init0      = .TRUE.
          poly%alloc_poly = .FALSE.

          nullify(poly%Rcoef)
          nullify(poly%Ccoef)
          nullify(poly%ind)

        END SUBROUTINE init0_poly

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE alloc_poly(poly,npoly,ndim,cplx)
           TYPE (para_poly), intent(inout)     :: poly
           integer, intent(in)                 :: npoly
           integer, intent(in), optional       :: ndim
           logical, intent(in), optional       :: cplx

           integer :: i,k,ip

           integer, pointer :: ind_exp(:)

           IF (.NOT. poly%init0) THEN
             write(out_unitp,*) ' ERROR in alloc_poly'
             write(out_unitp,*) ' poly has NOT been initiated with init0_poly'
             STOP
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
             poly%Ccoef(:)    = 0.d0
           ELSE
             allocate(poly%Rcoef(poly%nb_coef))
             poly%Rcoef(:)    = 0.d0
           END IF


           poly%alloc_poly = .TRUE.

        END SUBROUTINE alloc_poly

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE dealloc_poly(poly)
           TYPE (para_poly), intent(inout)     :: poly


           IF (.NOT. poly%init0) THEN
             write(out_unitp,*) ' ERROR in alloc_poly'
             write(out_unitp,*) ' poly has NOT been initiated with init0_poly'
             STOP
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

           CALL init0_poly(poly)

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


           IF (.NOT. poly%init0) THEN
             write(out_unitp,*) ' ERROR in build_ind_poly'
             write(out_unitp,*) ' poly has NOT been initiated with init0_poly'
             write(out_unitp,*) 'init0',poly%init0
             STOP
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

           IF (.NOT. poly%init0 .OR. .NOT. poly%alloc_poly) THEN
             write(out_unitp,*) ' ERROR in locate'
             write(out_unitp,*) ' poly has NOT been initiated with init0_poly'
             write(out_unitp,*) ' or is NOT allocated !'
             write(out_unitp,*) 'init0,alloc_poly',poly%init0,poly%alloc_poly
             STOP
           END IF

           locate = -1


           DO i=1,poly%nb_coef

             IF ( sum(abs(ind(:)-poly%ind(:,i))) == 0 ) THEN
                locate = i
                EXIT
             END IF

           END DO


!          write(out_unitp,*) 'ind,i',ind,i


        END FUNCTION locate
!       ==============================================================
!         write_poly : write a polynomial
!       ==============================================================
      !!@description: write_poly : write a polynomial
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE write_poly(poly)
           TYPE (para_poly), intent(in) :: poly

           integer :: i


           IF (.NOT. poly%init0) THEN
             write(out_unitp,*) ' ERROR in write_poly'
             write(out_unitp,*) ' poly has NOT been initiated with init0_poly'
             STOP
           END IF

           write(out_unitp,*) '--------------------------------------'
           write(out_unitp,*) 'ndim,npoly,nb_coef',                             &
                     poly%ndim,poly%npoly,poly%nb_coef
           write(out_unitp,*) 'init0,alloc_poly',poly%init0,poly%alloc_poly
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


        END SUBROUTINE write_poly

!       ==============================================================
!         p1PLUSp2TOp3 : p3 = p1 + p2
!       ==============================================================

      !!@description: p1PLUSp2TOp3 : p3 = p1 + p2
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE p1PLUSp2TOp3(p1,p2,p3)
           TYPE (para_poly), intent(in)     :: p1,p2
           TYPE (para_poly), intent(inout)    :: p3

!          working variables
           integer :: npoly3,ndim3,nb_coef3 ! parameters of p3
           logical :: cplx

!          CALL write_poly(p1)
!          CALL write_poly(p2)

           IF (.NOT. p1%init0 .OR. .NOT. p1%alloc_poly) THEN
             write(out_unitp,*) ' ERROR in p1PLUSp2TOp3'
             write(out_unitp,*) ' p1 has NOT been initiated with init0_poly'
             write(out_unitp,*) ' or is NOT allocated !'
             write(out_unitp,*) 'init0,alloc_poly',p1%init0,p1%alloc_poly
             STOP
           END IF
           IF (.NOT. p2%init0 .OR. .NOT. p2%alloc_poly) THEN
             write(out_unitp,*) ' ERROR in p1PLUSp2TOp3'
             write(out_unitp,*) ' p2 has NOT been initiated with init0_poly'
             write(out_unitp,*) ' or is NOT allocated !'
             write(out_unitp,*) 'init0,alloc_poly',p2%init0,p2%alloc_poly
             STOP
           END IF
           IF (.NOT. p3%init0) THEN
             write(out_unitp,*) ' ERROR in p1PLUSp2TOp3'
             write(out_unitp,*) ' p3 has NOT been initiated with init0_poly'
             STOP
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
             p3%Ccoef(:)  = 0.d0
             IF (p1%cplx) THEN
               p3%Ccoef(1:p1%nb_coef)  = p1%Ccoef(:)
             ELSE
               p3%Ccoef(1:p1%nb_coef)  = p1%Rcoef(:)
             END IF
             IF (p2%cplx) THEN
               p3%Ccoef(1:p2%nb_coef)  = p3%Ccoef(1:p2%nb_coef) +       &
                                         p2%Ccoef(:)
             ELSE
               p3%Ccoef(1:p2%nb_coef)  = p3%Ccoef(1:p2%nb_coef) +       &
                                         p2%Rcoef(:)
             END IF
           ELSE
             p3%Rcoef(:)  = 0.d0
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


        END SUBROUTINE p1PLUSp2TOp3
!       ==============================================================
!         p1TIMEp2TOp3 : p3 = p1 * p2
!       ==============================================================

      !!@description: p1TIMEp2TOp3 : p3 = p1 * p2
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE p1TIMEp2TOp3(p1,p2,p3)
           TYPE (para_poly), intent(in)     :: p1,p2
           TYPE (para_poly), intent(inout)    :: p3

!          working variables
           integer :: npoly3,ndim3,nb_coef3 ! parameters of p3
           integer :: ind3(p1%ndim)
           integer :: i1,i2,i3
           logical :: cplx

!          CALL write_poly(p1)
!          CALL write_poly(p2)

           IF (.NOT. p1%init0 .OR. .NOT. p1%alloc_poly) THEN
             write(out_unitp,*) ' ERROR in p1TIMEp2TOp3'
             write(out_unitp,*) ' p1 has NOT been initiated with init0_poly'
             write(out_unitp,*) ' or is NOT allocated !'
             write(out_unitp,*) 'init0,alloc_poly',p1%init0,p1%alloc_poly
             STOP
           END IF
           IF (.NOT. p2%init0 .OR. .NOT. p2%alloc_poly) THEN
             write(out_unitp,*) ' ERROR in p1TIMEp2TOp3'
             write(out_unitp,*) ' p2 has NOT been initiated with init0_poly'
             write(out_unitp,*) ' or is NOT allocated !'
             write(out_unitp,*) 'init0,alloc_poly',p2%init0,p2%alloc_poly
             STOP
           END IF
           IF (.NOT. p3%init0) THEN
             write(out_unitp,*) ' ERROR in p1TIMEp2TOp3'
             write(out_unitp,*) ' p3 has NOT been initiated with init0_poly'
             STOP
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
             p3%Ccoef(:)  = 0.d0
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

                 p3%Ccoef(i3) = p3%Ccoef(i3) + p1%Rcoef(i1)*p2%Ccoef(i2)
               END DO
               END DO
             ELSE IF (p1%cplx .AND. .NOT. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Ccoef(i3) = p3%Ccoef(i3) + p1%Ccoef(i1)*p2%Rcoef(i2)
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
             p3%Rcoef(:)  = 0.d0

             IF (p1%cplx .AND. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Rcoef(i3) = p3%Rcoef(i3) + p1%Ccoef(i1)*p2%Ccoef(i2)
               END DO
               END DO
             ELSE IF (.NOT. p1%cplx .AND. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Rcoef(i3) = p3%Rcoef(i3) + p1%Rcoef(i1)*p2%Ccoef(i2)
               END DO
               END DO
             ELSE IF (p1%cplx .AND. .NOT. p2%cplx) THEN
               DO i1=1,p1%nb_coef
               DO i2=1,p2%nb_coef

                 ind3(:) = p1%ind(:,i1) + p2%ind(:,i2)
                 i3 = locate(p3,ind3)

                 p3%Rcoef(i3) = p3%Rcoef(i3) + p1%Ccoef(i1)*p2%Rcoef(i2)
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


        END SUBROUTINE p1TIMEp2TOp3

!       ==============================================================
!         p2 = d1p1 (first derivative with respect to the variable id)
!       ==============================================================

      !!@description:  p2 = d1p1 (first derivative with respect to the variable
      !!               id)
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

!          CALL write_poly(p1)
!          CALL write_poly(p2)

           IF (.NOT. p1%init0 .OR. .NOT. p1%alloc_poly) THEN
             write(out_unitp,*) ' ERROR in d1p1TOp2'
             write(out_unitp,*) ' p1 has NOT been initiated with init0_poly'
             write(out_unitp,*) ' or is NOT allocated !'
             write(out_unitp,*) 'init0,alloc_poly',p1%init0,p1%alloc_poly
             STOP
           END IF
           IF (.NOT. p2%init0) THEN
             write(out_unitp,*) ' ERROR in d1p1TOp2'
             write(out_unitp,*) ' p2 has NOT been initiated with init0_poly'
             STOP
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

             p2%Ccoef(:) = 0.d0

             DO i1=2,p1%nb_coef
               IF (p1%ind(id,i1) == 0) CYCLE
               exp_q(:) = p1%ind(:,i1)
               exp_q(id) = exp_q(id) - 1

               i2 = locate(p2,exp_q)
               p2%Ccoef(i2) = p1%Ccoef(i1) * real( p1%ind(id,i1),kind=Rkind)

             END DO

           ELSE IF (p2%cplx .AND. .NOT. p1%cplx) THEN

             p2%Ccoef(:) = 0.d0

             DO i1=2,p1%nb_coef
               IF (p1%ind(id,i1) == 0) CYCLE
               exp_q(:) = p1%ind(:,i1)
               exp_q(id) = exp_q(id) - 1

               i2 = locate(p2,exp_q)
               p2%Ccoef(i2) = p1%Rcoef(i1) * real( p1%ind(id,i1),kind=Rkind)

             END DO

           ELSE ! p1 and p2 are real

             p2%Rcoef(:) = 0.d0

             DO i1=2,p1%nb_coef
               IF (p1%ind(id,i1) == 0) CYCLE
               exp_q(:) = p1%ind(:,i1)
               exp_q(id) = exp_q(id) - 1

               i2 = locate(p2,exp_q)
               p2%Rcoef(i2) = p1%Rcoef(i1) * real( p1%ind(id,i1),kind=Rkind)

             END DO

           END IF


        END SUBROUTINE d1p1TOp2
!       ==============================================================
!       p1_O_p2TOp3   : p3 = p1(+linear transfo) Rq linearTransfo : aX+b
!       ivar : X associated with the index ivar
!       ==============================================================

      !!@description:  p1_O_p2TOp3   : p3 = p1(+linear transfo) Rq linearTransfo
      !!               : aX+b ivar : X associated with the index ivar
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        SUBROUTINE p1_linearTransfoTOp3(p1,ivar,Ra,Ca,Rb,Cb,p3,cplx)
           TYPE (para_poly), intent(in)     :: p1
           real (kind=Rkind), intent(in)        :: Ra,Rb
           complex (kind=Rkind), intent(in)     :: Ca,Cb
           TYPE (para_poly), intent(inout)    :: p3
           logical, intent(in) :: cplx
           integer, intent(in) :: ivar


!          working variables
           integer :: npoly3,ndim3,nb_coef3 ! parameters of p3
           TYPE (para_poly)   :: pwork
           integer, pointer :: ind_exp(:)
           complex (kind=Rkind)    :: CCa,CCb
           integer :: exp_ivar,i1,i3,k

           real (kind=Rkind) :: binomial ! function

!          CALL write_poly(p1)

           IF (.NOT. p1%init0 .OR. .NOT. p1%alloc_poly) THEN
             write(out_unitp,*) ' ERROR in p1_linearTransfoTOp3'
             write(out_unitp,*) ' p1 has NOT been initiated with init0_poly'
             write(out_unitp,*) ' or is NOT allocated !'
             write(out_unitp,*) 'init0,alloc_poly',p1%init0,p1%alloc_poly
             STOP
           END IF
           IF (.NOT. p3%init0) THEN
             write(out_unitp,*) ' ERROR in p1_linearTransfoTOp3'
             write(out_unitp,*) ' p3 has NOT been initiated with init0_poly'
             STOP
           END IF

           IF (p3%alloc_poly) THEN
             IF (p3%ndim /= p1%ndim) THEN
               write(out_unitp,*) ' ERROR in p1_linearTransfoTOp3'
               write(out_unitp,*) ' p3 is allocated so'
               write(out_unitp,*) ' ndim of p3 and p1 MUST be identical'
               write(out_unitp,*) ' ndim of p1 p3:',p1%ndim,p3%ndim
               STOP
             END IF

             npoly3 = p1%npoly
             IF (p3%npoly < npoly3) THEN
               write(out_unitp,*) ' ERROR in p1_linearTransfoTOp3'
               write(out_unitp,*) ' p3 is allocated so'
               write(out_unitp,*) ' p3%npoly MUST be >= p1%npoly'
               write(out_unitp,*) ' npoly p1 p3:',p1%npoly,p3%npoly
               STOP
             END IF
           ELSE
             CALL dealloc_poly(p3)
             npoly3 = p1%npoly
             CALL alloc_poly(p3,ndim=p1%ndim,npoly=npoly3,              &
                                 cplx=(p1%cplx .OR. cplx))
           END IF

           IF (p1%cplx .OR. cplx) THEN
             IF (.NOT. p3%cplx) THEN
               write(out_unitp,*) ' ERROR in p1_linearTransfoTOp3'
               write(out_unitp,*) ' p3 MUST be complex'
               write(out_unitp,*) ' p3%cplx',p3%cplx
               STOP
             END IF
           END IF


           IF (ivar < 0 .OR. ivar > p1%ndim) THEN
             write(out_unitp,*) ' ERROR in p1_linearTransfoTOp3'
             write(out_unitp,*) ' ivar MUST be > 0 or < ndim'
             write(out_unitp,*) 'ivar,ndim',ivar,p1%ndim
             STOP
           END IF

           allocate(ind_exp(0:p1%ndim))

           IF (cplx .AND. p3%cplx) THEN
            CCa = Ca
            CCb = Cb
           ELSE IF (.NOT. cplx .AND. p3%cplx) THEN
            CCa = Ra
            CCb = Rb
           END IF

           IF (p3%cplx) THEN
             p3%Ccoef(:)  = 0.d0
             IF (p1%cplx) THEN
               DO i1=1,p1%nb_coef
                 ind_exp(0) = 0
                 ind_exp(1:p1%ndim) = p1%ind(:,i1)
                 exp_ivar = p1%ind(ivar,i1)

                 write(out_unitp,*) 'i1,exp_ivar,p1%ind',                       &
                             i1,exp_ivar,p1%ind(:,i1)
!                devlopment of coef*(aX+b)^exp_ivar
                 DO k=0,exp_ivar
                   ind_exp(ivar) = k
                   ind_exp(0) = sum(ind_exp(1:p1%ndim))
                   i3 = locate(p3,ind_exp(1:p1%ndim))
                   write(out_unitp,*) 'k,ind_exp,i3,p3%ind',                    &
                      k,ind_exp(1:p1%ndim),i3,p3%ind(:,i3)
                   p3%Ccoef(i3) = p3%Ccoef(i3) +                        &
              p1%Ccoef(i1)*binomial(exp_ivar,k)*CCa**k*CCb**(exp_ivar-k)
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
                   i3 = locate(p3,ind_exp(1:p1%ndim))
                   p3%Ccoef(i3) = p3%Ccoef(i3) +                        &
              p1%Rcoef(i1)*binomial(exp_ivar,k)*CCa**k*CCb**(exp_ivar-k)
                 END DO
               END DO
             END IF
           ELSE
             p3%Rcoef(:)  = 0.d0
             DO i1=1,p1%nb_coef
               ind_exp(0) = 0
               ind_exp(1:p1%ndim) = p1%ind(:,i1)
               exp_ivar = p1%ind(ivar,i1)
!              devlopment of coef*(aX+b)^exp_ivar
               DO k=0,exp_ivar
                 ind_exp(ivar) = k
                 ind_exp(0) = sum(ind_exp(1:p1%ndim))
                 i3 = locate(p3,ind_exp(1:p1%ndim))
                 p3%Rcoef(i3) = p3%Rcoef(i3) +                          &
              p1%Rcoef(i1)*binomial(exp_ivar,k)*Ra**k*Rb**(exp_ivar-k)
                 END DO
               END DO
           END IF

           deallocate(ind_exp)

           CALL write_poly(p3)
        END SUBROUTINE p1_linearTransfoTOp3

      END MODULE mod_poly

!=======================================================
      PROGRAM test
      USE mod_poly
      IMPLICIT NONE

      type (para_poly) :: p1,p2,p3

      integer :: ind(2),loc
      integer :: i,ip,iq,ind_exp(0:5)

      CALL init0_poly(p1)
      CALL init0_poly(p2)
      CALL init0_poly(p3)



      CALL alloc_poly(p1,npoly=1,cplx=.TRUE.)
      CALL alloc_poly(p2,npoly=2)
      CALL alloc_poly(p3,npoly=5)


      p1%Ccoef(:) = (/ (1.,1.),(2.,0.) /)
      p2%Rcoef(:) = (/ 0.,0.,1. /)

      CALL p1PLUSp2TOp3(p2,p1,p3)

      CALL write_poly(p1)
      CALL write_poly(p2)
      CALL write_poly(p3)

      CALL dealloc_poly(p1)
      CALL dealloc_poly(p2)
      CALL dealloc_poly(p3)


      write(out_unitp,*) '========================'
      write(out_unitp,*) '= plus ================='
      write(out_unitp,*) 'ndim = 2'
!     CALL alloc_poly(p1,npoly=1,ndim=2,cplx=.TRUE.)
      CALL alloc_poly(p1,npoly=1,ndim=2)
      CALL alloc_poly(p2,npoly=2,ndim=2)
      CALL alloc_poly(p3,npoly=5,ndim=2)


      p1%Rcoef(:) = (/ 0.,   1.0, -1.0/)
      p2%Rcoef(:) = (/ 0.,   0.,3. , 1.,0.,3./)

      CALL p1PLUSp2TOp3(p2,p1,p3)

      CALL write_poly(p1)
      CALL write_poly(p2)
      CALL write_poly(p3)

      CALL dealloc_poly(p1)
      CALL dealloc_poly(p2)
      CALL dealloc_poly(p3)

      write(out_unitp,*) '========================'
      write(out_unitp,*) '= time ================='
      write(out_unitp,*) 'ndim = 2'
      CALL alloc_poly(p1,npoly=1,ndim=2,cplx=.TRUE.)
      CALL alloc_poly(p2,npoly=2,ndim=2,cplx=.TRUE.)
      CALL alloc_poly(p3,npoly=5,ndim=2,cplx=.TRUE.)


      p1%Ccoef(:) = (/ (0.,0.),  ( 1.0,-1.),(-1.0,3.)/)

      p2%Ccoef(:) = (/ (3.,1.),  (1.,0.), (-1.,2.),                     &
                       (2.,2.), (4.,1.),(-3.,0.)/)

      CALL p1TIMEp2TOp3(p2,p1,p3)

      CALL write_poly(p1)
      CALL write_poly(p2)
      CALL write_poly(p3)

      CALL dealloc_poly(p1)
      CALL dealloc_poly(p2)
      CALL dealloc_poly(p3)

      write(out_unitp,*) '========================'
      write(out_unitp,*) '==== derivee ==========='
      write(out_unitp,*) 'ndim = 2'

      CALL alloc_poly(p1,npoly=2,ndim=2,cplx=.TRUE.)
      CALL alloc_poly(p2,npoly=1,ndim=2,cplx=.TRUE.)

      p1%Ccoef(:) = (/ (3.,1.),  (1.,0.), (-1.,2.),                     &
                       (2.,2.), (4.,1.),(-3.,0.)/)
      CALL write_poly(p1)
      CALL write_poly(p2)

      CALL d1p1TOp2(p1,p2,1)
      CALL write_poly(p2)
      CALL d1p1TOp2(p1,p2,2)
      CALL write_poly(p2)

      CALL dealloc_poly(p1)
      CALL dealloc_poly(p2)


      write(out_unitp,*) '========================'
      write(out_unitp,*) '==== linear transfo ===='
      write(out_unitp,*) 'ndim = 2'

      CALL alloc_poly(p1,npoly=2,ndim=2,cplx=.TRUE.)
      CALL alloc_poly(p2,npoly=2,ndim=2,cplx=.TRUE.)

      p1%Ccoef(:) = (/ (3.,1.),  (1.,0.), (-1.,2.),                     &
                       (2.,2.), (4.,1.),(-3.,0.)/)
      CALL write_poly(p1)

      CALL p1_linearTransfoTOp3(p1,ivar=1,Ra=0.d0,Ca=(2.d0,0d0),Rb=0.d0,&
                           Cb=(1.d0,0.d0),p3=p2,cplx=.true.)
      CALL write_poly(p2)

      CALL dealloc_poly(p1)
      CALL dealloc_poly(p2)

 20   CONTINUE
      write(out_unitp,*) '========================'
      write(out_unitp,*) '==== ind     ==========='
      write(out_unitp,*) 'ndim = 3'

      CALL alloc_poly(p1,npoly=4,ndim=3,cplx=.TRUE.)
      CALL write_poly(p1)
      CALL dealloc_poly(p1)
      STOP

      write(out_unitp,*) '========================'
      write(out_unitp,*) '========================'
      write(out_unitp,*) 'locate  '

      CALL alloc_poly(p1,npoly=2,ndim=2,cplx=.TRUE.)
      p1%Ccoef(:) = (/ (3.,1.),  (1.,0.), (-1.,2.),                     &
                       (2.,2.), (4.,1.),(-3.,0.)/)
      ind(:) = (/ 0,0 /)

 10   CONTINUE
!     read(in_unitp,*) ind


      loc = locate(p1,ind)
      write(out_unitp,*) 'loc',loc

!     GOTO 10

      CALL dealloc_poly(p1)

      END
         real(kind=Rkind) FUNCTION binomial(n,i)
         real(kind=Rkind) a
         integer i,k,n
         IF (n .LT. 0 .OR. i .GT. n .OR. i .LT. 0) THEN
           write(out_unitp,*) 'ERROR: binomial( n<0 i<0 i>n)',n,i
           STOP
         END IF
         a = 1.d0
         DO k=1,n
           a = a * real(k,kind=Rkind)
         END DO
         DO k=1,n-i
           a = a / real(k,kind=Rkind)
         END DO
         DO k=1,i
           a = a / real(k,kind=Rkind)
         END DO
         binomial = a

!        write(out_unitp,*) 'binomial',n,i,a
         END

