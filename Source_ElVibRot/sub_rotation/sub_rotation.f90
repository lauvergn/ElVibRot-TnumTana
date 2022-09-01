!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015 David Lauvergnat [1]
!      with contributions of
!        Josep Maria Luis (optimization) [2]
!        Ahai Chen (MPI) [1,4]
!        Lucien Dupuy (CRP) [5]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Institut de Química Computacional and Departament de Química,
!        Universitat de Girona, Catalonia, Spain
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!
!===========================================================================
!===========================================================================

      IMPLICIT NONE
      real (kind=8) :: gam_xyz(3),mu_xyz(3,3)
      complex (kind=8) :: gam_pmz(3),mu_pmz(3,3)
      integer :: JJ,K,M,Jout,Kout,Jout1,Kout1
      integer :: indK,indKout
      integer :: nOp1,nOp2
      Real (kind=8) :: Proj1,Proj2

      integer :: dim_rot
      complex (kind=8),pointer :: Hrot(:,:),Vecrot(:,:)

      real (kind=8) :: ZERO = 0.
      real (kind=8) :: ONE  = 1.
      real (kind=8) :: TWO  = 2.
      real (kind=8) :: HALF = 0.5
      complex (kind=8) :: EYE = cmplx(0.,1.,kind=8)


!     read gamma and mu with Jx, Jy, Jz operators
      read(in_unitp,*) gam_xyz(:)
      read(in_unitp,*) mu_xyz(:,:)
      read(in_unitp,*) JJ
      write(out_unitp,*) gam_xyz(:)
      write(out_unitp,*) mu_xyz(:,:)
      write(out_unitp,*) JJ

!     rewrite gamma and mu with J+, J-, Jz operators
      gam_pmz(1) = HALF*(EYE*gam_xyz(1)+gam_xyz(2))
      gam_pmz(2) = HALF*(EYE*gam_xyz(1)-gam_xyz(2))
      gam_pmz(3) = EYE*gam_xyz(3)


      mu_pmz(3,3) = mu_xyz(3,3)
      mu_pmz(1,3) = HALF*(EYE*mu_xyz(1,3)+mu_xyz(2,3))
      mu_pmz(3,1) = mu_pmz(1,3)
      mu_pmz(2,3) = HALF*(EYE*mu_xyz(1,3)-mu_xyz(2,3))
      mu_pmz(3,2) = mu_pmz(2,3)
      mu_pmz(1,1) = HALF**2 * ( mu_xyz(1,1) - mu_xyz(2,2) -             &
                                        TWO*EYE*mu_xyz(1,2) )
      mu_pmz(2,2) = HALF**2 * ( mu_xyz(1,1) - mu_xyz(2,2) +             &
                                        TWO*EYE*mu_xyz(1,2) )
      mu_pmz(1,2) = HALF**2 * ( mu_xyz(1,1) + mu_xyz(2,2) )
      mu_pmz(2,1) = mu_pmz(1,2)

      write(out_unitp,*) 'gam_pmz',gam_pmz
      write(out_unitp,*) 'mu_pmz',mu_pmz

      IF (JJ < 0) THEN
        write(out_unitp,*) ' ERROR in sub_rotation'
        write(out_unitp,*) ' JJ < 0',JJ
        STOP
      END IF

      dim_rot = 2*JJ+1  ! K in [-J.... 0 .... J]
      memory = product( [dim_rot,dim_rot] )
      allocate(Hrot(dim_rot,dim_rot),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"Hrot","main")
      memory = product( [dim_rot,dim_rot] )
      allocate(Vecrot(dim_rot,dim_rot),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"Vecrot","main")



      Hrot(:,:) = ZERO
      DO K=-JJ,JJ
!       - For gamma (J+,J-,Jz)
        DO nOp1=1,3
          CALL Rot_Op(K,JJ,Kout,Jout,Proj1,nOp1)
          indK    = K    + JJ+1
          indKout = Kout + JJ+1
          IF (abs(Kout) <= JJ)                                          &
                 Hrot(indK,indKout) = Hrot(indK,indKout) -              &
                                           gam_pmz(nOp1)*Proj1
        END DO

!       - For mu (J+^2,J-^2,Jz^2,J+J-....)
        DO nOp1=1,3
        DO nOp2=1,3
          CALL Rot_Op(K,JJ,Kout1,Jout1,Proj1,nOp1)
          CALL Rot_Op(Kout1,Jout1,Kout,Jout,Proj2,nOp2)
          indK    = K    + JJ+1
          indKout = Kout + JJ+1
          IF (abs(Kout) <= JJ)                                          &
                 Hrot(indK,indKout) = Hrot(indK,indKout) +              &
                               HALF*mu_pmz(nOp1,nOp2)*Proj1*Proj2
        END DO
        END DO
      END DO


      DO indK=1,dim_rot
      write(out_unitp,*) Hrot(:,indK)
      END DO

      memory = size(Hrot)
      deallocate(Hrot,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"Hrot","main")
      memory = size(Vecrot)
      deallocate(Vecrot,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"Vecrot","main")
      END
      SUBROUTINE Rot_Op(Kin,Jin,Kout,Jout,Proj,nOp)
      integer, intent(in) :: Kin,Jin
      integer, intent(in) :: nOp ! 1,2,3 => J+, J-, Jz
      integer, intent(inout) :: Kout,Jout
      real (kind=8), intent(inout) :: Proj


      Jout = Jin
      Kout = Kin
      Proj = 111111.
      IF (abs(Kin) > Jin) RETURN

      IF (nOp == 3) THEN !Jz
        Kout = Kin
        Proj = real(Kin,kind=8)
      ELSE IF (nOp == 1) THEN !J+
        Kout = Kin-1
        IF (Kout >= -Jin) Proj = sqrt(real(Jin*(Jin+1),kind=8)-         &
                                      real(Kin*(Kin-1),kind=8))
      ELSE IF (nOp == 2) THEN !J-
        Kout = Kin+1
        IF (Kout <= Jin) Proj = sqrt(real(Jin*(Jin+1),kind=8)-          &
                                     real(Kin*(Kin+1),kind=8))
      ELSE
        write(out_unitp,*) ' ERROR in Rot_Op'
        write(out_unitp,*) ' nOp MUST be : 1,2,3',nOp
        write(out_unitp,*) ' Remark : 1,2,3 => J+, J-, Jz'
        STOP
      END IF

!     write(out_unitp,*) 'Jin,Kin,nOp',Jin,Kin,nOp
!     write(out_unitp,*) 'Jout,Kout,Proj',Jout,Kout,Proj
      write(out_unitp,*) 'Kin,Kout,nOp,Proj',Kin,Kout,nOp,Proj

      end subroutine Rot_Op

