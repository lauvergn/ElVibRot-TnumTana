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
      implicit none

      integer :: npts,nptE
      real (kind=8) :: t0,tmax,dt
      real (kind=8) :: E,Emin,Emax,dE
      real (kind=8) :: conv
      real (kind=8), pointer :: t(:)
      complex (kind=8), pointer :: auto(:)
      complex (kind=8), pointer :: Expiwt(:)
      complex (kind=8), pointer :: funcE(:)
      character (len=50) :: file_auto

      real (kind=8) :: a,b
      integer :: i

      complex (kind=8), parameter :: EYE = (0,1)
      real (kind=8), parameter ::                                       &
       pi = 3.14159265358979323846264338327950288419716939937511d0


      namelist / param / Emin,Emax,conv,file_auto


      file_auto = 'file_auto'
      conv = 1.d0
      Emin = 0.d0
      Emax = -1.d0
      read(in_unitp,param)
      write(out_unitp,param)

!     read the time function (autocorrelation...)

       open(unit=10,file=file_auto)
       read(10,*) npts

          allocate(auto(npts))
          allocate(t(npts))

          DO i=1,npts
            read(10,*) t(i),a,b
            auto(i) = cmplx(a,b,kind=8)
          END DO
          dt = t(2)-t(1)
          tmax = t(npts)
          close(10)
          write(out_unitp,*) 'npts t0,tmax,dt: ',npts,t0,tmax,dt
!     END read the time function

!     Check the funcErgy grid...
          dE = 2.d0*Pi/tmax
          IF (Emax <0) Emax=2.d0*Pi/dt
          IF (Emax> 2.d0*Pi/dt) THEN
            write(out_unitp,*) 'Emax> 2.d0*Pi/dt',Emax,2.d0*Pi/tmax
            STOP
          END IF
          nptE = int((Emax-Emin)/dE)

          allocate(funcE(nptE))
          allocate(Expiwt(npts)
          DO i=1,nptE
            E=Emin+real(i-1,kind=8)*dE
            Expiwt(:) = exp(EYE*E*t(:))

            funcE(i) = 2.d0*dt*sum(Expiwt(:)*auto(:))
            write(out_unitp,*) E,funcE(i)
          END DO



          deallocate(funcE)
          deallocate(Expiwt)
          deallocate(t)
          deallocate(auto)
      END

