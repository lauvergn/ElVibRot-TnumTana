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

