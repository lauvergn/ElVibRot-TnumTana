!  quadprb.f  25 September 1998
!
      program quadprb
!
!***********************************************************************
!
!! QUADPRB calls a set of tests for the QUADRULE library.
!
      write ( *, * ) ' '
      write ( *, * ) 'QUADPRB'
      write ( *, * ) '  Sample problems for the QUADRULE library.'
      write ( *, * ) ' '

      call test01
      call test02
      call test03
      call test04
      call test05
      call test06
      call test07
      call test08
      call test09
      call test10
      call test11
      call test12
      call test13
      call test14
      call test15
      call test16
      call test17
      call test18
      call test19
      call test20
      call test21
      call test22
      call test23
      call test24

      stop
      end program quadprb
      subroutine test01
!
!***********************************************************************
!
!! TEST01 tests CHEBSET and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 9 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      a = 0.0 D+00
      b = 1.0 D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST01'
      write ( *, * ) '  CHEBSET sets up Chebyshev quadrature;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  The integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        if ( norder .ne. 8 ) then

          do i = 1, nfunc

            call funcset ( 'SET', i )

            call chebset ( norder, xtab, weight )

            call sumsub ( func, a, b, nsub, norder, xtab, weight,       &
              result(i) )

          end do

          write ( *, '(i2,2x,7f10.6)' )                                 &
            norder, ( result(i), i = 1, nfunc )

        end if

      end do

      return
      end subroutine test01
      subroutine test02
!
!***********************************************************************
!
!! TEST02 tests LEGCOM and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 10 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST02'
      write ( *, * ) '  LEGCOM sets up Gauss-Legendre quadrature;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  The integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call legcom ( norder, xtab, weight )

          call sumsub ( func, a, b, nsub, norder, xtab, weight,         &
            result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test02
      subroutine test03
!
!***********************************************************************
!
!! TEST03 tests LEGSET and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST03'
      write ( *, * ) '  LEGSET sets up Gauss-Legendre quadrature;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  The integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call legset ( norder, xtab, weight )

          call sumsub ( func, a, b, nsub, norder, xtab, weight,         &
            result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test03
      subroutine test04
!
!***********************************************************************
!
!! TEST04 tests LOBSET and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST04'
      write ( *, * ) '  LOBSET sets up Lobatto quadrature;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 2, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call lobset ( norder, xtab, weight )

          call sumsub ( func, a, b, nsub, norder, xtab, weight,         &
            result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test04
      subroutine test05
!
!***********************************************************************
!
!! TEST05 tests NCCSET and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST05'
      write ( *, * ) '  NCCSET sets up a closed Newton Cotes rule;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 2, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call nccset ( norder, xtab, weight )

          call sumsub ( func, a, b, nsub, norder, xtab, weight,         &
            result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test05
      subroutine test06
!
!***********************************************************************
!
!! TEST06 tests NCCCOM and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST06'
      write ( *, * ) '  NCCCOM computes a closed Newton Cotes rule;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 2, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call ncccom ( norder, xtab, weight )

          call sumsub ( func, a, b, nsub, norder, xtab, weight,         &
            result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test06
      subroutine test07
!
!***********************************************************************
!
!! TEST07 tests NCOSET and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 9 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST07'
      write ( *, * ) '  NCOSET sets up an open Newton-Cotes rule;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        if ( norder .le. 7 .or. norder .eq. 9 ) then

          do i = 1, nfunc

            call funcset ( 'SET', i )

            call ncoset ( norder, xtab, weight )

            call sumsub ( func, a, b, nsub, norder, xtab, weight,       &
              result(i) )

          end do

          write ( *, '(i2,2x,7f10.6)' )                                 &
            norder, ( result(i), i = 1, nfunc )

        end if

      end do

      return
      end subroutine test07
      subroutine test08
!
!***********************************************************************
!
!! TEST08 tests NCOCOM and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 9 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST08'
      write ( *, * ) '  NCOCOM sets up an open Newton-Cotes rule;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        if ( norder .le. 7 .or. norder .eq. 9 ) then

          do i = 1, nfunc

            call funcset ( 'SET', i )

            call ncocom ( norder, xtab, weight )

            call sumsub ( func, a, b, nsub, norder, xtab, weight,       &
              result(i) )

          end do

          write ( *, '(i2,2x,7f10.6)' )                                 &
            norder, ( result(i), i = 1, nfunc )

        end if

      end do

      return
      end subroutine test08
      subroutine test09
!
!***********************************************************************
!
!! TEST09 tests RADSET and SUMSUB.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 15 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func

      a = 0.0D+00
      b = 1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST09'
      write ( *, * ) '  RADSET sets up Radau quadrature;'
      write ( *, * ) '  SUMSUB carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  The integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call radset ( norder, xtab, weight )

          call sumsub ( func, a, b, nsub, norder, xtab, weight,         &
            result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test09
      subroutine test10
!
!***********************************************************************
!
!! TEST10 tests BASHSET and SUMMER.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 8 )
      parameter ( nfunc = 7 )
!
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST10'
      write ( *, * )                                                    &
        '  BASHSET sets up Adams-Bashforth quadrature;'
      write ( *, * ) '  SUMMER carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is [0,1].'
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call bashset ( norder, xtab, weight )

          call summer ( func, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test10
      subroutine test11
!
!***********************************************************************
!
!! TEST11 tests BDFPSET and BDFSUM.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 10 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) diftab(maxorder)
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST11'
      write ( *, * )                                                    &
        '  BDFPSET sets up Backward Difference Predictor quadrature;'
      write ( *, * ) '  BDFSUM carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is [0,1].'
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call bdfpset ( norder, weight, xtab )

          call bdfsum ( func, norder, weight, xtab, diftab, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test11
      subroutine test12
!
!***********************************************************************
!
!! TEST12 tests BDFCSET and BDFSUM.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 10 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) diftab(maxorder)
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)

      external func
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST12'
      write ( *, * )                                                    &
        '  BDFCSET sets up Backward Difference Corrector quadrature;'
      write ( *, * ) '  BDFSUM carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is [0,1].'
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call bdfcset ( norder, weight, xtab )

          call bdfsum ( func, norder, weight, xtab, diftab, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test12
      subroutine test13
!
!***********************************************************************
!
!! TEST13 tests LAGSET and LAGSUM.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      a = 1.0D+00

      write ( *, * ) ' '
      write ( *, * ) 'TEST13'
      write ( *, * ) '  LAGSET sets up Gauss-Laguerre quadrature;'
      write ( *, * ) '  LAGSUM carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  The integration interval is [ ',                &
        a, ', +Infinity ).'
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) '  The weight function is EXP ( - X ).'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call lagset ( norder, xtab, weight )

          call lagsum ( func, a, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test13
      subroutine test14
!
!***********************************************************************
!
!! TEST14 tests LAGCOM and LAGSUM.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) alpha
      real(kind=Rkind) b(maxorder)
      real(kind=Rkind) c(maxorder)
      real(kind=Rkind) csa
      real(kind=Rkind) csx
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) tsa
      real(kind=Rkind) tsx
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      a = 1.0D+00

      write ( *, * ) ' '
      write ( *, * ) 'TEST14'
      write ( *, * ) '  LAGCOM computes a Gauss-Laguerre rule;'
      write ( *, * ) '  LAGSUM carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  The integration interval is [ ',                &
        a, ', +Infinity ).'
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) '  The weight function is EXP ( - X ).'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      alpha = 0.0

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call lagcom ( norder, xtab, weight, alpha, b, c, csx, csa,    &
            tsx, tsa )

          call lagsum ( func, a, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test14
      subroutine test15
!
!***********************************************************************
!
!! TEST15 tests HERSET and SUMMER.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST15'
      write ( *, * ) '  HERSET sets up Gauss-Hermite quadrature;'
      write ( *, * ) '  SUMMER carries it out.'
      write ( *, * ) ' '
      write ( *, * )                                                    &
        '  The integration interval is ( -Infinity, +Infinity ).'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call herset ( norder, xtab, weight )

          call summer ( func, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
           norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test15
      subroutine test16
!
!***********************************************************************
!
!! TEST16 tests HERCOM and SUMMER.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 20 )
      parameter ( nfunc = 7 )
!
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST16'
      write ( *, * ) '  HERCOM computes a Gauss-Hermite rule;'
      write ( *, * ) '  SUMMER carries it out.'
      write ( *, * ) ' '
      write ( *, * )                                                    &
        '  The integration interval is ( -Infinity, +Infinity ).'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call hercom ( norder, xtab, weight )

          call summer ( func, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test16
      subroutine test17
!
!***********************************************************************
!
!! TEST17 tests LEGCOM and SUMSUB.
!
      integer norder
      parameter ( norder = 2 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      real(kind=Rkind) fx2sd1
      integer iexp
      integer nsub
      real(kind=Rkind) result
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      external fx2sd1
!
      a = - 1.0D+00
      b =   1.0D+00

      write ( *, * ) ' '
      write ( *, * ) 'TEST17'
      write ( *, * ) '  LEGCOM computes a Gauss-Legendre rule;'
      write ( *, * ) '  SUMSUB carries it out over subintervals.'
      write ( *, * ) ' '
      write ( *, * ) '  The integration interval is ', a, b
      write ( *, * ) '  Here, we use a fixed order NORDER = ', norder
      write ( *, * ) '  and use more and more subintervals.'
      write ( *, * ) ' '
      write ( *, * ) '  NSUB     Integral'
      write ( *, * ) ' '

      call legcom ( norder, xtab, weight )

      do iexp = 0, 9

        nsub = 2**iexp

        call sumsub ( fx2sd1, a, b, nsub, norder, xtab, weight, result )

        write ( *, '(i4,g16.8)' ) nsub, result

      end do

      return
      end subroutine test17
      subroutine test18
!
!***********************************************************************
!
!! TEST18 tests LEGSET, KRONSET and SUMMERGK.
!
      integer norderg
      integer norderk

      parameter ( norderg = 10 )
      parameter ( norderk = 2 * norderg + 1 )
!
      real(kind=Rkind) fx2sd1
      real(kind=Rkind) resultg
      real(kind=Rkind) resultk
      real(kind=Rkind) weightg(norderg)
      real(kind=Rkind) weightk(norderk)
      real(kind=Rkind) xtabg(norderg)
      real(kind=Rkind) xtabk(norderk)

      external fx2sd1
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST18'
      write ( *, * ) '  KRONSET sets up Kronrod quadrature;'
      write ( *, * ) '  LEGSET sets up Gauss-Legendre quadrature;'
      write ( *, * ) '  SUMMERGK carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is [-1, 1].'
      write ( *, * ) '  Integrand is X**2 / SQRT ( 1.1 - X**2 ).'
      write ( *, * ) ' '

      call legset ( norderg, xtabg, weightg )

      call kronset ( norderk, xtabk, weightk )

      call summergk ( fx2sd1, norderg, weightg, resultg,                &
        norderk, xtabk, weightk, resultk )

      write ( *, '(i2,2x,g16.8)' ) norderg, resultg
      write ( *, '(i2,2x,g16.8)' ) norderk, resultk
      write ( *, '(2x,2x,g16.8)' )          resultg-resultk

      return
      end subroutine test18
      subroutine test19
!
!***********************************************************************
!
!! TEST19 tests LEGSET, KRONSET and SUMSUBGK.
!
      integer norderg
      integer norderk

      parameter ( norderg = 7 )
      parameter ( norderk = 2 * norderg + 1 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      real(kind=Rkind) error
      real(kind=Rkind) fx2sd1
      integer nsub
      real(kind=Rkind) resultg
      real(kind=Rkind) resultk
      real(kind=Rkind) weightg(norderg)
      real(kind=Rkind) weightk(norderk)
      real(kind=Rkind) xtabg(norderg)
      real(kind=Rkind) xtabk(norderk)

      external fx2sd1
!
      a = - 1.0D+00
      b =   1.0D+00
      nsub = 5

      write ( *, * ) ' '
      write ( *, * ) 'TEST19'
      write ( *, * ) '  KRONSET sets up Kronrod quadrature;'
      write ( *, * ) '  LEGSET sets up Gauss-Legendre quadrature;'
      write ( *, * ) '  SUMSUBGK carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Integrand is X**2 / SQRT ( 1.1 - X**2 ).'
      write ( *, * ) ' '

      call legset ( norderg, xtabg, weightg )

      call kronset ( norderk, xtabk, weightk )

      call sumsubgk ( fx2sd1, a, b, nsub, norderg, weightg, resultg,    &
        norderk, xtabk, weightk, resultk, error )

      write ( *, '(i2,2x,g16.8)' ) norderg, resultg
      write ( *, '(i2,2x,g16.8)' ) norderk, resultk
      write ( *, '(2x,2x,g16.8)' )          error

      return
      end subroutine test19
      subroutine test20
!
!***********************************************************************
!
!! TEST20 tests MOULSET and SUMMER.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 10 )
      parameter ( nfunc = 7 )
!
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST20'
      write ( *, * ) '  MOULSET sets up Adams-Moulton quadrature;'
      write ( *, * ) '  SUMMER carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is [0,1].'
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call moulset ( norder, xtab, weight )

          call summer ( func, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test20
      subroutine test21
!
!***********************************************************************
!
!! TEST21 tests CHEBTOSET and SUMMER.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 6 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      a = - 1.0D+00
      b =   1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST21'
      write ( *, * ) '  CHEBTOSET sets up Gauss-Chebyshev quadrature,'
      write ( *, * ) '  ( open, first kind );'
      write ( *, * ) '  SUMMER carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) '  The weight function is 1 / sqrt ( 1 - X**2 )'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call chebtoset ( norder, xtab, weight )

          call summer ( func, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test21
      subroutine test22
!
!***********************************************************************
!
!! TEST22 tests CHEBTCSET and SUMMER.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 6 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      a = - 1.0D+00
      b =   1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST22'
      write ( *, * ) '  CHEBTCSET sets up Gauss-Chebyshev quadrature,'
      write ( *, * ) '  ( closed, first kind );'
      write ( *, * ) '  SUMMER carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) '  The weight function is 1 / sqrt ( 1 - X**2 )'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 2, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call chebtcset ( norder, xtab, weight )

          call summer ( func, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do
      return
      end subroutine test22
      subroutine test23
!
!***********************************************************************
!
!! TEST23 tests CHEBUSET and SUMMER.
!
      integer maxorder
      integer nfunc

      parameter ( maxorder = 4 )
      parameter ( nfunc = 7 )
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      character*6 fname
      real(kind=Rkind) func
      integer i
      integer norder
      integer nsub
      real(kind=Rkind) result(nfunc)
      real(kind=Rkind) weight(maxorder)
      real(kind=Rkind) xtab(maxorder)
!
      external func
!
      a = - 1.0D+00
      b =   1.0D+00
      nsub = 1

      write ( *, * ) ' '
      write ( *, * ) 'TEST23'
      write ( *, * ) '  CHEBUSET sets up Gauss-Chebyshev quadrature;'
      write ( *, * ) '  SUMMER carries it out.'
      write ( *, * ) ' '
      write ( *, * ) '  Integration interval is ', a, b
      write ( *, * ) '  Number of subintervals is ', nsub
      write ( *, * ) '  Quadrature order will vary.'
      write ( *, * ) '  Integrand will vary.'
      write ( *, * ) '  The weight function is sqrt ( 1 - X**2 )'
      write ( *, * ) ' '
      write ( *, '(4x, 7a10 ) ' ) ( fname(i), i = 1, nfunc )
      write ( *, * ) ' '

      do norder = 1, maxorder

        do i = 1, nfunc

          call funcset ( 'SET', i )

          call chebuset ( norder, xtab, weight )

          call summer ( func, norder, xtab, weight, result(i) )

        end do

        write ( *, '(i2,2x,7f10.6)' )                                   &
          norder, ( result(i), i = 1, nfunc )

      end do

      return
      end subroutine test23
      subroutine test24
!
!***********************************************************************
!
!! TEST24 compares LEGCOM and LEGSET.
!
      integer norder
      parameter ( norder = 19 )
!
      integer i
      integer iwdifmax
      integer ixdifmax
      real(kind=Rkind) wdifmax
      real(kind=Rkind) weight1(norder)
      real(kind=Rkind) weight2(norder)
      real(kind=Rkind) xdifmax
      real(kind=Rkind) xtab1(norder)
      real(kind=Rkind) xtab2(norder)
!
      write ( *, * ) ' '
      write ( *, * ) 'TEST24'
      write ( *, * ) '  LEGCOM computes Gauss-Legendre data;'
      write ( *, * ) '  LEGSET looks up the same data.'
      write ( *, * ) ' '
      write ( *, * ) '  Compare the data for NORDER = ', norder

      call legcom ( norder, xtab1, weight1 )
      call legset ( norder, xtab2, weight2 )

      xdifmax = 0.0D+00
      ixdifmax = 0

      wdifmax = 0.0D+00
      iwdifmax = 0

      do i = 1, norder

        if ( abs ( xtab1(i) - xtab2(i) ) .gt. xdifmax ) then
          xdifmax = abs ( xtab1(i) - xtab2(i) )
          ixdifmax = i
        end if

        if ( abs ( weight1(i) - weight2(i) ) .gt. wdifmax ) then
          wdifmax = abs ( weight1(i) - weight2(i) )
          iwdifmax = i
        end if

      end do

      if ( ixdifmax .ne. 0 ) then
        write ( *, * ) ' '
        write ( *, * ) '  Maximum abscissa difference is ', xdifmax
        write ( *, * ) '  for index I = ', ixdifmax
        write ( *, '(a,g25.18)' ) 'Computed:', xtab1(ixdifmax)
        write ( *, '(a,g25.18)' ) 'Stored:  ', xtab2(ixdifmax)
      else
        write ( *, * ) ' '
        write ( *, * ) '  Computed and stored abscissas are identical.'
      end if

      if ( iwdifmax .ne. 0 ) then
        write ( *, * ) ' '
        write ( *, * ) '  Maximum weight difference is   ', wdifmax
        write ( *, * ) '  for index I = ', iwdifmax
        write ( *, '(a,g25.18)' ) 'Computed:', weight1(iwdifmax)
        write ( *, '(a,g25.18)' ) 'Stored:  ', weight2(iwdifmax)
      else
        write ( *, * ) ' '
        write ( *, * ) '  Computed and stored weights are identical.'
      end if

      return
      end subroutine test24
      function f1sd1 ( x )
!
!***********************************************************************
!
      real(kind=Rkind) f1sd1
      real(kind=Rkind) x
!
      f1sd1 = 1.0D+00 / sqrt ( 1.1D+00 - x**2 )

      return
      end function f1sd1
      function fxsd1(x)
!
!***********************************************************************
!
      real(kind=Rkind) fxsd1
      real(kind=Rkind) x
!
      fxsd1 = x / sqrt ( 1.1D+00 - x**2 )

      return
      end function fxsd1
      function fx2sd1(x)
!
!***********************************************************************
!
      real(kind=Rkind) fx2sd1
      real(kind=Rkind) x
!
      fx2sd1 = x**2 / sqrt ( 1.1D+00 - x**2 )

      return
      end function fx2sd1
      function func ( x )
!
!***********************************************************************
!
!! FUNC evaluates a function of X, as chosen by the user.
!
      real(kind=Rkind) func
      integer ifunc
      real(kind=Rkind) x
!
      call funcset ( 'GET', ifunc )

      if ( ifunc .eq. 1 ) then
        func = 1.0D+00
      else if ( ifunc .eq. 2 ) then
        func = x
      else if ( ifunc .eq. 3 ) then
        func = x**2
      else if ( ifunc .eq. 4 ) then
        func = x**3
      else if ( ifunc .eq. 5 ) then
        func = x**4
      else if ( ifunc .eq. 6 ) then
        func = sin ( x )
      else if ( ifunc .eq. 7 ) then
        func = exp ( x )
      else
        func = 0.0D+00
      end if

      return
      end function func
      subroutine funcset ( action, i )
!
!***********************************************************************
!
!! FUNCSET sets the function to be returned by FUNC.
!
!
!  Parameters:
!
!    Input, character*3 ACTION, the action to be carried out.
!    'SET' means the call is made to set the function.
!    'GET' means the call is made to query the function.
!
!    Input/output, integer I.
!    If the input value of ACTION is 'SET', then I is also an input
!    quantity, namely, the index of the function to be chosen.
!    If the input value of ACTION is 'GET', then I is an output
!    quantity, equal to the index of the function that was last chosen.
!
      character*3 action
      integer i
      integer ival
!
      save ival
!
      data ival / 0 /
!
      if ( action .eq. 'SET' ) then
        ival = i
      else if ( action .eq. 'GET' ) then
        i = ival
      end if

      return
      end subroutine funcset
      function fname ( ifunc )
!
!***********************************************************************
!
!! FNAME returns the name of the function that will be evaluated in FUNC.
!
      character*6 fname
      integer ifunc
!
      if ( ifunc .eq. 1 ) then
        fname = '     1'
      else if ( ifunc .eq. 2 ) then
        fname = '     X'
      else if ( ifunc .eq. 3 ) then
        fname = '  X**2'
      else if ( ifunc .eq. 4 ) then
        fname = '  X**3'
      else if ( ifunc .eq. 5 ) then
        fname = '  X**4'
      else if ( ifunc .eq. 6 ) then
        fname = 'SIN(X)'
      else if ( ifunc .eq. 7 ) then
        fname = 'EXP(X)'
      else
        fname = '??????'
      end if

      return
      end function fname

