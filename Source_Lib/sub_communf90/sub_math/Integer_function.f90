
!================================================================
!    greatest common divisor
!================================================================
    FUNCTION gcd(a, b)
    integer :: gcd
    integer :: a,b

    integer :: aa,bb,t

    aa = a
    bb = b
    DO
      IF (bb == 0) EXIT
      t = bb
      bb = mod(aa,bb)
      aa = t
    END DO
    gcd = aa

    END FUNCTION gcd
!================================================================
!    fraction simplification a/b
!================================================================
    SUBROUTINE frac_simplification(a, b)
    integer :: gcd
    integer :: a,b

    integer :: aa,bb,t

    aa = a
    bb = b
    t = gcd(aa,bb)

    a = aa/t
    b = bb/t

    IF (b < 0) THEN
      b = -b
      a = -a
    END IF

    END SUBROUTINE frac_simplification
