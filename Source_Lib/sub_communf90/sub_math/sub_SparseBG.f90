function choose ( n, k )

!*****************************************************************************80
!
!! CHOOSE computes the binomial coefficient C(N,K).
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) N, K, are the values of N and K.
!
!    Output, integer ( kind=Ikind ) CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) choose
  integer ( kind=Ikind ) i
  integer ( kind=Ikind ) k
  integer ( kind=Ikind ) mn
  integer ( kind=Ikind ) mx
  integer ( kind=Ikind ) n
  integer ( kind=Ikind ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  choose = value

  return
end function choose
subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed 
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind=Ikind ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind=Ikind ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer ( kind=Ikind )  H, T, two internal parameters needed for the
!    computation.  The user should allocate space for these in the calling
!    program, include them in the calling sequence, but never alter them!
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) k

  integer ( kind=Ikind ) a(k)
  integer ( kind=Ikind ) h
  logical more
  integer ( kind=Ikind ) n
  integer ( kind=Ikind ) t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else

    if ( 1 < t ) then
      h = 0
    end if

    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end subroutine comp_next
subroutine dtable_close_write ( output_unit )

!*****************************************************************************80
!
!! DTABLE_CLOSE_WRITE closes a file used to write a DTABLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the output unit that was used.
!
  USE mod_system
  implicit none

  integer output_unit

  close ( unit = output_unit )

  return
end subroutine dtable_close_write
subroutine dtable_data_write ( output_unit, m, n, table )

!*****************************************************************************80
!
!! DTABLE_DATA_WRITE writes DTABLE data to a file.
!
!  Discussion:
!
!    This routine writes a single line of output for each point,
!    containing its spatial coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) OUTPUT_UNIT, the output unit.
!
!    Input, integer ( kind=Ikind ) M, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) N, the number of points.
!
!    Input, real ( kind=Rkind ) TABLE(M,N), the table data.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) m
  integer ( kind=Ikind ) n

  integer ( kind=Ikind ) output_unit
  integer ( kind=Ikind ) j
  character ( len = 40 ) string
  real     ( kind=Rkind ) table(m,n)
!
!  Create the format string.
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
  call s_blank_delete ( string )

  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do

  return
end subroutine dtable_data_write
subroutine dtable_header_write ( output_file_name, output_unit, m, n )

!*****************************************************************************80
!
!! DTABLE_HEADER_WRITE writes the header to a DTABLE file.
!
!  Discussion:
!
!    The file must already be open before this routine is called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer ( kind=Ikind ) OUTPUT_UNIT, the output unit.
!
!    Input, integer ( kind=Ikind ) M, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) N, the number of points.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) m
  integer ( kind=Ikind ) n
  character ( len = * ) output_file_name
  integer ( kind=Ikind ) output_unit
  character ( len = 40 ) string
  real    ( kind=Rkind ), parameter :: x = ONE

  call timestring ( string )

  write ( output_unit, '(a)'       ) '#  ' // trim ( output_file_name )
  write ( output_unit, '(a)'       ) '#  created by TABLE_IO.F90'
  write ( output_unit, '(a)'       ) '#  at ' // trim ( string )
  write ( output_unit, '(a)'       ) '#'
  write ( output_unit, '(a,i8)'    ) '#  Spatial dimension M = ', m
  write ( output_unit, '(a,i8)'    ) '#  Number of points N = ', n
  write ( output_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( x )
  write ( output_unit, '(a)'       ) '#'

  return
end subroutine dtable_header_write
subroutine dtable_open_write ( output_file_name, output_unit )

!*****************************************************************************80
!
!! DTABLE_OPEN_WRITE opens a file to write a DTABLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Output, integer OUTPUT_UNIT, the output unit to be used.
!
  USE mod_system
  implicit none

  character ( len = * ) output_file_name
  integer ( kind=Ikind ) output_status
  integer ( kind=Ikind ) output_unit

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_OPEN_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if

  return
end subroutine dtable_open_write
subroutine dtable_write ( output_file_name, m, n, table, header )

!*****************************************************************************80
!
!! DTABLE_WRITE writes a DTABLE to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer ( kind=Ikind ) M, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) N, the number of points.
!
!    Input, real ( kind=Rkind ) TABLE(M,N), the table data.
!
!    Input, logical HEADER, is TRUE if the header is to be included.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) m
  integer ( kind=Ikind ) n

  logical header
  character ( len = * ) output_file_name
  integer ( kind=Ikind ) output_unit
  real    ( kind=Rkind ) table(m,n)

  call dtable_open_write ( output_file_name, output_unit )

  if ( header ) then
    call dtable_header_write ( output_file_name, output_unit, m, n )
  end if

  call dtable_data_write ( output_unit, m, n, table )

  call dtable_close_write ( output_unit )

  return
end subroutine dtable_write
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind=Ikind ) IUNIT, the free unit number.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) i
  integer ( kind=Ikind ) ios
  integer ( kind=Ikind ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end subroutine get_unit
subroutine hermite_abscissa ( dim_num, point_num, grid_index, grid_base, &
  grid_point )

!*****************************************************************************80
!
!! HERMITE_ABSCISSA sets abscissas for multidimensional Gauss-Hermite quadrature.
!
!  Discussion:
!
!    The "nesting" as it occurs for Gauss-Hermite sparse grids simply
!    involves the use of a specified set of permissible orders for the
!    rule.  
!
!    The X array lists the (complete) Gauss-Hermite abscissas for rules 
!    of order 1, 3, 7, 15, 31, 63, and 127 in order. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) POINT_NUM, the number of points.
!
!    Input, integer ( kind=Ikind ) GRID_INDEX(DIM_NUM,POINT_NUM), for each
!    point and dimension, the index of the abscissa.
!
!    Input, integer ( kind=Ikind ) GRID_BASE(DIM_NUM), the "base" of the
!    rule being used in each dimension.
!
!    Output, real ( kind=Rkind ) GRID_POINT(DIM_NUM), the grid points of
!    abscissas.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  integer ( kind=Ikind ) point_num

  integer ( kind=Ikind ) dim
  integer ( kind=Ikind ) grid_base(dim_num)
  integer ( kind=Ikind ) grid_index(dim_num,point_num)
  real    ( kind=Rkind ) grid_point(dim_num,point_num)
  integer ( kind=Ikind ) i4_log_2
  integer ( kind=Ikind ) level
  integer ( kind=Ikind ) point
  integer ( kind=Ikind ) pointer
  integer ( kind=Ikind ), dimension ( 0:7 ) :: skip = [0, 1, 4, 11, 26, 57, 120, 247]
  real    ( kind=Rkind ), dimension ( 247 ) :: x = [&
    ZERO, &
   -0.122474487139158904909864203735_Rkind * TEN, &
    ZERO, &
    0.122474487139158904909864203735_Rkind * TEN, &
   -0.265196135683523349244708200652_Rkind * TEN, &
   -0.167355162876747144503180139830_Rkind * TEN, &
   -0.816287882858964663038710959027_Rkind, &
    ZERO, &
    0.816287882858964663038710959027_Rkind, &
    0.167355162876747144503180139830_Rkind * TEN, &
    0.265196135683523349244708200652_Rkind * TEN, &
   -0.449999070730939155366438053053_Rkind * TEN, &
   -0.366995037340445253472922383312_Rkind * TEN, &
   -0.296716692790560324848896036355_Rkind * TEN, &
   -0.232573248617385774545404479449_Rkind * TEN, &
   -0.171999257518648893241583152515_Rkind * TEN, &
   -0.113611558521092066631913490556_Rkind * TEN, &
   -0.565069583255575748526020337198_Rkind, &
    ZERO, &
    0.565069583255575748526020337198_Rkind, &
    0.113611558521092066631913490556_Rkind * TEN, &
    0.171999257518648893241583152515_Rkind * TEN, &
    0.232573248617385774545404479449_Rkind * TEN, &
    0.296716692790560324848896036355_Rkind * TEN, &
    0.366995037340445253472922383312_Rkind * TEN, &
    0.449999070730939155366438053053_Rkind * TEN, &
   -6.9956801237185402753248521473232_Rkind, &
   -6.2750787049428601427036567812530_Rkind, &
   -5.6739614446185883296332558789276_Rkind, &
   -5.1335955771123807045862968913996_Rkind, &
   -4.6315595063128599420667997654336_Rkind, &
   -4.1562717558181451724831352315314_Rkind, &
   -3.7007434032314694224497164589673_Rkind, &
   -3.2603207323135408104645401509648_Rkind, &
   -2.8316804533902054557015640151425_Rkind, &
   -2.4123177054804201051740184582119_Rkind, &
   -2.0002585489356389657975562598571_Rkind, &
   -1.5938858604721398261388419455550_Rkind, &
   -1.1918269983500464260821358649242_Rkind, &
   -0.79287697691530893968593032998830_Rkind, &
   -0.39594273647142311094670041663436_Rkind, &
    0.0000000000000000000000000000000_Rkind, &
    0.39594273647142311094670041663436_Rkind, &
    0.79287697691530893968593032998830_Rkind, &
    1.1918269983500464260821358649242_Rkind, &
    1.5938858604721398261388419455550_Rkind, &
    2.0002585489356389657975562598571_Rkind, &
    2.4123177054804201051740184582119_Rkind, &
    2.8316804533902054557015640151425_Rkind, &
    3.2603207323135408104645401509648_Rkind, &
    3.7007434032314694224497164589673_Rkind, &
    4.1562717558181451724831352315314_Rkind, &
    4.6315595063128599420667997654336_Rkind, &
    5.1335955771123807045862968913996_Rkind, &
    5.6739614446185883296332558789276_Rkind, &
    6.2750787049428601427036567812530_Rkind, &
    6.9956801237185402753248521473232_Rkind, &
   -10.435499877854168053468115427285_Rkind, &
   -9.8028759912974963635223935286507_Rkind, &
   -9.2792019543050391319404745506496_Rkind, &
   -8.8118581437284546442526628275570_Rkind, &
   -8.3807683451863219343010651043788_Rkind, &
   -7.9755950801420373181541806298501_Rkind, &
   -7.5901395198641066762479783194468_Rkind, &
   -7.2203167078889678461161324222529_Rkind, &
   -6.8632544331795368527353285876066_Rkind, &
   -6.5168348106821160605273395854042_Rkind, &
   -6.1794379922705969862418461787263_Rkind, &
   -5.8497884000810673462526582961482_Rkind, &
   -5.5268572526403031425047575122840_Rkind, &
   -5.2097979830408354861575136416263_Rkind, &
   -4.8979018644975742350745099214868_Rkind, &
   -4.5905665744435190229271294569091_Rkind, &
   -4.2872733352824404031727616199454_Rkind, &
   -3.9875699104197157485227052068068_Rkind, &
   -3.6910577000963465117322810559754_Rkind, &
   -3.3973817713303911852755941806287_Rkind, &
   -3.1062230279282566329138616746036_Rkind, &
   -2.8172919672837977750747135657355_Rkind, &
   -2.5303236304712010926855221718499_Rkind, &
   -2.2450734604812066298995918179330_Rkind, &
   -1.9613138583081485293922008411321_Rkind, &
   -1.6788312791720137520802800622638_Rkind, &
   -1.3974237486049625107570752063702_Rkind, &
   -1.1168987050996462690510970277840_Rkind, &
   -0.83707109558947615977737795461293_Rkind, &
   -0.55776166427908221668763665253822_Rkind, &
   -0.27879538567115223986687628627202_Rkind, &
    0.00000000000000000000000000000000_Rkind, &
    0.27879538567115223986687628627202_Rkind, &
    0.55776166427908221668763665253822_Rkind, &
    0.83707109558947615977737795461293_Rkind, &
    1.1168987050996462690510970277840_Rkind, &
    1.3974237486049625107570752063702_Rkind, &
    1.6788312791720137520802800622638_Rkind, &
    1.9613138583081485293922008411321_Rkind, &
    2.2450734604812066298995918179330_Rkind, &
    2.5303236304712010926855221718499_Rkind, &
    2.8172919672837977750747135657355_Rkind, &
    3.1062230279282566329138616746036_Rkind, &
    3.3973817713303911852755941806287_Rkind, &
    3.6910577000963465117322810559754_Rkind, &
    3.9875699104197157485227052068068_Rkind, &
    4.2872733352824404031727616199454_Rkind, &
    4.5905665744435190229271294569091_Rkind, &
    4.8979018644975742350745099214868_Rkind, &
    5.2097979830408354861575136416263_Rkind, &
    5.5268572526403031425047575122840_Rkind, &
    5.8497884000810673462526582961482_Rkind, &
    6.1794379922705969862418461787263_Rkind, &
    6.5168348106821160605273395854042_Rkind, &
    6.8632544331795368527353285876066_Rkind, &
    7.2203167078889678461161324222529_Rkind, &
    7.5901395198641066762479783194468_Rkind, &
    7.9755950801420373181541806298501_Rkind, &
    8.3807683451863219343010651043788_Rkind, &
    8.8118581437284546442526628275570_Rkind, &
    9.2792019543050391319404745506496_Rkind, &
    9.8028759912974963635223935286507_Rkind, &
    10.435499877854168053468115427285_Rkind, &
   -15.228338148167350978246954433464_Rkind, &
   -14.669595158833972632746354112896_Rkind, &
   -14.209085995284870755168244250887_Rkind, &
   -13.799722290211676634645246746673_Rkind, &
   -13.423518590070950062438258321855_Rkind, &
   -13.071208660474601901583995439649_Rkind, &
   -12.737235652415686338138003924072_Rkind, &
   -12.417939378869715805445879624069_Rkind, &
   -12.110749020947747600132123508132_Rkind, &
   -11.813772198267727195134584136191_Rkind, &
   -11.525565112572696599167888588564_Rkind, &
   -11.244994583785543445194384194300_Rkind, &
   -10.971150569840247423423040263881_Rkind, &
   -10.703288201027481347670940744690_Rkind, &
   -10.440787957772772867742591798027_Rkind, &
   -10.183127473450343888624126450357_Rkind, &
   -9.9298610495114250736847004273684_Rkind, &
   -9.6806044412474728038150712732737_Rkind, &
   -9.4350233389881650135019598506287_Rkind, &
   -9.1928244988460305715774195052527_Rkind, &
   -8.9537488108565404323807890169970_Rkind, &
   -8.7175658087076307363833999548548_Rkind, &
   -8.4840692689832473326097180339984_Rkind, &
   -8.2530736454457156579694124243888_Rkind, &
   -8.0244111514703375578594739796798_Rkind, &
   -7.7979293513870105420829120455591_Rkind, &
   -7.5734891556083454022834960763301_Rkind, &
   -7.3509631392269052701961258043733_Rkind, &
   -7.1302341220350710668064025713431_Rkind, &
   -6.9111939615465713197465633109366_Rkind, &
   -6.6937425208758294190074417381666_Rkind, &
   -6.4777867811645365448144903821487_Rkind, &
   -6.2632400742737354345609723857092_Rkind, &
   -6.0500214161419845694465474482388_Rkind, &
   -5.8380549248774187386601690807757_Rkind, &
   -5.6272693105464816659423455794909_Rkind, &
   -5.4175974259243240722848425872924_Rkind, &
   -5.2089758693153983587570258372239_Rkind, &
   -5.0013446320386360038520809107373_Rkind, &
   -4.7946467843764925009748509930857_Rkind, &
   -4.5888281947698372951606485031212_Rkind, &
   -4.3838372778464736294253744407459_Rkind, &
   -4.1796247675352031349421189892408_Rkind, &
   -3.9761435120673355916035814195920_Rkind, &
   -3.7733482881250526721004678400057_Rkind, &
   -3.5711956317782180447199756485249_Rkind, &
   -3.3696436841717397896643629240035_Rkind, &
   -3.1686520501953630191857798261495_Rkind, &
   -2.9681816685955910267761649521505_Rkind, &
   -2.7681946921824058801226545958892_Rkind, &
   -2.5686543769473501723144013022363_Rkind, &
   -2.3695249790490401080012474645702_Rkind, &
   -2.1707716587411506879498498083695_Rkind, &
   -1.9723603904195020079324743227565_Rkind, &
   -1.7742578780516791584676442103681_Rkind, &
   -1.5764314753267801315519597621879_Rkind, &
   -1.3788491099261778091441557053728_Rkind, &
   -1.1814792113700685848678583598423_Rkind, &
   -0.98429064194027277726568984213773_Rkind, &
   -0.78725263021825034151596831878971_Rkind, &
   -0.59033470680942102142230439346102_Rkind, &
   -0.39350664185130136568037826200185_Rkind, &
   -0.19673838392423251964272239737078_Rkind, &
    0.0000000000000000000000000000000_Rkind, &
    0.19673838392423251964272239737078_Rkind, &
    0.39350664185130136568037826200185_Rkind, &
    0.59033470680942102142230439346102_Rkind, &
    0.78725263021825034151596831878971_Rkind, &
    0.98429064194027277726568984213773_Rkind, &
    1.1814792113700685848678583598423_Rkind, &
    1.3788491099261778091441557053728_Rkind, &
    1.5764314753267801315519597621879_Rkind, &
    1.7742578780516791584676442103681_Rkind, &
    1.9723603904195020079324743227565_Rkind, &
    2.1707716587411506879498498083695_Rkind, &
    2.3695249790490401080012474645702_Rkind, &
    2.5686543769473501723144013022363_Rkind, &
    2.7681946921824058801226545958892_Rkind, &
    2.9681816685955910267761649521505_Rkind, &
    3.1686520501953630191857798261495_Rkind, &
    3.3696436841717397896643629240035_Rkind, &
    3.5711956317782180447199756485249_Rkind, &
    3.7733482881250526721004678400057_Rkind, &
    3.9761435120673355916035814195920_Rkind, &
    4.1796247675352031349421189892408_Rkind, &
    4.3838372778464736294253744407459_Rkind, &
    4.5888281947698372951606485031212_Rkind, &
    4.7946467843764925009748509930857_Rkind, &
    5.0013446320386360038520809107373_Rkind, &
    5.2089758693153983587570258372239_Rkind, &
    5.4175974259243240722848425872924_Rkind, &
    5.6272693105464816659423455794909_Rkind, &
    5.8380549248774187386601690807757_Rkind, &
    6.0500214161419845694465474482388_Rkind, &
    6.2632400742737354345609723857092_Rkind, &
    6.4777867811645365448144903821487_Rkind, &
    6.6937425208758294190074417381666_Rkind, &
    6.9111939615465713197465633109366_Rkind, &
    7.1302341220350710668064025713431_Rkind, &
    7.3509631392269052701961258043733_Rkind, &
    7.5734891556083454022834960763301_Rkind, &
    7.7979293513870105420829120455591_Rkind, &
    8.0244111514703375578594739796798_Rkind, &
    8.2530736454457156579694124243888_Rkind, &
    8.4840692689832473326097180339984_Rkind, &
    8.7175658087076307363833999548548_Rkind, &
    8.9537488108565404323807890169970_Rkind, &
    9.1928244988460305715774195052527_Rkind, &
    9.4350233389881650135019598506287_Rkind, &
    9.6806044412474728038150712732737_Rkind, &
    9.9298610495114250736847004273684_Rkind, &
    10.183127473450343888624126450357_Rkind, &
    10.440787957772772867742591798027_Rkind, &
    10.703288201027481347670940744690_Rkind, &
    10.971150569840247423423040263881_Rkind, &
    11.244994583785543445194384194300_Rkind, &
    11.525565112572696599167888588564_Rkind, &
    11.813772198267727195134584136191_Rkind, &
    12.110749020947747600132123508132_Rkind, &
    12.417939378869715805445879624069_Rkind, &
    12.737235652415686338138003924072_Rkind, &
    13.071208660474601901583995439649_Rkind, &
    13.423518590070950062438258321855_Rkind, &
    13.799722290211676634645246746673_Rkind, &
    14.209085995284870755168244250887_Rkind, &
    14.669595158833972632746354112896_Rkind, &
    15.228338148167350978246954433464_Rkind  &
        ]

  if ( any ( grid_base(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are less than 0.'
    stop
  end if

  if ( any ( 63 < grid_base(1:dim_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are greater than 63.'
    stop
  end if

  do point = 1, point_num
    do dim = 1, dim_num

      level = i4_log_2 ( grid_base(dim) + 1 )

      pointer = skip(level) + ( grid_index(dim,point) + grid_base(dim) + 1 )

      grid_point(dim,point) = x(pointer) 

    end do
  end do

  return
end subroutine hermite_abscissa
subroutine hermite_weights ( order, weight )

!*****************************************************************************80
!
!! HERMITE_WEIGHTS returns weights for certain Gauss-Hermite quadrature rules.
!
!  Discussion:
!
!    The allowed orders are 1, 3, 7, 15, 31, 63, and 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) ORDER, the order of the rule.
!    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
!
!    Output, real ( kind=Rkind ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric and should sum to SQRT(PI).
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) order

  real    ( kind=Rkind ) weight(order)

  if ( order == 1 ) then

    weight(1) = 1.77245385090551602729816748334_Rkind

  else if ( order == 3 ) then

    weight(1) = 0.295408975150919337883027913890_Rkind
    weight(2) = 0.118163590060367735153211165556_Rkind * TEN
    weight(3) = 0.295408975150919337883027913890_Rkind

  else if ( order == 7 ) then

    weight(1) = 0.971781245099519154149424255939_Rkind / TEN**3
    weight(2) = 0.545155828191270305921785688417_Rkind / TEN
    weight(3) = 0.425607252610127800520317466666_Rkind
    weight(4) = 0.810264617556807326764876563813_Rkind
    weight(5) = 0.425607252610127800520317466666_Rkind
    weight(6) = 0.545155828191270305921785688417_Rkind / TEN
    weight(7) = 0.971781245099519154149424255939_Rkind / TEN**3

  else if ( order == 15 ) then

    weight(1) =  0.152247580425351702016062666965_Rkind / TEN**8
    weight(2) =  0.105911554771106663577520791055_Rkind / TEN**5
    weight(3) =  0.100004441232499868127296736177_Rkind / TEN**3
    weight(4) =  0.277806884291277589607887049229_Rkind / TEN**2
    weight(5) =  0.307800338725460822286814158758_Rkind / TEN
    weight(6) =  0.158488915795935746883839384960_Rkind
    weight(7) =  0.412028687498898627025891079568_Rkind
    weight(8) =  0.564100308726417532852625797340_Rkind
    weight(9) =  0.412028687498898627025891079568_Rkind
    weight(10) = 0.158488915795935746883839384960_Rkind
    weight(11) = 0.307800338725460822286814158758_Rkind / TEN
    weight(12) = 0.277806884291277589607887049229_Rkind / TEN**2
    weight(13) = 0.100004441232499868127296736177_Rkind / TEN**3
    weight(14) = 0.105911554771106663577520791055_Rkind / TEN**5
    weight(15) = 0.152247580425351702016062666965_Rkind / TEN**8

  else if ( order == 31 ) then
  
    weight(  1) =   0.46189683944498305857470556847735_Rkind / TEN**21
    weight(  2) =   0.51106090079112519643027197715274_Rkind / TEN**17
    weight(  3) =   0.58995564987355133075257722133966_Rkind / TEN**14
    weight(  4) =   0.18603735214463569590294465062239_Rkind / TEN**11
    weight(  5) =   0.23524920032013205739850619940094_Rkind / TEN**9
    weight(  6) =   0.14611988344865057576066495091513_Rkind / TEN**7
    weight(  7) =   0.50437125589241034841778074689627_Rkind / TEN**6
    weight(  8) =   0.10498602757642934202945441341697_Rkind / TEN**4
    weight(  9) =   0.13952090395003623854995664958146_Rkind / TEN**3
    weight( 10) =   0.12336833073030489880608311394968_Rkind / TEN**2
    weight( 11) =   0.74827999140119116765002499116934_Rkind / TEN**2
    weight( 12) =   0.31847230731201222775249585776902_Rkind / TEN
    weight( 13) =   0.96717948160569462991143316029341_Rkind / TEN
    weight( 14) =   0.21213278866810461318136114862419_Rkind
    weight( 15) =   0.33877265789305344906000174083214_Rkind
    weight( 16) =   0.39577855609737786462923720809676_Rkind
    weight( 17) =   0.33877265789305344906000174083214_Rkind
    weight( 18) =   0.21213278866810461318136114862419_Rkind
    weight( 19) =   0.96717948160569462991143316029341_Rkind / TEN
    weight( 20) =   0.31847230731201222775249585776902_Rkind / TEN
    weight( 21) =   0.74827999140119116765002499116934_Rkind / TEN**2
    weight( 22) =   0.12336833073030489880608311394968_Rkind / TEN**2
    weight( 23) =   0.13952090395003623854995664958146_Rkind / TEN**3
    weight( 24) =   0.10498602757642934202945441341697_Rkind / TEN**4
    weight( 25) =   0.50437125589241034841778074689627_Rkind / TEN**06
    weight( 26) =   0.14611988344865057576066495091513_Rkind / TEN**7
    weight( 27) =   0.23524920032013205739850619940094_Rkind / TEN**9
    weight( 28) =   0.18603735214463569590294465062239_Rkind / TEN**11
    weight( 29) =   0.58995564987355133075257722133966_Rkind / TEN**14
    weight( 30) =   0.51106090079112519643027197715274_Rkind / TEN**17
    weight( 31) =   0.46189683944498305857470556847735_Rkind / TEN**21

  else if ( order == 63 ) then

    weight(  1) =   0.37099206434787551197827130470031_Rkind / TEN**47
    weight(  2) =   0.10400778615192299534481914814892_Rkind / TEN**41
    weight(  3) =   0.19796804708258311251124226474396_Rkind / TEN**37
    weight(  4) =   0.84687478191640015120141181138947_Rkind / TEN**34
    weight(  5) =   0.13071305930779945903630127634063_Rkind / TEN**30
    weight(  6) =   0.93437837175367456929765381518998_Rkind / TEN**28
    weight(  7) =   0.36027426635173044862245783257252_Rkind / TEN**25
    weight(  8) =   0.82963863115951789374753323156164_Rkind / TEN**23
    weight(  9) =   0.12266629909105281472971700203949_Rkind / TEN**20
    weight( 10) =   0.12288435628797061539461585325494_Rkind / TEN**18
    weight( 11) =   0.86925536958188009075932426691516_Rkind / TEN**17
    weight( 12) =   0.44857058689176221240330804981619_Rkind / TEN**15
    weight( 13) =   0.17335817955735154599902643794700_Rkind / TEN**13
    weight( 14) =   0.51265062385038307838565047455223_Rkind / TEN**12
    weight( 15) =   0.11808921844532942490513037158404_Rkind / TEN**10
    weight( 16) =   0.21508698297808025739828859845140_Rkind / TEN**9
    weight( 17) =   0.31371929535285447801497640621672_Rkind / TEN**8
    weight( 18) =   0.37041625984781705796752840204084_Rkind / TEN**7
    weight( 19) =   0.35734732949879669663960738150956_Rkind / TEN**6
    weight( 20) =   0.28393114498380927832990899215541_Rkind / TEN**5
    weight( 21) =   0.18709113003730498008961134765721_Rkind / TEN**4
    weight( 22) =   0.10284880800653635546698378640623_Rkind / TEN**3
    weight( 23) =   0.47411702610173128107201781718693_Rkind / TEN**3
    weight( 24) =   0.18409222622384813438539657470055_Rkind / TEN**2
    weight( 25) =   0.60436044551187631655712178246467_Rkind / TEN**2
    weight( 26) =   0.16829299199599730926458559757600_Rkind / TEN
    weight( 27) =   0.39858264027692992170237391875317_Rkind / TEN
    weight( 28) =   0.80467087993950415219587554532823_Rkind / TEN
    weight( 29) =   0.13871950817615293377792092082674_Rkind
    weight( 30) =   0.20448695346833761570957197160475_Rkind
    weight( 31) =   0.25799889943058042204920467417642_Rkind
    weight( 32) =   0.27876694884838411919175686949858_Rkind
    weight( 33) =   0.25799889943058042204920467417642_Rkind
    weight( 34) =   0.20448695346833761570957197160475_Rkind
    weight( 35) =   0.13871950817615293377792092082674_Rkind
    weight( 36) =   0.80467087993950415219587554532823_Rkind / TEN
    weight( 37) =   0.39858264027692992170237391875317_Rkind / TEN
    weight( 38) =   0.16829299199599730926458559757600_Rkind / TEN
    weight( 39) =   0.60436044551187631655712178246467_Rkind / TEN**2
    weight( 40) =   0.18409222622384813438539657470055_Rkind / TEN**2
    weight( 41) =   0.47411702610173128107201781718693_Rkind / TEN**3
    weight( 42) =   0.10284880800653635546698378640623_Rkind / TEN**3
    weight( 43) =   0.18709113003730498008961134765721_Rkind / TEN**4
    weight( 44) =   0.28393114498380927832990899215541_Rkind / TEN**5
    weight( 45) =   0.35734732949879669663960738150956_Rkind / TEN**6
    weight( 46) =   0.37041625984781705796752840204084_Rkind / TEN**7
    weight( 47) =   0.31371929535285447801497640621672_Rkind / TEN**8
    weight( 48) =   0.21508698297808025739828859845140_Rkind / TEN**9
    weight( 49) =   0.11808921844532942490513037158404_Rkind / TEN**10
    weight( 50) =   0.51265062385038307838565047455223_Rkind / TEN**12
    weight( 51) =   0.17335817955735154599902643794700_Rkind / TEN**13
    weight( 52) =   0.44857058689176221240330804981619_Rkind / TEN**15
    weight( 53) =   0.86925536958188009075932426691516_Rkind / TEN**17
    weight( 54) =   0.12288435628797061539461585325494_Rkind / TEN**18
    weight( 55) =   0.12266629909105281472971700203949_Rkind / TEN**20
    weight( 56) =   0.82963863115951789374753323156164_Rkind / TEN**23
    weight( 57) =   0.36027426635173044862245783257252_Rkind / TEN**25
    weight( 58) =   0.93437837175367456929765381518998_Rkind / TEN**28
    weight( 59) =   0.13071305930779945903630127634063_Rkind / TEN**30
    weight( 60) =   0.84687478191640015120141181138947_Rkind / TEN**34
    weight( 61) =   0.19796804708258311251124226474396_Rkind / TEN**37
    weight( 62) =   0.10400778615192299534481914814892_Rkind / TEN**41
    weight( 63) =   0.37099206434787551197827130470031_Rkind / TEN**47
 
  else if ( order == 127 ) then    

    weight(  1) =   0.12504497577050595552677230002883_Rkind / TEN**100
    weight(  2) =   0.17272798059419131415318615789672_Rkind / TEN**93
    weight(  3) =   0.89321681571986548608031150791499_Rkind / TEN**88
    weight(  4) =   0.77306185240893578449625186483810_Rkind / TEN**83
    weight(  5) =   0.20143957652648255497735460506196_Rkind / TEN**78
    weight(  6) =   0.21503714733610239701351039429345_Rkind / TEN**74
    weight(  7) =   0.11341924208594594813715533569504_Rkind / TEN**70
    weight(  8) =   0.33489139011795051950683388483136_Rkind / TEN**67
    weight(  9) =   0.60486548964016681064424451668405_Rkind / TEN**64
    weight( 10) =   0.71375092946352177824971347343892_Rkind / TEN**61
    weight( 11) =   0.57884563374885556636801095624030_Rkind / TEN**58
    weight( 12) =   0.33581166223858230300409326551248_Rkind / TEN**55
    weight( 13) =   0.14394641949253923568603163698953_Rkind / TEN**52
    weight( 14) =   0.46821808383216117724080263903889_Rkind / TEN**50
    weight( 15) =   0.11817054440684264071348471955361_Rkind / TEN**47
    weight( 16) =   0.23581659156008927203181682045005_Rkind / TEN**45
    weight( 17) =   0.37814427940797540210712758405540_Rkind / TEN**43
    weight( 18) =   0.49411031115771638145610738414006_Rkind / TEN**41
    weight( 19) =   0.53255303775425059266087298458297_Rkind / TEN**39
    weight( 20) =   0.47854390680131484999315199332765_Rkind / TEN**37
    weight( 21) =   0.36191883445952356128627543209554_Rkind / TEN**35
    weight( 22) =   0.23232083386343554805352497446119_Rkind / TEN**33
    weight( 23) =   0.12753331411008716683688974281454_Rkind / TEN**31
    weight( 24) =   0.60277753850758742112436095241270_Rkind / TEN**30
    weight( 25) =   0.24679773241777200207460855084439_Rkind / TEN**28
    weight( 26) =   0.88019567691698482573264198727415_Rkind / TEN**27
    weight( 27) =   0.27482489212040561315005725890593_Rkind / TEN**25
    weight( 28) =   0.75468218903085486125222816438456_Rkind / TEN**24
    weight( 29) =   0.18303134636280466270545996891835_Rkind / TEN**22
    weight( 30) =   0.39355990860860813085582448449811_Rkind / TEN**21
    weight( 31) =   0.75293161638581191068419292570042_Rkind / TEN**20
    weight( 32) =   0.12857997786722855037584105682618_Rkind / TEN**18
    weight( 33) =   0.19659326888445857792541925311450_Rkind / TEN**17
    weight( 34) =   0.26986511907214101894995783364250_Rkind / TEN**16
    weight( 35) =   0.33344414303198856330118301113874_Rkind / TEN**15
    weight( 36) =   0.37173303125150639885726463109574_Rkind / TEN**14
    weight( 37) =   0.37473954472839737091885387788983_Rkind / TEN**13
    weight( 38) =   0.34230094493397259538669512076007_Rkind / TEN**12
    weight( 39) =   0.28385303724993373166810860630552_Rkind / TEN**11
    weight( 40) =   0.21406920290454669208938772802828_Rkind / TEN**10
    weight( 41) =   0.14706331273431716244229273183839_Rkind / TEN**9
    weight( 42) =   0.92173940967434659264335883218167_Rkind / TEN**9
    weight( 43) =   0.52781663936972714041837056042506_Rkind / TEN**8
    weight( 44) =   0.27650497044951117835905283127679_Rkind / TEN**7
    weight( 45) =   0.13267855842539464770913063113371_Rkind / TEN**6
    weight( 46) =   0.58380944276113062188573331195042_Rkind / TEN**6
    weight( 47) =   0.23581561724775629112332165335800_Rkind / TEN**5
    weight( 48) =   0.87524468034280444703919485644809_Rkind / TEN**5
    weight( 49) =   0.29876790535909012274846532159647_Rkind / TEN**4
    weight( 50) =   0.93874435720072545206729594267039_Rkind / TEN**4
    weight( 51) =   0.27170762627931172053444716883938_Rkind / TEN**3
    weight( 52) =   0.72493929742498358979684249380921_Rkind / TEN**3
    weight( 53) =   0.17841208326763432884316727108264_Rkind / TEN**2
    weight( 54) =   0.40524855186046131499765636276283_Rkind / TEN**2
    weight( 55) =   0.85000263041544110385806705526917_Rkind / TEN**2
    weight( 56) =   0.16471142241609687824005585301760_Rkind / TEN
    weight( 57) =   0.29499296248213632269675010319119_Rkind / TEN
    weight( 58) =   0.48847387114300011006959603975676_Rkind / TEN
    weight( 59) =   0.74807989768583731416517226905270_Rkind / TEN
    weight( 60) =   0.10598520508090929403834368934301_Rkind
    weight( 61) =   0.13893945309051540832066283010510_Rkind
    weight( 62) =   0.16856236074207929740526975049765_Rkind
    weight( 63) =   0.18927849580120432177170145550076_Rkind
    weight( 64) =   0.19673340688823289786163676995151_Rkind
    weight( 65) =   0.18927849580120432177170145550076_Rkind
    weight( 66) =   0.16856236074207929740526975049765_Rkind
    weight( 67) =   0.13893945309051540832066283010510_Rkind
    weight( 68) =   0.10598520508090929403834368934301_Rkind
    weight( 69) =   0.74807989768583731416517226905270_Rkind / TEN
    weight( 70) =   0.48847387114300011006959603975676_Rkind / TEN
    weight( 71) =   0.29499296248213632269675010319119_Rkind / TEN
    weight( 72) =   0.16471142241609687824005585301760_Rkind / TEN
    weight( 73) =   0.85000263041544110385806705526917_Rkind / TEN**2
    weight( 74) =   0.40524855186046131499765636276283_Rkind / TEN**2
    weight( 75) =   0.17841208326763432884316727108264_Rkind / TEN**2
    weight( 76) =   0.72493929742498358979684249380921_Rkind / TEN**3
    weight( 77) =   0.27170762627931172053444716883938_Rkind / TEN**3
    weight( 78) =   0.93874435720072545206729594267039_Rkind / TEN**4
    weight( 79) =   0.29876790535909012274846532159647_Rkind / TEN**4
    weight( 80) =   0.87524468034280444703919485644809_Rkind / TEN**5
    weight( 81) =   0.23581561724775629112332165335800_Rkind / TEN**5
    weight( 82) =   0.58380944276113062188573331195042_Rkind / TEN**6
    weight( 83) =   0.13267855842539464770913063113371_Rkind / TEN**6
    weight( 84) =   0.27650497044951117835905283127679_Rkind / TEN**7
    weight( 85) =   0.52781663936972714041837056042506_Rkind / TEN**8
    weight( 86) =   0.92173940967434659264335883218167_Rkind / TEN**9
    weight( 87) =   0.14706331273431716244229273183839_Rkind / TEN**9
    weight( 88) =   0.21406920290454669208938772802828_Rkind / TEN**10
    weight( 89) =   0.28385303724993373166810860630552_Rkind / TEN**11
    weight( 90) =   0.34230094493397259538669512076007_Rkind / TEN**12
    weight( 91) =   0.37473954472839737091885387788983_Rkind / TEN**13
    weight( 92) =   0.37173303125150639885726463109574_Rkind / TEN**14
    weight( 93) =   0.33344414303198856330118301113874_Rkind / TEN**15
    weight( 94) =   0.26986511907214101894995783364250_Rkind / TEN**16
    weight( 95) =   0.19659326888445857792541925311450_Rkind / TEN**17
    weight( 96) =   0.12857997786722855037584105682618_Rkind / TEN**18
    weight( 97) =   0.75293161638581191068419292570042_Rkind / TEN**20
    weight( 98) =   0.39355990860860813085582448449811_Rkind / TEN**21
    weight( 99) =   0.18303134636280466270545996891835_Rkind / TEN**22
    weight(100) =   0.75468218903085486125222816438456_Rkind / TEN**24
    weight(101) =   0.27482489212040561315005725890593_Rkind / TEN**25
    weight(102) =   0.88019567691698482573264198727415_Rkind / TEN**27
    weight(103) =   0.24679773241777200207460855084439_Rkind / TEN**28
    weight(104) =   0.60277753850758742112436095241270_Rkind / TEN**30
    weight(105) =   0.12753331411008716683688974281454_Rkind / TEN**31
    weight(106) =   0.23232083386343554805352497446119_Rkind / TEN**33
    weight(107) =   0.36191883445952356128627543209554_Rkind / TEN**35
    weight(108) =   0.47854390680131484999315199332765_Rkind / TEN**37
    weight(109) =   0.53255303775425059266087298458297_Rkind / TEN**39
    weight(110) =   0.49411031115771638145610738414006_Rkind / TEN**41
    weight(111) =   0.37814427940797540210712758405540_Rkind / TEN**43
    weight(112) =   0.23581659156008927203181682045005_Rkind / TEN**45
    weight(113) =   0.11817054440684264071348471955361_Rkind / TEN**47
    weight(114) =   0.46821808383216117724080263903889_Rkind / TEN**50
    weight(115) =   0.14394641949253923568603163698953_Rkind / TEN**52
    weight(116) =   0.33581166223858230300409326551248_Rkind / TEN**55
    weight(117) =   0.57884563374885556636801095624030_Rkind / TEN**58
    weight(118) =   0.71375092946352177824971347343892_Rkind / TEN**61
    weight(119) =   0.60486548964016681064424451668405_Rkind / TEN**64
    weight(120) =   0.33489139011795051950683388483136_Rkind / TEN**67
    weight(121) =   0.11341924208594594813715533569504_Rkind / TEN**70
    weight(122) =   0.21503714733610239701351039429345_Rkind / TEN**74
    weight(123) =   0.20143957652648255497735460506196_Rkind / TEN**78
    weight(124) =   0.77306185240893578449625186483810_Rkind / TEN**83
    weight(125) =   0.89321681571986548608031150791499_Rkind / TEN**88
    weight(126) =   0.17272798059419131415318615789672_Rkind / TEN**93
    weight(127) =   0.12504497577050595552677230002883_Rkind / TEN**100
    
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HERMITE_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    write ( *, '(a)' ) '  Legal values are 1, 3, 7, 15, 31, 63 and 127.'
    stop

  end if

  return
end subroutine hermite_weights
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind=Ikind ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind=Ikind ) I4_HUGE, a "huge" I4.
!
  USE mod_system
  implicit none

  
  integer ( kind=Ikind ) i4_huge

  i4_huge = 2147483647

  return
end function i4_huge
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
!
!  Discussion:
!
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
!    An I4 is an integer ( kind=Ikind ) value.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) I, the number whose logarithm base 2
!    is desired.
!
!    Output, integer ( kind=Ikind ) I4_LOG_2, the integer part of the
!    logarithm base 2 of the absolute value of I.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) i
  integer ( kind=Ikind ) i_abs
  integer ( kind=Ikind ) i4_log_2
  integer ( kind=Ikind ) i4_huge

  if ( i == 0 ) then

    i4_log_2 = - i4_huge ( )

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end function i4_log_2
subroutine index_level_herm ( level, level_max, dim_num, point_num, grid_index, &
  grid_base, grid_level )

!*****************************************************************************80
!
!! INDEX_LEVEL_HERM: determine first level at which given index is generated.
!
!  Discussion:
!
!    We are constructing a sparse grid of Gauss-Hermite points.  The grid
!    is built up of product grids, with a characteristic LEVEL.  
!
!    We are concerned with identifying points in this product grid which
!    have actually been generated previously, on a lower value of LEVEL.
!
!    This routine determines the lowest value of LEVEL at which each of
!    the input points would be generated.
!
!    In 1D, given LEVEL, the number of points is ORDER = 2**(LEVEL+1) + 1,
!    (except that LEVEL = 0 implies ORDER = 1!), the BASE is (ORDER-1)/2, 
!    and the point INDEX values range from -BASE to +BASE.
!
!    The values of INDEX and BASE allow us to determine the abstract
!    properties of the point.  In particular, if INDEX is 0, the corresponding
!    Gauss-Hermite abscissa is 0, the special "nested" value we need
!    to take care of.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) LEVEL, the level at which these points were
!    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
!
!    Input, integer ( kind=Ikind ) LEVEL_MAX, the maximum level.
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) POINT_NUM, the number of points to be tested.
!
!    Input, integer ( kind=Ikind ) GRID_INDEX(DIM_NUM,POINT_NUM), the indices of
!    the points to be tested.
!
!    Input, integer ( kind=Ikind ) GRID_BASE(DIM_NUM), the "base", which is
!    essentially the denominator of the index.
!
!    Output, integer ( kind=Ikind ) GRID_LEVEL(POINT_NUM), the value of LEVEL at
!    which the point would first be generated.  This will be the same as
!    the input value of LEVEL, unless the point has an INDEX of 0 and
!    a corresponding BASE that is NOT zero.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  integer ( kind=Ikind ) point_num

  integer ( kind=Ikind ) dim
  integer ( kind=Ikind ) grid_base(dim_num)
  integer ( kind=Ikind ) grid_index(dim_num,point_num)
  integer ( kind=Ikind ) grid_level(point_num)
  integer ( kind=Ikind ) level
  integer ( kind=Ikind ) level_max
  integer ( kind=Ikind ) level_min
  integer ( kind=Ikind ) point
  
  level_min = max ( 0, level_max + 1 - dim_num )
!
!  If a point has a DIM-th component whose INDEX is 0, then the 
!  value of LEVEL at which this point would first be generated is
!  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
!
  do point = 1, point_num

    grid_level(point) = max ( level, level_min )

    do dim = 1, dim_num
      if ( grid_index(dim,point) == 0 ) then
        grid_level(point) = max ( grid_level(point) - grid_base(dim), level_min )
      end if
    end do

  end do

  return
end subroutine index_level_herm
subroutine level_to_order_open ( dim_num, level, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_OPEN converts a level to an order for open rules.
!
!  Discussion:
!
!    Sparse grids can naturally be nested.  A natural scheme is to use
!    a series of one-dimensional rules arranged in a series of "levels"
!    whose order roughly doubles with each step.
!
!    The arrangement described here works naturally for the Fejer Type 1,
!    Fejer Type 2, Gauss-Patterson, Newton Cotes Open, and
!    Newton Cotes Half Open rules.  It also can be used, partially, to describe
!    the growth of Gauss-Legendre and Gauss-Hermite rules.
!
!    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
!    point at the center, and for all values afterwards, we use the 
!    relationship
!
!      ORDER = 2**(LEVEL+1) - 1.
!
!    The following table shows how the growth will occur:
!
!    Level    Order
!
!    0          1
!    1          3 =  4 - 1
!    2          7 =  8 - 1
!    3         15 = 16 - 1
!    4         31 = 32 - 1
!    5         63 = 64 - 1
!
!    For the Fejer Type 1, Fejer Type 2, Gauss-Patterson, Newton Cotes Open, 
!    and Newton Cotes Open Half rules, the point growth is
!    nested.  If we have ORDER points on a particular LEVEL, the next level 
!    includes all these old points, plus ORDER+1 new points, formed in the 
!    gaps between successive pairs of old points plus an extra point at each 
!    end.
!
!    Level    Order = New + Old
!
!    0          1   =  1  +  0
!    1          3   =  2  +  1
!    2          7   =  4  +  3
!    3         15   =  8  +  7
!    4         31   = 16  + 15
!    5         63   = 32  + 31
!
!    If we use a series of Gauss-Legendre or Gauss-Hermite rules, then there
!    is almost no nesting, except that the central point is shared.  If we 
!    insist on producing a comparable series of such points, then the 
!    "nesting" behavior is as follows:
!
!    Level    Order = New + Old
!
!    0          1   =  1  +  0
!    1          3   =  2  +  1
!    2          7   =  6  +  1
!    3         15   = 14  +  1
!    4         31   = 30  +  1
!    5         63   = 62  +  1
!
!    Moreover, if we consider ALL the points used in such a set of "nested" 
!    Gauss-Hermite or Gauss-Legendre rules, then we must sum the "NEW" column,
!    and we see that we get roughly twice as many points as for the truly
!    nested rules.
!
!    In this routine, we assume that a vector of levels is given,
!    and the corresponding orders are desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) LEVEL(DIM_NUM), the nesting levels of the
!    1D rules.
!
!    Output, integer ( kind=Ikind ) ORDER(DIM_NUM), the order (number of points)
!    of the 1D rules.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num

  integer ( kind=Ikind ) dim
  integer ( kind=Ikind ) level(dim_num)
  integer ( kind=Ikind ) order(dim_num)

  do dim = 1, dim_num

    if ( level(dim) < 0 ) then
      order(dim) = -1
    else if ( level(dim) == 0 ) then
      order(dim) = 1
    else
      order(dim) = 2**( level(dim) + 1 ) - 1
    end if

  end do

  return
end subroutine level_to_order_open
subroutine monomial_integral_hermite ( dim_num, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_HERMITE integrates a Hermite monomial.
!
!  Discussion:
!
!    H(d,n) = Integral ( -Infinity < x < Infinity ) 
!      x1^n1 * x2^n2...*xd^nd * exp(-x1^2-x2^2...-xd^2 ) dx
!
!    H(d,n) is 0 if any n(i) odd.
!
!    H(d,n) = product ( 1 <= i <= d ) 
!      ( (n(i)-1)!! * sqrt(pi) / 2^(n(i)/2) for all n(i) even.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the dimension of the integral.
!
!    Input, integer ( kind=Ikind ) EXPON(DIM_NUM), the order of the integral.
!    0 <= EXPON(1:DIM_NUM).
!
!    Output, real ( kind=Rkind ) VALUE, the value of the integral.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  
  integer ( kind=Ikind ) dim
  integer ( kind=Ikind ), parameter :: i4_1 = 1
  integer ( kind=Ikind ), parameter :: i4_2 = 2
  integer ( kind=Ikind ) expon(dim_num)
  real    ( kind=Rkind ) r8_factorial2
  real    ( kind=Rkind ) r8_huge
  real    ( kind=Rkind ) value

  if ( any ( expon(1:dim_num) < 0 ) ) then

    value = - r8_huge ( )

  else if ( any ( mod ( expon(1:dim_num), i4_2 ) == 1 ) ) then

    value = ZERO

  else

    value = ONE
    do dim = 1, dim_num
      value = value * r8_factorial2 ( expon(dim) - i4_1 ) * sqrt ( pi ) &
        / 2.0_Rkind**( expon(dim) / 2 )
    end do
    
  end if

  return
end subroutine monomial_integral_hermite
subroutine monomial_quadrature_hermite ( dim_num, expon, point_num, weight, &
  x, quad_error )

!*****************************************************************************80
!
!! MONOMIAL_QUADRATURE_HERMITE applies a quadrature rule to a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) EXPON(DIM_NUM), the exponents.
!
!    Input, integer ( kind=Ikind ) POINT_NUM, the number of points in the rule.
!
!    Input, real ( kind=Rkind ) WEIGHT(POINT_NUM), the quadrature weights.
!
!    Input, real ( kind=Rkind ) X(DIM_NUM,POINT_NUM), the quadrature points.
!
!    Output, real ( kind=Rkind ) QUAD_ERROR, the quadrature error.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num

  real    ( kind=Rkind ) exact
  integer ( kind=Ikind ) expon(dim_num)
  integer ( kind=Ikind ) point_num
  real    ( kind=Rkind ) quad
  real    ( kind=Rkind ) quad_error
  real    ( kind=Rkind ) value(point_num)
  real    ( kind=Rkind ) weight(point_num)
  real    ( kind=Rkind ) x(dim_num,point_num)
!
!  Get the exact value of the integral of the unscaled monomial.
!
  call monomial_integral_hermite ( dim_num, expon, exact )
!
!  Evaluate the monomial at the quadrature points.
!
  call monomial_value ( dim_num, point_num, x, expon, value )
!
!  Compute the weighted sum.
!
  quad = dot_product ( weight, value )
!
!  If exact value is nonzero, use it to scale the data.
!
  if ( exact == ZERO ) then
    quad_error = abs ( quad )
  else
    quad_error = abs ( ( quad - exact ) / exact )
  end if
    
  return
end subroutine monomial_quadrature_hermite
subroutine monomial_value ( dim_num, point_num, x, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) POINT_NUM, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, real ( kind=Rkind ) X(DIM_NUM,POINT_NUM), the point coordinates.
!
!    Input, integer ( kind=Ikind ) EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind=Rkind ) VALUE(POINT_NUM), the value of the monomial.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  integer ( kind=Ikind ) point_num

  integer ( kind=Ikind ) dim
  integer ( kind=Ikind ) expon(dim_num)
  integer ( kind=Ikind ) point
  real    ( kind=Rkind ) value(point_num)
  real    ( kind=Rkind ) x(dim_num,point_num)

  value(1:point_num) = ONE

  do dim = 1, dim_num
    do point = 1, point_num
      if ( x(dim,point) /= ZERO ) then
        value(point) = value(point) * x(dim,point)**expon(dim)
      else if ( expon(dim) == 0 ) then
        value(point) = value(point)
      else
        value(point) = ZERO
      end if   
    end do
  end do

  return
end subroutine monomial_value
subroutine multigrid_index_z ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX_Z returns an indexed multidimensional grid.
!
!  Discussion:
!
!    For dimension DIM, the number of points is ORDER_1D(DIM).
!
!    We assume that ORDER_1D(DIM) is an odd number,
!      ORDER_1D(DIM) = N = 2 * M + 1
!    so that the points have coordinates
!      -M/M, -(M-1)/M, ..., -1/M, 0/M, 1/M, 2/M, 3/M, ..., (M-1)/M, M/M.
!    and we index them as
!      -M,   -(M-1),        -1,   0,   1,   2,   3,   ...,  M-1,    M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension of the points.
!
!    Input, integer ( kind=Ikind ) ORDER_1D(DIM_NUM), the order of the
!    rule in each dimension.
!
!    Input, integer ( kind=Ikind ) ORDER_ND, the product of the entries
!    of ORDER_1D.
!
!    Output, integer ( kind=Ikind ) INDX(DIM_NUM,ORDER_ND), the indices of
!    the points in the grid.  The second dimension of this array is equal
!    to the product of the entries of ORDER_1D.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  integer ( kind=Ikind ) order_nd

  integer ( kind=Ikind ) a(dim_num)
  logical more
  integer ( kind=Ikind ) order_1d(dim_num)
  integer ( kind=Ikind ) p
  integer ( kind=Ikind ) indx(dim_num,order_nd)

  more = .false.
  p = 0

  do

    call vec_colex_next2 ( dim_num, order_1d, a, more )

    if ( .not. more ) then
      exit
    end if

    p = p + 1
!
!  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1 = N - 1 = 2 * M.
!  Subtracting M sets the range to -M to +M, as we wish.
!
    indx(1:dim_num,p) = a(1:dim_num) - ( order_1d(1:dim_num) - 1 ) / 2

  end do

  return
end subroutine multigrid_index_z
subroutine product_weight_herm ( dim_num, order_1d, order_nd, w_nd )

!*****************************************************************************80
!
!! PRODUCT_WEIGHT_HERM: weights for a product Gauss-Hermite rule.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D Gauss-Hermite rules of varying order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) ORDER_1D(DIM_NUM), the order of the 1D rules.
!
!    Input, integer ( kind=Ikind ) ORDER_ND, the order of the product rule.
!
!    Output, real ( kind=Rkind ) W_ND(ORDER_ND), the product rule weights.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  integer ( kind=Ikind ) order_nd

  integer ( kind=Ikind ) dim
  integer ( kind=Ikind ) order_1d(dim_num)
  real    ( kind=Rkind ), allocatable, dimension ( : ) :: w_1d
  real    ( kind=Rkind ) w_nd(order_nd)

  w_nd(1:order_nd) = ONE

  do dim = 1, dim_num

    allocate ( w_1d(1:order_1d(dim)) )

    call hermite_weights ( order_1d(dim), w_1d )

    call r8vec_direct_product2 ( dim, order_1d(dim), w_1d, dim_num, &
      order_nd, w_nd )

    deallocate ( w_1d )
 
  end do

  return
end subroutine product_weight_herm
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Formula:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) N, the argument of the double factorial
!    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
!
!    Output, real ( kind=Rkind ) R8_FACTORIAL2, the value of the double
!    factorial of N.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) n
  real    ( kind=Rkind ) r8_factorial2
  real    ( kind=Rkind ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = ONE
    return
  end if

  r8_n = real ( n, kind=Rkind )
  r8_factorial2 = ONE

  do while ( ONE < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - TWO
  end do

  return
end function r8_factorial2
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind=Rkind ) R8_HUGE, a "huge" value.
!
  USE mod_system
  implicit none

  real    ( kind=Rkind ) r8_huge

  r8_huge = TEN**30

  return
end function r8_huge
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind=Ikind ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind=Rkind ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind=Ikind ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind=Ikind ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind=Rkind ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.  
!
!  Local Parameters:
!
!    Local, integer ( kind=Ikind ) START, the first location of a block of values
!    to set.
!
!    Local, integer ( kind=Ikind ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) factor_num
  integer ( kind=Ikind ) factor_order
  integer ( kind=Ikind ) point_num

  integer ( kind=Ikind ), save :: contig
  integer ( kind=Ikind ) factor_index
  real    ( kind=Rkind ) factor_value(factor_order)
  integer ( kind=Ikind ) j
  integer ( kind=Ikind ) k
  integer ( kind=Ikind ), save :: rep
  integer ( kind=Ikind ), save :: skip
  integer ( kind=Ikind ) start
  real    ( kind=Rkind ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = ONE
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end subroutine r8vec_direct_product2
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  USE mod_system
  implicit none

  character c
  integer ( kind=Ikind ) get
  integer ( kind=Ikind ) put
  integer ( kind=Ikind ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end subroutine s_blank_delete
subroutine sparse_grid_herm ( dim_num, level_max, point_num, grid_weight, &
  grid_point )
  
!****************************************************************************80 
! 
!!  SPARSE_GRID_HERM computes a sparse grid of Gauss-Hermite points. 
! 
!  Discussion: 
! 
!    The quadrature rule is associated with a sparse grid derived from 
!    a Smolyak construction using a 1D Gauss-Hermite quadrature rule.  
! 
!    The user specifies: 
!    * the spatial dimension of the quadrature region, 
!    * the level that defines the Smolyak grid. 
! 
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2008
! 
!  Author: 
! 
!    John Burkardt 
! 
!  Parameters: 
! 
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
! 
!    Input, integer ( kind=Ikind ) LEVEL_MAX, controls the size of the
!    sparse grid. 
! 
!    Input, integer ( kind=Ikind ) POINT_NUM, the number of points in the grid,
!    as determined by SPARSE_GRID_HERM_SIZE. 
! 
!    Output, real ( kind=Rkind ) GRID_WEIGHT(POINT_NUM), the weights.
! 
!    Output, real ( kind=Rkind ) GRID_POINT(DIM_NUM,POINT_NUM), the points.
! 
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  integer ( kind=Ikind ) point_num

  integer ( kind=Ikind ) choose
  integer ( kind=Ikind ) coeff
  
  
  integer ( kind=Ikind ), dimension ( dim_num ) :: grid_base2
  
  integer ( kind=Ikind ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind=Ikind ), allocatable, dimension ( : ) :: grid_level
  real    ( kind=Rkind ) grid_point(dim_num,point_num)
  real    ( kind=Rkind ) grid_point_temp(dim_num)
  real    ( kind=Rkind ) grid_weight(point_num)
  real    ( kind=Rkind ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind=Ikind ) h
  integer ( kind=Ikind ), parameter :: i4_1 = 1
  
  integer ( kind=Ikind ) level
  integer ( kind=Ikind ), dimension ( dim_num ) :: level_1d
  integer ( kind=Ikind ) level_max
  integer ( kind=Ikind ) level_min
  logical more
  integer ( kind=Ikind ), dimension ( dim_num ) :: order_1d
  integer ( kind=Ikind ) order_nd
  integer ( kind=Ikind ) point
  integer ( kind=Ikind ) point2
  integer ( kind=Ikind ) point3
  integer ( kind=Ikind ) point_num2
  integer ( kind=Ikind ) t

  grid_weight(1:point_num) = ZERO
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  point_num2 = 0

  level_min = max ( 0, level_max + 1 - dim_num )
  
  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!  The relationship is the same as for other OPEN rules.
!  The GL rule differs from the other OPEN rules only in the nesting behavior.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      grid_base2(1:dim_num) = ( order_1d(1:dim_num) - 1 ) / 2
!
!  The product of the 1D orders gives us the number of points in this subgrid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_level(1:order_nd) )
      allocate ( grid_weight2(1:order_nd) )
!
!  Compute the weights for this product grid.
!
      call product_weight_herm ( dim_num, order_1d, order_nd, grid_weight2 )
!
!  Now determine the coefficient of the weight.
!
      coeff = (-1)**( level_max - level ) &
        * choose ( dim_num - 1, level_max - level )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
!
      call multigrid_index_z ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Determine the first level of appearance of each of the points.
!  This allows us to flag certain points as being repeats of points
!  generated on a grid of lower level.  
!
!  This is SLIGHTLY tricky.
!
      call index_level_herm ( level, level_max, dim_num, order_nd, grid_index2, &
        grid_base2, grid_level )
!
!  Only keep those points which first appear on this level.
!
      do point = 1, order_nd
!
!  Either a "new" point (increase count, create point, create weight)
!
        if ( grid_level(point) == level ) then

          point_num2 = point_num2 + 1

          call hermite_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
            grid_base2(1:dim_num), grid_point(1:dim_num,point_num2) )

          grid_weight(point_num2) = real ( coeff, kind=Rkind ) * grid_weight2(point)
!
!  or an already existing point (create point temporarily, find match,
!  add weight to matched point's weight).
!
        else
 
          call hermite_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
            grid_base2(1:dim_num), grid_point_temp(1:dim_num) )
            
          point3 = -1
          
          do point2 = 1, point_num2
            if ( all ( grid_point(1:dim_num,point2) == grid_point_temp(1:dim_num) ) ) then
              point3 = point2
              exit
            end if
          end do
          
          if ( point3 == -1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPARSE_GRID_HERMITE - Fatal error!'
            write ( *, '(a)' ) '  Could not match point.'
            stop
          end if
          
          grid_weight(point3) = grid_weight(point3) + &
            real ( coeff, kind=Rkind ) * grid_weight2(point)
                   
        end if

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_level )
      deallocate ( grid_weight2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  return
end subroutine sparse_grid_herm
subroutine sparse_grid_herm_index ( dim_num, level_max, point_num, grid_index, &
  grid_base )

!*****************************************************************************80
!
!! SPARSE_GRID_HERM_INDEX indexes points in a sparse Gauss-Hermite grid.
!
!  Discussion:
!
!    The sparse grid is assumed to be formed from 1D Gauss-Hermite rules
!    of ODD order, which have the property that only the central abscissa,
!    X = 0.0, is "nested".
!
!    The necessary dimensions of GRID_INDEX can be determined by 
!    calling SPARSE_GRID_HERM_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind=Ikind ) POINT_NUM, the total number of points in
!    the grids.
!
!    Output, integer ( kind=Ikind ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level 
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, integer ( kind=Ikind ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the Gauss-Hermite rules associated with each point 
!    and dimension.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num
  integer ( kind=Ikind ) point_num

  
  
  integer ( kind=Ikind ) grid_base(dim_num,point_num)
  integer ( kind=Ikind ), dimension ( dim_num ) :: grid_base2
  integer ( kind=Ikind ) grid_index(dim_num,point_num)
  integer ( kind=Ikind ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind=Ikind ), allocatable, dimension ( : ) :: grid_level
  integer ( kind=Ikind ) h
  
  integer ( kind=Ikind ) level
  integer ( kind=Ikind ), dimension ( dim_num ) :: level_1d
  integer ( kind=Ikind ) level_max
  integer ( kind=Ikind ) level_min
  logical more
  integer ( kind=Ikind ), dimension ( dim_num ) :: order_1d
  integer ( kind=Ikind ) order_nd
  integer ( kind=Ikind ) point
  integer ( kind=Ikind ) point_num2
  integer ( kind=Ikind ) t
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  point_num2 = 0

  level_min = max ( 0, level_max + 1 - dim_num )
  
  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!  The relationship is the same as for other OPEN rules.
!  The GL rule differs from the other OPEN rules only in the nesting behavior.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      grid_base2(1:dim_num) = ( order_1d(1:dim_num) - 1 ) / 2
!
!  The product of the 1D orders gives us the number of points in this subgrid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_level(1:order_nd) )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
!
      call multigrid_index_z ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Determine the first level of appearance of each of the points.
!  This allows us to flag certain points as being repeats of points
!  generated on a grid of lower level.  
!
!  This is SLIGHTLY tricky.
!
      call index_level_herm ( level, level_max, dim_num, order_nd, grid_index2, &
        grid_base2, grid_level )
!
!  Only keep those points which first appear on this level.
!
      do point = 1, order_nd

        if ( grid_level(point) == level ) then

          point_num2 = point_num2 + 1

          grid_index(1:dim_num,point_num2) = grid_index2(1:dim_num,point)
          grid_base(1:dim_num,point_num2) = grid_base2(1:dim_num)

        end if

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_level )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  return
end subroutine sparse_grid_herm_index
subroutine sparse_grid_herm_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_HERM_SIZE sizes a sparse grid of Gauss-Hermite points.
!
!  Discussion:
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      LEVEL_MIN <= LEVEL <= LEVEL_MAX.
!
!    where LEVEL_MAX is user specified, and 
!
!      LEVEL_MIN = max ( 0, LEVEL_MAX + 1 - DIM_NUM ).
!
!    The grids are only very weakly nested, since Gauss-Hermite rules
!    only have the origin in common.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind=Ikind ) POINT_NUM, the number of points in the grid.
!
  USE mod_system
  implicit none
  
  integer ( kind=Ikind ) dim_num
  
  integer ( kind=Ikind ) dim
  integer ( kind=Ikind ) h
  integer ( kind=Ikind ) level
  integer ( kind=Ikind ) level_1d(dim_num)
  integer ( kind=Ikind ) level_max
  integer ( kind=Ikind ) level_min
  logical more
  
  integer ( kind=Ikind ) order_1d(dim_num)
  
  integer ( kind=Ikind ) point_num
  integer ( kind=Ikind ) t
!
!  Special case.
!
  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  point_num = 0

  level_min = max ( 0, level_max + 1 - dim_num )
  
  do level = level_min, level_max
!
!  The middle loop generates the next partition that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      do dim = 1, dim_num
!
!  If we can reduce the level in this dimension by 1 and
!  still not go below LEVEL_MIN.
!
        if ( level_min < level .and. 1 < order_1d(dim) ) then
          order_1d(dim) = order_1d(dim) - 1
        end if
		
      end do

      point_num = point_num + product ( order_1d(1:dim_num) )

      if ( .not. more ) then
        exit
      end if
	  
    end do
  end do

  return
end subroutine sparse_grid_herm_size
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  USE mod_system
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = [&
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ']
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  USE mod_system
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = [&
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ']
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestring
subroutine vec_colex_next2 ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_COLEX_NEXT2 generates vectors in colex order.
!
!  Discussion:
!
!    The vectors are produced in colexical order, starting with
!
!    (0,        0,        ...,0),
!    (1,        0,        ...,0),
!     ...
!    (BASE(1)-1,0,        ...,0)
!
!    (0,        1,        ...,0)
!    (1,        1,        ...,0)
!    ...
!    (BASE(1)-1,1,        ...,0)
!
!    (0,        2,        ...,0)
!    (1,        2,        ...,0)
!    ...
!    (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
!
!  Examples:
!
!    DIM_NUM = 2,
!    BASE = ( 3, 3 )
!
!    0   0
!    1   0
!    2   0
!    0   1
!    1   1
!    2   1
!    0   2
!    1   2
!    2   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind=Ikind ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind=Ikind ) BASE(DIM_NUM), the bases to be used in each dimension.
!    In dimension I, entries will range from 0 to BASE(I)-1.
!
!    Input/output, integer ( kind=Ikind ) A(DIM_NUM).  On each return, A
!    will contain entries in the range 0 to N-1.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output 
!    vector and stop calling the routine.
!
  USE mod_system
  implicit none

  integer ( kind=Ikind ) dim_num

  integer ( kind=Ikind ) a(dim_num)
  integer ( kind=Ikind ) base(dim_num)
  integer ( kind=Ikind ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 0
    more = .true.

  else

    do i = 1, dim_num

      a(i) = a(i) + 1

      if ( a(i) < base(i) ) then
        return
      end if

      a(i) = 0

    end do

    more = .false.

  end if

  return
end subroutine vec_colex_next2
