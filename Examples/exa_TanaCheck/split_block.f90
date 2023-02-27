 MODULE mod_rec_vec
 IMPLICIT NONE

 integer :: type_vec     = 0  ! 0: in star, all coupled, 1: in line, 2: jacobi
 integer :: max_layer    = 4
 integer :: max_v        = 2
 integer :: max_b        = 2
 integer :: nrho         = 2
 logical :: cart         = .FALSE.
 logical :: SpheConv_xzy = .FALSE.
 logical :: cos_th       = .TRUE.
 logical :: frame0vec    = .FALSE. ! possibility to use a frame with 0 vector (nb_vect=0)
 logical :: type100      = .FALSE. ! possibility to use constraints
 logical :: zmat_order   = .TRUE.  ! possibility to use zmat_order

 integer,parameter :: Line_len=256

 character(len=*), parameter :: vec_name_poly  = " &vector /"
 character(len=*), parameter :: vec_name_cart  = " &vector cart=t /"
 character(len=*), parameter :: vec_name_xzy   = " &vector Spherical_convention='x-zy' /"
 character(len=*), parameter :: vec_name_frame = " &vector frame=t nb_vect=0 /"

 TYPE recBF
   character(len=:), allocatable :: vec_name
   integer                       :: layer = 0
   TYPE(recBF), pointer          :: tab_BF(:) => null()
 END TYPE recBF

CONTAINS

RECURSIVE SUBROUTINE dealloc_recBF(BF)
  IMPLICIT NONE
  TYPE(recBF), intent(inout) :: BF

  integer :: i

  write(6,*) ' dealloc BF alloc vec_name',allocated(BF%vec_name)

  IF (allocated(BF%vec_name)) THEN
    deallocate(BF%vec_name)
  END IF
  BF%layer = 0
  IF (associated(BF%tab_BF)) THEN
    DO i=1,size(BF%tab_BF)
      CALL dealloc_recBF(BF%tab_BF(i))
    END DO
    deallocate(BF%tab_BF)
    nullify(BF%tab_BF)
  END IF

END SUBROUTINE dealloc_recBF
RECURSIVE SUBROUTINE write_recBF(BF,nio,rec)
  IMPLICIT NONE
  TYPE(recBF), intent(in) :: BF
  integer,      intent(in) :: nio
  integer, optional,      intent(in) :: rec

  integer :: i,rec_loc
  character(len=:), allocatable :: vec_name

  IF (present(rec)) THEN
    rec_loc = rec
  ELSE
    rec_loc = 0
  END IF

  IF (allocated(BF%vec_name)) THEN
    vec_name = trim(BF%vec_name)
    DO i=1,rec_loc
      vec_name = trim("   " // vec_name)
    END DO

    write(nio,'(a)') vec_name ; flush(nio)
    deallocate(vec_name)
  END IF
  !write(6,*) 'asso tab_BF',associated(BF%tab_BF)

  IF (associated(BF%tab_BF)) THEN
    rec_loc = rec_loc + 1
    DO i=1,size(BF%tab_BF)
      CALL write_recBF(BF%tab_BF(i),nio,rec_loc)
    END DO
  END IF

END SUBROUTINE write_recBF
RECURSIVE SUBROUTINE recBF2_TO_recBF1(BF1,BF2)
  IMPLICIT NONE
  TYPE(recBF), intent(in) :: BF2
  TYPE(recBF), intent(inout) :: BF1

  integer :: i

   BF1%vec_name = trim(BF2%vec_name)
   BF1%layer    = BF2%layer

   IF (associated(BF2%tab_BF)) THEN
     allocate(BF1%tab_BF(size(BF2%tab_BF)))
     DO i=1,size(BF1%tab_BF)
       CALL recBF2_TO_recBF1(BF1%tab_BF(i),BF2%tab_BF(i))
     END DO
   END IF


END SUBROUTINE recBF2_TO_recBF1
RECURSIVE FUNCTION get_layer_OF_recBF(BF) RESULT(layer)
  IMPLICIT NONE
  TYPE(recBF), intent(in) :: BF
  integer :: layer

  integer, allocatable :: tab_layer(:)
  integer :: i



   IF (associated(BF%tab_BF)) THEN
     allocate(tab_layer(size(BF%tab_BF)))

     DO i=1,size(BF%tab_BF)
       tab_layer(i) = get_layer_OF_recBF(BF%tab_BF(i))
     END DO
     layer = 1 + maxval(tab_layer)

     deallocate(tab_layer)

   ELSE
     layer = 1
   END IF


END FUNCTION get_layer_OF_recBF
FUNCTION trim(string,ltrim)
  character(len=:), allocatable     :: trim
  character(len=*), intent(in)      :: string
  logical, optional,intent(in)      :: ltrim

  logical :: ltrim_loc

  IF (allocated(trim)) deallocate(trim)

  IF (present(ltrim)) THEN
    ltrim_loc = ltrim
  ELSE
    ltrim_loc = .TRUE.
  END IF

  IF (ltrim_loc) THEN
    allocate(character(len=len_trim(string)) :: trim)
    trim = trim(string)
  ELSE
    allocate(character(len=len(string)) :: trim)
    trim = string
  END IF

END FUNCTION trim

FUNCTION TO_string(i)
  character (len=:), allocatable  :: TO_string
  integer, intent(in) :: i

  character(len=Line_len) :: name_int
  integer :: clen

  write(name_int,*) i

  clen = len_trim(adjustl(name_int))
  allocate(character(len=clen) :: TO_string)

  TO_string = trim(adjustl(name_int))

END FUNCTION TO_string

SUBROUTINE get_next_tab_v(tab_v,nv,nb,end_tab)
  integer, intent(in) :: nv,nb
  integer, allocatable, intent(inout) :: tab_v(:)
  logical, intent(inout) :: end_tab
  integer :: i

  end_tab = .FALSE.

   !write(6,*) 'IN get_next_tab_v'
   !write(6,*) ' nv,nb',nv,nb  ; flush(6)
   IF (.NOT. allocated(tab_v)) THEN

     allocate(tab_v(nb))

     tab_v(:) = 1
     IF (nb == 1) THEN
       tab_v(1) = nv-1
     ELSE
       tab_v(1) = 0
     END IF
   END IF

   IF (nb == 1) THEN
     tab_v(1) = tab_v(1) + 1
     end_tab  = (tab_v(1) > nv)
   ELSE
     DO i=1,nb-1
       tab_v(i)  = tab_v(i) + 1
       tab_v(nb) = nv - sum(tab_v(1:nb-1))

       IF (tab_v(nb) > 0) EXIT ! it means that the table values are OK.

       tab_v(i) = 1
     END DO

     end_tab = (tab_v(nb) <= 0) ! it means that no other correct table values are possible.
   END IF

   IF (end_tab) deallocate(tab_v)

   !IF (allocated(tab_v)) write(6,*) 'tab_v',tab_v
   !write(6,*) 'OUT get_next_tab_v' ; flush(6)

END SUBROUTINE get_next_tab_v
SUBROUTINE get_next_DP_ind(tab_ind,tab_n,end_tab)
  integer, allocatable, intent(in)    :: tab_n(:)
  integer, allocatable, intent(inout) :: tab_ind(:)
  logical, intent(inout) :: end_tab

  integer :: i,n

  end_tab = .FALSE.
  n = size(tab_n)

  IF (.NOT. allocated(tab_ind)) THEN
    tab_ind    = tab_n
    tab_ind(:) = 1
    tab_ind(1) = 0
  END IF

  tab_ind(1) = tab_ind(1) + 1

  DO i=1,n-1
    IF (tab_ind(i) > tab_n(i)) THEN
      tab_ind(i)   = 1
      tab_ind(i+1) = tab_ind(i+1) + 1
    END IF
  END DO
  IF (tab_ind(n) > tab_n(n)) THEN
    end_tab = .TRUE.
    deallocate(tab_ind)
  END IF

END SUBROUTINE get_next_DP_ind

RECURSIVE SUBROUTINE Split_blocks_new(nv)
  IMPLICIT NONE
  integer, intent(in)         :: nv


  TYPE(recBF) :: BF_vec_poly
  TYPE(recBF) :: BF_vec_cart
  TYPE(recBF) :: BF_vec_frame
  TYPE(recBF) :: BF_vec_xzy

  TYPE(recBF) :: BF_vec(nv)


  integer :: nb,i,nio,iv,iiv,i_BF,nb_BF,nb_BF_tmp
  integer,allocatable :: tab_v(:),tab_ind(:),tab_ndim(:)
  logical :: end_tab


  character(len=:), allocatable     :: cart_loc
  character(len=:), allocatable     :: data_name


  write(6,*) 'BEGINNING Split_blocks_new'
  write(6,*) 'nv',nv


  ! initialization for nv=1 (just one BF), but 3 possibilities...

  BF_vec_poly%vec_name  = vec_name_poly
  BF_vec_poly%layer     = 1
  BF_vec_cart%vec_name  = vec_name_cart
  BF_vec_cart%layer     = 1
  BF_vec_xzy%vec_name   = vec_name_xzy
  BF_vec_xzy%layer      = 1
  BF_vec_frame%vec_name = vec_name_frame
  BF_vec_frame%layer    = 1

  allocate(BF_vec(1)%tab_BF(1)) ! number of i_BF for v=1
  BF_vec(1)%tab_BF(1)%vec_name = trim(" &vector /")  ! it is not going to be used
  BF_vec(1)%tab_BF(1)%layer    = 1


  DO iv=2,max_v
    ! first the number of terms, for the allocation of BF_vec(iv)%tab_BF(:)
    i_BF = 0
    DO nb=1,max_b
    DO

      CALL get_next_tab_v(tab_v,iv-1,nb,end_tab)
      IF (end_tab) EXIT
      write(6,*) 'iv,nb',iv,nb,'tab_v: ',tab_v,' sum:',sum(tab_v) ; flush(6)

      allocate(tab_ndim(nb))
      DO i=1,nb
        iiv = tab_v(i)
        tab_ndim(i) = size(BF_vec(iiv)%tab_BF)
      END DO
      DO ! because in BF_vec(iiv) there are several BF (in tab_BF)
        CALL get_next_DP_ind(tab_ind,tab_ndim,end_tab)
        IF (end_tab) EXIT

        i_BF = i_BF + 1

      END DO
      deallocate(tab_ndim)


    END DO
    END DO
    nb_BF = i_BF
    write(6,*) 'iv,nb_BF',iv,nb_BF ; flush(6)

    allocate(BF_vec(iv)%tab_BF(nb_BF))

    i_BF = 0
    DO nb=1,max_b
    DO

      CALL get_next_tab_v(tab_v,iv-1,nb,end_tab)
      IF (end_tab) EXIT


      allocate(tab_ndim(nb))
      DO i=1,nb
        iiv = tab_v(i)
        tab_ndim(i) = size(BF_vec(iiv)%tab_BF)
      END DO
      DO ! because in BF_vec(iiv) there are several BF (in tab_BF)
        CALL get_next_DP_ind(tab_ind,tab_ndim,end_tab)
        IF (end_tab) EXIT

        i_BF = i_BF + 1

        allocate(BF_vec(iv)%tab_BF(i_BF)%tab_BF(nb))
        IF (zmat_order) THEN
          BF_vec(iv)%tab_BF(i_BF)%vec_name = trim( &
            " &Vector frame=t nb_vect=" // TO_string(nb) // " zmat_order=t /")
        ELSE
          BF_vec(iv)%tab_BF(i_BF)%vec_name = trim( &
            " &Vector frame=t nb_vect=" // TO_string(nb) // " zmat_order=f /")
        END IF

        DO i=1,nb
          iiv = tab_v(i)
          IF (iiv == 1) THEN ! special case with 1 vector (because of the 3 possibilities)
            IF (cart .AND. i > 1) THEN
              CALL recBF2_TO_recBF1(BF_vec(iv)%tab_BF(i_BF)%tab_BF(i), &
                                    BF_vec_cart)
            ELSE IF (SpheConv_xzy .AND. i > 1) THEN
              CALL recBF2_TO_recBF1(BF_vec(iv)%tab_BF(i_BF)%tab_BF(i), &
                                    BF_vec_xzy)
            ELSE IF (frame0vec) THEN
              CALL recBF2_TO_recBF1(BF_vec(iv)%tab_BF(i_BF)%tab_BF(i), &
                                    BF_vec_frame)
            ELSE
              CALL recBF2_TO_recBF1(BF_vec(iv)%tab_BF(i_BF)%tab_BF(i), &
                                    BF_vec_poly)
            END IF
          ELSE
            CALL recBF2_TO_recBF1(BF_vec(iv)%tab_BF(i_BF)%tab_BF(i), &
                                  BF_vec(iiv)%tab_BF(tab_ind(i)))
          END IF
          !BF_vec(iv)%tab_BF(i_BF)%tab_BF(i) = BF_vec(iiv)%tab_BF(tab_ind(i))
        END DO

      END DO
      deallocate(tab_ndim)


    END DO
    END DO

  END DO

  !write(6,*) ' All BF from 2 up to max_v: ',max_v
  !DO iv=2,max_v
  iv = max_v
  DO i_BF=1,size(BF_vec(iv)%tab_BF)
    write(6,*)
    write(6,*) '============================='
    write(6,*) ' iv,i_BF: ',iv,i_BF
    write(6,*) ' iv,i_BF, layer: ',iv,i_BF,get_layer_OF_recBF(BF_vec(iv)%tab_BF(i_BF))
    !CALL write_recBF(BF_vec(iv)%tab_BF(i_BF),6)

    IF (get_layer_OF_recBF(BF_vec(iv)%tab_BF(i_BF)) > max_layer) CYCLE


    data_name = trim( "dat_nv" // TO_string(iv) // "-num" // TO_string(i_BF) )

    nio=10
    open(unit=nio,file=data_name)

    CALL write_header(iv,nio,type_vec)
    write(nio,*)
    CALL write_recBF(BF_vec(iv)%tab_BF(i_BF),nio)

    CALL write_footer(iv,nio)

    close(nio)


  END DO
  !END DO
  write(6,*) '============================='

  write(6,*) 'END Split_blocks_new'

END SUBROUTINE Split_blocks_new

SUBROUTINE write_header(nv,nio,type_vec)
  integer, intent(in) ::  nv,nio
  integer ::  type_vec

  integer ::  i,j,nb_x
  real (kind=8) :: mass(nv)

  nb_X = 0
  IF (type_vec == 2) nb_X=nv-1

  CALL random_number(mass)
  mass = 1.d0+mass


  write(nio,*) "&variables"
  write(nio,*) "    ", 'nrho=',nrho
  write(nio,*) "    ", 'Old_Qtransfo=f'
  write(nio,*) "    ", 'Tana=t Tana_Init_Only=t VSCFform=t LaTeXForm=f '
  write(nio,*) "    ", 'nb_Qtransfo=3'
  write(nio,*) '/'
  write(nio,*) " &Coord_transfo name_transfo='bunch' nb_vect=",nv," nb_X=",nb_X," inTOout=f /"
  write(nio,*) 10.,(mass(j),j=1,nv)
  SELECT CASE (type_vec)
  CASE (0) ! vectors in star
    do j=1, nv
      write(nio,*) 1,j+1
    end do

  CASE (1) ! vectors in line
    do j=1, nv
      write(nio,*) j,j+1
    end do

  CASE (2) ! vectors in line
    DO j=1,nb_X
      write(nio,*) "&dummyX tab_At_TO_X=",(i,i=1,j+1), "type_dummyX='COM' /"
    END DO
    write(nio,*) 1,2
    do j=2, nv
      write(nio,*) j+nb_X+1,j+1
    end do

  CASE default ! vectors in star
    do j=1, nv
      write(nio,*) 1,j+1
    end do
  END SELECT

  write(nio,*) "&Coord_transfo name_transfo='poly' cos_th=",cos_th," / "
  write(nio,*)

END SUBROUTINE write_header
SUBROUTINE write_footer(nv,nio)
  integer, intent(in) ::  nv,nio
  integer ::  nb_var,k
  real(kind=8) :: Qrand
  integer, allocatable :: type_var(:)

  nb_var = max(1,3*nv-3)
  allocate(type_var(nb_var))
  DO k=1,nb_var
   CALL random_number(Qrand)
   IF (Qrand < 0.3) THEN
     type_var(k) = 100
   ELSE
     type_var(k) = 1
   END IF
  END DO
  IF (.NOT. type100) type_var(:) = 1
  write(nio,*)
  write(nio,*) " &Coord_transfo name_transfo='active' /"
  write(nio,*) type_var
  write(nio,*) "&minimum Read_nameQ=f /"

  DO k=1, nb_var
    CALL random_number(Qrand)
    write(nio,*) 0.5d0 + 0.3*Qrand
  END DO
  write(nio,*)
  deallocate(type_var)

END SUBROUTINE write_footer

END MODULE mod_rec_vec

SUBROUTINE read_arg()
  USE mod_rec_vec
  IMPLICIT NONE
  integer :: nv


  character(len=:), allocatable :: arg,arg2
  integer :: long,i,err_read

  namelist / rec_vec / type_vec,max_layer,max_v,max_b,nrho,cart,   &
                       cos_th,frame0vec,type100,zmat_order,SpheConv_xzy


  DO i=1, COMMAND_ARGUMENT_COUNT(),2
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=long )
    allocate( character(len=long) :: arg )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=long )
    allocate( character(len=long) :: arg2 )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

    SELECT CASE(arg)
    CASE("-type_vec","-tv")
      read(arg2,*) type_vec
    CASE("-v","-vec")
      read(arg2,*) max_v
    CASE("-b","-block")
      read(arg2,*) max_b
    CASE("-l","-layer")
      read(arg2,*) max_layer
    CASE("-cart")
      read(arg2,*) cart
    CASE("-cos_th")
      read(arg2,*) cos_th
    CASE("-frame0vec")
      read(arg2,*) frame0vec
    CASE("-nrho")
      read(arg2,*) nrho
    CASE("-SpheConv_xzy","-SCxzy")
      read(arg2,*) SpheConv_xzy
    END SELECT

    print *,"Argument de rang ", i, " ==> ", arg , ' arg_val: ',arg2

    deallocate(arg)
    deallocate(arg2)
  END DO



  read(5,rec_vec,IOSTAT=err_read)
  IF (err_read /= 0) THEN
    write(6,*) ' WARNING: no namelist or errors while reading the namelist "rec_vec"'
  END IF
  !write(6,rec_vec)
  write(6,*) 'type_vec  ',type_vec
  write(6,*) 'max_v     ',max_v
  write(6,*) 'max_b     ',max_b
  write(6,*) 'max_layer ',max_layer
  write(6,*) 'cart,cos_th,frame0vec,SpheConv_xzy ',cart,cos_th,frame0vec,SpheConv_xzy
  write(6,*) 'nrho      ',nrho



  nv = max_v
  write(6,*) 'split ',nv ; flush(6)

  CALL Split_blocks_new(nv)

  write(6,*) '=================================='
  write(6,*) '=================================='

END SUBROUTINE read_arg

PROGRAM test
  USE mod_rec_vec
  IMPLICIT NONE
  integer :: nv


  character(len=:), allocatable :: arg,arg2
  integer :: long,i

  write(6,*) '=================================='
  write(6,*) 'AUTOMATIC GENERATION of Tana data'
  write(6,*) '=================================='
  write(6,*) 'Code written by David Lauvergnat [1]'
  write(6,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
  write(6,*) '=================================='

  CALL read_arg()

  nv = max_v
  write(6,*) 'split ',nv ; flush(6)

  CALL Split_blocks_new(nv)

  write(6,*) '=================================='
  write(6,*) '=================================='

END PROGRAM test
