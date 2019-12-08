!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================
      MODULE mod_Qtransfo
      use mod_system
      USE mod_dnSVM
      use mod_constant, only: table_atom

      USE mod_CartesianTransfo
      USE mod_QTOXanaTransfo
      USE mod_BunchPolyTransfo
      USE mod_ZmatTransfo
      USE mod_RectilinearNM_Transfo
      USE mod_OneDTransfo
      USE mod_ThreeDTransfo
      USE mod_Rot2CoordTransfo
      USE mod_FlexibleTransfo
      USE mod_GeneTransfo
      USE mod_HyperSpheTransfo
      USE mod_LinearNMTransfo
      USE mod_RPHTransfo
      USE mod_ActiveTransfo

      IMPLICIT NONE

      PRIVATE

        TYPE, PUBLIC :: Type_Qtransfo
          logical                           :: print_done      = .FALSE.
          character (len=Name_len)          :: name_transfo    = "identity"
          integer                           :: num_transfo     = 0
          logical                           :: inTOout         = .TRUE.
          integer                           :: nb_var          = 0
          integer                           :: nb_act          = 0
          integer                           :: nb_transfo      = 0
          integer                           :: opt_transfo     = 0 ! option for the transformation
          logical                           :: skip_transfo    = .FALSE.
          integer                           :: opt_param       = 0
          logical                           :: Primitive_coord = .FALSE.


          TYPE (Type_CartesianTransfo)      :: CartesianTransfo
          TYPE (Type_QTOXanaTransfo)        :: QTOXanaTransfo
          TYPE (Type_ZmatTransfo)           :: ZmatTransfo
          TYPE (Type_RectilinearNM_Transfo) :: RectilinearNM_Transfo
          TYPE (Type_BunchTransfo),pointer  :: BunchTransfo       => null()
          TYPE (Type_BFTransfo)             :: BFTransfo

          TYPE (Type_LinearTransfo)         :: LinearTransfo
          TYPE (Type_FlexibleTransfo)       :: FlexibleTransfo
          TYPE (Type_GeneTransfo)           :: GeneTransfo

          TYPE (Type_oneDTransfo),pointer      :: oneDTransfo(:)      => null()
          TYPE (Type_ThreeDTransfo),pointer    :: ThreeDTransfo       => null()
          TYPE (Type_Rot2CoordTransfo),pointer :: Rot2CoordTransfo(:) => null()
          TYPE (Type_HyperSpheTransfo)         :: HyperSpheTransfo
          integer, pointer                     :: list_Qin_TO_Qout(:) => null() ! "order" transfo

          TYPE (Type_NMTransfo), pointer    :: NMTransfo           => null()
          TYPE (Type_RPHTransfo), pointer   :: RPHTransfo          => null()
          TYPE (Type_ActiveTransfo), pointer:: ActiveTransfo       => null()

          integer                           :: nb_Qin       = 0  ! size the input coordinates
          integer                           :: nb_Qout      = 0 ! size the output coordinates
          integer, pointer                  :: type_Qin(:)  => null() ! size nb_Qin
          integer, pointer                  :: type_Qout(:) => null() ! true pointer (will point to the previous type_Qin)
                                                                      ! except for the first transfo
          character (len=Name_len), pointer :: name_Qin(:)  => null()
          character (len=Name_len), pointer :: name_Qout(:) => null() ! true pointer (will point to the previous name_Qin)
                                                                      ! except for the first transfo

        END TYPE Type_Qtransfo

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_Qtransfodim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_Qtransfodim1
      END INTERFACE

      PUBLIC alloc_array,dealloc_array,dealloc_Qtransfo
      PUBLIC read_Qtransfo,Write_Qtransfo,Sub_Check_LinearTransfo,sub_Type_Name_OF_Qin
      PUBLIC Qtransfo1TOQtransfo2,calc_Qtransfo

      CONTAINS

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE read_Qtransfo(Qtransfo,nb_Qin,mendeleev)

        TYPE (Type_Qtransfo), intent(inout) :: Qtransfo
        integer, intent(inout)              :: nb_Qin
        TYPE (table_atom), intent(in)       :: mendeleev

        character (len=Name_len) :: name_transfo,name_dum
        integer :: nat,nb_vect,nbcol,nb_flex_act,nb_transfo,nb_G,nb_X
        integer :: opt_transfo
        logical :: skip_transfo
        logical :: cos_th,purify_hess,eq_hess,k_Half,inTOout
        logical :: hessian_old,hessian_onthefly,hessian_cart,d0c_read
        character (len=line_len)      :: file_hessian
        logical :: hessian_read,k_read,with_vectors,not_all
        logical :: check_LinearTransfo
        integer :: i,it,i_Q,iF_inout,iat,iQin,iQout,nb_read
        real (kind=Rkind) ::  at
        real (kind=Rkind), pointer ::  M_mass(:,:)
        character (len=Name_len), pointer :: name_at(:)


        namelist /Coord_transfo/ name_transfo,nat,nb_vect,cos_th,       &
                                 nb_G,nb_X,opt_transfo,skip_transfo,    &
                                 inTOout,with_vectors,                  &
                                 nb_transfo,purify_hess,eq_hess,k_Half, &
                              hessian_old,hessian_onthefly,file_hessian,&
                               hessian_cart,hessian_read,k_read,nb_read,&
                                 d0c_read,not_all,check_LinearTransfo
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory,err_io
      character (len=*), parameter :: name_sub = "read_Qtransfo"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
!-----------------------------------------------------------
       nullify(M_mass)

        name_transfo = "identity"
        IF (Qtransfo%num_transfo > 1 .AND. nb_Qin < 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_Qin < 1',nb_Qin
          write(out_unitp,*) 'and it is NOT the initial transformation',&
                 Qtransfo%num_transfo
          STOP
        END IF
        it               = Qtransfo%num_transfo
        opt_transfo      = 0
        skip_transfo     = .FALSE.
        inTOout          = .TRUE.
        nat              = 0
        nb_vect          = 0
        nb_G             = 0
        nb_X             = 0
        cos_th           = .TRUE.
        purify_hess      = .FALSE.
        eq_hess          = .FALSE.
        k_Half           = .FALSE.
        with_vectors     = .TRUE.
        hessian_old      = .TRUE.
        hessian_cart     = .TRUE.
        hessian_onthefly = .FALSE.
        file_hessian     = 'xx_freq.fchk'
        hessian_read     = .FALSE.
        k_read           = .FALSE.
        d0c_read         = .FALSE.
        nb_read          = 0
        nb_transfo       = 1
        not_all          = .FALSE.
        check_LinearTransfo = .TRUE.

        read(in_unitp,Coord_transfo,IOSTAT=err_io)
        IF (err_io < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the namelist "Coord_transfo"'
          write(out_unitp,*) ' end of file or end of record'
          write(out_unitp,*) ' Probably, nb_transfo is to large in the namelist "variables"'
          write(out_unitp,*) '   or you have forgotten a coordinate tranformation ...'
          write(out_unitp,*) '   or you have forgotten the "Cartesian transfo"'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        IF (err_io > 0) THEN
          write(out_unitp,Coord_transfo)
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the namelist "Coord_transfo"'
          write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        write(out_unitp,*) '=========================================='
        write(out_unitp,*) '=========================================='

        IF (debug) write(out_unitp,Coord_transfo)

        Qtransfo%name_transfo = name_transfo
        Qtransfo%inTOout      = inTOout
        Qtransfo%opt_transfo  = opt_transfo
        Qtransfo%skip_transfo = skip_transfo
        write(out_unitp,'(a,a)' ) ' transfo:               ',Qtransfo%name_transfo
        write(out_unitp,'(a,i0)') ' Option of the transfo: ',Qtransfo%opt_transfo
        write(out_unitp,'(a,l)' ) ' Skip the transfo:      ',Qtransfo%skip_transfo
        write(out_unitp,'(a,i0)') ' num_transfo:           ',Qtransfo%num_transfo
        write(out_unitp,'(a,l)' ) ' inTOout:               ',Qtransfo%inTOout
        write(out_unitp,'(a)'   ) '------------------------------------------'
        CALL flush_perso(out_unitp)


        name_transfo = Qtransfo%name_transfo
        CALL string_uppercase_TO_lowercase(name_transfo)

        SELECT CASE (name_transfo)
        CASE ('identity')
          Qtransfo%nb_Qin  = nb_Qin
          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qid")
          Qtransfo%type_Qin(:) = Qtransfo%type_Qout(:)
          CONTINUE ! nothing !

         CASE ('order')
          Qtransfo%nb_Qin  = nb_Qin
          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qorder")

          CALL alloc_array(Qtransfo%list_Qin_TO_Qout,(/Qtransfo%nb_Qin/),&
                          "Qtransfo%list_Qin_TO_Qout",name_sub)
          read(in_unitp,*,IOSTAT=err_io) Qtransfo%list_Qin_TO_Qout(1:nb_Qin)
          IF (err_io /= 0) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) '  while reading "list_Qin_TO_Qout"'
             write(out_unitp,*) ' end of file or end of record'
             write(out_unitp,*) ' Check your data !!'
             STOP
          END IF
          IF (print_level > 0) write(out_unitp,*) 'list_Qin_TO_Qout',Qtransfo%list_Qin_TO_Qout(:)
          Qtransfo%type_Qin(:) = 0
          DO i=1,Qtransfo%nb_Qin
            i_Q = Qtransfo%list_Qin_TO_Qout(i)
            IF (i_Q < 0 .OR. i_Q > Qtransfo%nb_Qin) THEN
              Qtransfo%type_Qin(i) = 0
            ELSE
             Qtransfo%type_Qin(i) = Qtransfo%type_Qout(i_Q)
            END IF
          END DO
!          IF (count(Qtransfo%type_Qin(:) == 0) > 0) THEN
!             write(out_unitp,*) ' ERROR in ',name_sub
!             write(out_unitp,*) '  type_Qin "type_Qin"',Qtransfo%type_Qin(:)
!
!             write(out_unitp,*) '  Wrong "list_Qin_TO_Qout"',           &
!                                        Qtransfo%list_Qin_TO_Qout(:)
!             write(out_unitp,*) ' Check your data !!'
!             STOP
!          END IF
        CASE ('linear')
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .FALSE.
          Qtransfo%LinearTransfo%transp = .FALSE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear")
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('linear_transp')
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv    = .FALSE.
          Qtransfo%LinearTransfo%transp = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear")
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('linear_inv')
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear")
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('linear_inv_transp','linear_transp_inv')
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .TRUE.
          Qtransfo%LinearTransfo%transp = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear")
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('lc_projection_inv')
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = .FALSE.
          IF (nb_transfo < 1) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) '  Wrong number of transformation:',nb_transfo
             write(out_unitp,*) '  for the LC_projection_inv transformation'
             write(out_unitp,*) ' Check your data !!'
             STOP
          END IF
          Qtransfo%nb_transfo = nb_transfo
          CALL Read_LC_projectionTransfo(Qtransfo%LinearTransfo,        &
                                  nb_transfo,opt_transfo,not_all,nb_Qin)
          CALL sub_Type_Name_OF_Qin(Qtransfo,"QLCproj")
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('nm')
          Qtransfo%nb_Qin                     = nb_Qin
          Qtransfo%LinearTransfo%inv          = .FALSE.
          Qtransfo%LinearTransfo%check_LinearTransfo = .FALSE.

          CALL alloc_array(Qtransfo%NMTransfo,'Qtransfo%NMTransfo',name_sub)

          Qtransfo%NMTransfo%purify_hess      = purify_hess
          Qtransfo%NMTransfo%eq_hess          = eq_hess
          Qtransfo%NMTransfo%k_Half           = k_Half
          Qtransfo%NMTransfo%hessian_old      = hessian_old
          Qtransfo%NMTransfo%hessian_onthefly = hessian_onthefly
          Qtransfo%NMTransfo%hessian_cart     = hessian_cart
          IF ((hessian_read .OR. k_read) .AND. nb_read < 1) nb_read = 1
          Qtransfo%NMTransfo%hessian_read     = hessian_read
          Qtransfo%NMTransfo%k_read           = k_read
          Qtransfo%NMTransfo%d0c_read         = d0c_read
          Qtransfo%NMTransfo%nb_read          = nb_read

          Qtransfo%NMTransfo%file_hessian%name      = trim(file_hessian)
          Qtransfo%NMTransfo%file_hessian%unit      = 0
          Qtransfo%NMTransfo%file_hessian%formatted = .TRUE.
          Qtransfo%NMTransfo%file_hessian%append    = .FALSE.
          Qtransfo%NMTransfo%file_hessian%old       = hessian_old

          CALL Read_NMTransfo(Qtransfo%NMTransfo,nb_Qin)
          CALL sub_Type_Name_OF_Qin(Qtransfo,"QNM")

          CALL alloc_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)
          CALL mat_id(Qtransfo%LinearTransfo%mat,nb_Qin,nb_Qin)
          CALL mat_id(Qtransfo%LinearTransfo%mat_inv,nb_Qin,nb_Qin)
          CALL Sub_Check_LinearTransfo(Qtransfo)


        CASE ('rph')
          Qtransfo%nb_Qin  = nb_Qin
          CALL alloc_array(Qtransfo%RPHTransfo,'Qtransfo%RPHTransfo',name_sub)
          IF (Qtransfo%opt_transfo /= 0) THEN
             CALL Read_RPHTransfo(Qtransfo%RPHTransfo,nb_Qin,Qtransfo%opt_transfo)
          END IF

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QRPH")
          Qtransfo%type_Qin(:) = 0

        CASE ('hyperspherical')
          Qtransfo%nb_Qin  = nb_Qin
          CALL Read_HyperSpheTransfo(Qtransfo%HyperSpheTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QhyperSphe")
          Qtransfo%type_Qin(:) = Qtransfo%type_Qout(:)
          DO i=2,size(Qtransfo%HyperSpheTransfo%list_HyperSphe)
            Qtransfo%type_Qin(Qtransfo%HyperSpheTransfo%list_HyperSphe(i)) = 4 ! an angle
          END DO
          !Qtransfo%type_Qin(:) = 0

        CASE ('oned')
          Qtransfo%nb_Qin  = nb_Qin
          IF (nb_transfo < 1) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) '  Wrong number of transformation:',nb_transfo
             write(out_unitp,*) '  for the oneD transformation'
             write(out_unitp,*) ' Check your data !!'
             STOP
          END IF
          Qtransfo%nb_transfo = nb_transfo
          CALL Read_oneDTransfo(Qtransfo%oneDTransfo,nb_transfo,nb_Qin)

          IF (.NOT. skip_transfo) THEN
            Qtransfo%opt_param = 0
            DO i=1,nb_transfo
              Qtransfo%opt_param = Qtransfo%opt_param +                 &
                              count(Qtransfo%oneDTransfo(i)%opt_cte > 0)
            END DO
          END IF

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QoneD")
          Qtransfo%type_Qin(:) = 0

        CASE ('threed')
          Qtransfo%nb_Qin  = nb_Qin

          CALL Read_ThreeDTransfo(Qtransfo%ThreeDTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Q3D")
          Qtransfo%type_Qin(:) =0

        CASE ('rot2coord')
          Qtransfo%nb_Qin  = nb_Qin

          CALL Read_Rot2CoordTransfo(Qtransfo%Rot2CoordTransfo,nb_transfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qrot2Coord")
          Qtransfo%type_Qin(:) = 0


        CASE ('flexible')
          Qtransfo%nb_Qin  = nb_Qin
          !CALL Read_FlexibleTransfo(Qtransfo%FlexibleTransfo,nb_Qin)
          CALL Qtransfo%FlexibleTransfo%QtransfoRead(nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qflex")
          Qtransfo%type_Qin(:) = Qtransfo%type_Qout(:)

        CASE ('gene')
          Qtransfo%nb_Qin  = nb_Qin
          CALL Read_GeneTransfo(Qtransfo%GeneTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qgene")
          Qtransfo%type_Qin(:) = 0

        CASE ('active') ! the last read transformation
          Qtransfo%nb_Qin  = nb_Qin
          CALL alloc_array(Qtransfo%ActiveTransfo,'Qtransfo%ActiveTransfo',name_sub)
          IF (Qtransfo%opt_transfo == 1) THEN
            CALL Read2_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
          ELSE
            CALL Read_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
          END IF
          ! the set of type_Qin and name_Qin are done in type_var_analysis

        CASE ('zmat') ! It should be one of the first transfo read
          IF (nat < 2) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) ' nat < 2',nat
              write(out_unitp,*) ' Check your data !!'
              STOP
          END IF
          Qtransfo%Primitive_Coord    = .TRUE.

          Qtransfo%ZmatTransfo%cos_th = cos_th
          Qtransfo%ZmatTransfo%nat0   = nat
          Qtransfo%ZmatTransfo%nat    = nat + 1
          Qtransfo%ZmatTransfo%nb_var = max(1,3*nat-6)
          Qtransfo%ZmatTransfo%ncart  = 3*(nat+1)
          Qtransfo%nb_Qin             = max(1,3*nat-6)
          Qtransfo%nb_Qout            = 3*(nat+1)
          IF (debug) write(out_unitp,*) 'nat0,nat,nb_var,ncart',        &
                     Qtransfo%ZmatTransfo%nat0,Qtransfo%ZmatTransfo%nat,&
                  Qtransfo%ZmatTransfo%nb_var,Qtransfo%ZmatTransfo%ncart

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qzmat")
          Qtransfo%ZmatTransfo%type_Qin => Qtransfo%type_Qin
          Qtransfo%ZmatTransfo%name_Qin => Qtransfo%name_Qin

          CALL Read_ZmatTransfo(Qtransfo%ZmatTransfo,mendeleev)

        CASE ('bunch','bunch_poly') ! It should one of the first transfo

         IF (.NOT. associated(Qtransfo%BunchTransfo)) THEN
           allocate(Qtransfo%BunchTransfo,stat=err_mem)
           memory = 1
           CALL error_memo_allo(err_mem,memory,'Qtransfo%BunchTransfo',name_sub,'Type_BunchTransfo')
         END IF

         IF (nb_vect < 1 .AND. nat > 1) nb_vect = nat-1

          IF (nb_vect < 1) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' nb_vect < 1',nb_vect
             write(out_unitp,*) ' Check your data !!'
             STOP
          END IF
          IF (name_transfo == 'bunch_poly') with_vectors = .FALSE.
          IF (.NOT. with_vectors)       Qtransfo%inTOout = .FALSE.
          Qtransfo%BunchTransfo%nb_vect= nb_vect

          IF (Qtransfo%inTOout) THEN

            nat = 2*nb_vect+1
            Qtransfo%BunchTransfo%nb_G   = nb_G
            Qtransfo%BunchTransfo%nb_X   = nb_X
            Qtransfo%BunchTransfo%nat0   = nat
            Qtransfo%BunchTransfo%nat    = nat + 1
            Qtransfo%BunchTransfo%nb_var = max(1,3*nb_vect-3)
            Qtransfo%BunchTransfo%ncart  = 3*(nat+1)
            Qtransfo%nb_Qin              = 3*nb_vect ! or 3*nat !!
            Qtransfo%nb_Qout             = 3*(nat+1)

            CALL sub_Type_Name_OF_Qin(Qtransfo,"QBunch")
            Qtransfo%BunchTransfo%type_Qin => Qtransfo%type_Qin
            Qtransfo%BunchTransfo%name_Qin => Qtransfo%name_Qin

            IF (debug) write(out_unitp,*) 'nat0,nat,nb_vect,ncart,nb_G',&
                  Qtransfo%BunchTransfo%nat0,Qtransfo%BunchTransfo%nat, &
              Qtransfo%BunchTransfo%nb_vect,Qtransfo%BunchTransfo%ncart,&
                              Qtransfo%BunchTransfo%nb_G

            CALL Read_BunchTransfo(Qtransfo%BunchTransfo,mendeleev)

          ELSE

            Qtransfo%BunchTransfo%nat_act = nb_vect + 1
            nat = Qtransfo%BunchTransfo%nat_act + nb_G + nb_X
            Qtransfo%BunchTransfo%nb_G   = nb_G
            Qtransfo%BunchTransfo%nb_X   = nb_X
            Qtransfo%BunchTransfo%nat0   = nat
            Qtransfo%BunchTransfo%nat    = nat + 1
            Qtransfo%BunchTransfo%nb_var = max(1,3*nb_vect-3)
            Qtransfo%BunchTransfo%ncart  = 3*(nat+1)
            Qtransfo%nb_Qin              = 3*nb_vect
            Qtransfo%nb_Qout             = 3*(nat+1)

            CALL sub_Type_Name_OF_Qin(Qtransfo,"QBunch")
            Qtransfo%BunchTransfo%type_Qin => Qtransfo%type_Qin
            Qtransfo%BunchTransfo%name_Qin => Qtransfo%name_Qin

            IF (debug) write(out_unitp,*) 'nat0,nat,nb_vect,ncart,nb_G',&
                  Qtransfo%BunchTransfo%nat0,Qtransfo%BunchTransfo%nat, &
              Qtransfo%BunchTransfo%nb_vect,Qtransfo%BunchTransfo%ncart,&
                              Qtransfo%BunchTransfo%nb_G

            CALL Read2_BunchTransfo(Qtransfo%BunchTransfo,mendeleev,with_vectors)

            IF (with_vectors) THEN
              CALL M_Tana_FROM_Bunch2Transfo(Qtransfo%BunchTransfo)
            END IF
          END IF

        CASE ('poly')
          IF ( .NOT. associated(Qtransfo%BunchTransfo)) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) 'For Poly transfo, ... '
            write(out_unitp,*) ' Qtransfo%BunchTransfo MUST be associoted TO'
            write(out_unitp,*) ' mole%tab_Qtransfo(1)%BunchTransfo.'

            write(out_unitp,*) ' Check the fortran !!'
            STOP
          END IF

          Qtransfo%Primitive_Coord      = .TRUE.
          nb_vect                       = nb_Qin/3 + 1
          Qtransfo%BFTransfo%nb_var     = nb_Qin
          Qtransfo%BFTransfo%Def_cos_th = cos_th
          Qtransfo%nb_Qin               = nb_Qin
          Qtransfo%nb_Qout              = 3*nb_vect
          CALL alloc_array(Qtransfo%type_Qin,(/Qtransfo%nb_Qin/),       &
                          "Qtransfo%type_Qin",name_sub)
          CALL alloc_array(Qtransfo%name_Qin,(/Qtransfo%nb_Qin/),Name_len,&
                          "Qtransfo%name_Qin",name_sub)


          Qtransfo%BFTransfo%type_Qin => Qtransfo%type_Qin
          Qtransfo%BFTransfo%name_Qin => Qtransfo%name_Qin
          Qtransfo%BFTransfo%Frame    = .TRUE.
          Qtransfo%BFTransfo%euler(:) = .FALSE.
          iF_inout = 0
          CALL RecRead_BFTransfo(Qtransfo%BFTransfo,                    &
                                 Qtransfo%BunchTransfo,iF_inout)

          IF (debug) THEN
            write(out_unitp,*) ' Type and name of polyspherical coordinates'
            DO i=1,Qtransfo%nb_Qin
              write(out_unitp,*) 'i,type,name',i,Qtransfo%type_Qin(i),Qtransfo%name_Qin(i)
            END DO

            CALL RecWrite_BFTransfo(Qtransfo%BFTransfo,.TRUE.)
          END IF

          ! Calculation of M_Tana if needed
          IF (count(Qtransfo%BunchTransfo%M_Tana /= ZERO) == 0) THEN
            CALL M_Tana_FROM_Bunch2Transfo(Qtransfo%BunchTransfo)
          END IF
          nullify(Qtransfo%BunchTransfo)

        CASE ('qtox_ana')
          Qtransfo%Primitive_Coord    = .TRUE.
          Qtransfo%nb_Qin             = max(1,3*nat-6)
          Qtransfo%nb_Qout            = 3*nat+3
          CALL alloc_array(Qtransfo%type_Qin,(/Qtransfo%nb_Qin/),       &
                          "Qtransfo%type_Qin",name_sub)
          CALL alloc_array(Qtransfo%name_Qin,(/Qtransfo%nb_Qin/),Name_len,&
                          "Qtransfo%name_Qin",name_sub)
          DO i=1,Qtransfo%nb_Qin
            CALL make_nameQ(Qtransfo%name_Qin(i),"Qana",i,it)
          END DO
          Qtransfo%QTOXanaTransfo%type_Qin => Qtransfo%type_Qin

          Qtransfo%QTOXanaTransfo%nat0      = nat
          Qtransfo%QTOXanaTransfo%nat       = nat + 1
          Qtransfo%QTOXanaTransfo%nat_act   = nat
          Qtransfo%QTOXanaTransfo%nb_var    = max(1,3*nat-6)
          Qtransfo%QTOXanaTransfo%ncart     = 3*(nat+1)
          Qtransfo%QTOXanaTransfo%ncart_act = 3*nat

          IF (debug) write(out_unitp,*) 'nat0,nat,nb_var,ncart',        &
                                         Qtransfo%QTOXanaTransfo%nat0,  &
                                         Qtransfo%QTOXanaTransfo%nat,   &
                                         Qtransfo%QTOXanaTransfo%nb_var,&
                                         Qtransfo%QTOXanaTransfo%ncart

          CALL Read_QTOXanaTransfo(Qtransfo%QTOXanaTransfo,mendeleev)

        CASE ('cartesian') ! It should be one of the first transfo read
          Qtransfo%nb_Qin             = nb_Qin ! ncart_act
          Qtransfo%nb_Qout            = nb_Qin ! ncart_act
          CALL alloc_array(Qtransfo%type_Qin,(/Qtransfo%nb_Qin/),       &
                          "Qtransfo%type_Qin",name_sub)
          CALL alloc_array(Qtransfo%name_Qin,(/Qtransfo%nb_Qin/),Name_len,&
                          "Qtransfo%name_Qin",name_sub)
          DO i=1,Qtransfo%nb_Qin
            CALL make_nameQ(Qtransfo%name_Qin(i),"Qxyz_transfo",i,it)
          END DO

          CALL Read_CartesianTransfo(Qtransfo%CartesianTransfo)

        CASE default ! ERROR: wrong transformation !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The transformation is UNKNOWN: ',        &
                                                    trim(name_transfo)
          CALL Write_list_Qtransfo(out_unitp)
          STOP

        END SELECT


        ! for Qout type, name ....
        IF (Qtransfo%num_transfo == 1) THEN ! cartessian coordinates (Qout)
          CALL alloc_array(Qtransfo%type_Qout,(/Qtransfo%nb_Qout/),     &
                          "Qtransfo%type_Qout",name_sub)
          CALL alloc_array(Qtransfo%name_Qout,(/Qtransfo%nb_Qout/),Name_len,&
                          "Qtransfo%name_Qout",name_sub)
          Qtransfo%type_Qout(:) = 1 ! cartesian type
          DO i=1,Qtransfo%nb_Qout
            iat = (i-1)/3 +1
            IF (mod(i,3) == 1) CALL make_nameQ(Qtransfo%name_Qout(i),"X",iat,0)
            IF (mod(i,3) == 2) CALL make_nameQ(Qtransfo%name_Qout(i),"Y",iat,0)
            IF (mod(i,3) == 0) CALL make_nameQ(Qtransfo%name_Qout(i),"Z",iat,0)
          END DO
        END IF


        IF (debug) CALL Write_QTransfo(Qtransfo)
        write(out_unitp,*) '=========================================='
        write(out_unitp,*) '=========================================='

      END SUBROUTINE read_Qtransfo

      SUBROUTINE Write_list_Qtransfo(nio)
        integer, intent(in) :: nio

        write(nio,*) ' The possible coordinate transformations are:'
        write(nio,*)
        write(nio,*) '"zmat"'
        write(nio,*) '"Rec_NM"'
        write(nio,*) '"QTOX_ana"'
        write(nio,*) '"bunch_poly"'
        write(nio,*) '"bunch"'

        write(nio,*)
        write(nio,*) '"poly"'

        write(nio,*)
        write(nio,*) '"identity"'

        write(nio,*) '"linear"'
        write(nio,*) '"linear_transp"'
        write(nio,*) '"linear_inv"'
        write(nio,*) '"linear_inv_transp" or "linear_transp_inv"'
        write(nio,*) '"LC_projection_inv"'

        write(nio,*) '"hyperspherical"'
        write(nio,*) '"flexible"'
        write(nio,*) '"gene"'
        write(nio,*) '"order"'
        write(nio,*) '"oneD"'
        write(nio,*) '"ThreeD"'
        write(nio,*) '"Rot2Coord"'

        write(nio,*)
        write(nio,*) '"NM"'
        write(nio,*) '"RPH"'
        write(nio,*)
        write(nio,*) '"active"'

        write(nio,*) ' Special transformation:'
        write(nio,*) '"Cartesian"'


       !write(nio,*) '""'
        CALL flush_perso(nio)
      END SUBROUTINE Write_list_Qtransfo

      SUBROUTINE dealloc_Qtransfo(Qtransfo)
        TYPE (Type_Qtransfo), intent(inout) :: Qtransfo

        character (len=Name_len) :: name_transfo
        character (len=*),parameter :: name_sub='dealloc_Qtransfo'
        !logical, parameter :: debug = .TRUE.
        logical, parameter :: debug = .FALSE.

        name_transfo = Qtransfo%name_transfo
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING : ',name_sub,' : ',name_transfo
          CALL flush_perso(out_unitp)
        END IF

        Qtransfo%print_done      = .FALSE.
        Qtransfo%name_transfo    = "identity"
        Qtransfo%num_transfo     = 0
        Qtransfo%opt_transfo     = 0
        Qtransfo%nb_var          = 0
        Qtransfo%nb_act          = 0
        Qtransfo%nb_transfo      = 0
        Qtransfo%skip_transfo    = .FALSE.
        Qtransfo%opt_param       = 0
        Qtransfo%Primitive_Coord = .FALSE.



        ! ==== LinearTransfo ========================
        CALL dealloc_LinearTransfo(Qtransfo%LinearTransfo)

        ! ==== NMTransfo ========================
        IF (associated(Qtransfo%NMTransfo)) THEN
          CALL dealloc_NMTransfo(Qtransfo%NMTransfo)
          CALL dealloc_array(Qtransfo%NMTransfo,'Qtransfo%NMTransfo',name_sub)
        END IF

        ! ==== FlexibleTransfo ========================
        CALL dealloc_FlexibleTransfo(Qtransfo%FlexibleTransfo)

        ! ==== RPHTransfo ========================
        IF (associated(Qtransfo%RPHTransfo)) THEN
          CALL dealloc_RPHTransfo(Qtransfo%RPHTransfo)
          CALL dealloc_array(Qtransfo%RPHTransfo,'Qtransfo%RPHTransfo',name_sub)
        END IF

        ! ==== geneTransfo ========================
        CALL dealloc_GeneTransfo(Qtransfo%GeneTransfo)

        ! ==== HyperSpheTransfo ========================
        Qtransfo%HyperSpheTransfo%nb_HyperSphe = 0
        IF (associated(Qtransfo%HyperSpheTransfo%list_HyperSphe) ) THEN
          CALL dealloc_array(Qtransfo%HyperSpheTransfo%list_HyperSphe,  &
                            "Qtransfo%HyperSpheTransfo%list_HyperSphe",name_sub)
        END IF

        ! ==== oneDTransfo ========================
        CALL dealloc_oneDTransfo(Qtransfo%oneDtransfo)

        ! ==== ThreeDTransfo ========================
        CALL dealloc_ThreeDTransfo(Qtransfo%ThreeDTransfo)

        ! ==== Rot2CoordTransfo ========================
        CALL dealloc_Rot2CoordTransfo(Qtransfo%Rot2CoordTransfo)

        ! ==== ZmatTransfo ========================
        CALL dealloc_ZmatTransfo(Qtransfo%ZmatTransfo)
        ! ==== RectilinearNM_Transfo ================
        CALL dealloc_RectilinearNM_Transfo(Qtransfo%RectilinearNM_Transfo)
        ! ==== BunchTransfo ========================
        CALL dealloc_BunchTransfo(Qtransfo%BunchTransfo)
        ! ==== BFTransfo ===========================
        CALL dealloc_BFTransfo(Qtransfo%BFTransfo)
        ! ==== QTOXanaTransfoTransfo ===============
        CALL dealloc_QTOXanaTransfo(Qtransfo%QTOXanaTransfo)

        ! ==== orderTransfo ===========================
        IF (associated(Qtransfo%list_Qin_TO_Qout) ) THEN
          CALL dealloc_array(Qtransfo%list_Qin_TO_Qout,                 &
                            "Qtransfo%list_Qin_TO_Qout",name_sub)
        END IF

        ! ==== activeTransfo ===========================
        IF (associated(Qtransfo%ActiveTransfo)) THEN
          CALL dealloc_ActiveTransfo(Qtransfo%ActiveTransfo)
          CALL dealloc_array(Qtransfo%ActiveTransfo,'Qtransfo%ActiveTransfo',name_sub)
        END IF

        ! ==== CartesianTransfo ========================
        CALL dealloc_CartesianTransfo(Qtransfo%CartesianTransfo)


        ! ==== Coordinates ========================
        Qtransfo%nb_Qin  = 0
        Qtransfo%nb_Qout = 0

        IF (associated(Qtransfo%type_Qin) ) THEN
          CALL dealloc_array(Qtransfo%type_Qin,                         &
                            "Qtransfo%type_Qin",name_sub)
        END IF
        nullify(Qtransfo%type_Qout) ! because it is a true pointer

        IF (associated(Qtransfo%name_Qin) ) THEN
          CALL dealloc_array(Qtransfo%name_Qin,                         &
                            "Qtransfo%name_Qin",name_sub)
        END IF
        nullify(Qtransfo%name_Qout)  ! because it is a true pointer

        IF (debug) THEN
          write(out_unitp,*) 'END : ',name_sub,' : ',name_transfo
          CALL flush_perso(out_unitp)
        END IF
      END SUBROUTINE dealloc_Qtransfo

      SUBROUTINE alloc_array_OF_Qtransfodim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_Qtransfo), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Qtransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_Qtransfo')

      END SUBROUTINE alloc_array_OF_Qtransfodim1
      SUBROUTINE dealloc_array_OF_Qtransfodim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_Qtransfo), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Qtransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_Qtransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Qtransfodim1
      SUBROUTINE Qtransfo1TOQtransfo2(Qtransfo1,Qtransfo2)
        TYPE (Type_Qtransfo), intent(in)    :: Qtransfo1
        TYPE (Type_Qtransfo), intent(inout) :: Qtransfo2
        integer :: it,n
        character (len=Name_len) :: name_transfo
!     -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Qtransfo1TOQtransfo2'
!     -----------------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'name_transfo: ',Qtransfo1%name_transfo
        CALL Write_Qtransfo(Qtransfo1)
        CALL flush_perso(out_unitp)
      END IF
!     -----------------------------------------------------------------
      Qtransfo2%print_done      = .FALSE.

      Qtransfo2%name_transfo    = Qtransfo1%name_transfo
      Qtransfo2%inTOout         = Qtransfo1%inTOout

      Qtransfo2%nb_var          = Qtransfo1%nb_var
      Qtransfo2%nb_act          = Qtransfo1%nb_act
      Qtransfo2%nb_transfo      = Qtransfo1%nb_transfo

      Qtransfo2%nb_Qin          = Qtransfo1%nb_Qin
      Qtransfo2%nb_Qout         = Qtransfo1%nb_Qout

      Qtransfo2%num_transfo     = Qtransfo1%num_transfo

      Qtransfo2%opt_transfo     = Qtransfo1%opt_transfo
      Qtransfo2%skip_transfo    = Qtransfo1%skip_transfo

      Qtransfo2%opt_param       = Qtransfo1%opt_param
      Qtransfo2%Primitive_Coord = Qtransfo1%Primitive_Coord


      CALL alloc_array(Qtransfo2%type_Qin,shape(Qtransfo1%type_Qin),    &
                      "Qtransfo2%type_Qin",name_sub)
      Qtransfo2%type_Qin(:) = Qtransfo1%type_Qin(:)

      CALL alloc_array(Qtransfo2%name_Qin,shape(Qtransfo1%name_Qin),Name_len,&
                      "Qtransfo2%name_Qin",name_sub)
      Qtransfo2%name_Qin(:) = Qtransfo1%name_Qin(:)

      ! for type_Qout and name_Qout, it will be done after (from another type_Qin, name_Qin)
      ! except for num_transfo=0
      IF (Qtransfo2%num_transfo == 0) THEN
        CALL alloc_array(Qtransfo2%type_Qout,shape(Qtransfo1%type_Qout),  &
                        "Qtransfo2%type_Qout",name_sub)
        Qtransfo2%type_Qout(:) = Qtransfo1%type_Qout(:)

        CALL alloc_array(Qtransfo2%name_Qout,shape(Qtransfo1%name_Qout),Name_len,&
                        "Qtransfo2%name_Qout",name_sub)
        Qtransfo2%name_Qout(:) = Qtransfo1%name_Qout(:)
      END IF

      name_transfo = Qtransfo2%name_transfo
      CALL string_uppercase_TO_lowercase(name_transfo)

      SELECT CASE (name_transfo)
      CASE ('identity')
        CONTINUE ! nothing !

      CASE ('order')
        n = size(Qtransfo1%list_Qin_TO_Qout)
        CALL alloc_array(Qtransfo2%list_Qin_TO_Qout,                    &
                                     shape(Qtransfo1%list_Qin_TO_Qout), &
                        "Qtransfo2%list_Qin_TO_Qout",name_sub)
        Qtransfo2%list_Qin_TO_Qout(:) = Qtransfo1%list_Qin_TO_Qout(:)

      CASE ('linear','linear_inv','lc_projection_inv',                  &
            'linear_transp','linear_transp_inv','linear_inv_transp')
        n = size(Qtransfo1%LinearTransfo%mat,dim=1)
        CALL alloc_LinearTransfo(Qtransfo2%LinearTransfo,n)
        Qtransfo2%LinearTransfo%mat     = Qtransfo1%LinearTransfo%mat
        Qtransfo2%LinearTransfo%mat_inv = Qtransfo1%LinearTransfo%mat_inv
        Qtransfo2%LinearTransfo%inv     = Qtransfo1%LinearTransfo%inv
        Qtransfo2%LinearTransfo%transp  = Qtransfo1%LinearTransfo%transp
        Qtransfo2%LinearTransfo%check_LinearTransfo = Qtransfo1%LinearTransfo%check_LinearTransfo

      CASE ('nm')
        IF (associated(Qtransfo1%LinearTransfo%mat)) THEN
          n = size(Qtransfo1%LinearTransfo%mat,dim=1)
          CALL alloc_LinearTransfo(Qtransfo2%LinearTransfo,n)
          Qtransfo2%LinearTransfo%mat = Qtransfo1%LinearTransfo%mat
          Qtransfo2%LinearTransfo%mat_inv = Qtransfo1%LinearTransfo%mat_inv
          Qtransfo2%LinearTransfo%inv = Qtransfo1%LinearTransfo%inv
          Qtransfo2%LinearTransfo%transp  = Qtransfo1%LinearTransfo%transp
          Qtransfo2%LinearTransfo%check_LinearTransfo = Qtransfo1%LinearTransfo%check_LinearTransfo
        END IF
        IF (associated(Qtransfo1%NMTransfo)) THEN
          CALL alloc_array(Qtransfo2%NMTransfo,                         &
                          'Qtransfo2%NMTransfo',name_sub)
          CALL NMTransfo1TONMTransfo2(Qtransfo1%NMTransfo,Qtransfo2%NMTransfo)
        END IF

      CASE ('rph')
        IF (associated(Qtransfo1%RPHTransfo)) THEN
          CALL alloc_array(Qtransfo2%RPHTransfo,                        &
                          'Qtransfo2%RPHTransfo',name_sub)
          CALL RPHTransfo1TORPHTransfo2(Qtransfo1%RPHTransfo,           &
                                        Qtransfo2%RPHTransfo)
        END IF

      CASE ('hyperspherical')
        Qtransfo2%HyperSpheTransfo%nb_HyperSphe =                       &
                         Qtransfo1%HyperSpheTransfo%nb_HyperSphe

        CALL alloc_array(Qtransfo2%HyperSpheTransfo%list_HyperSphe,     &
                       shape(Qtransfo1%HyperSpheTransfo%list_HyperSphe),&
                        "Qtransfo2%HyperSpheTransfo%list_HyperSphe",name_sub)
        Qtransfo2%HyperSpheTransfo%list_HyperSphe =                     &
                              Qtransfo1%HyperSpheTransfo%list_HyperSphe

      CASE ('oned')
        IF (Qtransfo2%nb_transfo < 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  Wrong number of transformation:',       &
                                           Qtransfo2%nb_transfo
          write(out_unitp,*) '  for the oneD transformation'
          write(out_unitp,*) ' Check the fortran source !!'
          STOP
        END IF

        CALL oneDTransfo1TOoneDTransfo2(Qtransfo1%oneDTransfo,          &
                                        Qtransfo2%oneDTransfo)

      CASE ('threed')
        CALL ThreeDTransfo1TOThreeDTransfo2(Qtransfo1%ThreeDTransfo,    &
                                            Qtransfo2%ThreeDTransfo)

      CASE ('rot2coord')
        CALL Rot2CoordTransfo1TORot2CoordTransfo2(                      &
                  Qtransfo1%Rot2CoordTransfo,Qtransfo2%Rot2CoordTransfo)

      CASE ('flexible')
        CALL FlexibleTransfo1TOFlexibleTransfo2(                        &
                    Qtransfo1%FlexibleTransfo,Qtransfo2%FlexibleTransfo)

      CASE ('gene')
        CALL GeneTransfo1TOGeneTransfo2(Qtransfo1%GeneTransfo,          &
                                        Qtransfo2%GeneTransfo)

      CASE ('active')
        IF (associated(Qtransfo1%ActiveTransfo)) THEN
          CALL alloc_array(Qtransfo2%ActiveTransfo,                     &
                          "Qtransfo2%ActiveTransfo",name_sub)
          CALL ActiveTransfo1TOActiveTransfo2(Qtransfo1%ActiveTransfo,  &
                                              Qtransfo2%ActiveTransfo)
        END IF

      CASE ('zmat')
        CALL ZmatTransfo1TOZmatTransfo2(Qtransfo1%ZmatTransfo,          &
                                        Qtransfo2%ZmatTransfo)
        Qtransfo2%ZmatTransfo%type_Qin => Qtransfo2%type_Qin
        Qtransfo2%ZmatTransfo%name_Qin => Qtransfo2%name_Qin

      CASE ('rec_nm')
        CALL RectilinearNM_Transfo1TORectilinearNM_Transfo2(            &
                                       Qtransfo1%RectilinearNM_Transfo, &
                                       Qtransfo2%RectilinearNM_Transfo)
        Qtransfo2%ZmatTransfo%type_Qin => Qtransfo2%type_Qin
        Qtransfo2%ZmatTransfo%name_Qin => Qtransfo2%name_Qin

      CASE ('bunch','bunch_poly')
        CALL BunchTransfo1TOBunchTransfo2(Qtransfo1%BunchTransfo,       &
                                          Qtransfo2%BunchTransfo)
        Qtransfo2%BunchTransfo%type_Qin => Qtransfo2%type_Qin
        Qtransfo2%BunchTransfo%name_Qin => Qtransfo2%name_Qin

      CASE ('poly')
        CALL Rec_BFTransfo1TOBFTransfo2(Qtransfo1%BFTransfo,            &
                                        Qtransfo2%BFTransfo)
        Qtransfo2%BFTransfo%type_Qin => Qtransfo2%type_Qin
        Qtransfo2%BFTransfo%name_Qin => Qtransfo2%name_Qin

      CASE ('qtox_ana')
        CALL QTOXanaTransfo1TOQTOXanaTransfo2(Qtransfo1%QTOXanaTransfo, &
                                              Qtransfo2%QTOXanaTransfo)
        Qtransfo2%QTOXanaTransfo%type_Qin => Qtransfo2%type_Qin

      CASE ('cartesian')
        CALL CartesianTransfo1TOCartesianTransfo2(                      &
                  Qtransfo1%CartesianTransfo,Qtransfo2%CartesianTransfo)

      CASE default
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The transformation is UNKNOWN: ',          &
                                     trim(Qtransfo1%name_transfo)
        CALL Write_list_Qtransfo(out_unitp)
        write(out_unitp,*) ' Check the source!'
        STOP
      END SELECT

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!     -----------------------------------------------------------------
      END SUBROUTINE Qtransfo1TOQtransfo2
      SUBROUTINE calc_Qtransfo(dnQin,dnQout,Qtransfo,nderiv,inTOout)

        TYPE (Type_dnVec),    intent(inout)        :: dnQin,dnQout
        TYPE (Type_Qtransfo), intent(in)           :: Qtransfo
        integer,              intent(in)           :: nderiv
        logical,              intent(in), optional :: inTOout

        logical           :: inTOout_loc
        TYPE (Type_dnS)   :: dnR
        integer           :: iv,it,i,iQ,iQin,iQout
        TYPE (Type_dnVec), pointer :: tab_dnXVect(:)   ! dim: nb_vect_tot
        character (len=Name_len) :: name_transfo


!     -----------------------------------------------------------------
      integer :: nderiv_debug = 1
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_Qtransfo'
!     -----------------------------------------------------------------

      IF (present(inTOout)) THEN
        inTOout_loc = inTOout
      ELSE
        inTOout_loc = .TRUE.
      END IF

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'New Qtransfo',it,Qtransfo%name_transfo
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'Qtransfo%nb_act',Qtransfo%nb_act
        write(out_unitp,*) 'inTOout',inTOout_loc
        write(out_unitp,*) 'nb_Qin,nb_Qout',Qtransfo%nb_Qin,Qtransfo%nb_Qout

        IF (inTOout_loc) THEN
          write(out_unitp,*) 'dnOin :'
          CALL Write_dnSVM(dnQin,nderiv_debug)
        ELSE
          write(out_unitp,*) 'dnOout :'
          CALL Write_dnSVM(dnQout,nderiv_debug)
        END IF
      END IF
!     -----------------------------------------------------------------
      IF (inTOout_loc) THEN
        CALL alloc_dnSVM(dnQout,Qtransfo%nb_Qout,Qtransfo%nb_act,nderiv)
      ELSE
        CALL alloc_dnSVM(dnQin,Qtransfo%nb_Qin,Qtransfo%nb_act,nderiv)
      END IF

      name_transfo = Qtransfo%name_transfo
      CALL string_uppercase_TO_lowercase(name_transfo)
      IF (Qtransfo%skip_transfo) name_transfo = 'identity' ! it should be done before

      SELECT CASE (name_transfo)
      CASE ('identity')
        IF (inTOout_loc) THEN
          CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
        ELSE
          CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
        END IF

      CASE ('order')
        CALL alloc_dnSVM(dnR,dnQout%nb_var_deriv,nderiv)
        IF (inTOout_loc) THEN
          DO iQin=1,Qtransfo%nb_Qin
            CALL sub_dnVec_TO_dnS(dnQin,dnR,iQin,nderiv)
            iQout = Qtransfo%list_Qin_TO_Qout(iQin)
            CALL sub_dnS_TO_dnVec(dnR,dnQout,iQout,nderiv)
          END DO
        ELSE
          DO iQin=1,Qtransfo%nb_Qin
            iQout = Qtransfo%list_Qin_TO_Qout(iQin)
            CALL sub_dnVec_TO_dnS(dnQout,dnR,iQout,nderiv)
            CALL sub_dnS_TO_dnVec(dnR,dnQin,iQin,nderiv)
          END DO
        END IF
        CALL dealloc_dnSVM(dnR)

      CASE ('linear','linear_inv','lc_projection_inv',                  &
            'linear_transp','linear_transp_inv','linear_inv_transp','nm')
        CALL calc_LinearTransfo(dnQin,dnQout,Qtransfo%LinearTransfo,    &
                                                     nderiv,inTOout_loc)

      CASE ('rph')
        IF (associated(Qtransfo%RPHTransfo)) THEN
            CALL calc_RPHTransfo(dnQin,dnQout,Qtransfo%RPHTransfo,      &
                                                     nderiv,inTOout_loc)
        ELSE
          IF (inTOout_loc) THEN
            CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
          ELSE
            CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
          END IF
        END IF

      CASE ('hyperspherical')
        CALL calc_HyperSpheTransfo(dnQin,dnQout,Qtransfo%HyperSpheTransfo,&
                                   nderiv,inTOout_loc)

      CASE ('oned')
        CALL calc_oneDTransfo(dnQin,dnQout,Qtransfo%oneDTransfo,        &
                              nderiv,inTOout_loc)

      CASE ('threed')
        CALL calc_ThreeDTransfo(dnQin,dnQout,Qtransfo%ThreeDTransfo,    &
                                nderiv,inTOout_loc)

      CASE ('rot2coord')
        CALL calc_Rot2CoordTransfo(dnQin,dnQout,                        &
                                             Qtransfo%Rot2CoordTransfo, &
                                                     nderiv,inTOout_loc)

      CASE ('flexible')
        CALL calc_FlexibleTransfo(dnQin,dnQout,Qtransfo%FlexibleTransfo,&
                                  nderiv,inTOout_loc)

      CASE ('gene')
        CALL calc_GeneTransfo(dnQin,dnQout,Qtransfo%GeneTransfo,        &
                                  nderiv,inTOout_loc)

      CASE ('active') ! it has to be the first one, but the last one read
        CALL calc_ActiveTransfo(dnQin,dnQout,Qtransfo%ActiveTransfo,    &
                                                     nderiv,inTOout_loc)

      CASE ('zmat') ! it can be one of the last one
        IF (inTOout_loc) THEN
          CALL calc_ZmatTransfo(dnQin,dnQout,Qtransfo%ZmatTransfo,nderiv)
        ELSE
          CALL calc_ZmatTransfo_outTOin(dnQin,dnQout,Qtransfo%ZmatTransfo,nderiv)
        END IF

      CASE ('rec_nm')
        CALL calc_RectilinearNM_Transfo(dnQin,dnQout,                   &
                                        Qtransfo%RectilinearNM_Transfo, &
                                                     nderiv,inTOout_loc)

      CASE ('bunch','bunch_poly') ! it has to be one of the last one
          CALL calc_BunchTransfo(dnQin,dnQout,Qtransfo%BunchTransfo,    &
                                 nderiv,inTOout_loc)

      CASE ('poly')
        IF (inTOout_loc) THEN

          ! initialization : allocation....
          nullify(tab_dnXVect)
          CALL alloc_array(tab_dnXVect,(/Qtransfo%BFTransfo%nb_vect_tot/),&
                          "tab_dnXVect",name_sub)
          DO iv=1,Qtransfo%BFTransfo%nb_vect_tot
            CALL alloc_dnSVM(tab_dnXVect(iv),3,dnQin%nb_var_deriv,nderiv)
          END DO

          iQin = 0
          CALL calc_PolyTransfo(dnQin,iQin,dnQout,tab_dnXVect,0,        &
                                Qtransfo%BFTransfo,nderiv)


          DO iv=1,Qtransfo%BFTransfo%nb_vect_tot
            CALL dealloc_dnSVM(tab_dnXVect(iv))
          END DO
          CALL dealloc_array(tab_dnXVect,"tab_dnXVect",name_sub)

        ELSE
          CALL calc_PolyTransfo_outTOin(dnQin,dnQout,                   &
                                        Qtransfo%BFTransfo,nderiv)
        END IF

      CASE ('qtox_ana') ! it has to be one of the last one
        IF (nderiv > 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' nderiv MUST = 0',nderiv
           write(out_unitp,*) ' USED, num_x=t and num_g=t in'
           write(out_unitp,*) ' the namelist "variables" or "geom"'
           STOP
        END IF
        CALL Q_TO_X_ana(dnQin%d0, size(dnQin%d0),                       &
                        dnQout%d0,size(dnQout%d0),inTOout_loc)

      CASE ('cartesian') ! it has to be one of the last one
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' Do NOT use this subroutine'
         write(out_unitp,*) ' CALL directly "calc_CartesianTransfo_new"'
         STOP
        !CALL calc_CartesianTransfo(dnQin,dnQout,Qtransfo%CartesianTransfo,&
        !                           nderiv,inTOout_loc)

      CASE default
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The transformation is UNKNOWN: ',          &
                                   trim(Qtransfo%name_transfo)
        CALL Write_list_Qtransfo(out_unitp)
        write(out_unitp,*) ' Check the source!'
        STOP
      END SELECT

!     -----------------------------------------------------------------
      IF (debug) THEN
        IF (inTOout_loc) THEN
          write(out_unitp,*) 'dnOout :'
          CALL Write_dnSVM(dnQout,nderiv_debug)
        ELSE
          write(out_unitp,*) 'dnOin :'
          CALL Write_dnSVM(dnQin,nderiv_debug)
        END IF
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE calc_Qtransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_Qtransfo(Qtransfo,force_print)
        TYPE (Type_Qtransfo) :: Qtransfo
        logical, optional    :: force_print

        character (len=Name_len) :: name_transfo,name_dum
        integer :: nat,nb_var,nb_vect,nbcol,nb_flex_act
        integer :: err
        integer :: i,it,i_Q
        logical :: force_print_loc

        character (len=*), parameter :: name_sub = "Write_Qtransfo"


        IF (present(force_print)) THEN
          force_print_loc = force_print
        ELSE
          force_print_loc = .FALSE.
        END IF

        IF (Qtransfo%print_done .AND. .NOT. force_print_loc) RETURN
        write(out_unitp,*) 'BEGINNING ',name_sub

        Qtransfo%print_done = .TRUE.

        write(out_unitp,*) 'name_transfo,num_transfo: ',     &
                       trim(Qtransfo%name_transfo),Qtransfo%num_transfo

        write(out_unitp,*) 'Primitive_Coord: ',Qtransfo%Primitive_Coord

        write(out_unitp,*) ' Option of the transfo: ',Qtransfo%opt_transfo
        write(out_unitp,*) ' Skip the transfo: ',Qtransfo%skip_transfo

        write(out_unitp,*) ' Parameter(s) to be optimized?: ',Qtransfo%opt_param

        write(out_unitp,*) 'nb_var,nb_act',                             &
                         Qtransfo%nb_var,Qtransfo%nb_act
        write(out_unitp,*) 'nb_Qin,nb_Qout',                            &
                         Qtransfo%nb_Qin,Qtransfo%nb_Qout

        CALL flush_perso(out_unitp)
        write(out_unitp,*) '---------------------------------------'
        IF (associated(Qtransfo%name_Qout) .AND. associated(Qtransfo%type_Qout)) THEN
          DO i_Q=1,Qtransfo%nb_Qout
            write(out_unitp,*) 'i_Q,name_Qout,type_Qout',i_Q," ",       &
                   trim(Qtransfo%name_Qout(i_Q)),                       &
                   Qtransfo%type_Qout(i_Q)
            CALL flush_perso(out_unitp)

          END DO
        ELSE
          write(out_unitp,*) 'asso name_Qout and type_Qout',            &
           associated(Qtransfo%name_Qout),associated(Qtransfo%type_Qout)
        END IF

        IF (associated(Qtransfo%name_Qin) .AND. associated(Qtransfo%type_Qin)) THEN
          write(out_unitp,*) '---------------------------------------'
          DO i_Q=1,Qtransfo%nb_Qin
            write(out_unitp,*) 'i_Q,name_Qin,type_Qin',i_Q," ",         &
                   trim(Qtransfo%name_Qin(i_Q)),                        &
                   Qtransfo%type_Qin(i_Q)
          END DO
        ELSE
          write(out_unitp,*) 'asso name_Qin and type_Qin',              &
           associated(Qtransfo%name_Qin),associated(Qtransfo%type_Qin)
        END IF
        write(out_unitp,*) '---------------------------------------'

        name_transfo = Qtransfo%name_transfo
        CALL string_uppercase_TO_lowercase(name_transfo)

        SELECT CASE (name_transfo)
        CASE ('identity')
          CONTINUE ! nothing !

        CASE ('order')
           write(out_unitp,*) 'list_Qin_TO_Qout',Qtransfo%list_Qin_TO_Qout(:)

        CASE ('linear','linear_inv','lc_projection_inv',                &
            'linear_transp','linear_transp_inv','linear_inv_transp')
          write(out_unitp,*)  'Mat of LinearTransfo: '
          CALL Write_Mat(Qtransfo%LinearTransfo%mat,out_unitp,4)

          write(out_unitp,*)  'Mat_inv of LinearTransfo: '
          CALL Write_Mat(Qtransfo%LinearTransfo%mat_inv,out_unitp,4)

        CASE ('nm')
          IF (associated(Qtransfo%NMTransfo)) THEN
            CALL Write_NMTransfo(Qtransfo%NMTransfo)
          END IF
          IF (associated(Qtransfo%LinearTransfo%mat)) THEN
            write(out_unitp,*)  'Mat of LinearTransfo (NM): '
            CALL Write_Mat(Qtransfo%LinearTransfo%mat,out_unitp,4)
          END IF
          IF (associated(Qtransfo%LinearTransfo%mat_inv)) THEN
            write(out_unitp,*)  'Mat_inv of LinearTransfo (NM): '
            CALL Write_Mat(Qtransfo%LinearTransfo%mat_inv,out_unitp,4)
          END IF

        CASE ('rph')
          IF (associated(Qtransfo%RPHTransfo)) THEN
            CALL Write_RPHTransfo(Qtransfo%RPHTransfo)
          END IF

        CASE ('hyperspherical')
          write(out_unitp,*) 'nb_HyperSphe: ',                          &
                 Qtransfo%HyperSpheTransfo%nb_HyperSphe
          write(out_unitp,*) 'list_HyperSphe: ',                        &
                 Qtransfo%HyperSpheTransfo%list_HyperSphe(:)

        CASE ('oned')
          write(out_unitp,*) 'oneD transfo'
          IF (Qtransfo%nb_transfo < 1) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  Wrong number of transformation:',   &
                                                   Qtransfo%nb_transfo
              write(out_unitp,*) '  for the oneD transformation'
              write(out_unitp,*) ' Check the fortran source !!'
              STOP
          END IF
          CALL Write_oneDTransfo(Qtransfo%oneDTransfo)

        CASE ('threed')
          CALL Write_ThreeDTransfo(Qtransfo%ThreeDTransfo)

        CASE ('rot2coord')
          CALL Write_Rot2CoordTransfo(Qtransfo%Rot2CoordTransfo)

        CASE ('flexible')
          nb_flex_act = Qtransfo%FlexibleTransfo%nb_flex_act
          write(out_unitp,*) 'nb_flex_act',nb_flex_act,':',             &
                 Qtransfo%FlexibleTransfo%list_act(1:nb_flex_act)
          write(out_unitp,*) 'flex: ',                                  &
                               Qtransfo%FlexibleTransfo%list_flex(:)

        CASE ('gene')
          CALL Write_GeneTransfo(Qtransfo%GeneTransfo)

        CASE ('active')
          CALL Write_ActiveTransfo(Qtransfo%ActiveTransfo)

        CASE ('zmat')
          CALL Write_ZmatTransfo(Qtransfo%ZmatTransfo)

        CASE ('bunch','bunch_poly') ! It should one of the first transfo
          CALL Write_BunchTransfo(Qtransfo%BunchTransfo)

        CASE ('poly')
          CALL RecWrite_BFTransfo(Qtransfo%BFTransfo)

        CASE ('qtox_ana')
          CALL Write_QTOXanaTransfo(Qtransfo%QTOXanaTransfo)

        CASE ('cartesian')
          CALL Write_CartesianTransfo(Qtransfo%CartesianTransfo)

        CASE default ! ERROR: wrong transformation !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The transformation is UNKNOWN: ',        &
                                         trim(Qtransfo%name_transfo)
          CALL Write_list_Qtransfo(out_unitp)
          write(out_unitp,*) ' Check the source!'
          STOP
        END SELECT
        write(out_unitp,*) 'END ',name_sub

        CALL flush_perso(out_unitp)
      END SUBROUTINE Write_Qtransfo

      SUBROUTINE sub_Type_Name_OF_Qin(Qtransfo,name_coord)
        IMPLICIT NONE
        TYPE(type_qtransfo), intent(inout) :: Qtransfo
        character (len=*),   intent(in)    :: name_coord

        integer :: i
        integer :: it
        character (len=*), parameter :: name_sub = 'sub_Type_Name_OF_Qin'

        IF (.NOT. associated(Qtransfo%type_Qin)) THEN
          CALL alloc_array(Qtransfo%type_Qin,(/Qtransfo%nb_Qin/),       &
                          "Qtransfo%type_Qin",name_sub)
        END IF

        IF (.NOT. associated(Qtransfo%name_Qin)) THEN
          CALL alloc_array(Qtransfo%name_Qin,(/Qtransfo%nb_Qin/),Name_len,&
                          "Qtransfo%name_Qin",name_sub)
        END IF

        it = Qtransfo%num_transfo
        DO i=1,Qtransfo%nb_Qin
          CALL make_nameQ(Qtransfo%name_Qin(i),trim(adjustl(name_coord)),i,it)
          Qtransfo%type_Qin(i) = 0
        END DO

      END SUBROUTINE sub_Type_Name_OF_Qin

      SUBROUTINE Sub_Check_LinearTransfo(Qtransfo)
      TYPE (Type_Qtransfo), intent(inout) :: Qtransfo


       integer        :: i_Qout,i_Qin,typ_Q

!      -----------------------------------------------------------------
      character (len=*), parameter :: name_sub='Sub_Check_LinearTransfo'
!      logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
!      -----------------------------------------------------------------
       IF (.NOT. Qtransfo%LinearTransfo%check_LinearTransfo) RETURN
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*)
         CALL Write_Qtransfo(Qtransfo)
         write(out_unitp,*)
       END IF
!      -----------------------------------------------------------------
      IF (.NOT. associated(Qtransfo%type_Qout) ) THEN
        write(out_unitp,*) ' ERROR in name_sub'
        write(out_unitp,*) ' Qtransfo%type_Qout is not associated'
        write(out_unitp,*) ' CHECK the fortran !'
        STOP
      END IF
      Qtransfo%type_Qin(:) = 0

      DO i_Qin=1,Qtransfo%nb_Qin
        typ_Q =-1
        DO i_Qout=1,Qtransfo%nb_Qin
          IF (Qtransfo%LinearTransfo%mat(i_Qout,i_Qin) /= ZERO) THEN
            IF (typ_Q == -1) THEN
              typ_Q = Qtransfo%type_Qout(i_Qout)
              Qtransfo%type_Qin(i_Qin) = typ_Q
            ELSE
            IF (typ_Q /= Qtransfo%type_Qout(i_Qout) ) THEN
                write(out_unitp,*) '==================================='
                write(out_unitp,*) '==================================='
                write(out_unitp,*) '==================================='
                CALL Write_Qtransfo(Qtransfo)
                write(out_unitp,*) '==================================='
                write(out_unitp,*) '==================================='
                write(out_unitp,*) '==================================='
                write(out_unitp,*) 'ERROR: in ',name_sub
                write(out_unitp,*) 'i_Qout,i_Qin',i_Qout,i_Qin
                write(out_unitp,*) 'type_Qout and type_Qin',            &
                                       Qtransfo%type_Qout(i_Qout),typ_Q
                write(out_unitp,*)
                write(out_unitp,*) '==================================='
                write(out_unitp,*) '==================================='
                write(out_unitp,*) '==================================='
                STOP
             END IF
            END IF
          END IF
        END DO
      END DO

!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'type_Qin  : ',Qtransfo%type_Qin
         write(out_unitp,*) 'type_Qout : ',Qtransfo%type_Qout
         write(out_unitp,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

      END SUBROUTINE Sub_Check_LinearTransfo

      END MODULE mod_Qtransfo
