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
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
      MODULE mod_GeneTransfo
      use mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_GeneTransfo
        integer          :: nb_Coord_transfo  = 0
        integer, pointer :: list_gene(:) => null()
        integer, pointer :: list_Coord_transfo(:) => null()

      END TYPE Type_GeneTransfo

      PUBLIC :: Type_GeneTransfo, Read_GeneTransfo, Write_GeneTransfo, calc_GeneTransfo
      PUBLIC :: alloc_GeneTransfo, dealloc_GeneTransfo, GeneTransfo1TOGeneTransfo2

      CONTAINS

!=======================================================================
!     Felxible transfo
!=======================================================================
      !!@description: Felxible transf
      !!@param: TODO
      SUBROUTINE Read_GeneTransfo(GeneTransfo,nb_Qin)

      TYPE (Type_GeneTransfo), intent(inout) :: GeneTransfo
      integer, intent(in) :: nb_Qin

      integer :: i,it,iQ,err

      character (len=*), parameter :: name_sub='Read_GeneTransfo'


      CALL alloc_array(GeneTransfo%list_gene,[nb_Qin],                &
                      "GeneTransfo%list_gene",name_sub)
      GeneTransfo%list_gene(:) = 0

      read(in_unitp,*,IOSTAT=err) GeneTransfo%list_gene(:)
      IF (err /= 0) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) '  while reading "list_gene"'
         write(out_unitp,*) '  end of file or end of record'
         write(out_unitp,*) ' Check your data !!'
         STOP
      END IF

      GeneTransfo%nb_Coord_transfo = count(GeneTransfo%list_gene == 1)
      IF (GeneTransfo%nb_Coord_transfo < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The number of coordinates in the general transfo is 0'
        write(out_unitp,*) ' Check your data !!'
        write(out_unitp,*) 'list_gene',GeneTransfo%list_gene(:)
        STOP
      END IF

      CALL alloc_array(GeneTransfo%list_Coord_transfo,                  &
                                    [GeneTransfo%nb_Coord_transfo],   &
                      "GeneTransfo%list_Coord_transfo",name_sub)
      GeneTransfo%list_Coord_transfo(:) = 0

      iQ = 0
      DO i=1,nb_Qin
        IF (GeneTransfo%list_gene(i) == 1) THEN
          iQ = iQ+1
          GeneTransfo%list_Coord_transfo(iQ) = i
        END IF
      END DO

      CALL Write_GeneTransfo(GeneTransfo)

      END SUBROUTINE Read_GeneTransfo

!=======================================================================
!     Felxible transfo
!=======================================================================
      !!@description: Felxible transf
      !!@param: TODO
      SUBROUTINE Write_GeneTransfo(GeneTransfo)

      TYPE (Type_GeneTransfo), intent(in) :: GeneTransfo


      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_GeneTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub

      write(out_unitp,*) 'list_gene',GeneTransfo%list_gene(:)
      write(out_unitp,*) 'nb_Coord_transfo',GeneTransfo%nb_Coord_transfo
      write(out_unitp,*) 'list_Coord_transfo',GeneTransfo%list_Coord_transfo(:)

      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_GeneTransfo

      SUBROUTINE alloc_GeneTransfo(GeneTransfo,nb_Qin,nb_Coord_transfo)

      TYPE (Type_GeneTransfo), intent(inout) :: GeneTransfo
      integer                                :: nb_Qin,nb_Coord_transfo

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_GeneTransfo'


      CALL dealloc_GeneTransfo(GeneTransfo)

      CALL alloc_array(GeneTransfo%list_gene,[nb_Qin],                &
                      "GeneTransfo%list_gene",name_sub)
      GeneTransfo%list_gene(:) = 0

      GeneTransfo%nb_Coord_transfo = nb_Coord_transfo

      CALL alloc_array(GeneTransfo%list_Coord_transfo,                  &
                                    [GeneTransfo%nb_Coord_transfo],   &
                      "GeneTransfo%list_Coord_transfo",name_sub)
      GeneTransfo%list_Coord_transfo(:) = 0


      END SUBROUTINE alloc_GeneTransfo

      SUBROUTINE dealloc_GeneTransfo(GeneTransfo)

      TYPE (Type_GeneTransfo), intent(inout) :: GeneTransfo

      character (len=*), parameter :: name_sub='dealloc_GeneTransfo'

      IF (associated(GeneTransfo%list_gene)) THEN
        CALL dealloc_array(GeneTransfo%list_gene,                       &
                          "GeneTransfo%list_gene",name_sub)
      END IF

      IF (associated(GeneTransfo%list_Coord_transfo)) THEN
        CALL dealloc_array(GeneTransfo%list_Coord_transfo,                  &
                          "GeneTransfo%list_Coord_transfo",name_sub)
      END IF

      GeneTransfo%nb_Coord_transfo = 0

      END SUBROUTINE dealloc_GeneTransfo

      SUBROUTINE GeneTransfo1TOGeneTransfo2(GeneTransfo1,GeneTransfo2)

!      for the Activerix and Tnum --------------------------------------
       TYPE (Type_GeneTransfo), intent(in)    :: GeneTransfo1
       TYPE (Type_GeneTransfo), intent(inout) :: GeneTransfo2

!----- for debuging ----------------------------------
      character (len=*), parameter ::                                   &
                                  name_sub = 'GeneTransfo1TOGeneTransfo2'
       logical, parameter :: debug=.FALSE.
!       logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'GeneTransfo1'
        CALL Write_GeneTransfo(GeneTransfo1)
        CALL flush_perso(out_unitp)
      END IF
!---------------------------------------------------------------------

      CALL dealloc_GeneTransfo(GeneTransfo2)

      GeneTransfo2%nb_Coord_transfo = GeneTransfo1%nb_Coord_transfo

      CALL alloc_GeneTransfo(GeneTransfo2,size(GeneTransfo1%list_gene), &
                                          GeneTransfo1%nb_Coord_transfo)

      GeneTransfo2%list_gene          = GeneTransfo1%list_gene
      GeneTransfo2%list_Coord_transfo = GeneTransfo1%list_Coord_transfo

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'GeneTransfo2'
        CALL Write_GeneTransfo(GeneTransfo2)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE GeneTransfo1TOGeneTransfo2


      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_GeneTransfo(dnQin,dnQout,GeneTransfo,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)        :: dnQin,dnQout
        TYPE (Type_GeneTransfo), intent(in)     :: GeneTransfo
        integer, intent(in)                     :: nderiv
        logical                                 :: inTOout



        TYPE (Type_dnS)   :: dnQgene

        TYPE (Type_dnS)   :: dnQ


        integer :: id,jd,kd

        integer :: iQ,iq_gene,it
        integer :: pd_Qin,pd_Qene
        integer :: qd_Qin,qd_Qene
        integer :: rd_Qin,rd_Qene

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_GeneTransfo'
       logical, parameter :: debug=.FALSE.
!       logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL Write_GeneTransfo(GeneTransfo)
        CALL flush_perso(out_unitp)
      END IF
!---------------------------------------------------------------------

       CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
       CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

       IF (inTOout) THEN
         iq_gene = 0
         CALL alloc_dnS(dnQgene,GeneTransfo%nb_Coord_transfo,nderiv)
         CALL alloc_dnS(dnQ,dnQin%nb_var_deriv,nderiv)

         DO iQ=1,dnQin%nb_var_vec
           IF (GeneTransfo%list_gene(iQ) == 0) THEN
             CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQin,dnQout,iQ,nderiv)
           ELSE
             iq_gene = iq_gene + 1
             write(out_unitp,*) 'STOP in ',name_sub
             write(out_unitp,*) ' You should:'
             write(out_unitp,*) ' - Uncomment the "CALL calc_dnQgene" at'
             write(out_unitp,*) ' line 255 of the Tnum/Qtransfo/GeneTransfo.f90 file'
             write(out_unitp,*) ' - Define the "calc_dnQgene" subroutine in '
             write(out_unitp,*) ' the sub_system.f file'

             STOP
!             CALL calc_dnQgene(iq_gene,dnQgene,                         &
!                               dnQin%d0(GeneTransfo%list_Coord_transfo),&
!                         GeneTransfo%nb_Coord_transfo,nderiv,it,inTOout)

             CALL sub_ZERO_TO_dnS(dnQ,nderiv)

             IF (nderiv >= 0) THEN
               dnQ%d0 = dnQgene%d0
             END IF

             ! first derivative
             ! dQO_i/dQact(.) = Sum_p (dQgene_i/dQ_p * dQN_p/dQact(.) )
             IF (nderiv >= 1) THEN
               dnQ%d1(:) = ZERO
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 dnQ%d1(:) = dnQ%d1(:) +                                &
                                dnQgene%d1(pd_Qene) * dnQin%d1(pd_Qin,:)
               END DO
             END IF

             ! 2d derivative (two contributions)
             IF (nderiv >= 2) THEN
               dnQ%d2(:,:) = ZERO
              ! d2QO_i/dQact(.)dQact(.) = Sum_p (dQgene_i/dQ_p * d2QN_p/dQact(.)dQact(.) )
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 dnQ%d2(:,:) = dnQ%d2(:,:) +                            &
                              dnQgene%d1(pd_Qene) * dnQin%d2(pd_Qin,:,:)
               END DO

              ! d2QO_i/dQact(id)dQact(jd) = Sum_pq (d2Qgene_i/dQ_pdQ_q * dQN_p/dQact(id) * dQN_q/dQact(id) )
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO qd_Qene=1,GeneTransfo%nb_Coord_transfo

                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 qd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)

                 DO id=1,dnQ%nb_var_deriv
                 DO jd=1,dnQ%nb_var_deriv
                   dnQ%d2(id,jd) = dnQ%d2(id,jd) +                      &
                              dnQgene%d2(pd_Qene,qd_Qene) *             &
                              dnQin%d1(pd_Qin,id) * dnQin%d1(qd_Qin,jd)
                 END DO
                 END DO

               END DO
               END DO

             END IF

             ! 3d derivative (3 contributions)
             IF (nderiv >= 3) THEN
               dnQ%d3(:,:,:) = ZERO

               ! first contribution
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 dnQ%d3(:,:,:) = dnQ%d3(:,:,:) +                        &
                            dnQgene%d1(pd_Qene) * dnQin%d3(pd_Qin,:,:,:)
               END DO

               ! 2d contribution
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO qd_Qene=1,GeneTransfo%nb_Coord_transfo

                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 qd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)

                 DO id=1,dnQ%nb_var_deriv
                 DO jd=1,dnQ%nb_var_deriv
                 DO kd=1,dnQ%nb_var_deriv

                   dnQ%d3(id,jd,kd) = dnQ%d3(id,jd,kd) +                &
                              dnQgene%d2(pd_Qene,qd_Qene) * (           &
                         dnQin%d2(pd_Qin,id,kd) * dnQin%d1(qd_Qin,jd) + &
                         dnQin%d1(pd_Qin,id) * dnQin%d2(qd_Qin,jd,kd) + &
                         dnQin%d2(pd_Qin,id,jd) * dnQin%d1(qd_Qin,kd) )
                 END DO
                 END DO
                 END DO

               END DO
               END DO

               ! 3d contribution
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO qd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO rd_Qene=1,GeneTransfo%nb_Coord_transfo

                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 qd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)
                 rd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)

                 DO id=1,dnQ%nb_var_deriv
                 DO jd=1,dnQ%nb_var_deriv
                 DO kd=1,dnQ%nb_var_deriv

                   dnQ%d3(id,jd,kd) = dnQ%d3(id,jd,kd) +                &
                                  dnQgene%d3(pd_Qene,qd_Qene,rd_Qene) * &
                                                  dnQin%d1(pd_Qin,id) * &
                                                  dnQin%d1(qd_Qin,jd) * &
                                                  dnQin%d1(rd_Qin,kd)
                 END DO
                 END DO
                 END DO

               END DO
               END DO
               END DO

             END IF

             IF (nderiv >= 4) THEN
               write(out_unitp,*) 'ERROR in ',name_sub
               write(out_unitp,*) 'nderiv is too large (>3) ',nderiv
               STOP
             END IF


             CALL sub_dnS_TO_dnVec(dnQ,dnQout,iQ,nderiv)
           END IF
         END DO
         CALL dealloc_dnS(dnQgene)
         CALL dealloc_dnS(dnQ)

       ELSE

         iq_gene = 0
         CALL alloc_dnS(dnQgene,GeneTransfo%nb_Coord_transfo,nderiv)
         CALL alloc_dnS(dnQ,dnQin%nb_var_deriv,nderiv)

         DO iQ=1,dnQin%nb_var_vec
           IF (GeneTransfo%list_gene(iQ) == 0) THEN
             CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQout,dnQin,iQ,nderiv)
           ELSE
             iq_gene = iq_gene + 1
!             CALL calc_dnQgene(iq_gene,dnQgene,                         &
!                              dnQout%d0(GeneTransfo%list_Coord_transfo),&
!                         GeneTransfo%nb_Coord_transfo,nderiv,it,inTOout)

             CALL sub_ZERO_TO_dnS(dnQ,nderiv)

             IF (nderiv >= 0) THEN
               dnQ%d0 = dnQgene%d0
             END IF

             ! first derivative
             ! dQO_i/dQact(.) = Sum_p (dQgene_i/dQ_p * dQN_p/dQact(.) )
             IF (nderiv >= 1) THEN
               dnQ%d1(:) = ZERO
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 dnQ%d1(:) = dnQ%d1(:) +                                &
                               dnQgene%d1(pd_Qene) * dnQout%d1(pd_Qin,:)
               END DO
             END IF

             ! 2d derivative (two contributions)
             IF (nderiv >= 2) THEN
               dnQ%d2(:,:) = ZERO
              ! d2QO_i/dQact(.)dQact(.) = Sum_p (dQgene_i/dQ_p * d2QN_p/dQact(.)dQact(.) )
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 dnQ%d2(:,:) = dnQ%d2(:,:) +                            &
                             dnQgene%d1(pd_Qene) * dnQout%d2(pd_Qin,:,:)
               END DO

              ! d2QO_i/dQact(id)dQact(jd) = Sum_pq (d2Qgene_i/dQ_pdQ_q * dQN_p/dQact(id) * dQN_q/dQact(id) )
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO qd_Qene=1,GeneTransfo%nb_Coord_transfo

                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 qd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)

                 DO id=1,dnQ%nb_var_deriv
                 DO jd=1,dnQ%nb_var_deriv
                   dnQ%d2(id,jd) = dnQ%d2(id,jd) +                      &
                              dnQgene%d2(pd_Qene,qd_Qene) *             &
                             dnQout%d1(pd_Qin,id) * dnQout%d1(qd_Qin,jd)
                 END DO
                 END DO

               END DO
               END DO

             END IF

             ! 3d derivative (3 contributions)
             IF (nderiv >= 3) THEN
               dnQ%d3(:,:,:) = ZERO

               ! first contribution
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 dnQ%d3(:,:,:) = dnQ%d3(:,:,:) +                        &
                           dnQgene%d1(pd_Qene) * dnQout%d3(pd_Qin,:,:,:)
               END DO

               ! 2d contribution
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO qd_Qene=1,GeneTransfo%nb_Coord_transfo

                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 qd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)

                 DO id=1,dnQ%nb_var_deriv
                 DO jd=1,dnQ%nb_var_deriv
                 DO kd=1,dnQ%nb_var_deriv

                   dnQ%d3(id,jd,kd) = dnQ%d3(id,jd,kd) +                &
                              dnQgene%d2(pd_Qene,qd_Qene) * (           &
                       dnQout%d2(pd_Qin,id,kd) * dnQout%d1(qd_Qin,jd) + &
                       dnQout%d1(pd_Qin,id) * dnQout%d2(qd_Qin,jd,kd) + &
                        dnQout%d2(pd_Qin,id,jd) * dnQout%d1(qd_Qin,kd) )
                 END DO
                 END DO
                 END DO

               END DO
               END DO

               ! 3d contribution
               DO pd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO qd_Qene=1,GeneTransfo%nb_Coord_transfo
               DO rd_Qene=1,GeneTransfo%nb_Coord_transfo

                 pd_Qin = GeneTransfo%list_Coord_transfo(pd_Qene)
                 qd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)
                 rd_Qin = GeneTransfo%list_Coord_transfo(qd_Qene)

                 DO id=1,dnQ%nb_var_deriv
                 DO jd=1,dnQ%nb_var_deriv
                 DO kd=1,dnQ%nb_var_deriv

                   dnQ%d3(id,jd,kd) = dnQ%d3(id,jd,kd) +                &
                                  dnQgene%d3(pd_Qene,qd_Qene,rd_Qene) * &
                                                 dnQout%d1(pd_Qin,id) * &
                                                 dnQout%d1(qd_Qin,jd) * &
                                                 dnQout%d1(rd_Qin,kd)
                 END DO
                 END DO
                 END DO

               END DO
               END DO
               END DO

             END IF

             IF (nderiv >= 4) THEN
               write(out_unitp,*) 'ERROR in ',name_sub
               write(out_unitp,*) 'nderiv is too large (>3) ',nderiv
               STOP
             END IF


             CALL sub_dnS_TO_dnVec(dnQ,dnQin,iQ,nderiv)
           END IF
         END DO
         CALL dealloc_dnS(dnQgene)
         CALL dealloc_dnS(dnQ)

       END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_GeneTransfo

      END MODULE mod_GeneTransfo
