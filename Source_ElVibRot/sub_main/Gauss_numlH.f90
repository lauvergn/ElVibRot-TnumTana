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
      Program Gauss_num
      USE mod_system
      USE mod_Coord_KEO
      USE mod_poly
      USE mod_GWP
      USE mod_propa
      IMPLICIT NONE

!.....DECLARE ARRAY Tnum
       TYPE (CoordType) :: mole
       TYPE (Tnum)    :: para_Tnum

!     for Qsym Qact ....
      TYPE (param_Q) :: para_Q


      integer :: ndimA
      real (kind=Rkind), pointer :: d0g(:,:),d1g(:,:,:),d2g(:,:,:,:)
      real (kind=Rkind), pointer :: d0GG(:,:),d1GG(:,:,:)
      real (kind=Rkind), pointer :: d2GG(:,:,:,:)
      real (kind=Rkind)              :: mu0(3)

!     - for the pot on the fly -------------------------
      logical :: onthefly,lo_calc_EG,newQsym,Adiago,Asym,ReNorm,test
      integer :: calc_EG
      integer :: nderiv
!     ------------------------------------------------------

      real (kind=Rkind), pointer :: DeltaQsym(:),PQsym(:)

!.....DECLARE ARRAY Gauss
      TYPE (para_GWP) :: GWP,GWP0
      TYPE (para_poly) :: poly,poly0,conj_poly
      TYPE (para_poly) :: w1_poly,w2_poly
      TYPE (para_LHA)  :: param_LHA
      complex (kind=Rkind) :: Ca,Cb,auto_c
      real (kind=Rkind)    :: None_symA
      integer, pointer :: ind_exp(:)
      integer          :: id0,id1,id2
      real (kind=Rkind),pointer :: freq(:),d0c(:,:),A(:,:)
      logical :: traj

!----- variables for the propagation ----------------------------
      TYPE (param_propa) :: para_propa
      integer            :: no
      real (kind=Rkind)      :: TFmaxE

       character(len=Name_len) pos,restarti,restartf
       character(len=Name_len) imp,sig
       character(len=Name_len) testenergy
       integer NDUMP,NWP
       integer NRS
       integer :: NCOUNT,NFIN,NCO
       integer :: ND,N,IREP,NGAMMA,NAK
       integer :: INDEX,INDEXI,INDEXF
       real(kind=Rkind) :: X,XPAS,XEND,eneRGY,GAMMA
       real(kind=Rkind) :: AK,X0,XFS,E
       real(kind=Rkind) :: BIDON,BID
       real(kind=Rkind),pointer :: Y(:),FU(:),YOUT(:)
       real(kind=Rkind),pointer :: arg(:),VTEMP(:)
       real(kind=Rkind),pointer :: deltaQ(:),sigma(:)
       complex(kind=Rkind),pointer :: ZI(:,:),ZS(:,:)
       complex(kind=Rkind),pointer :: trav(:)
       INTEGER,pointer :: inverse_index(:)
       complex(kind=Rkind) :: EYE,phase
       complex(kind=Rkind) :: TRACE

!.....DECLARE DIVERS ------------------------------------------
       integer :: i,j,k,iact,NN
       real(kind=Rkind)             :: pot0

!.....NAMELIST Tnum
      NAMELIST /niveau/pot0,onthefly,calc_EG,newQsym,Adiago,Asym,ReNorm,&
                       test

!.....NAMELIST Gauss
      NAMELIST/file/pos,imp,sig,restartf,testenergy,restarti
      NAMELIST/gauss/X0,XEND,XPAS,GAMMA,NAK,AK,NRS,NGAMMA,NWP,TFmaxE,   &
                     traj

!.....FILES

      pos='Qmoy'
      imp='Pmoy'
      sig='sig'
      testenergy='gaussiene_out'
      restarti='rsi'
      restartf='rsf'
      read(in_unitp,file)

      open(1,file='hessienne')
      open(8,file=restarti)
      open(9,file=restartf)
      open(10,file=testenergy)

      open(31,file=pos)
      open(32,file='matA')
      write(32,*)'x,xfs,A'
      open(41,file=sig)
      open(71,file=imp)

      open(60,file='linear-PZ')
      write(out_unitp0,*)'x,PZ'
      open(61,file='linear-Z')
      write(out_unitp1,*)'x,Z'
      open(62,file='phase')
      write(out_unitp2,*)'x,Phase'
      open(64,file='vecY')
  3   FORMAT(22(E12.5,1X))


!.....Tnum

!-------------------------------------------------
!     - read the CoordType -------------------------
!     --------------------------------------------
      CALL lect0_zmat(mole,para_Tnum)
      para_Tnum%nrho    = 1
!     --------------------------------------------
!     - initialization of para_Q -----------------
      para_Q%nb_var    = mole%nb_var
      CALL alloc_param_Q(para_Q)
!     --------------------------------------------
!-------------------------------------------------


      pot0  = -huge(1.d0)
      onthefly = .TRUE.
      calc_EG  = 4
      Adiago   = .TRUE.
      Asym     = .TRUE.
      ReNorm   = .TRUE.
      newQsym  = .FALSE.
      test     = .FALSE.

      read(in_unitp,niveau)
      write(out_unitp,niveau)

      memory = product( (/ mole%nb_var /) )
      allocate(DeltaQsym(mole%nb_var),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"DeltaQsym","main")
      memory = product( (/ mole%nb_var /) )
      allocate(PQsym(mole%nb_var),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"PQsym","main")
      DeltaQsym(:) = 0.
      PQsym(:)  = 0.
      IF (newQsym) THEN
        write(out_unitp,*) 'variables,DeltaQ,PQ :'
        DO i=1,mole%nb_var
         read(in_unitp,*) para_Q%Qsym0(i),DeltaQsym(i),PQsym(i)
         write(out_unitp,*) 'Qsym0,DeltaQsym,PQsym(',i,')=',                    &
                      para_Q%Qsym0(i),DeltaQsym(i),PQsym(i)
        END DO
      ELSE
        write(out_unitp,*) 'variables :'
        DO i=1,mole%nb_var
         read(in_unitp,*) para_Q%Qsym0(i)
         write(out_unitp,*) 'Qsym0(',i,')=',para_Q%Qsym0(i)
        END DO
      END IF
      CALL Qsym_TO_Qact(para_Q%Qsym0,para_Q%Qact0,para_Q%nb_var,        &
                        mole%liste_QactTOQdyn)
      para_Q%Qsym(:) = para_Q%Qsym0(:)
      CALL Qsym_TO_Qact(para_Q%Qsym,para_Q%Qact,para_Q%nb_var,          &
                        mole%liste_QactTOQdyn)

!......Allocation Tnum

        ndimA = mole%nb_act+3

        memory = product( (/ ndimA,ndimA /) )
        allocate(d0g(ndimA,ndimA),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"d0g","main")
        memory = product( (/ ndimA,ndimA,mole%nb_act /) )
        allocate(d1g(ndimA,ndimA,mole%nb_act),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"d1g","main")
        memory = product( (/ ndimA,ndimA,mole%nb_act,mole%nb_act /) )
        allocate(d2g(ndimA,ndimA,mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"d2g","main")
        memory = product( (/ ndimA,ndimA /) )
        allocate(d0GG(ndimA,ndimA),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"d0GG","main")
        memory = product( (/ ndimA,ndimA,mole%nb_act /) )
        allocate(d1GG(ndimA,ndimA,mole%nb_act),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"d1GG","main")
        memory = product( (/ ndimA,ndimA,mole%nb_act,mole%nb_act /) )
        allocate(d2GG(ndimA,ndimA,mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"d2GG","main")

!......INPUT Gauss
!     NAMELIST/gauss/X0,XEND,XPAS,GAMMA,NAK,AK,NRS,NGAMMA,NWP,TFmaxE
       ngamma = 0
       gamma  = 1.d0
       nak    = 0
       ak     = 0.d0
       x0     = 0.d0
       xend   = 100.d0
       xpas   = 10.d0
       nwp    = 5
       nrs    = 0
       TFmaxE = 5000.d0 ! maximal enrgy in cm-1
       traj   = .FALSE.
       read(in_unitp,gauss)

!......Allocation Gauss
       CALL init0_GWP(GWP)
       CALL init0_GWP(GWP0)
       CALL init0_poly(poly)
       CALL init0_poly(conj_poly)
       CALL init0_poly(poly0)
       CALL init0_poly(w1_poly)
       CALL init0_poly(w2_poly)
       CALL alloc_GWP(GWP,ndim=mole%nb_act,                             &
                      cplx=.TRUE.,                                      &
                      linearization=.TRUE.,trajectory=traj)
       memory = product( (/ mole%nb_act /) )
       allocate(ind_exp(mole%nb_act),stat=err_mem) ! change alloc done
       CALL error_memo_allo(err_mem,memory,"ind_exp","main")

       CALL init0_LHA(param_LHA)
       CALL alloc_LHA(param_LHA,ndim=mole%nb_act)

       ND=mole%nb_act
       N = 2*ND ! Qt and Pt
       IF (.NOT. GWP%trajectory) THEN
         N = N + 2 ! phase and norm
         IF (GWP%linearization) THEN
            N = N + 4*ND**2 ! Z and PZ
         ELSE
            N = N + 2*ND**2 ! A
         END IF
       END IF
       write(out_unitp,*) 'ND,N',ND,N

       memory = product( (/ N /) )
       allocate(Y(N),FU(N),YOUT(N),stat=err_mem) ! change alloc done
       CALL error_memo_allo(err_mem,memory,"Y(N),FU(N),YOUT","main")
       memory = product( (/ ND /) )
       allocate(ARG(ND),stat=err_mem) ! change alloc done
       CALL error_memo_allo(err_mem,memory,"ARG","main")
       memory = product( (/ ND /) )
       allocate(deltaQ(ND),sigma(ND),stat=err_mem) ! change alloc done
       CALL error_memo_allo(err_mem,memory,"deltaQ(ND),sigma","main")
       memory = product( (/ ND /) )
       allocate(VTEMP(ND),stat=err_mem) ! change alloc done
       CALL error_memo_allo(err_mem,memory,"VTEMP","main")
       memory = product( (/ ND /) )
       allocate(ZI(ND,ND),ZS(ND,ND),TRAV(ND),INVERSE_INDEX(ND),         &
                                                           stat=err_mem) ! change alloc done
       CALL error_memo_allo(err_mem,memory,"ZI(ND,ND),ZS(ND,ND),TRAV(ND) &
                                                 ,INVERSE_INDEX","main")

!.....CONSTANT
       eye=cmplx(0.d0,1.d0,kind=Rkind)
       irep=0
       ndump=int((xend-X0)/xpas)/nwp     !AU PLUS 5 WP NWP = 5
!==================================================
!     For para_propa and the autocorr function
        para_propa%WPTmax       = XEND
        para_propa%WPdeltaT     = XPAS
        para_propa%nb_micro     = 1

        para_propa%file_autocorr%name = 'file_auto'
        para_propa%file_spectrum%name = 'file_spectrum'

        para_propa%TFnexp2      = 4+int(log(XEND/XPAS)/log(2.d0))
        para_propa%TFmaxE       = TFmaxE / get_ConvUnit('E','cm-1')

      CALL file_open(para_propa%file_autocorr,no)
!==================================================

!....MICH
!....Lee and Helle JCP 76 3035 1982

!......INPUT Gauss
!......INITIAL CONDITIONS

       Y(:)=0.d0
!      - RESTART if NRS > 0 ---------------------
       IF(NRS == 0) then

!        variables P and Q
         DO i=1,ND
           iact = mole%liste_QactTOQdyn(i)
           GWP%Qmean(i) = para_Q%Qsym(iact)
           GWP%Pmean(i) = PQsym(iact)
         END DO

!        variables P (old fashion)
         IF (nak > 0 .AND. nak <= ND) THEN
           GWP%Pmean(:)=0.d0
           GWP%Pmean(nak)=ak
         END IF
         write(out_unitp,*) 'Pmean',nak,GWP%Pmean

!        WIDTH AND CORRELATION MATRIX
!        ZERO FOR REAL PART OF THE A MATRIX
!        INPUT OF THE IMAGINARY PART OF THE A MATRIX

!        COMPUTATION OF  deltaQ

         CALL time_perso('init')
         nderiv = 2
         IF (GWP%trajectory) nderiv = 1
         CALL dnOp_grid(para_Q%Qact,                                         &
                        param_LHA%d0Ene,param_LHA%d1Ene,param_LHA%d2Ene,&
                        nderiv,                                         &
                        param_LHA%d0mu,param_LHA%d1mu,mole,             &
                        para_Tnum,onthefly)
!        CALL gr_hes(para_Q%Qsym,
!    *               param_LHA%d0Ene,param_LHA%d1Ene,param_LHA%d2Ene,ND)
         param_LHA%pot0 = pot0
         mu0(:) = param_LHA%d0mu
         IF ( pot0 == -huge(1.0d0) ) param_LHA%pot0 = param_LHA%d0Ene
         param_LHA%d0Ene = param_LHA%d0Ene - param_LHA%pot0
         write(out_unitp,*) 'ene0,pot0 PREMIER APPEL',                          &
                         param_LHA%d0Ene,param_LHA%pot0
         write(out_unitp,*) 'mu0',mu0

         nderiv = 2
         IF (GWP%trajectory) nderiv = 1
         CALL calc2_d0d1d2g_G_bis(para_Q,                               &
                                  d0g,d1g,d2g,                          &
                                  d0GG,d1GG,d2GG,nderiv,                &
                                  para_Tnum,mole)

         param_LHA%d0G(:,:) = d0gg(1:nd,1:nd)
         CALL time_perso('init')


         IF (.NOT. GWP%trajectory) THEN
           memory = product( (/ nd /) )
           allocate(freq(nd),stat=err_mem) ! change alloc done
           CALL error_memo_allo(err_mem,memory,"freq","main")
           memory = product( (/ nd,nd /) )
           allocate(d0c(nd,nd),stat=err_mem) ! change alloc done
           CALL error_memo_allo(err_mem,memory,"d0c","main")
           memory = product( (/ nd,nd /) )
           allocate(A(nd,nd),stat=err_mem) ! change alloc done
           CALL error_memo_allo(err_mem,memory,"A","main")
           CALL calc_freq_width(nd,A,d0c,freq,                          &
                                param_LHA%d2Ene,param_LHA%d0G)
           write(out_unitp,*) 'freq cm-1',freq * get_ConvUnit('E','cm-1')
           write(out_unitp,*) 'width matrix'
           CALL Write_Mat(A,out_unitp,5)
           memory = size(freq)
           deallocate(freq,stat=err_mem) ! change dealloc done
           CALL error_memo_allo(err_mem,-memory,"freq","main")
           memory = size(d0c)
           deallocate(d0c,stat=err_mem) ! change dealloc done
           CALL error_memo_allo(err_mem,-memory,"d0c","main")
           write(out_unitp,*)


!          k*mu
           DO i=1,ND
             deltaQ(i)=param_LHA%d2Ene(i,i)/param_LHA%d0G(i,i)
           END DO

!          deltaQ=sqrt(1/(2*sqrt(k*mu)))
           deltaQ(:) = sqrt(deltaQ(:))*2.d0
           deltaQ(:) = 1.d0/sqrt(deltaQ(:))


           write(out_unitp,*)'   '
           write(out_unitp,*)'hes(i,i), masse,  delta '
           write(out_unitp,*)'   '
           DO i=1,nd
             write(out_unitp,21)i,param_LHA%d2Ene(i,i),1.d0/param_LHA%d0G(i,i), &
                        deltaQ(i)
 21          format(i3,4(1x,f18.6))
           END DO

!          HERE VARIATION OF THE DELTA OF THE CHOSEN VARIABLE
           IF (ngamma > 0 .AND. ngamma <= ND) THEN
             deltaQ(ngamma) = deltaQ(ngamma)*gamma
           END IF


!          COMPUTATION OF ALPHA   alpha=1/(sig)**2   sig = 2 deltaQ
           IF (Adiago) THEN
             GWP%CAmean(:,:) = 0.
             bid = 1.d0
             DO i=1,ND
               GWP%CAmean(i,i) =                                        &
                  cmplx(0.d0,0.25d0/deltaQ(i)**2,kind=Rkind)
               bid = bid * (2.d0*imag(GWP%CAmean(i,i))/pi)**0.25d0
             END DO
!            INITIALISATION OF THE PHASE
             GWP%Cphase = -log(bid)*eye
           ELSE
             GWP%CAmean(:,:) = -eye*A(:,:)
           END IF
           memory = size(A)
           deallocate(A,stat=err_mem) ! change dealloc done
           CALL error_memo_allo(err_mem,-memory,"A","main")

!          INITIALISATION OF PZ_0 = 2 A_0
           GWP%CPZ(:,:) = 2.d0 * GWP%CAmean(:,:)
!          INITIALISATION OF Z_0 = Identity
           GWP%CZ(:,:) = 0.d0
           DO i=1,ND
             GWP%CZ(i,i) = 1.d0
           END DO
         END IF


!        variables Q shifted with DeltaQsym
         para_Q%Qsym(:) = para_Q%Qsym0(:) + DeltaQsym(:)

!        ----------------------------------------------
!        Qmean, Pmean, CZ, CPZ, gamma => Y(:)
         CALL GWP_TO_Y(GWP,Y,N)
         write(out_unitp,*) Y(:)
         write(8,*) Y(:)
!-----------------------------------------------------

      ELSE

        read(8,*) Y(:)
        CALL Y_TO_GWP(Y,GWP,N)

      ENDIF


!......INITIALIZATION


      IF (GWP%trajectory) THEN
        CALL GWP1_TO_GWP2(GWP,GWP0)
      ELSE
        CALL GWP_TO_poly(GWP,poly)
!       Normalization
        CALL p1TOp2(p1=poly,p2=conj_poly,conjug=.TRUE.)
        CALL calc_auto_cor(auto_c,conj_poly,poly)
        write(out_unitp,*) 'norm of GWP:',real(auto_c,kind=Rkind)
!       Renormalization of GWP
        GWP%Cphase = GWP%Cphase +eye*0.5d0*log(real(auto_c,kind=Rkind))
        CALL GWP_TO_poly(GWP,poly)
        CALL p1TOp2(p1=poly,p2=conj_poly,conjug=.TRUE.)
        CALL calc_auto_cor(auto_c,conj_poly,poly)
        write(out_unitp,*) 'New norm of GWP:',real(auto_c,kind=Rkind)
        CALL GWP_TO_Y(GWP,Y,N)

        CALL GWP1_TO_GWP2(GWP,GWP0)
        CALL GWP_TO_poly(GWP0,poly)
!       the conjugate of poly => poly0  (for <GWP0 I GWP>)
        CALL p1TOp2(p1=poly,p2=poly0,conjug=.TRUE.)
      END IF


      IF (test) STOP

!------------------------------------------------------
!     Time loop
!     -------------------------------------------------
      CALL time_perso('GWP propa')
      DO  NCOUNT=0,int(XEND/XPAS)
        X = X0 + real(NCOUNT,kind=Rkind)*XPAS

        lo_calc_EG = .NOT. GWP%trajectory
        CALL FCN(X,Y,FU,ND,N,NdimA,para_Q,mole,                         &
                 onthefly,lo_calc_EG,GWP,param_LHA,                     &
                 para_Tnum)

!       -------------------------------------------------
!       Classical energy calculation
        CALL calc_ene(GWP,param_LHA%d0G,energy,param_LHA%d0Ene)

!       autocorrelation function
        IF (GWP%trajectory) THEN
!         ===============================
!         for trajectory only ......
!
          write(out_unitp,*) 'GWP%Qmean mu',GWP%Qmean,param_LHA%d0mu(:)
!         auto_c = dot_product(GWP%Qmean,GWP0%Qmean)
          auto_c = dot_product(param_LHA%d0mu,mu0)
!
!
!         ===============================
        ELSE
          CALL GWP_TO_poly(GWP,poly)
          CALL calc_auto_cor(auto_c,poly0,poly)
        END IF
        write(out_unitp,*) 'auto_c',auto_c
        CALL Write_AutoCorr(no,X,auto_c)
!       -------------------------------------------------

!       -------------------------------------------------
!       - analysis of the GWP + writing -----------------
        IF (.NOT. GWP%trajectory) THEN
!         extract sigma
          DO i=1,ND
            sigma(i)=1.d0/sqrt(imag(GWP%CAmean(i,i)))
          END DO
        END IF

!       output dynamical information
        xfs=x*0.0241d0

        write(out_unitp,31)x,xfs,param_LHA%d0Ene,energy
        write(10,31)x,xfs,param_LHA%d0Ene,energy
 31     format('energy:',4(1x,f18.6))

        write(31,3)x,xfs,GWP%Qmean
        write(71,3)x,xfs,GWP%Pmean
        IF (.NOT. GWP%trajectory) THEN
          write(32,3)x,xfs,GWP%CAmean
          write(41,3)x,xfs,sigma

          write(out_unitp0,*) x,GWP%CPZ
          write(out_unitp1,*) x,GWP%CZ
          write(out_unitp2,*) x,GWP%Cphase
        END IF
        write(out_unitp4,*) x,xfs,Y
        CALL flush_perso(out_unitp)
        CALL flush_perso(10)
        CALL flush_perso(31)
        CALL flush_perso(32)
        CALL flush_perso(41)
        CALL flush_perso(71)
        CALL flush_perso(60)
        CALL flush_perso(61)
        CALL flush_perso(62)
        CALL flush_perso(64)
!       -------------------------------------------------

        IF (X > XEND) EXIT
!       bonnes variables : Y GWP (it0)
        CALL RK4(Y,FU,ND,N,X,XPAS,YOUT,para_Q,mole,                     &
                 NdimA,onthefly,calc_EG,GWP,param_LHA,                  &
                 para_Tnum)
        Y(:) = YOUT(:)
!       bonnes variables : Y (it0+Delta)
!       Y => Q P A PZ Z
        CALL Y_TO_GWP(Y,GWP,N)
!       bonnes variables : Y GWP (it0+Delta)
!       CALL write_GWP(GWP)


        IF (.NOT. GWP%trajectory) THEN
!         Symmetrization de A
          None_symA = 0.
          DO i=1,ND
          DO j=i+1,ND
           None_symA=max(None_symA,abs(GWP%CAmean(i,j)-GWP%CAmean(j,i)))
            IF (Asym) THEN
              GWP%CAmean(i,j) = 0.5d0*(GWP%CAmean(i,j)+GWP%CAmean(j,i))
              GWP%CAmean(j,i) = GWP%CAmean(i,j)
            END IF
          END DO
          END DO
          write(out_unitp,*) ' None_symA ',None_symA
!         bonnes variables : GWP (it0+Delta)
!         INITIALISATION OF PZ_0 = 2 A_0
          GWP%CPZ(:,:) = 2.d0 * GWP%CAmean(:,:)
!         INITIALISATION OF Z_0 = Identity
          GWP%CZ(:,:) = 0.d0
          DO i=1,ND
            GWP%CZ(i,i) = 1.d0
          END DO
!         CALL write_GWP(GWP)

!         Normalization
!         transfer GWP to a complex polynomial (degree=2)
          CALL GWP_TO_poly(GWP,poly)
          CALL p1TOp2(p1=poly,p2=conj_poly,conjug=.TRUE.)
          CALL calc_auto_cor(auto_c,conj_poly,poly)
          write(out_unitp,*) 'norm of GWP:',real(auto_c,kind=Rkind)
          IF ( abs(auto_c-1.d0) > 0.2) THEN
            write(out_unitp,*) ' ERROR in propagation GWP'
            write(out_unitp,*) ' The norm is NOT conserved',                    &
                                   real(auto_c,kind=Rkind)
            STOP
          END IF
!         Renormalization of GWP
          IF (ReNorm) THEN
            GWP%Cphase = GWP%Cphase +eye*0.5d0 *                        &
                              log(real(auto_c,kind=Rkind))
            write(out_unitp,*) 'X,renorm phase',X,GWP%Cphase
          END IF
        END IF

        CALL GWP_TO_Y(GWP,Y,N)
!       bonnes variables : Y GWP (it0+Delta)
      END DO
!------------------------------------------------------
      CALL time_perso('GWP propa')

      write(9,*) Y(:)

      close(no)
      CALL time_perso('TF_autocorr')
      CALL TF_autocorr(para_propa)
      CALL time_perso('TF_autocorr')


      memory = size(d0g)
      deallocate(d0g,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d0g","main")
      memory = size(d1g)
      deallocate(d1g,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d1g","main")
      memory = size(d2g)
      deallocate(d2g,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d2g","main")
      memory = size(d0GG)
      deallocate(d0GG,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d0GG","main")
      memory = size(d1GG)
      deallocate(d1GG,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d1GG","main")
      memory = size(d2GG)
      deallocate(d2GG,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d2GG","main")

      CALL dealloc_CoordType(mole)
      CALL dealloc_param_Q(para_Q)
      memory = size(DeltaQsym)
      deallocate(DeltaQsym,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"DeltaQsym","main")
      memory = size(PQsym)
      deallocate(PQsym,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"PQsym","main")
      CALL dealloc_GWP(GWP)
      CALL dealloc_LHA(param_LHA)

      close(10)
      close(31)
      close(32)
      close(41)
      close(71)
      close(60)
      close(61)
      close(62)
      close(64)

      end program Gauss_num
      SUBROUTINE FCN(X,Y,FU,ND,N,NdimA,para_Q,mole,                     &
                     onthefly,lo_calc_EG,GWP,param_LHA,                 &
                     para_Tnum)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_GWP
      IMPLICIT NONE

      TYPE (para_GWP) :: GWP
      TYPE (para_LHA)  :: param_LHA
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
!     for Qsym Qact ....
      TYPE (param_Q) :: para_Q

      INTEGER N,ND,NDFIN
      INTEGER K,INDEX,I,j
      INTEGER INDEXI,INDEXF
      integer :: ncount
      INTEGER NDERIV,NDIMA
      logical :: onthefly,lo_calc_EG
      REAL(kind=Rkind) VEP
      REAL(kind=Rkind) X,TRHQP
      REAL(kind=Rkind) E
      REAL(kind=Rkind) Y(N),FU(N)
      REAL(kind=Rkind) d0g(ndimA,ndimA)
      REAL(kind=Rkind) d1g(ndimA,ndimA,mole%nb_act)
      REAL(kind=Rkind) d2g(ndimA,ndimA,mole%nb_act,mole%nb_act)
      REAL(kind=Rkind) d0GG(ndimA,ndimA)
      REAL(kind=Rkind) d1GG(ndimA,ndimA,mole%nb_act)
      REAL(kind=Rkind) d2GG(ndimA,ndimA,mole%nb_act,mole%nb_act)

      REAL(kind=Rkind) VTEMP(ND)
      REAL(kind=Rkind) HQ(ND),HP(ND)
      REAL(kind=Rkind) HQQ(ND,ND),HPP(ND,ND)
      REAL(kind=Rkind) HQP(ND,ND),HPQ(ND,ND)
      complex(kind=Rkind) TRACE,EYE
      complex(kind=Rkind) ZP(ND,ND)
      complex(kind=Rkind) PZP(ND,ND)
      complex(kind=Rkind) phaseP
      real(kind=Rkind) :: Lag

      REAL(kind=Rkind) HESS(ND,ND)


      EYE = cmplx(0.0,1.0,kind=Rkind)
!     write(out_unitp,*) 'FCN',Y
!     ------------------------------------------------
!     Y => Q P PZ Z (A)
      CALL Y_TO_GWP(Y,GWP,N)

      DO i=1,ND
        para_Q%Qsym(mole%liste_QactTOQdyn(i))=GWP%Qmean(i)
      END DO

      nderiv = 1
      IF (lo_calc_EG) nderiv = 2
      CALL dnOp_grid(para_Q%Qact,                                            &
                     param_LHA%d0Ene,param_LHA%d1Ene,HESS,              &
                     nderiv,                                            &
                     param_LHA%d0mu,param_LHA%d1mu,mole,                &
                     para_Tnum,onthefly)
!     CALL gr_hes(para_Q%Qsym,
!    *            param_LHA%d0Ene,param_LHA%d1Ene,HESS,ND)
      IF (lo_calc_EG) param_LHA%d2Ene(:,:) = HESS(:,:)
      param_LHA%d0Ene = param_LHA%d0Ene - param_LHA%pot0

      CALL calc2_d0d1d2g_G_bis(para_Q,                                  &
                               d0g,d1g,d2g,                             &
                               d0GG,d1GG,d2GG,nderiv,                   &
                               para_Tnum,mole)

      param_LHA%d0G(:,:)     = d0gg(1:nd,1:nd)
      param_LHA%d1G(:,:,:)   = d1GG(1:nd,1:nd,:)
      IF (lo_calc_EG) param_LHA%d2G(:,:,:,:) = d2GG(1:nd,1:nd,:,:)

!*****CALCUL DE HQ

      DO k=1,ND
        VTEMP(:) = MATMUL(param_LHA%d1G(:,:,k),GWP%Pmean)
        HQ(k)    =  DOT_PRODUCT(GWP%Pmean,VTEMP)/2.d0+param_LHA%d1Ene(k)
      END DO

!*****CALCUL DE HP

      HP(:) = MATMUL(param_LHA%d0G,GWP%Pmean)


      IF (.NOT. GWP%trajectory) THEN
!       CALCUL DE HQQ

        DO j=1,ND
        DO k=1,ND
         vtemp(:) = MATMUL(param_LHA%d2G(:,:,j,k),GWP%Pmean)
         HQQ(j,k) = DOT_PRODUCT(GWP%Pmean,VTEMP)/2.d0 +                 &
                     param_LHA%d2Ene(j,k)
        END DO
        END DO

!       CALCUL DE HPP

        HPP(:,:) = param_LHA%d0G(:,:)

!       CALCUL DE HQP and HPQ

        DO j=1,ND
        DO k=1,ND
          HQP(j,k) = DOT_PRODUCT(param_LHA%d1G(k,:,j),GWP%Pmean)
          HPQ(k,j) = HQP(j,k)
        END DO
        END DO
      END IF

!!#### CALCUL DES VITESSES


!*****CALCUL DE Q POINT

      FU(1:ND) = HP(:)

!*****CALCUL DE P POINT

      FU(ND+1:2*ND) = -HQ(:)

      IF (GWP%trajectory) RETURN

!*****CALCUL DE PZ POINT

      PZP(:,:) = -MATMUL(HQQ,GWP%CZ)-MATMUL(HQP,GWP%CPZ)

!*****CALCUL DE Z POINT

      ZP(:,:) = MATMUL(HPP,GWP%CPZ)+MATMUL(HPQ,GWP%CZ)

!     Real part of PZP
      ncount = 2*ND
      FU(ncount+1:ncount+ND**2)   =                                     &
             real(reshape(transpose(PZP(:,:)),(/ ND**2 /) ),kind=Rkind)
!      Imaginary part of PZP
      ncount = ncount + ND**2
      FU(ncount+1:ncount+ND**2)   =                                     &
             imag(reshape(transpose(PZP(:,:)),(/ ND**2 /) ))

!     Real part of PZP
      ncount = ncount + ND**2
      FU(ncount+1:ncount+ND**2)   =                                     &
             real(reshape(transpose(ZP(:,:)),(/ ND**2 /) ),kind=Rkind)
!     Imaginary part of PZP
      ncount = ncount + ND**2
      FU(ncount+1:ncount+ND**2)   =                                     &
             imag(reshape(transpose(ZP(:,:)),(/ ND**2 /) ))

!*****CALCUL DE phase POINT : i.tr(HPP.A)+i/2.tr(HPQ)-L
      IF (N == 2*(ND+2*ND**2+1)) THEN
         VTEMP(:) = MATMUL(param_LHA%d0G,GWP%Pmean)
         Lag=0.5d0*DOT_PRODUCT(GWP%Pmean,VTEMP)-param_LHA%d0Ene

         ZP(:,:)  = matmul(HPP,GWP%CAmean)
         phaseP = 0.d0
         DO i=1,ND
           phaseP = phaseP + ZP(i,i) + HPQ(i,i)*0.5d0
         END DO
         phaseP = eye * phaseP - Lag
         FU(N-1) = real(phaseP,kind=Rkind)
         FU(N) = imag(phaseP)
      END IF

!     write(out_unitp,*) 'FCN',FU

      end subroutine FCN

      SUBROUTINE RK4(Y,DYDX,ND,N,X,H,YOUT,para_Q,mole,                  &
                     NdimA,onthefly,calc_EG,GWP,param_LHA,              &
                     para_Tnum)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_GWP
      IMPLICIT NONE

      TYPE (para_GWP) :: GWP
      TYPE (para_LHA)  :: param_LHA

      TYPE (CoordType) mole
      TYPE (Tnum)    :: para_Tnum

!     for Qsym Qact ....
      TYPE (param_Q) :: para_Q

      REAL(kind=Rkind) pot0,ene
      logical :: onthefly,lo_calc_EG
      integer :: calc_EG

      INTEGER N,ND,I,ndimA
      REAL(kind=Rkind) H,HH,H6,X,XH,DQ1,dQ2
      REAL(kind=Rkind) Y(N),DYDX(N),YOUT(N),YT(N),DYT(N),DYM(N)

      HH = H*0.5D0
      H6 = H/6.d0
      XH = X + HH

      lo_calc_EG = .NOT. GWP%trajectory .AND. (calc_EG == 4)

      YT(:) = Y(:) + DYDX(:)*HH

      CALL FCN(XH,YT,DYT,ND,N,NdimA,para_Q,mole,                        &
               onthefly,lo_calc_EG,GWP,param_LHA,                       &
               para_Tnum)


      YT(:) = Y(:) + DYT(:)*HH

      CALL FCN(XH,YT,DYM,ND,N,NdimA,para_Q,mole,                        &
               onthefly,lo_calc_EG,GWP,param_LHA,                       &
               para_Tnum)

      YT(:)  = Y(:)   + H*DYM(:)
      DYM(:) = DYT(:) + DYM(:)

      CALL FCN(X+H,YT,DYT,ND,N,NdimA,para_Q,mole,                       &
               onthefly,lo_calc_EG,GWP,param_LHA,                       &
               para_Tnum)

      YOUT(:) = Y(:) + H6*( DYDX(:) + DYT(:) + 2.d0*DYM(:) )

      end subroutine RK4
      SUBROUTINE GWP_TO_Y(GWP,Y,N)
      USE mod_system
      USE mod_GWP
      IMPLICIT NONE

       TYPE (para_GWP) :: GWP
       integer         :: N
       real(kind=Rkind)    :: Y(N)
       integer         :: ND,ncount

       IF (.NOT. GWP%cplx .OR. .NOT. GWP%linearization) THEN
         write(out_unitp,*) ' ERROR in GWP_TO_Y'
         write(out_unitp,*) ' Real or without lineariztion impossible !'
         write(out_unitp,*) ' cplx,linearization',GWP%cplx,GWP%linearization
         STOP
       END IF

!      write(out_unitp,*) ' BEGINNING GWP_TO_Y'

       ND = GWP%ndim
       Y(:) = 0.d0

       ncount = 0
       Y(ncount+1:ncount+ND)      = GWP%Qmean(:)
       ncount = ncount + ND
       Y(ncount+1:ncount+ND)      = GWP%Pmean(:)

       IF (GWP%trajectory) RETURN

!      Real part of PZ
       ncount = ncount + ND
       Y(ncount+1:ncount+ND**2)   =                                     &
          real(reshape(transpose(GWP%CPZ(:,:)),(/ ND**2 /) ),kind=Rkind)
!      Imaginary part of PZ
       ncount = ncount + ND**2
       Y(ncount+1:ncount+ND**2)   =                                     &
             imag(reshape(transpose(GWP%CPZ(:,:)),(/ ND**2 /) ))

!      Real part of Z
       ncount = ncount + ND**2
       Y(ncount+1:ncount+ND**2)   =                                     &
           real(reshape(transpose(GWP%CZ(:,:)),(/ ND**2 /) ),kind=Rkind)
!      Imaginary part of Z
       ncount = ncount + ND**2
       Y(ncount+1:ncount+ND**2)   =                                     &
             imag(reshape(transpose(GWP%CZ(:,:)),(/ ND**2 /) ))

       IF (N == 2*(ND+2*ND**2+1)) THEN
!        phase
         ncount = ncount + ND**2
         Y(ncount+1) = real(GWP%Cphase,kind=Rkind)
         Y(ncount+2) = imag(GWP%Cphase)
       END IF

!      write(out_unitp,*) ' END GWP_TO_Y'

      end subroutine GWP_TO_Y
      SUBROUTINE Y_TO_GWP(Y,GWP,N)
      USE mod_system
      USE mod_GWP
      IMPLICIT NONE

       TYPE (para_GWP) :: GWP
       integer         :: N
       real(kind=Rkind)    :: Y(N)



       integer             :: ND,ncount,i
       real(kind=Rkind)    :: Rmat(GWP%ndim,GWP%ndim)
       real(kind=Rkind)    :: Imat(GWP%ndim,GWP%ndim)

       complex(kind=Rkind) :: ZI(GWP%ndim,GWP%ndim)
       complex(kind=Rkind) :: ZS(GWP%ndim,GWP%ndim)
       complex(kind=Rkind) :: trav(GWP%ndim)
       integer             :: inverse_index(GWP%ndim)

       IF (.NOT. GWP%cplx .OR. .NOT. GWP%linearization) THEN
         write(out_unitp,*) ' ERROR in Y_TO_GWP'
         write(out_unitp,*) ' Real or without lineariztion impossible !'
         write(out_unitp,*) ' cplx,linearization',GWP%cplx,GWP%linearization
         STOP
       END IF

!      write(out_unitp,*) ' BEGINNING Y_TO_GWP'

       ND = GWP%ndim



       ncount = 0
       GWP%Qmean(:)       = Y(ncount+1:ncount+ND)
       ncount = ncount + ND
       GWP%Pmean(:)       = Y(ncount+1:ncount+ND)

       IF (GWP%trajectory) RETURN

!      Real and imaginary parts of PZ
       ncount = ncount + ND
       Rmat(:,:) = reshape(Y(ncount+1:ncount+ND**2),(/ ND,ND /) )
       ncount = ncount + ND**2
       Imat(:,:) = reshape(Y(ncount+1:ncount+ND**2),(/ ND,ND /) )
       GWP%CPZ(:,:) = transpose(cmplx(Rmat,Imat,kind=Rkind))

!      Real and imaginary parts of Z
       ncount = ncount + ND**2
       Rmat(:,:) = reshape(Y(ncount+1:ncount+ND**2),(/ ND,ND /) )
       ncount = ncount + ND**2
       Imat(:,:) = reshape(Y(ncount+1:ncount+ND**2),(/ ND,ND /) )
       GWP%CZ(:,:) = transpose(cmplx(Rmat,Imat,kind=Rkind))


!      Calcul de A

        ZS(:,:) = GWP%CZ(:,:)
        CALL inversion_cplx(ZI,ZS,trav,inverse_index,ND)
!!      checking: is ZS below=Id?
        ZS(:,:) = matmul(GWP%CZ,ZI)
        DO i=1,ND
         ZS(i,i) = ZS(i,i) - 1.d0
        END DO
!       write(out_unitp,*) 'sum(CZ.ZI-Id):',sum(abs(ZS))

        GWP%CAmean(:,:) = 0.5d0 * matmul(GWP%CPZ,ZI)


       IF (N == 2*(ND+2*ND**2+1)) THEN
!        phase
         ncount = ncount + ND**2
         GWP%Cphase = cmplx(Y(ncount+1),Y(ncount+2),kind=Rkind)
       END IF

!      write(out_unitp,*) ' END Y_TO_GWP'

      end subroutine Y_TO_GWP
!==================================================
!
!     P2(Q) = eye*( DQt*A*DQ + P*DQ + phase) avec DQ=Q-Qmean
!
!==================================================
      SUBROUTINE GWP_TO_poly(GWP,poly)
      USE mod_system
      USE mod_poly
      USE mod_GWP
      IMPLICIT NONE

       TYPE (para_GWP)  :: GWP
       TYPE (para_poly) :: poly

       integer :: ind_exp(GWP%ndim)
       integer :: i,j,i_coef
       TYPE (para_poly) :: w1_poly

       complex (kind=Rkind) :: Ca,Cb
       real (kind=Rkind)    :: Ra,Rb

       complex (kind=Rkind) :: eye

       eye = cmplx(0.d0,1.d0,kind=Rkind)

       IF (.NOT. GWP%init0) THEN
         write(out_unitp,*) ' ERROR in GWP_TO_poly'
         write(out_unitp,*) ' GWP has NOT been initiated with init0_GWP'
         STOP
       END IF
       IF (GWP%trajectory) THEN
         write(out_unitp,*) ' ERROR in GWP_TO_poly'
         write(out_unitp,*) ' GWP%trajectory IS .TRUE.'
         write(out_unitp,*) ' You should NOT use this subroutine'
         STOP
       END IF
       IF (.NOT. poly%init0) THEN
         write(out_unitp,*) ' ERROR in GWP_TO_poly'
         write(out_unitp,*) ' poly has NOT been initiated with init0_poly'
         STOP
       END IF

!      write(out_unitp,*) ' BEGINNING GWP_TO_poly'
!      CALL write_GWP(GWP)


       CALL dealloc_poly(poly)
       CALL alloc_poly(poly,npoly=2,ndim=GWP%ndim,cplx=GWP%cplx)
       CALL init0_poly(w1_poly)

!      - First step GWP => poly en DQ=Q-Qmean
       IF (GWP%cplx) THEN
!        phase
         ind_exp(:) = 0
         i_coef = locate(poly,ind_exp)
         poly%Ccoef(i_coef) = eye * GWP%Cphase

!        Pmean*DQ
         DO i=1,GWP%ndim
           ind_exp(:) = 0
           ind_exp(i) = 1
           i_coef = locate(poly,ind_exp)
           poly%Ccoef(i_coef) = cmplx(0.d0,GWP%Pmean(i),kind=Rkind)
         END DO

!        DQ*Amean*DQ
         DO i=1,GWP%ndim
         DO j=1,GWP%ndim
           ind_exp(:) = 0
           ind_exp(i) = 1
           ind_exp(j) = ind_exp(j) + 1
           i_coef = locate(poly,ind_exp)
           poly%Ccoef(i_coef) = poly%Ccoef(i_coef) + eye*GWP%CAmean(i,j)
         END DO
         END DO
       ELSE
         STOP
       END IF

!      CALL write_poly(poly)


       DO i=1,GWP%ndim
         Ra = 1.d0
         Rb = -GWP%Qmean(i)
         CALL p1_linearTransfoTOp2(poly,i,Ra,Ca,Rb,Cb,                  &
                                   w1_poly,.FALSE.)
         CALL p1TOp2(p1=w1_poly,p2=poly)
       END DO
       CALL dealloc_poly(w1_poly)

!      CALL write_poly(poly)
!      write(out_unitp,*) ' END GWP_TO_poly'

      END SUBROUTINE GWP_TO_poly
!==================================================
!
!     Classical energy
!
!==================================================
      SUBROUTINE calc_ene(GWP,G,energyClas,ene0)
      USE mod_system
      USE mod_GWP
      IMPLICIT NONE

       TYPE (para_GWP), intent(inout)     :: GWP
       real (kind=Rkind), intent(inout)   :: energyClas
       real (kind=Rkind), intent(in)      :: G(GWP%ndim,GWP%ndim)
       real (kind=Rkind), intent(in)      :: ene0


       complex (kind=Rkind) :: vtemp(GWP%ndim)

       vtemp(:) = matmul(G,GWP%Pmean)
       energyClas = 0.5d0*dot_product(GWP%Pmean,vtemp) + ene0

      end subroutine calc_ene
!==================================================
!
!    Analytical autocorrelation
!    v1 work only if the coef of qi.qj are zero
!
!==================================================
      SUBROUTINE calc1_auto_cor(auto_c,poly0,poly)
      USE mod_system
      USE mod_poly
      IMPLICIT NONE

!.....DECLARE ARRAY Gauss
      TYPE (para_poly) :: poly,poly0
      TYPE (para_poly) :: w1_poly,w2_poly
      complex (kind=Rkind) :: Ca,Cb,auto_c
      complex (kind=Rkind) :: lnauto
      real (kind=Rkind)    :: Ra,Rb
      integer          :: ind_exp(poly0%ndim)
      integer          :: id0,id1,id2

!.....DECLARE DIVERS ------------------------------------------
      integer :: i,j,K
      real (kind=Rkind), parameter ::                                   &
       pi = 3.14159265358979323846264338327950288419716939937511d0


      CALL init0_poly(w1_poly)
      CALL init0_poly(w2_poly)

      CALL write_poly(poly0)
      CALL write_poly(poly)
      CALL p1PLUSp2TOp3(p1=poly0,p2=poly,p3=w2_poly)
      CALL write_poly(w2_poly)

      auto_c = 1.0
      lnauto = 0.0
      DO i=1,poly0%ndim
        Ca = 1.
        ind_exp(:) = 0
        ind_exp(i) = 1
        id1 = locate(w2_poly,ind_exp)
        ind_exp(:) = 0
        ind_exp(i) = 2
        id2 = locate(w2_poly,ind_exp)
!       Cb =  -0.5d0*w2_poly%Ccoef(id1)/w2_poly%Ccoef(id2)
!       CALL p1_linearTransfoTOp2(w2_poly,i,1.d0,Ca,0.d0,Cb,
!    *                            w1_poly,.TRUE.)
!       CALL p1TOp2(p1=w1_poly,p2=w2_poly)
        lnauto = lnauto-0.25d0*w2_poly%Ccoef(id1)**2/w2_poly%Ccoef(id2)
        auto_c = auto_c * sqrt(-pi/w2_poly%Ccoef(id2))
      END DO
      ind_exp(:) = 0
      id0 = locate(w2_poly,ind_exp)
      lnauto = lnauto + w2_poly%Ccoef(id0)
      auto_c = auto_c * exp(lnauto)
      write(out_unitp,*) 'auto_c',auto_c
      CALL write_poly(w2_poly)

      CALL dealloc_poly(w1_poly)
      CALL dealloc_poly(w2_poly)

      end subroutine calc1_auto_cor
!==================================================
!
!    Analytical autocorrelation
!    v2
!
!==================================================
      SUBROUTINE calc_auto_cor(auto_c,poly0,poly)
      USE mod_system
      USE mod_poly
      IMPLICIT NONE

!.....DECLARE ARRAY Gauss
      TYPE (para_poly) :: poly,poly0
      TYPE (para_poly) :: LTpoly,NewPoly,dNewPoly
      TYPE (para_poly) :: w1_poly
      complex (kind=Rkind) :: Ca,auto_c
      integer          :: ind_exp(poly0%ndim)
      integer          :: id0,id1,id2

!.....DECLARE DIVERS ------------------------------------------
      integer :: i,j,K
      real (kind=Rkind), parameter ::                                   &
       pi = 3.14159265358979323846264338327950288419716939937511d0


      CALL init0_poly(NewPoly)
      CALL init0_poly(dNewPoly)
      CALL init0_poly(LTpoly)
      CALL init0_poly(w1_poly)

      CALL p1PLUSp2TOp3(p1=poly0,p2=poly,p3=NewPoly)
!     write(out_unitp,*) ' p0PLUSp'
!     CALL write_poly(NewPoly)


      auto_c = 1.0
      DO i=1,poly0%ndim
!       - LTpoly = a qi^2 + qi dNewPoly + Remainder (degree =2)

!       - partial LTpoly
!       - LTpoly = a qi^2 + Remainder (degree =2)
!       => remove all terms with qi^1
        CALL p1TOp2(p1=NewPoly,p2=LTpoly)
        DO k=1,LTpoly%nb_coef
          IF (LTpoly%ind(i,k) == 1) LTpoly%Ccoef(k) = 0.d0
        END DO
!       write(out_unitp,*) ' partial LTpoly',i
!       CALL write_poly(LTpoly)

!       - linear variable (i) transformation with dp0PLUSp
        CALL d1p1TOp2(p1=NewPoly,p2=dNewPoly,id=i)
        ind_exp(:) = 0
        ind_exp(i) = 1
        id1 = locate(dNewPoly,ind_exp)
        dNewPoly%Ccoef(id1) = 0.d0
!       write(out_unitp,*) ' dNewPoly',i
!       CALL write_poly(dNewPoly)


!       - calc : dNewPoly^2/(-4 a)
        CALL p1TIMEp2TOp3(p1=dNewPoly,p2=dNewPoly,p3=w1_poly)
        ind_exp(:) = 0
        ind_exp(i) = 2
        id2 = locate(NewPoly,ind_exp)
        Ca =  -0.25d0/NewPoly%Ccoef(id2)
        w1_poly%Ccoef(:) = w1_poly%Ccoef(:) * Ca
!       write(out_unitp,*) ' dNewPoly^2/(-4 a)'
!       CALL write_poly(w1_poly)

!       - full LTpoly
!       - LTpoly = a ui^2 + Remainder (degree =2)
        CALL p1PLUSp2TOp3(p1=LTpoly,p2=w1_poly,p3=NewPoly)
!       write(out_unitp,*) ' full LTpoly',i
!       CALL write_poly(NewPoly)

!       - integration along ui
        auto_c = auto_c * sqrt(-pi/NewPoly%Ccoef(id2))
      END DO

      ind_exp(:) = 0
      id0 = locate(NewPoly,ind_exp)
      auto_c = auto_c * exp(NewPoly%Ccoef(id0))
!     write(out_unitp,*) 'auto_c',auto_c

      CALL dealloc_poly(w1_poly)
      CALL dealloc_poly(LTpoly)
      CALL dealloc_poly(NewPoly)
      CALL dealloc_poly(dNewPoly)

      end subroutine calc_auto_cor
!------------------------------------
!     V = v0 +1/2Sum hes(i,j) (Qi-Qeqi)(Qj-Qeqj)
!
!------------------------------------
      SUBROUTINE gr_hes(Q,ene,gr,hes,n)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind) :: Q(n),gr(n),hes(n,n)
      real (kind=Rkind) :: Qeq(n),DQ(n)
      real (kind=Rkind) :: ene
      integer           :: n

      integer :: i,k


      hes(:,:) = 0.
      DO i=1,n
        hes(i,i) = 0.5
        Qeq(i) = 2.1
      END DO
      DQ(:) = Q(:)-Qeq(:)

      DO k=1,n
        gr(k) = dot_product(hes(:,k),DQ)
      END DO

      ene = -10.d0 + 0.5d0*dot_product(gr,DQ)

      end subroutine gr_hes

