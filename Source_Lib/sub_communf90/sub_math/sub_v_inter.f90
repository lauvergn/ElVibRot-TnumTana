!================================================================
!    fonction v_inter(x)
!
!    v6  1D Huffaker
!    v7  1D BO
!    v8  1D (x-Req)/(c(1)*x+c(2)*Req)
!     c(1)=0 c(2)=1 : Dunham
!     c(1)=1 c(2)=0 : SPF
!     c(1)=0 c(2)=1 : Ogilvie
!    v9  1D x**i*(x-1)**k
!
!    v12 1D poly legendre          en cos(x)
!    v13 1D poly legendre paire    en cos(x)
!    v14 1D poly legendre impaire  en cos(x)
!
!    v15 1D x^n
!    v16 1D x^n paire
!    v17 1D x^n impaire
!    v18 1D x^n paire sans n=0
!    v19 1D x^n sans n=0
!
!    v22 1D poly legendre          en x
!    v23 1D poly legendre paire    en x
!    v24 1D poly legendre impaire  en x
!
!    v26 1D serie de fourier en x
!================================================================

       real*8 FUNCTION v_inter2(x,ndim,ntyp,nom,nsurf)
       USE mod_system
       IMPLICIT NONE

       integer max_dim,nsurf,newsurf
       parameter (max_dim=6)
       real*8 x(ndim)
       integer n(0:max_dim)

!      - function ---------------------------------------------

       real*8 v1_inter,v2_inter
       real*8 vgene_inter,vgene_inter2

       real*8   v1,v2,v3,v4,v_poly,v6,v7,v8,v9
!      external v1,v2,v3,v4,v_poly,v6,v7,v8,v9
       real*8   v10,v11,v12,v13,v14,v15,v16,v17,v18,v19
!      external v10,v11,v12,v13,v14,v15,v16,v17,v18,v19
       real*8   v20,v21,v22,v23,v24,v25,v26,v27,v28,v29
!      external v20,v21,v22,v23,v24,v25,v26,v27,v28,v29
       real*8   v30,v31,v32,v33,v34,v35,v36,v37,v38,v39
!      external v30,v31,v32,v33,v34,v35,v36,v37,v38,v39
!      ---------------------------------------------------------

       integer max_fit
       parameter (max_fit=2000)
       real*8 F(max_fit)

       integer nn,ntyp,ndim

       character*14 nom
       logical begin,exist
       data begin/.true./
       SAVE begin,n,F,newsurf

!---------------------------------------------------------------
!      initialisation la premiere fois

       IF (begin .OR. nsurf .NE. newsurf) THEN
         newsurf = nsurf
         CALL read_para1d(F,nn,n,ndim,max_fit,nom,exist,nsurf)
         write(out_unitp,*) nom,nn,n,ndim,nsurf
         IF ( .NOT. exist) STOP
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois
!---------------------------------------------------------------

       nn = n(0)
!      write(out_unitp,*) 'BEGINING v_inter2',ntyp,ndim,nn,n

       SELECT CASE (ntyp)
         CASE (1)
           v_inter2 = vgene_inter(x,ndim,v1,F,nn,n)
         CASE (2)
           v_inter2 = vgene_inter(x,ndim,v2,F,nn,n)
         CASE (3)
           v_inter2 = vgene_inter(x,ndim,v3,F,nn,n)
         CASE (4)
          !v_inter2=vgene_inter(x,ndim,v4,F,nn,n)
           v_inter2=vgene_inter2(x,ndim,F,nn,n,ntyp)
         CASE (5)
          !v_inter2=vgene_inter(x,ndim,v_poly,F,nn,n)
          v_inter2=vgene_inter2(x,ndim,F,nn,n,ntyp)
         CASE (6)
           v_inter2 = vgene_inter(x,ndim,v6,F,nn,n)
         CASE (7)
           v_inter2 = vgene_inter(x,ndim,v7,F,nn,n)
         CASE (8)
           v_inter2 = vgene_inter(x,ndim,v8,F,nn,n)
         CASE (9)
           v_inter2 = vgene_inter(x,ndim,v9,F,nn,n)
         CASE (10)
           v_inter2 = vgene_inter(x,ndim,v10,F,nn,n)
         CASE (11)
           v_inter2 = vgene_inter(x,ndim,v11,F,nn,n)

         CASE (12)
           v_inter2=vgene_inter2(x,ndim,F,nn,n,ntyp)
         CASE (13)
           v_inter2 = vgene_inter(x,ndim,v13,F,nn,n)
         CASE (14)
           v_inter2 = vgene_inter(x,ndim,v14,F,nn,n)
         CASE (15)
           v_inter2 = vgene_inter(x,ndim,v15,F,nn,n)
         CASE (16)
           v_inter2 = vgene_inter(x,ndim,v16,F,nn,n)
         CASE (17)
           v_inter2 = vgene_inter(x,ndim,v17,F,nn,n)
         CASE (18)
           v_inter2 = vgene_inter(x,ndim,v18,F,nn,n)
         CASE (19)
           v_inter2 = vgene_inter(x,ndim,v19,F,nn,n)


         CASE (22)
           v_inter2 = vgene_inter(x,ndim,v22,F,nn,n)
         CASE (23)
           v_inter2 = vgene_inter(x,ndim,v23,F,nn,n)
         CASE (24)
           v_inter2 = vgene_inter(x,ndim,v24,F,nn,n)


         CASE (26)
           v_inter2 = vgene_inter(x,ndim,v26,F,nn,n)
         CASE (27)
           v_inter2 = vgene_inter(x,ndim,v27,F,nn,n)
         CASE (28)
           v_inter2 = vgene_inter(x,ndim,v28,F,nn,n)
         CASE (29)
           v_inter2 = vgene_inter(x,ndim,v29,F,nn,n)
         CASE (30)
           v_inter2 = vgene_inter(x,ndim,v30,F,nn,n)

         CASE (32)
           v_inter2 = vgene_inter(x,ndim,v32,F,nn,n)
         CASE (33)
           v_inter2 = vgene_inter(x,ndim,v33,F,nn,n)
         CASE (34)
           v_inter2 = vgene_inter(x,ndim,v34,F,nn,n)

         CASE Default
           v_inter2=ZERO
       END SELECT



!      write(out_unitp,*) 'END v_inter2',x,v_inter2

       end function v_inter2
!================================================================
!    fonction v1_inter(x,ndim)
!================================================================

       FUNCTION v1_inter(x,ndim)
       USE mod_system
       IMPLICIT NONE

       real*8 :: v1_inter


       integer ndim
       real*8 x(ndim)
       integer nn,n(0:ndim)

       real*8   v1
       external v1

       integer max_fit
       parameter (max_fit=2000)
       real*8 F(max_fit)
       real*8 vgene_inter

       character*14 nom
       logical begin,exist
       data begin/.true./
       SAVE


!---------------------------------------------------------------
!      initialisation la premiere fois
       IF (begin) THEN
         nom='inter-t1'
         CALL read_para1d(F,nn,n,ndim,max_fit,nom,exist,1)
         IF ( .NOT. exist) THEN
           nn=1
           F(1)=ONE
         END IF
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois
!---------------------------------------------------------------
!

       v1_inter = vgene_inter(x,ndim,v1,F,nn,n)

       RETURN
       end function v1_inter
!================================================================
!    fonction v2_inter(x)
!================================================================

       FUNCTION v2_inter(x,ndim)
       USE mod_system
       IMPLICIT NONE

       real*8 :: v2_inter

       integer ndim
       real*8 x(ndim)
       integer n(0:ndim),nn

       real*8   v2
       external v2

       integer max_fit
       parameter (max_fit=2000)
       real*8 F(max_fit)

       real*8 vgene_inter

       character*14 nom
       logical begin,exist
       data begin/.true./
       SAVE


!---------------------------------------------------------------
!      initialisation la premiere fois
       IF (begin) THEN
         nom='inter-t2'
         CALL read_para1d(F,nn,n,ndim,max_fit,nom,exist,1)
         IF ( .NOT. exist) THEN
           nn=0
         END IF
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois
!---------------------------------------------------------------

       v2_inter = vgene_inter(x,ndim,v2,F,nn,n)


       RETURN
       end function v2_inter
!================================================================
!    fonction vgene_inter(x,ndim) 1 D
!================================================================
       FUNCTION vgene_inter(x,ndim,v_typ,F,nn,n)
       USE mod_system
       IMPLICIT NONE

       real*8 :: vgene_inter


       integer   ndim
       real*8    x(ndim)
       integer   nn,n(0:ndim)
       real*8    z
       real*8    F(nn)
       integer   kl

!      - function ----------------------------------------------
       real*8 v_typ
!      ---------------------------------------------------------

!      write(out_unitp,*) 'BEGINING vgene_inter',ndim,nn,n


       z=ZERO
       DO kl=1,nn
         z = z + F(kl) * v_typ(x,ndim,kl,n)
!        write(out_unitp,*) z,F(kl),ndim,nn,n
       END DO

       vgene_inter = z

!      write(out_unitp,*) 'END vgene_inter',x,z,ndim,nn,n

       end function vgene_inter
!================================================================
!    fonction vgene_inter(x,ndim) 1 D
!================================================================
       FUNCTION vgene_inter2(x,ndim,F,nn,n,ntyp)
       USE mod_system
       IMPLICIT NONE

       real(kind=Rkind) :: vgene_inter


       integer   ndim
       real(kind=Rkind) ::    x(ndim)
       integer   ntyp,nn,n(0:ndim)
       real(kind=Rkind) ::    F(nn)

       integer   kl
       real(kind=Rkind) ::    z

!      - function ----------------------------------------------
       real(kind=Rkind) :: v
!      ---------------------------------------------------------

!      write(out_unitp,*) 'BEGINING vgene_inter2',ndim,nn,n
!      write(out_unitp,*) 'F(:)',F
!      write(out_unitp,*) 'x',x


       z=ZERO
       DO kl=1,nn
         z = z + F(kl) * v(x,ndim,kl,n,ntyp)
         !write(out_unitp,*) z,F(kl),ndim,nn,n
         !flush(6)
       END DO

       vgene_inter2 = z

!      write(out_unitp,*) 'END vgene_inter2',x,z,ndim,nn,n

       end function vgene_inter2

