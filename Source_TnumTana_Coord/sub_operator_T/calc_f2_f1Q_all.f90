!
!======================================================================
!
!      T(2D) de NH3
!
!======================================================================
!
!     CALL calc_f2_f1Q_NH3(Q,Tdef2,Tdef1,vep,rho,nb_act,nb_var)
      SUBROUTINE calc_f2_f1Q_NH3(Q,f2Q,f1Q,vep,rho,nq,nb_var)
      USE mod_system
      IMPLICIT NONE


       integer nb_var,nq
       real(kind=Rkind) Q(nb_var)
       real(kind=Rkind) vep,rho
       real(kind=Rkind) f1Q(nq)
       real(kind=Rkind) f2Q(nq,nq)


       real(kind=Rkind) mN,mH,qq
       parameter (mN=25525.4342723654045_Rkind)
       parameter (mH=1837.10882307318116_Rkind)
       parameter (qq=mN/mH)


       real(kind=Rkind) x,R,x2,x3,x4,x6,R2

       integer i,j


      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING calc_f2_f1Q'
      write(out_unitp,*) 'Q',Q
      END IF
!---------------------------------------------------------------------

       IF (nq .NE.  2) THEN
         write(out_unitp,*) ' ERROR in calc_f2_f1Q_NH3'
         write(out_unitp,*) ' nq should be equal to 2',nq
         STOP
       END IF

       x  = Q(2)
       R  = Q(1)

       x2 = x*x
       x3 = x*x2
       x4 = x2*x2
       x6 = x4*x2

       R2 = R*R

!      veritables fonctions f2Q et f1Q
       f1Q(1)   = x*(15.d0+2.d0*qq-15.d0*x2)/(6.d0*mN*R2)
       f1Q(2)   = (3.d0*x2-1.d0)/(2.d0*mN*R)

       f2Q(1,1) = (3.d0+qq-3.d0*x2)*(x2-1.d0)/(6.d0*mN*R2)
       f2Q(2,2) = -(qq+3.d0*x2)/(6.d0*mN)
       f2Q(1,2) = -(x-x3)/(mN*R)
       f2Q(2,1) = f2Q(1,2)


       vep = 3.d0*(1.d0+qq)*(3.d0+qq)**2+(qq-15.d0)*(3.d0+qq)**2 * x2-  &
             3.d0*(qq-3.d0)*(21.d0+5.d0*qq) * x4-                       &
             9.d0*(qq-3.d0)**2 * x6
       vep = vep/( 6.d0*mN*(3.d0+qq+(qq-3.d0)*x2)**2 * R2 )
       rho = 1.d0

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'x R',x,R
      write(out_unitp,*) 'f1Q',f1Q
      write(out_unitp,*) 'f2Q',f2Q
      write(out_unitp,*) 'vep,rho',vep,rho
      write(out_unitp,*) 'END calc_f2_f1Q'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc_f2_f1Q_NH3
!
!======================================================================
!
!      Evalue pour un jeu de coordonnes d0Q(nb_var) les valeurs
!      des fonctions f2Q et f1Q.
!      Ces fonctions sont "devant" les elements :
!           f2Q(Q1,Q2) d2./dQ1dQ2
!      et   f1Q(Q1) d./dQ1
!
!      -------------------------------------------------------
!      -------------------------------------------------------
!      pour une molecule triatomique en Jacobi :
!      avec un element de volume dT = dx dR dr
!          x = cos(T)
!
!
!                       T   A2
!                  A1_____ / r
!                     R   /
!                        A3
!      ordre des variables x(T), R, r eq : Q(1), Q(2) Q(3)
!      avec x = Q(1) = cos(T)
!      masses dans l'ordre des atome A1, A2, A3
!
!      1D
!
!======================================================================
!
      SUBROUTINE calc_f2_f1Q_jac1(Q,f2Q,f1Q,vep,rho,nq,nb_var)
      USE mod_system
      IMPLICIT NONE


       integer nb_var,nq
       real(kind=Rkind) Q(nb_var)
       real(kind=Rkind) vep,rho
       real(kind=Rkind) f1Q(nb_var)
       real(kind=Rkind) f2Q(nb_var,nb_var)


       real(kind=Rkind) MM_inv,m_inv,inertie_inv
       real(kind=Rkind) x,RR,r

       integer i,j


       integer nb_at
       parameter (nb_at=3)
       real(kind=Rkind) masses(nb_at)
       data masses/1837.10882307318116d0,21874.1407256995735d0,         &
                   25525.4342723654045d0/

       SAVE masses

      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING calc_f2_f1Q'
      write(out_unitp,*) 'Q',Q
      write(out_unitp,*) 'masses',masses
      END IF
!---------------------------------------------------------------------

       IF (nq .NE. 1) THEN
         write(out_unitp,*) ' ERROR in calc_f2_f1Q_cart'
         write(out_unitp,*) ' nq should be equal to 1',nq
         STOP
       END IF

       x  = Q(1)
       RR = Q(2)
       r  = Q(3)

!      variable temporaire 1/M, 1/m et (1/MR2 + 1/mr2)
       MM_inv      = 1.d0/masses(1) + 1.d0/(masses(2)+masses(3))
       m_inv       = 1.d0/masses(2) + 1.d0/masses(3)
       inertie_inv = MM_inv/(RR*RR) + m_inv/(r*r)

!      veritables fonctions f2Q et f1Q
       f1Q(1)   =          inertie_inv * x
       f2Q(1,1) = -0.5d0 * inertie_inv * (1.d0-x*x)


       vep = 0.d0
       rho = 1.d0

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'x RR r',x,RR,r
      write(out_unitp,*) 'M RR',1.d0/MM_inv,RR
      write(out_unitp,*) 'm r',1.d0/m_inv,r
      write(out_unitp,*) 'inertie_inv',inertie_inv
      write(out_unitp,*) 'f1Q',f1Q
      write(out_unitp,*) 'f2Q',f2Q
      write(out_unitp,*) 'END calc_f2_f1Q'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc_f2_f1Q_jac1
!
!======================================================================
!
!      Evalue pour un jeu de coordonnes d0Q(nb_var) les valeurs
!      des fonctions f2Q et f1Q.
!      Ces fonctions sont "devant" les elements :
!           f2Q(Q1,Q2) d2./dQ1dQ2
!      et   f1Q(Q1) d./dQ1
!
!      -------------------------------------------------------
!      -------------------------------------------------------
!      pour une molecule triatomique en Jacobi :
!      avec un element de volume dT = dx dR dr
!          x = cos(T)
!
!
!                       T   A2
!                  A1_____ / r
!                     R   /
!                        A3
!      ordre des variables x(T), R, r eq : Q(1), Q(2) Q(3)
!      avec x = Q(1) = cos(T)
!      masses dans l'ordre des atome A1, A2, A3
!
!======================================================================
!
      SUBROUTINE calc_f2_f1Q_jac3(Q,f2Q,f1Q,vep,rho,nq,nb_var)
      USE mod_system
      IMPLICIT NONE


       integer nb_var,nq
       real(kind=Rkind) Q(nb_var)
       real(kind=Rkind) vep,rho
       real(kind=Rkind) f1Q(nb_var)
       real(kind=Rkind) f2Q(nb_var,nb_var)


       real(kind=Rkind) MM_inv,m_inv,inertie_inv
       real(kind=Rkind) x,RR,r

       integer i,j


       integer nb_at
       parameter (nb_at=3)
       real(kind=Rkind) masses(nb_at)
       data masses/1800.d0,1800.d0,1800.d0/

       SAVE masses

      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING calc_f2_f1Q'
      write(out_unitp,*) 'Q',Q
      write(out_unitp,*) 'masses',masses
      END IF
!---------------------------------------------------------------------

       IF (nb_var .NE. nq) THEN
         write(out_unitp,*) ' ERROR in calc_f2_f1Q_cart'
         write(out_unitp,*) ' nq should be equal to nb_var',nq,nb_var
         STOP
       END IF

       IF (nq .NE. nb_at) THEN
         write(out_unitp,*) ' ERROR in calc_f2_f1Q_cart'
         write(out_unitp,*) ' nq should be equal to nb_at',nq,nb_at
         STOP
       END IF
!      -------------------------------------------------------
!      inititialisation
       DO i=1,nb_var
         f1Q(i) = 0.d0
         DO j=1,nb_var
           f2Q(i,j) = 0.d0
         END DO
       END DO
!      -------------------------------------------------------

       x  = Q(1)
       RR = Q(2)
       r  = Q(3)

!      variable temporaire 1/M, 1/m et (1/MR2 + 1/mr2)
       MM_inv      = 1.d0/masses(1) + 1.d0/(masses(2)+masses(3))
       m_inv       = 1.d0/masses(2) + 1.d0/masses(3)
       inertie_inv = MM_inv/(RR*RR) + m_inv/(r*r)

!      veritables fonctions f2Q et f1Q
       f1Q(1)   =          inertie_inv * x
       f2Q(1,1) = -0.5d0 * inertie_inv * (1.d0-x*x)
       f2Q(2,2) = -0.5d0 * MM_inv
       f2Q(3,3) = -0.5d0 * m_inv


!      pour un test avec f2Q non diago
!      f2Q(2,3) = 0.2d0 * (f2Q(2,2)+f2Q(3,3))
!      f2Q(3,2) = f2Q(2,3)


       vep = 0.d0
       rho = 1.d0

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'x RR r',x,RR,r
      write(out_unitp,*) 'M RR',1.d0/MM_inv,RR
      write(out_unitp,*) 'm r',1.d0/m_inv,r
      write(out_unitp,*) 'inertie_inv',inertie_inv
      write(out_unitp,*) 'f1Q',f1Q
      write(out_unitp,*) 'f2Q',f2Q
      write(out_unitp,*) 'END calc_f2_f1Q'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc_f2_f1Q_jac3
!======================================================================
!
!      Evalue pour un jeu de coordonnes d0Q(nb_var) les valeurs
!      des fonctions f2Q et f1Q.
!      Ces fonctions sont "devant" les elements :
!           f2Q(Q1,Q2) d2./dQ1dQ2
!      et   f1Q(Q1) d./dQ1
!
!      -------------------------------------------------------
!      -------------------------------------------------------
!      pour une molecule triatomique en Jacobi :
!      avec un element de volume dT = sin(a)da dR dr
!
!
!                           C
!                          / r
!                       a /
!                  B____A/
!                     RR
!
!      ordre des variables RR, a, r eq : Q(1), Q(2) Q(3)
!      masses dans l'ordre des atome H C et N
!
!      Rq: la variable r est constante
!
!======================================================================
!
      SUBROUTINE calc_f2_f1Q_val(Q,f2Q,f1Q,vep,rho,nq,nb_var)
      USE mod_system
      IMPLICIT NONE

       integer nq,nb_var

       real(kind=Rkind) Q(3)
       real(kind=Rkind) vep,rho
       real(kind=Rkind) f1Q(2)
       real(kind=Rkind) f2Q(2,2)

       real(kind=Rkind) mA,mB,mC,mAB,mAC,mBC,mABC
       parameter (mA=1837.10882307318116)
       parameter (mB=21874.1407256995735)
       parameter (mC=25525.4342723654045)
       parameter (mAB=mA+mB,mBC=mB+mC,mAC=mA+mC)
       parameter (mABC=mA+mB+mC)


       real(kind=Rkind) a,RR,r
       real(kind=Rkind) sa,ca,c2a

       integer i,j

!---------------------------------------------------------------------
      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING calc_f2_f1Q'
      write(out_unitp,*) 'Q',Q
      END IF
!---------------------------------------------------------------------


       RR = Q(1)
       a  = Q(2)
       r  = Q(3)

!     write(out_unitp,*) 'RR a r',RR,a,r
!     write(out_unitp,*) 'RR a r',Q(1),Q(2),Q(3)
!     STOP

       sa = sin(a)
       ca = cos(a)
       c2a = cos(a+a)

       rho=sa




       f2Q(1,1) = 0.5d0*(-1.d0/mA-1.d0/mB+mC/(mA*mAC)*ca*ca )
       f2Q(1,2) = -sa*(mC*r*ca-mAC*RR) / (mA*(mAC)*r*RR)
       f2Q(2,1) = f2Q(1,2)
       f2Q(2,2) = mAC*(mAB*mC*r*r+mB*RR*(-2.d0*mC*r*ca+mAC*RR))
       f2Q(2,2) = f2Q(2,2) - mB*mC*mC*r*r*sa*sa
       f2Q(2,2) = -f2Q(2,2)/(2.d0*mA*mB*mC*mAC*r*r*RR*RR)

       f1Q(1) = ca/(r*mA) - mC*(1.d0+3.d0*c2a)/(4.d0*mAC*mA*RR)
       f1Q(2) = mAC*ca/(mC*r*r*sa)
       f1Q(2) = f1Q(2) + (mA*mA-mB*mC+mA*mBC + 2.d0*mB*mC*c2a)*         &
                            ca/(mB*mAC*RR*RR*sa)
       f1Q(2) = f1Q(2) - 2.d0*c2a/(sa*RR*r)
       f1Q(2) = -f1Q(2)/(mA+mA)


       vep = (mC*r+3.d0*mC*r*c2a-4.d0*mAC*ca*RR)/(4.d0*mA*mAC*r*RR*RR)


!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'Q',Q
      write(out_unitp,*) 'rho vep',rho,vep
      write(out_unitp,*) 'f1Q',f1Q
      write(out_unitp,*) 'f2Q',f2Q
      write(out_unitp,*) 'END calc_f2_f1Q'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine calc_f2_f1Q_val

