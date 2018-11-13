!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong, Josep Maria Luis
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!===========================================================================
!===========================================================================

      SUBROUTINE TF_autocorr(para_propa)
      USE mod_system
      USE mod_propa
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      real (kind=Rkind)  :: microDeltaT


      complex (kind=Rkind),pointer :: C(:),W(:)
      real (kind=Rkind) :: E,DE,T,Tmax
      real (kind=Rkind) :: A,B,cmc


      integer :: no,ni

      integer :: NPT,NPT2,I,max_iE
      integer :: necr

      character (len=Name_len) :: name
      real (kind=Rkind)        :: ca,cb,cc,auTOcm_inv
      integer :: err_mem,memory

      CALL file_open(para_propa%file_autocorr,ni)
      CALL file_open(para_propa%file_spectrum,no)
      write(out_unitp,*) 'ni,no',ni,no

      microDeltaT = para_propa%WPdeltaT/para_propa%nb_micro
      auTOcm_inv  = get_Conv_au_TO_unit('E','cm-1')

      NPT  = para_propa%WPTmax/microDeltaT
      NPT2 = 2**para_propa%TFnexp2
      IF (NPT2 < NPT) THEN
        write(out_unitp,*) ' ERROR in TF_autocorr'
        write(out_unitp,*) ' NPT2(=2^TFnexp2) is inferior to npt',NPT2,NPT
        write(out_unitp,*) ' Increase "TFnexp2" in the "&propa" namlist'
        STOP
      END IF

      nullify(C)
      nullify(W)
      CALL alloc_array(C,(/NPT2/),"C","TF_autocorr")
      CALL alloc_array(W,(/NPT2/),"C","TF_autocorr")

      write(out_unitp,*)  'NPT,NPT2,microWPdeltaT',NPT,NPT2,microDeltaT
      CALL flush_perso(out_unitp)

      DO I=1,NPT
        CALL Read_AutoCorr(ni,t,C(I))
        !write(out_unitp,*) t,C(I) ; CALL flush_perso(out_unitp)
      END DO
      C(1) = C(1) * HALF
      close(ni)
      C(NPT+1:NPT2) = cmplx(ZERO,ZERO,kind=Rkind)


      TMAX = NPT2*microDeltaT
      DE   = TWO*pi/TMAX * auTOcm_inv

      write(out_unitp,*) 'DeltaE(TF)  (cm-1)',DE
      write(out_unitp,*) 'TFmaxE  (cm-1)', para_propa%TFmaxE*auTOcm_inv
      write(out_unitp,*) 'TFminE  (cm-1)', para_propa%TFminE*auTOcm_inv

      CALL PREFFT(NPT2,1,para_propa%TFnexp2,W)
      CALL FFT(NPT2,1,microDeltaT,para_propa%TFnexp2,W,C)

      E = para_propa%TFminE*auTOcm_inv
      max_iE = min(NPT2,int(para_propa%TFmaxE * auTOcm_inv/DE))
      IF (max_iE == 0) max_iE = NPT2
      DO i=1,max_iE
        write(no,33) E,real(C(i),kind=Rkind),aimag(C(i)),abs(C(i))
        E = E + DE
      END DO
      close(no)

  33  FORMAT(5(E13.6,' '))

      CALL dealloc_array(C,"C","TF_autocorr")
      CALL dealloc_array(W,"C","TF_autocorr")

      END SUBROUTINE TF_autocorr

