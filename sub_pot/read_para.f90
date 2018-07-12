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

!
!================================================================
!    read the parameters of the fits
!================================================================
      SUBROUTINE read_para0d(F,nn,max_points,nom1,exist)
      USE mod_system
      USE mod_file
      IMPLICIT NONE

       integer :: nn
       character (len=*) :: nom1
       integer :: max_points
       real (kind=Rkind) :: F(max_points)
       logical :: exist

       integer :: no,ios,kl

       write(out_unitp,*) nom1
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)

       IF (ios .EQ. 0) THEN
         read(no,*) nn
         IF (nn .GT. max_points) THEN
           write(out_unitp,*) ' ERROR : nb de points du fit (',nn,') >'
           write(out_unitp,*) '         a max_points (',max_points,')'
           write(out_unitp,*) '         STOP in read_para0d'
           STOP
         END IF
         DO kl=1,nn
           read(no,*) F(kl)
!          write(out_unitp,*) F(kl)
         END DO
         CLOSE(no)
         exist = .TRUE.
       ELSE
         write(out_unitp,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


!      write(out_unitp,*) 'nom1,exist ',nom1,exist
!      write(out_unitp,*) 'nn,max_points',nn,max_points

       RETURN
       end subroutine read_para0d
!
!================================================================
!    read the parameters of the fits
!================================================================
      SUBROUTINE read_para1d(F,nn,max_points,nb_fit,nom1,exist)
      USE mod_system
      USE mod_file
      IMPLICIT NONE

       integer :: nn,nb_fit
       character (len=*) :: nom1
       integer :: max_points
       real (kind=Rkind) :: F(max_points)
       logical :: exist

       integer :: no,ios,kl,i

!      write(out_unitp,*) nom1
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)

       IF (ios .EQ. 0) THEN
         DO i=1,nb_fit
           read(no,*) nn
           IF (nn .GT. max_points) THEN
             write(out_unitp,*) ' ERROR : nb de points du fit (',nn,') >'
             write(out_unitp,*) '         a max_points (',max_points,')'
             write(out_unitp,*) '         STOP in read_para0d'
             STOP
           END IF
           DO kl=1,nn
            read(no,*) F(kl)
!           write(out_unitp,*) F(kl)
           END DO
         END DO
         CLOSE(no)
         exist = .TRUE.
       ELSE
         write(out_unitp,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


!      write(out_unitp,*) 'nom1,exist ',nom1,exist
!      write(out_unitp,*) 'nn,max_points',nn,max_points

       RETURN
       end subroutine read_para1d
!
!================================================================
!    read parameters for the fit
!    (several fits)
!================================================================
      SUBROUTINE read_para2d(F,nn,nb_fit,max_fit,max_points,nom1,exist)
      USE mod_system
      USE mod_file
      IMPLICIT NONE

       integer :: max_points,max_fit,nb_fit
       integer :: nn(max_fit)
       character (len=*) :: nom1
       real (kind=Rkind) :: F(max_points,max_fit)
       logical :: exist

       integer :: no,ios,kl,i

       write(out_unitp,*) nom1,nb_fit,max_fit,max_points
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)
       IF (ios .EQ. 0) THEN

         IF (nb_fit == 0) read(no,*) nb_fit

         IF (nb_fit > max_fit) THEN
           write(out_unitp,*) ' ERROR in read_para2d'
           write(out_unitp,*) ' The index nfit (',nb_fit,                &
                     ') is greater than max_fit',max_fit
           STOP
         END IF
         write(out_unitp,*) nom1,nb_fit
         IF (nb_fit <= 0) THEN
           write(out_unitp,*) ' ERROR in read_para2d'
           write(out_unitp,*) ' nb_fit is < 1 !!',nb_fit,' in the file :',nom1
           STOP
         END IF
         DO i=1,nb_fit
           read(no,*) nn(i)
           write(out_unitp,*) 'nom1,nb_fit,i,nn ',nom1,nb_fit,i,nn(i)
           IF (nn(i) .GT. max_points) THEN
             write(out_unitp,*) ' ERROR : nb de points du fit (',nn(i),') >'
             write(out_unitp,*) '         a max_points (',max_points,')'
             write(out_unitp,*) '         STOP in read_para2d'
             STOP
           END IF
           DO kl=1,nn(i)
            read(no,*) F(kl,i)
!           write(out_unitp,*) F(kl,i)
           END DO
         END DO
         CLOSE(no)
         exist = .TRUE.
       ELSE
         write(out_unitp,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


!      write(out_unitp,*) 'nom1,exist ',nom1,exist
!      write(out_unitp,*) 'nn,max_points',nn,max_points

       end subroutine read_para2d
!
!================================================================
!    read the parameters of the fits
!================================================================
      SUBROUTINE read_para3d(F,n,ndim,nb_fit,max_fit,max_points,        &
                              nom1,exist)
      USE mod_system
      USE mod_file
      IMPLICIT NONE

       integer :: max_points,max_fit,ndim
       integer :: n(0:ndim,max_fit),nb_fit
       character (len=*) :: nom1
       real (kind=Rkind) :: F(max_points,max_fit)
       logical :: exist

       integer :: no,ios,kl,i

       write(out_unitp,*) 'read_para3d: nom1,nb_fit,max_fit,max_points: ',&
                                      nom1,nb_fit,max_fit,max_points
       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)

       IF (ios == 0) THEN

         IF (nb_fit == 0) read(no,*) nb_fit

         IF (nb_fit > max_fit) THEN
           write(out_unitp,*) ' ERROR in read_para3d'
           write(out_unitp,*) ' The index nfit (',nb_fit,               &
                     ') is greater than max_fit',max_fit
           STOP
         END IF
         write(out_unitp,*) 'nom1,nb_fit,ndim: ',nom1,nb_fit,ndim
         IF (nb_fit <= 0) THEN
           write(out_unitp,*) ' ERROR in read_para3d'
           write(out_unitp,*) ' nb_fit is < 1 !!',nb_fit,' in the file :',nom1
           STOP
         END IF
         DO i=1,nb_fit
           read(no,*) n(0:ndim,i)
           write(out_unitp,*) 'nom1,nb_fit,i,n ',nom1,nb_fit,i,n(0:ndim,i)
           IF (n(0,i) > max_points) THEN
             write(out_unitp,*) ' ERROR : nb de points du fit (',n(0,i),') >'
             write(out_unitp,*) '         a max_points (',max_points,')'
             write(out_unitp,*) '         STOP in read_para3d'
             STOP
           END IF
           DO kl=1,n(0,i)
            read(no,*) F(kl,i)
!           write(out_unitp,*) F(kl,i)
           END DO
         END DO
         CLOSE(no)
         exist = .TRUE.
       ELSE
         write(out_unitp,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF

       end subroutine read_para3d
!================================================================
!    read the parameters of the fits
!    + the type of the functions
!================================================================
      SUBROUTINE read_para4d(F,n,ndim,nt,max_points,nom1,exist)
      USE mod_system
      USE mod_file
      IMPLICIT NONE

       integer :: max_points,ndim,nt
       integer :: n(0:ndim)
       character (len=*) :: nom1
       real (kind=Rkind) :: F(max_points)
       logical :: exist

       integer :: no,ios,kl,i

       write(out_unitp,*) 'read_para4d: nom1,max_points: ',nom1,max_points

       CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,       &
                       old=.TRUE.,err_file=ios)
       IF (ios == 0) THEN

         read(no,*) nt
         read(no,*) i ! for nb_fit (not used)

         write(out_unitp,*) 'nom1,nt,ndim: ',nom1,nt,ndim
         read(no,*) n(0:ndim)
         write(out_unitp,*) 'nom1,n ',nom1,n(0:ndim)
         IF (n(0) > max_points) THEN
             write(out_unitp,*) ' ERROR : nb de points du fit (',n(0),') >'
             write(out_unitp,*) '         a max_points (',max_points,')'
             write(out_unitp,*) '         STOP in read_para4d'
             STOP
           END IF
           DO kl=1,n(0)
            read(no,*) F(kl)
!           write(out_unitp,*) F(kl)
           END DO
         CLOSE(no)
         exist = .TRUE.
       ELSE
         write(out_unitp,*) 'The file (',nom1,') does not exist !!'
         exist = .FALSE.
       END IF


       end subroutine read_para4d

