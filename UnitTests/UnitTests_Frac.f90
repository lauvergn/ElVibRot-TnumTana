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
PROGRAM UnitTest_Frac
  USE mod_NumParameters
  USE mod_Frac
  IMPLICIT NONE


   TYPE(Frac_t)              :: frac1,frac2,frac3
   TYPE(Frac_t), allocatable :: tab_frac1(:)
   TYPE(Frac_t), allocatable :: tab_frac2(:)
   integer                   :: i


   frac3 = '1/3'
   write(out_unitp,*) 'frac3="1/3"   1/3 ?: ',frac_TO_string(frac3)
   frac3 = '-1/3'
   write(out_unitp,*) 'frac3="-1/3" -1/3 ?: ',frac_TO_string(frac3)
   frac3 = '-1/3'
   write(out_unitp,*) 'frac3="1/-3" -1/3 ?: ',frac_TO_string(frac3)
   frac3 = '-1/-3'
   write(out_unitp,*) 'frac3="-1/-3" 1/3 ?: ',frac_TO_string(frac3)

   frac1 = Frac_t(1,2)
   frac2 = 2
   write(out_unitp,*) 'frac1=        1/2 ?: ',frac_TO_string(frac1)
   write(out_unitp,*) 'frac2=        2   ?: ',frac_TO_string(frac2)
   write(out_unitp,'(1x,a,f4.1)') 'frac1=        0.5 ?: ',frac_TO_real(frac1)
   write(out_unitp,*)
   write(out_unitp,*) '+frac1=       1/2 ?: ',frac_TO_string(+frac1)
   write(out_unitp,*) 'frac1+frac2=  5/2 ?: ',frac_TO_string(frac1+frac2)
   write(out_unitp,*) 'frac1+2=      5/2 ?: ',frac_TO_string(frac1+2)
   write(out_unitp,*) '2+frac1=      5/2 ?: ',frac_TO_string(2+frac1)
   write(out_unitp,*)
   write(out_unitp,*) '-frac1=      -1/2 ?: ',frac_TO_string(-frac1)
   write(out_unitp,*) 'frac1-frac2= -3/2 ?: ',frac_TO_string(frac1-frac2)
   write(out_unitp,*) 'frac1-2=     -3/2 ?: ',frac_TO_string(frac1-2)
   write(out_unitp,*) '2-frac1=      3/2 ?: ',frac_TO_string(2-frac1)
   write(out_unitp,*)
   write(out_unitp,*) 'frac1*frac2=  1   ?: ',frac_TO_string(frac1*frac2)
   write(out_unitp,*) 'frac1*2=      1   ?: ',frac_TO_string(frac1*2)
   write(out_unitp,*) '2*frac1=      1   ?: ',frac_TO_string(2*frac1)
   write(out_unitp,*)
   write(out_unitp,*) 'frac1/frac2=  1/4 ?: ',frac_TO_string(frac1/frac2)
   write(out_unitp,*) 'frac1/2=      1/4 ?: ',frac_TO_string(frac1/2)
   write(out_unitp,*) '2/frac1=      4   ?: ',frac_TO_string(2/frac1)
   write(out_unitp,*)

   write(out_unitp,*) 'frac1==frac2    F ?: ',(frac1 == frac2)
   write(out_unitp,*) '    1==frac2    F ?: ',(1 == frac2)
   write(out_unitp,*) 'frac1==1        F ?: ',(frac1 == 1)
   write(out_unitp,*) 'frac1==frac1    T ?: ',(frac1 == frac1)
   write(out_unitp,*) 'frac2==2        T ?: ',(frac2 == 2)
   write(out_unitp,*) '    2==frac2    T ?: ',(2 == frac2)
   write(out_unitp,*)

   write(out_unitp,*) 'frac1/=frac2    T ?: ',(frac1 /= frac2)
   write(out_unitp,*) '    1/=frac2    T ?: ',(1 /= frac2)
   write(out_unitp,*) 'frac1/=1        T ?: ',(frac1 /= 1)
   write(out_unitp,*) 'frac1/=frac1    F ?: ',(frac1 /= frac1)
   write(out_unitp,*) 'frac2/=2        F ?: ',(frac2 /= 2)
   write(out_unitp,*) '    2/=frac2    F ?: ',(2 /= frac2)
   write(out_unitp,*)

   write(out_unitp,*) 'frac1<frac2     T ?: ',(frac1 < frac2)
   write(out_unitp,*) '    1<frac2     T ?: ',(1 < frac2)
   write(out_unitp,*) 'frac1<1         T ?: ',(frac1 < 1)
   write(out_unitp,*) 'frac1<frac1     F ?: ',(frac1 < frac1)
   write(out_unitp,*) 'frac2<2         F ?: ',(frac2 < 2)
   write(out_unitp,*) '    2<frac2     F ?: ',(2 < frac2)
   write(out_unitp,*)

   write(out_unitp,*) 'frac1>frac2     F ?: ',(frac1 > frac2)
   write(out_unitp,*) '    1>frac2     F ?: ',(1 > frac2)
   write(out_unitp,*) 'frac1>1         F ?: ',(frac1 > 1)
   write(out_unitp,*) 'frac1>frac1     F ?: ',(frac1 > frac1)
   write(out_unitp,*) 'frac2>2         F ?: ',(frac2 > 2)
   write(out_unitp,*) '    2>frac2     F ?: ',(2 > frac2)
   write(out_unitp,*)

   write(out_unitp,*) 'frac1>=frac2    F ?: ',(frac1 >= frac2)
   write(out_unitp,*) '    1>=frac2    F ?: ',(1 >= frac2)
   write(out_unitp,*) 'frac1>=1        F ?: ',(frac1 >= 1)
   write(out_unitp,*) 'frac1>=frac1    T ?: ',(frac1 >= frac1)
   write(out_unitp,*) 'frac2>=2        T ?: ',(frac2 >= 2)
   write(out_unitp,*) '    2>=frac2    T ?: ',(2 >= frac2)
   write(out_unitp,*)

   write(out_unitp,*) 'frac1<=frac2    T ?: ',(frac1 <= frac2)
   write(out_unitp,*) '    1<=frac2    T ?: ',(1 <= frac2)
   write(out_unitp,*) 'frac1<=1        T ?: ',(frac1 <= 1)
   write(out_unitp,*) 'frac1<=frac1    T ?: ',(frac1 <= frac1)
   write(out_unitp,*) 'frac2<=2        T ?: ',(frac2 <= 2)
   write(out_unitp,*) '    2<=frac2    T ?: ',(2 <= frac2)
   write(out_unitp,*)

   write(out_unitp,*)
   write(out_unitp,*)
   frac1 = Frac_t(1,6)
   frac2 = Frac_t(18,12)
   write(out_unitp,*) 'frac1         1/6 ?: ',frac_TO_string(frac1)
   write(out_unitp,*) 'frac2=18/12:  3/2 ?: ',frac_TO_string(frac2)
   write(out_unitp,*) 'frac1+frac2:  5/3 ?: ',frac_TO_string(frac1+frac2)

   STOP

   write(out_unitp,*) 'table of Frac:'
   tab_frac1 = [ Frac_t(1,2), Frac_t(1,3), Frac_t(1,4) ]
   write(out_unitp,*) 'tab_frac1            ',                          &
               (  (frac_TO_string(tab_frac1(i)) // ' '),i=1,size(tab_frac1))


   tab_frac2 = tab_frac1
   write(out_unitp,*) 'tab_frac2=tab_frac1  ',tab_frac2
   write(out_unitp,*)

   tab_frac2 = 5*tab_frac1
   write(out_unitp,*) '5*tab_frac1          ',tab_frac2
   tab_frac2 = tab_frac1*5
   write(out_unitp,*) 'tab_frac1*5          ',tab_frac2
   tab_frac2 = frac1*tab_frac1
   write(out_unitp,*) 'frac1*tab_frac1      ',tab_frac2
   tab_frac2 = tab_frac1*frac1
   write(out_unitp,*) 'tab_frac1*frac1      ',tab_frac2
   tab_frac2 = tab_frac1*tab_frac1
   write(out_unitp,*) 'tab_frac1*tab_frac1  ',tab_frac2
   write(out_unitp,*)

   tab_frac2 = 5/tab_frac1
   write(out_unitp,*) '5/tab_frac1          ',tab_frac2
   tab_frac2 = tab_frac1/5
   write(out_unitp,*) 'tab_frac1/5          ',tab_frac2
   tab_frac2 = tab_frac1/tab_frac1
   write(out_unitp,*) 'tab_frac1/tab_frac1  ',tab_frac2
   tab_frac2 = tab_frac1/0
   write(out_unitp,*) 'tab_frac1/0          ',tab_frac2
   write(out_unitp,*)

   tab_frac2 = tab_frac1 + tab_frac1
   write(out_unitp,*) 'tab_frac1+tab_frac1  ',tab_frac2
   tab_frac2 = tab_frac1 + 5
   write(out_unitp,*) 'tab_frac1+5          ',tab_frac2
   tab_frac2 = 5+tab_frac1
   write(out_unitp,*) '5+tab_frac1          ',tab_frac2
   write(out_unitp,*)

   tab_frac2 = tab_frac1 - tab_frac1
   write(out_unitp,*) 'tab_frac1-tab_frac1  ',tab_frac2
   tab_frac2 = tab_frac1 - 5
   write(out_unitp,*) 'tab_frac1-5          ',tab_frac2
   tab_frac2 = 5-tab_frac1
   write(out_unitp,*) '5-tab_frac1          ',tab_frac2
   write(out_unitp,*)

END PROGRAM UnitTest_Frac
