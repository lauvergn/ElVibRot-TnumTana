/*==========================================================================
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
!===========================================================================*/
#include <stdio.h>
#include <math.h>


void Qact_TO_Qcart_TnumTanaDriver_FOR_c(double *, int *, double *, int *);
void Qcart_TO_Qact_TnumTanaDriver_FOR_c(double *, int *, double *, int *);
void Init_TnumTana_FOR_Driver_FOR_c(int *, int *, int *);

int main(void)
{
  int        init,nb_act,nb_cart;
  double     Qact[6];
  double     Qcart[12];

  int        j;


  Qact[0]=2.2791110577179996;
  Qact[1]=2.0824950811500056;
  Qact[2]=2.1237280405061139;
  Qact[3]=2.0824950811500056;
  Qact[4]=2.1237280405061139;
  Qact[5]=3.1415926535897931;

  init = 0;
  Init_TnumTana_FOR_Driver_FOR_c(&nb_act,&nb_cart,&init);

  printf("nb_act = %i\n",nb_act);
  printf("nb_cart = %i\n",nb_cart);

  

  for (j=0;j<nb_act;j++)
    printf(" Qact[:] = %30.22le \n",Qact[j]);

   Qact_TO_Qcart_TnumTanaDriver_FOR_c(Qact,&nb_act,Qcart,&nb_cart);

  for (j=0;j<nb_cart;j=j+3)
    printf(" X %30.22le %30.22le %30.22le \n",Qcart[j+0],Qcart[j+2],Qcart[j+2]);

   Qcart_TO_Qact_TnumTanaDriver_FOR_c(Qact,&nb_act,Qcart,&nb_cart);

  for (j=0;j<nb_act;j++)
    printf(" Qact[:] = %30.22le \n",Qact[j]);
}
