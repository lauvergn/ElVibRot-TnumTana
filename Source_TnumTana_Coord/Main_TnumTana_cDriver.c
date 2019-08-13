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
