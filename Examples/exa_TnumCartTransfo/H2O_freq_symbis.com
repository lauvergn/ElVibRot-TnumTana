 %chk=H2O_freq_symbis
 #n freq=noraman hf sto-3g opt=(z-matrix,verytight) nosym
 # scf=tight FormCheck=All
 # integral=(grid=ultrafine) 
 # iop(1/18=40,2/9=111,2/11=1,7/33=1,1/33=1)
 # iop(3/27=30,11/27=30)
  
  xx 
  
           0             1
 O
 X  O 1.
 H1 O R1 X  A
 H2 O R2 X  A H1 180.

R1 0.9895
R2 0.9895
A  50.

