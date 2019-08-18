 %chk=xx
 #n freq=noraman hf sto-3g opt=(z-matrix,verytight)
 # scf=tight FormCheck=All
 # integral=(grid=ultrafine) 
 # iop(1/18=40,2/9=111,2/11=1,7/33=1,1/33=1)
 # iop(3/27=30,11/27=30)
  
  xx 
  
           0             1
 O
 H1 O R1
 H2 O R2 H1 A

R1 0.9895
R2 0.9895
A  100.015

