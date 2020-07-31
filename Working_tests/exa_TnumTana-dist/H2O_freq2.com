%chk=h2o_freq
#hf sto-3g  opt=(tight,z-matrix) freq nosym


0 1
O
X  O 1.
H1 O R1 X  A
H2 O R2 X  A  H1 180.

R1 = 1.
R2 = 1.
A = 50.

