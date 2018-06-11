%nprocshared=20
%mem=5GB
#p m062x/6-311+g(2df,2p) opt=(verytight,gdiis) freq IOP(2/16=3) 

Gaussian input prepared by ASE

0 2
O                 1.0100000000       -0.1639000000        0.0000000000
O                -0.1631000000        0.4452000000       -0.0000000000
H                -0.8469000000       -0.2813000000       -0.0000000000



