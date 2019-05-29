%nprocshared=20
%mem=5GB
#p m062x/cc-pVTZ opt=(calcfc,maxcycles=900) freq IOP(7/33=1,2/16=3) scf=(maxcycle=900) 

Gaussian input prepared by ASE

0 1
O                 0.5845000000       -0.3481000000        0.2757000000
O                -0.6151000000       -0.3874000000       -0.2543000000
H                 1.1402000000        0.3472000000       -0.1517000000
H                -1.1097000000        0.3882000000        0.1303000000



