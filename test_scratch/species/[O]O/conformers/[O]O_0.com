%nprocshared=20
%mem=5GB
#p m062x/cc-pVTZ opt=(calcfc,maxcycles=900) freq IOP(7/33=1,2/16=3) scf=(maxcycle=900) 

Gaussian input prepared by ASE

0 2
O                -0.1608000000        0.4351000000       -0.0000000000
O                 1.0186000000       -0.1617000000        0.0000000000
H                -0.8579000000       -0.2734000000       -0.0000000000



