%nprocshared=20
%mem=5GB
#p m062x/cc-pVTZ opt=(calcfc,maxcycles=900) freq IOP(7/33=1,2/16=3) scf=(maxcycle=900) 

Gaussian input prepared by ASE

0 2
O                -0.1524406321        0.4576922330       -0.0000000000
O                 0.9907171681       -0.1737913882        0.0000000000
H                -0.8251838069       -0.2709811306       -0.0000000000



