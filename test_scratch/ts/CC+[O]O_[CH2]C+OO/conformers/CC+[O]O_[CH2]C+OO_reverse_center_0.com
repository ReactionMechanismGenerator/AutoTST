%nprocshared=20
%mem=5GB
#p m062x/cc-pVTZ Opt=(ts,calcfc,noeigentest,ModRedun,maxcycles=900) scf=(maxcycle=900) 

Gaussian input prepared by ASE

0 2
O                 1.9644030000        0.5339730000       -0.0113000000
O                 1.2945950000       -0.5670360000        0.3634990000
C                -2.1876300000        0.2676530000        0.2560940000
C                -0.9733280000       -0.2754610000       -0.4751240000
H                -3.0617620000        0.2990440000       -0.3964650000
H                -2.0118910000        1.2783330000        0.6214980000
H                -2.4415100000       -0.3563950000        1.1120630000
H                -0.6757720000        0.3567000000       -1.3124230000
H                -1.1225780000       -1.2919300000       -0.8359850000
H                 2.2494620000        0.3468190000       -0.9178050000
H                -0.0421870000       -0.3212180000        0.2257050000

1 3 F
1 5 F
1 6 F
1 7 F
1 8 F
1 9 F
1 10 F
3 5 F
3 6 F
3 7 F
3 8 F
3 9 F
3 10 F
5 6 F
5 7 F
5 8 F
5 9 F
5 10 F
6 7 F
6 8 F
6 9 F
6 10 F
7 8 F
7 9 F
7 10 F
8 9 F
8 10 F
9 10 F
