%nprocshared=20
%mem=5GB
#p m062x/cc-pVTZ opt=(ts,calcfc,noeigentest,maxcycles=900) freq scf=(maxcycle=900) IOP(7/33=1,2/16=3) 

Gaussian input prepared by ASE

0 2
O                 1.4858710000       -0.5402760000        0.3271820000
O                 1.1634040000        0.4474590000       -0.5973120000
C                -1.0838410000        0.6989150000        0.3635410000
C                -1.5654170000       -0.6042930000       -0.1981150000
H                 0.0984540000        0.8039650000       -0.2024560000
H                -0.7935350000        0.6861150000        1.4124650000
H                -1.6330640000        1.5879140000        0.0748700000
H                -2.5153000000       -0.8939930000        0.2603530000
H                -0.8438100000       -1.3947620000       -0.0025840000
H                -1.7214650000       -0.5366700000       -1.2732970000
H                 2.1100650000       -0.0777650000        0.8991350000


