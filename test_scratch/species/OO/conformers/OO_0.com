%nprocshared=20
%mem=5GB
#p m062x/cc-pVTZ opt=(calcfc,maxcycles=900) freq IOP(7/33=1,2/16=3) scf=(maxcycle=900) 

Gaussian input prepared by ASE

0 1
O                -0.4815752002        0.4486028265       -0.0161843387
O                 0.5621927091       -0.3364087883       -0.0190349565
H                -1.2494802337       -0.1745376583        0.0578327924
H                 0.8409854112       -0.3660868683        0.9339237839



