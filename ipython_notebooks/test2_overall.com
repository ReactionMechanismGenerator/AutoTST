%nprocshared=20
%mem=5GB
#p hf/6-31g* opt=(ts,calcfc,noeigentest) freq 

Gaussian input prepared by ASE

0 2
O                -3.1668700073        0.1708626607       -0.7505687466
O                -3.8863727486        1.3133346049       -0.5229368425
C                 0.2295197688        0.3652120341       -1.0900205158
C                 0.9175774348       -0.1336025439       -0.0290478628
C                -0.9058984558        1.2760405992       -1.0524531635
C                 2.0635661868       -1.0908888612       -0.2375742588
C                 0.6218654899        0.2007282106        1.4106268666
H                -4.2820245070        1.4992158577       -1.3664928788
H                 0.5166914990        0.0197901499       -2.0697586483
H                -2.0139370098        0.5940594053       -0.9395815095
H                -1.0881120871        1.8131459093       -1.9724383102
H                -1.0134932156        1.9096170101       -0.1866023692
H                 1.8699217735       -2.0394374996        0.2583484008
H                 2.2397150449       -1.2903241553       -1.2880592009
H                 2.9821584048       -0.6940194740        0.1891328087
H                -0.2336688304        0.8501772285        1.5318843679
H                 0.4257824959       -0.7086887584        1.9732513903
H                 1.4806195537        0.6804865035        1.8751691685


