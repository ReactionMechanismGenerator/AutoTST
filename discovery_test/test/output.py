# Coordinates for [CH]=CC=C in Input Orientation (angstroms):
#   C   -0.2686   -0.5255   -0.2591
#   C    0.6900    0.3360    0.4316
#   C   -1.4757   -0.1170   -0.6365
#   C    1.8858   -0.0397    0.8150
#   H    0.0536   -1.5412   -0.4574
#   H    0.3638    1.3591    0.6305
#   H   -1.8104    0.8950   -0.4448
#   H   -2.1623   -0.7780   -1.1457
#   H    2.7239    0.4113    1.3192
conformer(
    label = '[CH]=CC=C',
    E0 = (367.461, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(53.0391, 'amu')),
        NonlinearRotor(
            inertia = ([10.2779, 109.027, 119.305], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([167.281, 297.023, 505.835, 577.221, 753.934, 797.527, 883.089, 938.592, 959.471, 1025.3, 1163.87, 1228.6, 1294.86, 1423.27, 1621.77, 1681.46, 3021.42, 3102.28, 3119.84, 3193.84, 3214.56], 'cm^-1'),
        ),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
)

# Coordinates for [O]O in Input Orientation (angstroms):
#   O   -0.1632    0.4357    0.0000
#   O    0.9953   -0.1661    0.0000
#   H   -0.8320   -0.2696    0.0000
conformer(
    label = '[O]O',
    E0 = (18.538, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(32.9977, 'amu')),
        NonlinearRotor(
            inertia = ([0.800947, 14.5112, 15.3121], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(frequencies=([1245.9, 1447.19, 3645.15], 'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
)

# Coordinates for OO in Input Orientation (angstroms):
#   O   -0.5826    0.4330    0.2094
#   O    0.5758   -0.3828    0.3045
#   H   -1.1901   -0.1479   -0.2608
#   H    1.1969    0.0976   -0.2531
conformer(
    label = 'OO',
    E0 = (-122.037, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(34.0055, 'amu')),
        NonlinearRotor(
            inertia = ([1.63801, 18.4022, 19.029], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([378.1, 1029.14, 1344.99, 1458.5, 3786.54, 3787.27], 'cm^-1'),
        ),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

# Coordinates for [CH]=C[C]=C in Input Orientation (angstroms):
#   C    0.7829   -0.3844    0.0320
#   C   -1.6863    0.3638   -0.4935
#   C   -0.4701   -0.0005   -0.2385
#   C    1.9284    0.1962   -0.4800
#   H    0.9001   -1.2264    0.7157
#   H   -2.1929    1.1099    0.1105
#   H   -2.2416   -0.0601   -1.3241
#   H    2.9796    0.0015   -0.3544
conformer(
    label = '[CH]=C[C]=C',
    E0 = (551.027, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(52.0313, 'amu')),
        NonlinearRotor(
            inertia = ([9.75325, 110.672, 116.97], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([212.287, 239.629, 471.402, 522.168, 558.199, 823.409, 844.226, 891.327, 964.761, 1026.62, 1149.04, 1317.57, 1423.17, 1906.4, 3024.56, 3044.53, 3114.43, 3228.23], 'cm^-1'),
        ),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
)

# Coordinates for TS in Input Orientation (angstroms):
#   C   -0.7032   -0.5905    0.1258
#   C   -0.9294    0.7955    0.4329
#   C   -1.4641   -1.6387   -0.0933
#   C   -1.4831    1.6481   -0.4119
#   H    0.5977   -0.7885    0.0366
#   H   -0.5549    1.1596    1.3909
#   H   -2.5422   -1.5758    0.0210
#   H   -1.0468   -2.5959   -0.3769
#   H   -1.7055    2.7031   -0.4053
#   O    2.0537    0.5984    0.0717
#   O    1.7817   -0.7404   -0.0916
#   H    2.0468    0.9464   -0.8288
conformer(
    label = 'TS',
    E0 = (474.097, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(86.0368, 'amu')),
        NonlinearRotor(
            inertia = ([118.281, 209.82, 312.632], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([63.6799, 73.8973, 121.156, 187.791, 236.769, 340.268, 374.544, 461.43, 538.576, 653.924, 713.542, 829.586, 884.529, 919.366, 949.909, 984.144, 1076.18, 1090.07, 1251.29, 1409.56, 1410.14, 1467.63, 1557.42, 1713.59, 3053.35, 3063.86, 3150.56, 3207.86, 3734.94], 'cm^-1'),
        ),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    frequency = (-1861.27, 'cm^-1'),
)

#   ======= =========== =========== =========== ===============
#   Temp.   k (TST)     Tunneling   k (TST+T)   Units
#   ======= =========== =========== =========== ===============
#     300 K   5.244e-06     80.6089   4.227e-04 cm^3/(mol*s)
#     400 K   5.876e-02     9.09393   5.344e-01 cm^3/(mol*s)
#     500 K   1.892e+01     3.80743   7.206e+01 cm^3/(mol*s)
#     600 K   1.008e+03     2.47844   2.499e+03 cm^3/(mol*s)
#     800 K   1.823e+05     1.66109   3.028e+05 cm^3/(mol*s)
#    1000 K   4.998e+06      1.3895   6.944e+06 cm^3/(mol*s)
#    1500 K   6.147e+08     1.16614   7.169e+08 cm^3/(mol*s)
#    2000 K   9.179e+09     1.09546   1.006e+10 cm^3/(mol*s)
#   ======= =========== =========== =========== ===============


#   ======= ============ =========== ============ ============= =========
#   Temp.    Kc (eq)        Units     krev (TST)   krev (TST+T)   Units
#   ======= ============ =========== ============ ============= =========
#     300 K   5.450e-08              9.623e+01     7.757e+03      cm^3/(mol*s)
#     400 K   4.842e-06              1.213e+04     1.104e+05      cm^3/(mol*s)
#     500 K   7.546e-05              2.508e+05     9.549e+05      cm^3/(mol*s)
#     600 K   4.855e-04              2.077e+06     5.147e+06      cm^3/(mol*s)
#     800 K   5.184e-03              3.517e+07     5.841e+07      cm^3/(mol*s)
#    1000 K   2.198e-02              2.274e+08     3.159e+08      cm^3/(mol*s)
#    1500 K   1.543e-01              3.985e+09     4.647e+09      cm^3/(mol*s)
#    2000 K   4.107e-01              2.235e+10     2.448e+10      cm^3/(mol*s)
#   ======= ============ =========== ============ ============= =========


# krev (TST) = Arrhenius(A=(5.79565,'cm^3/(mol*s)'), n=3.20622, Ea=(38.5481,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Fitted to 8 data points; dA = *|/ 1.54098, dn = +|- 0.0564657, dEa = +|- 0.31666 kJ/mol""") 
# krev (TST+T) = Arrhenius(A=(7.76192e-08,'cm^3/(mol*s)'), n=5.43369, Ea=(14.5306,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Fitted to 8 data points; dA = *|/ 16.6696, dn = +|- 0.367402, dEa = +|- 2.06039 kJ/mol""") 

kinetics(
    label = '[CH]=CC=C+[O]O_OO+[CH]=C[C]=C',
    kinetics = Arrhenius(
        A = (5.20225e-08, 'cm^3/(mol*s)'),
        n = 5.72802,
        Ea = (59.9526, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (303.03, 'K'),
        Tmax = (2500, 'K'),
        comment = 'Fitted to 59 data points; dA = *|/ 2.7068, dn = +|- 0.130686, dEa = +|- 0.718902 kJ/mol',
    ),
)

