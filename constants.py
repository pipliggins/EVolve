"""
Stores universally useful constants.

Stores constants for molecular mass, gas equation constants
and partition coefficients.

"""

# Universal gas constant
R = 8.3144598

# Avogadro's constant  mol-1
N_A = 6.0221409e+23

# Boltzmann's constant in erg K-1
K = 1.3806504e-16

m = {}  # atomic and molecular masses  g/mol
m['h'] = 1.00794
m['o'] = 15.9994
m['c'] = 12.011
m['s'] = 32.066
m['n'] = 14.0067
m['h2'] = 2.01588
m['o2'] = 31.9988
m['h2o'] = 18.01528
m['co2'] = 44.0098
m['co'] = 28.0104
m['ch4'] = 16.04276
m['so2'] = 64.0648
m['h2s'] = 34.08188
m['hcn'] = 27.02564
m['he'] = 4.002602
m['s2'] = 64.132
m['ocs'] = 60.0764
m['c2h2'] = 26.03788
m['nh3'] = 17.03052
m['n2'] = 28.0134
m['si'] = 28.0855
m['ti'] = 47.867
m['al'] = 26.981539
m['mn'] = 54.938044
m['mg'] = 24.305
m['ca'] = 40.078
m['na'] = 22.990
m['k'] = 39.0983
m['p'] = 30.97376
m['li'] = 6.941
m['fe'] = 55.845
m['sio2'] = 60.0843
m['tio2'] = 79.8658
m['al2o3'] = 101.961278
m['fe2o3'] = 159.6882
m['feo'] = 71.8444
m['feot'] = 71.8444
m['mno'] = 70.937445
m['mgo'] = 40.3044
m['cao'] = 56.0774
m['na2o'] = 61.978938
m['k2o'] = 94.196
m['p2o5'] = 141.944524
m['li2o'] = 29.8814

D = {}  # Partition coefficients for volatile partitioning during melting
D['h2o'] = 0.01         # Aubaud 2004
D['n_nno'] = 0.00062    # Li 2013 at NiNiO buffer, for 60% Oliv 40% pyrox
D['n_fefeo'] = 0.0073   # Li 2013 at FeFeO buffer, for 60% Oliv 40% pyrox
D['s'] = 0.01           # Callegaro et al 2020 for reduced systems
D['co2'] = 5.5e-4       # Rosenthal 2015
D['h2'] = 0.01          # Aubaud 2004 - we assume this is comparable to water.