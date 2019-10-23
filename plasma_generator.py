import numpy as np
import matplotlib.pyplot as plt
from plasma import Plasmas, get_params

# Major Radius (m)
Rmajor = np.array([2.5, 2.6, 2.7, 2.8])

# Aspect ration (R/a)
aspect_ratio = np.array([1.5, 1.6, 1.7])

# Minor radius
aminor = np.outer(Rmajor, 1/aspect_ratio)

# Elongation
elong = ([2.7, 2.8, 2.9])

# triangularity
tri = ([0.0, 0.1, 0.2, 0.3])

# Core temp in keV
T0 = ([20.0, 22.0, 25.0])

# Temperature pedestal height fraction
Tped_height = ([0.2, 0.3, 0.4])

# Density in core (#10^19)
n0 = ([16.0, 18.0, 20.0])

# Density pedestal height fraction
nped_height = ([0.7, 0.8, 0.9])

# Pedestal radial position
ped_rad = ([0.85, 0.9])

# No. of radial zones
nrhos = 100

# No. of poloidal points
ntheta = 1000

# Create Plasmas class with these values
pl = Plasmas(Rmajor, aspect_ratio, elong, tri, T0,
             Tped_height, n0, nped_height, ped_rad, nrhos, ntheta)

# Calculate fusion powers of different scenarios
# pl.calc_fusion_power(fuspow_type=1)
pl.calc_fusion_power()

# Parameter to filter
key = 'pfus'

# Point to filter around
value = 1.0

# Tolerance of parameters
tol = 0.2

# Filter around given parameter
filter_da = pl.filter_params(key='pfus', value=1.0, tol=0.2)

filter_params = get_params(filter_da)

nparams = len(filter_params)
names = filter_da.coords.to_index().names


for i in range(len(filter_params)):
    print(filter_params[i])
print('There are {} permutations that are are with {} of {} for {} '.format(
    nparams, tol, value, key))
print(names)
