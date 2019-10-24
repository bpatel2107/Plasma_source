import numpy as np
import matplotlib.pyplot as plt
from plasma import Plasmas, get_params

# Major Radius (m)
Rmajor = np.array([2.5, 2.6])

# Aspect ration (R/a)
aspect_ratio = np.array([1.5, 1.6])

# Minor radius
aminor = np.outer(Rmajor, 1/aspect_ratio)

# Elongation
elong = ([2.8, 2.9])

# triangularity
tri = ([0.0, 0.1])

# Core temp in keV
T0 = ([20.0, 22.0, 25.0])

# Temperature pedestal height fraction
Tped_height = ([0.2, 0.3, 0.4])

# Separatrix temeperature
Tsep_height = [0.01, 0.02]

# Temperature peaking factor
alpha_T = ([1.0, 1.1])

# Density in core (#10^19)
n0 = ([14.0, 16.0])

# Density pedestal height fraction
nped_height = ([0.8, 0.9])

# Density pedestal height fraction
nsep_height = [0.01, 0.02]

# Density peaking factor
alpha_n = ([1.0, 1.1])

# Pedestal radial position
ped_rad = ([0.85, 0.9])

# No. of radial zones
nrhos = 100

# No. of poloidal points
ntheta = 1000

# Create Plasmas class with these values
pl = Plasmas(Rmajor, aspect_ratio, elong, tri, T0,
             Tped_height, Tsep_height,
             n0, nped_height, nsep_height, ped_rad,
             alpha_n, alpha_T, nrhos, ntheta)

# Parameter to filter
key = 'pfus'

# Point to filter around
value = 1.0

# Tolerance of parameters
tol = 0.2

#
fustype = ['Fausser', 'Wesson', 'SCENE']

# Loop through different ways of calculating fusion power
for i in range(3):
    print('!!! {} Method !!!'.format(fustype[i]))

    # Calculate fusion powers of different scenarios
    pl.calc_fusion_power(fuspow_type=i)

    # Filter around given parameter
    filter_da = pl.filter_params(key='pfus', value=1.0, tol=0.2)

    # Get list of suitable parameters
    filter_params = get_params(filter_da)

    nparams = len(filter_params)
    names = filter_da.coords.to_index().names

    # for i in range(len(filter_params)):
    #    print(filter_params[i])
    print('There are {} permutations that are are with {} of {} for {} '.format(
        nparams, tol, value, key))

    # Get first value in form to be read by datarray
    d = {name: par for name, par in zip(names, filter_params[0])}

    for key, val in zip(d.keys(), d.values()):
        print('{} = {}'.format(key, val))

    print('{} method for fusion power {:.3f} GW\n'.format(
        fustype[i], filter_da.sel(d).data))
