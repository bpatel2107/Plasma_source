import numpy as np
import matplotlib.pyplot as plt
from plasma import Plasmas, get_params

pi = np.pi

Rmajor = np.array([2.5, 2.6])
aspect_ratio = ([1.61, 1.7])

aminor = Rmajor/aspect_ratio

elong = ([2.8, 2.9])

tri = ([0.0, 0.1])

# Core temp in keV
T0 = ([25.0, 22.0])

Tped_height = ([0.2, 0.3])

# Density in core (#10^19)
n0 = ([16.0, 18.0])

nped_height = ([0.8, 0.9])

ped_rad = ([0.85, 0.9])

nrhos = 100

ntheta = 1000

pl = Plasmas(Rmajor, aspect_ratio, elong, tri, T0,
             Tped_height, n0, nped_height, ped_rad, nrhos, ntheta)


pl.calc_fusion_power(fuspow_type=1)
#print('SCENE calc = {} GW'.format(pl.ds['pfus'].data))

pl.calc_fusion_power()
#print('Wesson calc = {} GW'.format(pl.ds['pfus'].data))

filter_da = pl.filter_params(key='pfus', value=1.0, tol=0.2)

filter_list = get_params(filter_da)
