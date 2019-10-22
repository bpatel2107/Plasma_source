import numpy as np
import matplotlib.pyplot as plt
from plasma import Plasmas

pi = np.pi

Rmajor = np.array([2.5, 2.6])
aspect_ratio = 1.61

aminor = Rmajor/aspect_ratio

elong = 2.85

tri = 0.0

# Core temp in keV
T0 = 25

Tped_height = 0.2

Tped = T0 * Tped_height

# Density in core (#10^19)
n0 = 16.0

nped_height = 0.8

ped_rad = 0.9

nped = n0 * nped_height

p0 = n0*T0 * 1.6

nrhos = 20
rho = np.linspace(0, 1, nrhos)
rho_grid = np.linspace(0, 1, nrhos-1)
rho = rho[:, np.newaxis]

ntheta = 1000

pl = Plasmas(Rmajor, aspect_ratio, elong, tri, T0,
             Tped_height, n0, nped_height, nrhos, ntheta)


T = T0 - (T0-Tped)/ped_rad * rho_grid

n = n0 - (n0-nped)/ped_rad * rho_grid

theta = np.linspace(0, 2*pi, 1000)
dtheta = theta[1] - theta[0]

# Generate plasma boundary
Rbnd = Rmajor + aminor * \
    np.outer(rho, np.cos(theta + np.arcsin(tri)*np.sin(theta)))
Zbnd = elong * aminor * np.outer(rho,  np.sin(theta))

for i in range(nrhos):
    plt.plot(Rbnd[i, :], Zbnd[i, :])
plt.show()
area = np.sum((Rbnd-Rmajor)**2 + Zbnd**2, axis=1) * dtheta / 2

darea = np.diff(area)

vol = 2 * pi * Rmajor * area
dvol = 2 * pi * Rmajor * darea

# arg = -0.476 * np.abs(np.log(1.45e-5*Tavg*1e3))**2.25

# pden = 1.27e4 * navg**2 * np.exp(arg)

# pfus = pden * vol

# print('SCENE calc says {} GW'.format(pfus*1e-9))

pfus_wes = 1.1e-24 * T**2 * \
    (n*10**19)**2/4 * 17.5e6 * 1.6e-19 * dvol * 1e-9

pfus = np.sum(pfus_wes)

print('Wesson calc says {} GW'.format(pfus))
