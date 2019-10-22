from numpy import pi, linspace, diff, newaxis, meshgrid, shape
import numpy as np
import xarray as xr


class Plasmas:
    """
    Generates plasma states with the given parameters

    """

    def __init__(self, Rmajor=2.5, aspect_ratio=1.67, elong=2.8, tri=0.0,
                 T0=25.0, Tped_height=0.2, n0=16.0, nped_height=0.8,
                 nrhos=20, ntheta=1000):

        # Radial zones
        self.nrhos = nrhos
        self.rhos = linspace(0, 1, nrhos)

        # Poloidal range
        self.ntheta = ntheta
        self.theta = linspace(0, 2*pi, ntheta)
        self.dtheta = self.theta[1] - self.theta[0]

        # Shaping parameters
        self.Rmajor = Rmajor
        self.aspect_ratio = aspect_ratio
        self.elong = elong
        self.tri = tri

        # Kinetic profiles
        self.T0 = T0
        self.Tped_height = Tped_height
        self.n0 = n0
        self.nped_height = nped_height

        # Dataset containing all the grid
        self.ds = self._create_grid(Rmajor=Rmajor, aspect_ratio=aspect_ratio, elong=elong, tri=tri,
                                    T0=T0, Tped_height=Tped_height, n0=n0, nped_height=nped_height, rhos=self.rhos, theta=self.theta)

        # Calculate area/volume of each plasma
        self.calc_plasma_vol()

    def _create_grid(self, **kwargs):

        ds = xr.Dataset(coords=kwargs)
        return ds

    def calc_plasma_vol(self, nrhos=20, ntheta=1000):

        ds = self.ds
        aminor = ds['Rmajor']/ds['aspect_ratio']

        Rbnd = ds['Rmajor'] + aminor * ds['rhos'] * \
            np.cos(ds['theta'] + np.arcsin(ds['tri']) * np.sin(ds['theta']))

        Zbnd = ds['elong'] * aminor * ds['rhos'] * np.sin(ds['theta'])

        rsqu = (Rbnd - ds['Rmajor'])**2 + Zbnd**2
        area = rsqu.sum(dim='theta') * self.dtheta/2

        print(area)
