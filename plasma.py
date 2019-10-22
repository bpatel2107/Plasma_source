from numpy import pi, linspace, diff, newaxis, meshgrid, shape
import numpy as np
import xarray as xr


class Plasmas:
    """
    Generates plasma states with the given parameters

    """

    def __init__(self, Rmajor=2.5, aspect_ratio=1.67, elong=2.8, tri=0.0,
                 T0=25.0, Tped_height=0.2, n0=16.0, nped_height=0.8, ped_rad=0.9,
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

        self.ped_rad = ped_rad

        # Dataset containing all the grid
        self.ds = self._create_grid(Rmajor=Rmajor, aspect_ratio=aspect_ratio,
                                    elong=elong, tri=tri, T0=T0,
                                    Tped_height=Tped_height, n0=n0,
                                    nped_height=nped_height, ped_rad=ped_rad,
                                    rhos=self.rhos, theta=self.theta)

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

        darea = area.diff('rhos')

        vol = 2 * pi * ds['Rmajor'] * area
        dvol = 2 * pi * ds['Rmajor'] * darea

        self.area = area
        self.darea = darea
        self.vol = vol
        self.dvol = dvol

    def gen_kinetic_profiles(self):

        ds = self.ds

        T = ds['T0'] - ds['T0'] * \
            (1 - ds['Tped_height'])/ds['ped_rad'] * ds['rhos']

        n = ds['n0'] - ds['n0'] * \
            (1 - ds['nped_height'])/ds['ped_rad'] * ds['rhos']

        self.ds['T'] = T
        self.ds['n'] = n

    def calc_fusion_power(self, fuspow_type=0):

        ds = self.ds
        self.fuspow_type = fuspow_type

        self.gen_kinetic_profiles()

        if fuspow_type == 0:

            pden = 1.1e-24 * ds['T']**2 * ds['n']**2 / \
                4 * 17.5e6 * 1.6e-19 * self.dvol * 1e-9 * 1e38

        else:

            arg = -0.476 * np.abs(np.log(1.45e-5*ds['T']*1e3))**2.25

            pden = 1.27e4 * ds['n']**2 * np.exp(arg) * self.vol * 1e-9

        self.ds['pfus'] = pden.sum('rhos')

    def filter_params(self, key=None, value=None, tol=None):

        if key is None:
            raise ValueError('Key needed to filter_params')

        if value is None:
            raise ValueError('Filter valued needed')

        if tol is None:
            raise ValueError('Tolerance to value must be set')

        ds = self.ds[key]

        filter_ds = ds.where(abs(ds.data - value) < tol)

        return filter_ds


def get_params(da=None):

    if da is None:
        raise ValueError('DataArray needed as input')

    keys = da.to_dict()['coords'].keys()

    filter_list = []
    for pfus in da.data:

        print(pfus)
        # if not np.isnan(da.data):
        # print(da)
        # filter_list.append(da.coords)

        # return filter_list
