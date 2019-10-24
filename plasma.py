from numpy import pi, linspace, diff, newaxis, meshgrid, shape
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


class Plasmas:
    """
    Generates plasma states with the given parameters

    """

    def __init__(self, Rmajor=2.5, aspect_ratio=1.67, elong=2.8, tri=0.0,
                 T0=25.0, Tped_height=0.2, Tsep_height=0.01, n0=16.0,
                 nped_height=0.8, nsep_height=0.01, ped_rad=0.9,
                 alpha_T=1.0, alpha_n=1.0, nrhos=20, ntheta=1000):

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
        self.Tsep_height = Tsep_height
        self.alpha_T = alpha_T

        self.n0 = n0
        self.nped_height = nped_height
        self.nsep_height = nsep_height
        self.alpha_n = alpha_n

        self.ped_rad = ped_rad

        # Dataset containing all the grid
        self.ds = self._create_grid(Rmajor=Rmajor, aspect_ratio=aspect_ratio,
                                    elong=elong, tri=tri, T0=T0,
                                    Tped_height=Tped_height,
                                    Tsep_height=Tsep_height,
                                    n0=n0, nped_height=nped_height,
                                    nsep_height=nsep_height, ped_rad=ped_rad,
                                    alpha_T=alpha_T, alpha_n=alpha_n,
                                    rhos=self.rhos, theta=self.theta)

        # Calculate area/volume of each plasma
        self.calc_plasma_vol(self.nrhos, self.ntheta)

    def _create_grid(self, **kwargs):

        ds = xr.Dataset(coords=kwargs)
        return ds

    def calc_plasma_vol(self, nrhos=20, ntheta=1000):
        # Use Green theorem to find area of plasma boundary

        ds = self.ds
        aminor = ds['Rmajor']/ds['aspect_ratio']

        Rbnd = ds['Rmajor'] + aminor * ds['rhos'] * \
            np.cos(ds['theta'] + np.arcsin(ds['tri']) * np.sin(ds['theta']))

        Zbnd = ds['elong'] * aminor * ds['rhos'] * np.sin(ds['theta'])

        dR = Rbnd.diff('theta', label='lower')
        dZ = Zbnd.diff('theta', label='lower')

        greens = (Rbnd.isel(theta=slice(0, -1))*dZ -
                  Zbnd.isel(theta=slice(0, -1))*dR)

        area = greens.sum(dim='theta') * 0.5

        darea = area.diff('rhos')

        vol = 2 * pi * ds['Rmajor'] * area
        dvol = 2 * pi * ds['Rmajor'] * darea

        self.area = area
        self.darea = darea
        self.vol = vol
        self.dvol = dvol

    def gen_kinetic_profiles(self):
        # Generates kinetic profiles using profile form in
        # Fus Eng Des 87 (2012) 787 - 792

        ds = self.ds

        Tped = ds['T0']*ds['Tped_height']
        nped = ds['n0']*ds['nped_height']

        Tsep = ds['T0']*ds['Tsep_height']
        nsep = ds['n0']*ds['nsep_height']

        # Core profiles
        Tcore = Tped + (ds['T0'] - Tped) * \
            (1 - (ds['rhos']/ds['ped_rad'])**2) ** ds['alpha_T']
        ncore = nped + (ds['n0'] - nped) * \
            (1 - (ds['rhos']/ds['ped_rad'])**2) ** ds['alpha_n']

        # Edge profiles
        Tedge = Tsep + (Tped - Tsep) * (1 - ds['rhos'])/(1 - ds['ped_rad'])
        nedge = nsep + (nped - nsep) * (1 - ds['rhos'])/(1 - ds['ped_rad'])

        # Use core profile inside ped radius, else use edge profile
        T = Tcore.where(Tcore.rhos < ds['ped_rad'], Tedge)
        n = ncore.where(ncore.rhos < ds['ped_rad'], nedge)

        self.ds['T'] = T
        self.ds['n'] = n

    def calc_fusion_power(self, fuspow_type=0):
        # Calculates fusion power using 3 differents rule
        # fuspow_type = 0: Fus Eng Des 87 (2012) 787 - 792
        # fuspow_type = 1: Parametric form in Wesson Tokamaks 4th edition
        # fuspow_type = 2: Equation used in SCENE for power

        ds = self.ds
        self.fuspow_type = fuspow_type

        self.gen_kinetic_profiles()

        T = ds['T']
        n = ds['n']

        if fuspow_type == 0:
            # Sadler-Van Belle formula

            C1 = 2.5663271e-18
            C2 = 19.983026
            C3 = 2.5077133e-2
            C4 = 2.5773408e-3
            C5 = 6.1880463e-5
            C6 = 6.6024089e-2
            C7 = 8.1215505e-3

            U = 1 - T * (C3 + T*(C4 - C5*T))/(1 + T * (C6 + C7*T))

            rr = C1/(U**(5./6.) * T**(2/3)) * np.exp(-C2 * (U/T)**(1/3))

            pden = n**2 * 1.0e38 / 4 * rr * 17.5e6 * 1.6e-19 *  \
                self.dvol * 1e-9

        elif fuspow_type == 1:
            # Wesson fourth edition

            pden = 1.1e-24 * T**2 * n**2 / \
                4 * 17.5e6 * 1.6e-19 * self.dvol * 1e-9 * 1e38

        else:
            # SCENE's method

            arg = -0.476 * np.abs(np.log(1.45e-5*T*1e3))**2.25

            pden = 5 * 1.27e4 * n**2 * np.exp(arg) * self.dvol * 1e-9

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

    keys = da.coords.to_index()
    values = da.data.flatten()
    filter_params = []

    for i in range(len(values)):
        if not np.isnan(values[i]):
            filter_params.append(keys[i])

    return filter_params
