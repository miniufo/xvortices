# -*- coding: utf-8 -*-
'''
Created on 2021.11.07

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
'''
import numpy as np
import xarray as xr
from numpy import deg2rad, sin, cos, arcsin, arccos, arctan2, hypot


'''
Here defines the core function of the data interpolation
'''
def load_cylind(ds, olon, olat, azimNum=36, radiNum=11, radMax=10,
                lonname='lon', latname='lat'):
    """
    Load scalar data from a lat/lon grid to a cylindrical grid translating
    with a vortex.

    Parameters
    ----------
    ds: xarray.DataArray or a xarray.Dataset or (list of) DataArray
        A given lat/lon grid variable or dataset to be interpolated
    olon: (list of) float, numpy.array, or xarray.DataArray
        Central longitude of the cylindrical coordinate, in degree
    olat: (list of) float, numpy.array, or xarray.DataArray
        Central latitude of the cylindrical coordinate, in degree
    azimNum: int
        Number of azimuthal grid points
    radiNum: int
        Number of radial grid points
    radMax: float
        Maximum radius in degree
    lonname: str
        Name of longitude in ds
    latname: str
        Name of latitude in ds
    """
    azim = xr.DataArray(np.linspace(0, 360-360/azimNum, azimNum),
                        dims='azim',
                        coords={'azim':np.linspace(0, 360-360/azimNum, azimNum)})
    radi = xr.DataArray(np.linspace(0, radMax, radiNum),
                        dims='radi',
                        coords={'radi':np.linspace(0, radMax, radiNum)})
    
    olon_r = deg2rad(olon)
    olat_r = deg2rad(olat)
    azim_r = deg2rad(azim)
    radi_r = deg2rad(radi)
    
    lats_r = arcsin(sin(olat_r)*cos(radi_r) + cos(olat_r)*sin(radi_r)*cos(azim_r))
    dlam_r = 1.0/cos(lats_r) * arcsin(sin(radi_r)*sin(azim_r))
    lons_r = olon_r - dlam_r
    etas_r = arccos(sin(olat_r)*sin(dlam_r)*sin(azim_r) - cos(dlam_r)*cos(azim_r))
    
    etas_r = xr.where(etas_r.azim<180, -etas_r+np.pi, etas_r+np.pi)
    
    lats = np.rad2deg(lats_r)
    lons = np.rad2deg(lons_r)
    # etas = np.rad2deg(etas_r)
    
    if type(ds) in [list, np.ndarray, np.array]:
        vs_interp = [v.interp(coords={lonname:lons, latname:lats}
                              ).drop_vars([latname,lonname]) for v in ds]
    elif type(ds) in [xr.Dataset]:
        vs_interp = [ds[v].interp(coords={lonname:lons, latname:lats}
                              ).drop_vars([latname,lonname]) for v in ds.data_vars]
    else:
        vs_interp = ds.interp(coords={lonname:lons, latname:lats}
                              ).drop_vars([latname,lonname])
    
    return vs_interp, lons, lats, etas_r


def project_to_cylind(u, v, etas):
    """
    Re-project zonal/meridional (u/v) velocity components onto
    azimuthal/radial (uaz/vra) components.

    Parameters
    ----------
    u: xarray.DataArray
        Zonal velocity component
    v: xarray.DataArray
        Meridional velocity component
    etas: xarray.DataArray
        Local angle between radial direction and local north.
    """
    uaz = -u*cos(etas) - v*sin(etas) # azimuth component
    vra = -u*sin(etas) + v*cos(etas) # radial  component
    
    return uaz.rename('ut'), vra.rename('vr')


def storm_relative(uc, vc, uaz, vra):
    """
    Calculate storm-relative velocity given the translating
    velocity of the center.

    Parameters
    ----------
    uc: xarray.DataArray
        Zonal velocity of the center.
    vc: xarray.DataArray
        Meridional velocity of the center.
    uaz: xarray.DataArray
        Azimuth velocity component
    vra: xarray.DataArray
        radial velocity component
    """
    cd = arctan2(vc, uc)
    cs = hypot(uc, vc)
    
    ang = cd - deg2rad(uaz.azim) - np.pi/2.0
    
    cAzim = sin(ang) * cs
    cRadi = cos(ang) * cs
    
    uaz_rel = uaz - cAzim
    vra_rel = vra - cRadi
    
    return uaz_rel, vra_rel

