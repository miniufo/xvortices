# -*- coding: utf-8 -*-
'''
Created on 2022.06.16

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
'''
import cartopy.crs as ccrs
import cartopy.feature
import itertools
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.patch import geos_to_path
from matplotlib.collections import PolyCollection, LineCollection


_R_earth = 6371200.0
_deg2m = np.pi / 180.0


'''
Here defines some util functions for this package
'''
def plot3D(lons, lats, da, azimngle=-60, elevangle=30, title=None,
           reverseZ=False, lonR=8, latR=5, fontsize=18, reso='50m',
           figsize=(13, 8), cmap='jet', vmin=None, vmax=None, alpha=0.7):
    """Plot 3D structure of a variable

    Parameters
    ----------
    lons: xarray.DataArray
        Longitudes of the data (3D spatial variable)
    lats: xarray.DataArray
        Latitude of the data (3D spatial variable)
    da: xarray.DataArray
        Data variable (3D spatial variable)
    azimngle: float
        Azimuthal angle for the view
    elevangle: float
        Elevation angle for the view
    title: str
        Title of the plot
    reverseZ: bool
        Reverse the vertical axis.  Useful in pressure coordinate of atmosphere.
    lonR: float
        Radial range of the plot along longitude (degree)
    latR: float
        Radial range of the plot along latitude (degree)
    fontsize: int
        Font size
    reso: str
        Resolution of the cartopy map, one of ['110m', '50m', '10m']
    figsize: tuple
        Figure size
    cmap: str
        Colormap for the plot
    vmin: float
        Minimum value corresponds to the cmap
    vmax: float
        Maximum value corresponds to the cmap
    alpha: float
        Transparency within 0 ~ 1
    """
    dims = da.dims

    azishift = len(da[dims[-1]]) // 4 # rotate a quarter from north to east

    var  = da  .roll({dims[-1]: azishift}).pad({dims[-1]: (0,1)}, 'wrap')
    lon  = lons.roll({dims[-1]: azishift}).pad({dims[-1]: (0,1)}, 'wrap')
    lat  = lats.roll({dims[-1]: azishift}).pad({dims[-1]: (0,1)}, 'wrap')

    if len(dims) != 3:
        raise Exception('only 3D spatial data can be plotted')

    # prepare some coordinates, and attach rgb values to each
    z, y, x = xr.broadcast(var[dims[0]], lat, lon)

    lonmax, lonmin = lon.max().values, lon.min().values
    latmax, latmin = lat.max().values, lat.min().values

    xrange = [lonmin - lonR, lonmax + lonR]
    yrange = [latmin - latR, latmax + latR]
    zrange = [z.min(), z.max()]

    if reverseZ:
        zrange = zrange[::-1]

    if title == None:
        title = '3D structure' + ('' if var.name == None else ' of ' + var.name)

    # combine the color components
    vmin = var.min().values if vmin == None else vmin
    vmax = var.max().values if vmax == None else vmax

    norm   = colors.Normalize(vmin=vmin, vmax=vmax, clip=False)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    fcolor = mapper.to_rgba(var.values.ravel()).reshape(var.shape+(4,))
    fcolor = fcolor[:-1, :-1, :-1, :]
    fcolor[..., 3] = alpha  # change transparency

    ################# plot 3D cylind #################
    fig = plt.figure(figsize=figsize)
    # fig.subplots_adjust(left=0, right=0.5, bottom=0.3, top=1)
    ax3d = fig.add_axes([0, 0, 1, 1], projection='3d',
                        xlim=xrange, ylim=yrange, zlim=zrange)
    ax3d.set_box_aspect((2.0*lonR/latR*0.86, 2, 1))

    cond = np.ones_like(x[:-1, :-1, :-1])
    cond[..., azishift*3:] = 0
    ax3d.voxels(x, y, z, filled=cond, facecolors=fcolor, edgecolors=fcolor,
                linewidth=0, zorder=200)
    ax3d.set_title(title, fontsize=fontsize)
    ax3d.view_init(elevangle, azimngle)
    ax3d.set_xlabel('longitude', fontsize=fontsize-2)
    ax3d.set_ylabel('latitude' , fontsize=fontsize-2)
    ax3d.set_zlabel(var[dims[0]].name, fontsize=fontsize-2)

    ############## get the extent as a shapely geometry and clip ##############
    ax2d = plt.figure().add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree(),
                                 xlim=xrange, ylim=yrange)
    ax2d.contourf(x[0], y[0], z[0], alpha=0, transform=ccrs.PlateCarree())
    clip_geom = ax2d._get_extent_geom().buffer(0)
    plt.close(ax2d.figure)

    LAND = NaturalEarthFeature('physical', 'land', reso, edgecolor='face',
                facecolor=np.array((240, 240, 220)) / 256., zorder=-1)

    OCEAN = NaturalEarthFeature('physical', 'ocean', reso, edgecolor='face',
                facecolor=np.array((152, 183, 226)) / 256., zorder=-1)

    add_feature3d(ax3d, OCEAN, clip_geom, zs=zrange[0])
    add_feature3d(ax3d, LAND , clip_geom, zs=zrange[0])
    # add_feature3d(ax3d, cartopy.feature.COASTLINE, zs=zrange[0])

    plt.show()


"""
Below are the private helper methods

Some code snippet from:
https://stackoverflow.com/questions/23785408/3d-cartopy-similar-to-matplotlib-basemap
"""
def add_feature3d(ax3d, feature, clip_geom=None, zs=None):
    """Add cartopy feature to a given 3D axes

    Parameters
    ----------
    ax3d: mpl_toolkits.mplot3d.Axes3D
        A given 3D axes
    feature: cartopy.feature
        A cartopy feature of land or ocean or coastline
    clip_geom: numpy.ndarray
        An array to clip the feature
    zs: float
        Which z level to add the feature
    """
    concat = lambda iterable: list(itertools.chain.from_iterable(iterable))

    target_projection = ccrs.PlateCarree()
    geoms = list(feature.geometries())

    if target_projection != feature.crs:
        # Transform the geometries from the feature's CRS into the
        # desired projection.
        geoms = [target_projection.project_geometry(geom, feature.crs)
                 for geom in geoms]

    if clip_geom:
        # Clip the geometries based on the extent of the map
        # (because mpl3d can't do it for us)
        geoms = [geom.intersection(clip_geom) for geom in geoms if geom.is_valid]

    # Convert the geometries to paths so we can use them in matplotlib.
    paths = concat(geos_to_path(geom) for geom in geoms if geom.is_valid)

    # Bug: mpl3d can't handle edgecolor='face'
    kwargs = feature.kwargs
    if kwargs.get('edgecolor') == 'face':
        kwargs['edgecolor'] = kwargs['facecolor']

    polys = concat(path.to_polygons(closed_only=False) for path in paths)

    fcolor = kwargs.get('facecolor', 'none')

    if isinstance(fcolor, str) and fcolor == 'none':
        lc = LineCollection(polys, **kwargs)
    else:
        lc = PolyCollection(polys, closed=False, **kwargs)

    ax3d.add_collection3d(lc, zs=zs)

