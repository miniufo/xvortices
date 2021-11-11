# -*- coding: utf-8 -*-
"""
Created on 2021.11.07

@author: MiniUFO
Copyright 2018. All rights reserved. Use is subject to license terms.
"""

#%% cal lat/lon in cylindrical
import xarray as xr
from besttracks.besttracks import parse_TCs

condTC = lambda tc:tc.year==2004 and tc.name=='Haima'
condRc = None

haima = parse_TCs('d:/Data/Typhoons/CMA/original/*.txt',
                  rec_cond=condRc, tc_cond=condTC,
                  agency='CMA')[0]
dset = xr.open_dataset('E:/OneDrive/Python/MyPack/xvortices/tests/Haima2004.nc')


# align times for best-track and gridded dataset
tstr = (haima.records.TIME.iloc[ 0] if haima.records.TIME.iloc[ 0]>dset.time[ 0]
            else dset.time[ 0]).values
tend = (haima.records.TIME.iloc[-1] if haima.records.TIME.iloc[-1]<dset.time[-1]
            else dset.time[-1]).values

haima = haima.sel(lambda df:(tstr<=df.TIME) & (df.TIME<=tend)).translate_velocity()
dset  = dset.sel({'time':slice(tstr, tend)})


#%%
import xarray as xr
from xvortices.xvortices import load_cylind, project_to_cylind, storm_relative

azimNum, radiNum, radMax = 72, 31, 6

olon = haima.get_as_xarray('LON')
olat = haima.get_as_xarray('LAT')

[u, v, w, h], lons, lats, etas = load_cylind(dset,
                                             olon=olon, olat=olat,
                                             azimNum=azimNum, radiNum=radiNum,
                                             radMax=radMax)
uaz, vra = project_to_cylind(u, v, etas)

uaz_rel, vra_rel = storm_relative(haima.get_as_xarray('uo'),
                                  haima.get_as_xarray('vo'), uaz, vra)

#%%
uazm = uaz.mean('azim')
uazrm = uaz_rel.mean('azim')
vram = vra.mean('azim')
vrarm = vra_rel.mean('azim')

#%% plot azimuthal storm-relative
import proplot as pplt
import numpy as np

tlev = 4

fig, axes = pplt.subplots(nrows=2, ncols=2, figsize=(12,10), sharex=3, sharey=3)

fontsize = 16
levels = np.linspace(-10, 16, 27) # u

ax = axes[0,0]
m=ax.contourf(uazm[tlev], levels=levels, cmap='jet')
ax.set_title('azimuth-mean rotational wind', fontsize=fontsize)
ax.invert_yaxis()

ax = axes[0,1]
m=ax.contourf(uazrm[tlev], levels=levels, cmap='jet')
ax.set_title('azimuth-mean rotational wind (relative)', fontsize=fontsize)
ax.colorbar(m, loc='r', length=1)

ax = axes[1,0]
m=ax.contourf(vram[tlev], levels=np.linspace(-5, 6, 23), cmap='jet')
ax.set_title('azimuth-mean radial wind', fontsize=fontsize)
ax.invert_yaxis()

ax = axes[1,1]
m=ax.contourf(vrarm[tlev], levels=np.linspace(-5, 6, 23), cmap='jet')
ax.set_title('azimuth-mean radial wind (relative)', fontsize=fontsize)
ax.colorbar(m, loc='r', length=1)


axes.format(abc='(a)', xlabel='radius (degree)', ylabel='pressure (hPa)',
            yscale='log', ylocator=np.linspace(100, 1000, 10))


#%% plot horizontal interpolation
import proplot as pplt
import numpy as np

zlev = 0
tlev = 16

fig, axes = pplt.subplots(nrows=2, ncols=2, figsize=(12,9), sharex=3, sharey=3,
                          proj=pplt.Proj('cyl', lon_0=0))

fontsize = 16
levels1 = np.linspace(-5, 12, 11) # ut
levels2 = np.linspace(-12, 8, 11) # vr

ax = axes[0,0]
m=ax.tricontourf(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
                 uaz[tlev,zlev].data.ravel(), levels=levels1, cmap='jet')
ax.set_title('azimuth velocity', fontsize=fontsize)

ax = axes[0,1]
m=ax.tricontourf(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
                  uaz_rel[tlev,zlev].data.ravel(), levels=levels1, cmap='jet')
ax.set_title('azimuth velocity (relative)', fontsize=fontsize)
ax.colorbar(m, loc='r', length=0.9)

ax = axes[1,0]
m=ax.tricontourf(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
                  vra[tlev,zlev].data.ravel(), levels=levels2, cmap='jet')
ax.set_title('radial velocity', fontsize=fontsize)

ax = axes[1,1]
m=ax.tricontourf(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
                  vra_rel[tlev,zlev].data.ravel(), levels=levels2, cmap='jet')
ax.set_title('radial velocity (relative)', fontsize=fontsize)
ax.colorbar(m, loc='r', length=0.9)

axes.format(abc='(a)', lonlabels='b', latlabels='l',
            latlim=[10, 60], lonlim=[85, 145],
            coast=True, lonlines=5, latlines=5)


#%% plot horizontal interpolation
import proplot as pplt

zlev = 0
tlev = 16
varN = ut_rel
varO = u

fig, axes = pplt.subplots(nrows=1, ncols=3, figsize=(12,4.2), sharex=3, sharey=3,
                          proj=pplt.Proj('cyl', lon_0=0))

fontsize = 16
levels = np.linspace(-8, 12, 11) # ut
# levels = np.linspace(-0.5, 0.5, 11) # w
# levels = np.linspace(287, 301, 15) # t

ax = axes[0]
m=ax.scatter(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
             c=varN[tlev,zlev].data.ravel(),
             s=4, levels=levels, cmap='jet')
ax.set_title('var in cylindrical coord.', fontsize=fontsize)
ax.set_ylim([29, 51])
ax.set_xlim([103, 133])

ax = axes[1]
m=ax.tricontourf(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
                 varN[tlev,zlev].data.ravel(), levels=levels, cmap='jet')
ax.set_title('var in cylindrical coord.', fontsize=fontsize)
ax.set_ylim([29, 51])
ax.set_xlim([103, 133])

ax = axes[2]
m=ax.contourf(varO[tlev,zlev], levels=levels, cmap='jet')
ax.set_title('var in spherical coord.', fontsize=fontsize)
ax.set_ylim([29, 51])
ax.set_xlim([103, 133])

fig.colorbar(m, loc='b', length=1)

axes.format(abc=True, abcloc='l', abcstyle='(a)', lonlabels='b', latlabels='l',
            coast=True, lonlines=5, latlines=3)



#%% plot horizontal velocity
import proplot as pplt

zlev = 3
tlev = 3

fig, axes = pplt.subplots(nrows=2, ncols=2, figsize=(12,10), sharex=3, sharey=3,
                          proj=pplt.Proj('cyl', lon_0=0))

fontsize = 16

ax = axes[0,0]
ax.quiver(lons[tlev].data, lats[tlev].data,
          vs[0][tlev, zlev].data, vs[1][tlev, zlev].data,
          width=0.0016, headwidth=12., headlength=12., scale=200)
ax.set_title('velocity in cylindrical coord.', fontsize=fontsize)
ax.set_ylim([14, 36])
ax.set_xlim([108, 135])

ax = axes[0,1]
ax.quiver(u.lon.data, u.lat.data,
          u[tlev, zlev].data, v[tlev, zlev].data,
          width=0.0016, headwidth=12., headlength=12., scale=200)
ax.set_title('velocity in cylindrical coord.', fontsize=fontsize)
ax.set_ylim([14, 36])
ax.set_xlim([108, 135])

ax = axes[1,0]
m=ax.tricontourf(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
                 ut[tlev,zlev].data.ravel(), cmap='jet')
ax.colorbar(m, loc='b', length=1)
ax.quiver(lons[tlev].data, lats[tlev].data,
          vs[0][tlev, zlev].data, vs[1][tlev, zlev].data,
          width=0.0016, headwidth=12., headlength=12., scale=200)
ax.set_title('azimuth component', fontsize=fontsize)
ax.set_ylim([14, 36])
ax.set_xlim([108, 135])

ax = axes[1,1]
m=ax.tricontourf(lons[tlev].data.ravel(), lats[tlev].data.ravel(),
                 vr[tlev,zlev].data.ravel(), cmap='jet')
ax.colorbar(m, loc='b', length=1)
ax.quiver(lons[tlev].data, lats[tlev].data,
          vs[0][tlev, zlev].data, vs[1][tlev, zlev].data,
          width=0.0016, headwidth=12., headlength=12., scale=200)
ax.set_title('radial component', fontsize=fontsize)
ax.set_ylim([14, 36])
ax.set_xlim([108, 135])


axes.format(abc='(a)', lonlabels='b', latlabels='l',
            coast=True, lonlines=5, latlines=3, grid=False)



