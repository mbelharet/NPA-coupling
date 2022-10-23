# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
zmax = 200
varname = 'PHY'

# ## Reading mesh

mesh = xr.open_dataset('one-way/mesh_mask.nc').isel(x=1, y=1)
mesh

depth = mesh['gdept_1d'].values[0]
depth.shape

iok = np.nonzero(depth <= zmax)[0]
iok

depth = depth[iok]
depth

# ## Two way forcings

data = xr.open_mfdataset('one-way/DYFAMED*1d*.nc').isel(x=1, y=1)
data

oneway = data[varname].values
oneway.shape

# ## Two-way forcings

data = xr.open_mfdataset('two-ways/DYFAMED*1d*nc').isel(y_grid_T=1, x_grid_T=1, y=1, x=1)
data

twoway = data[varname].values
twoway.shape

# ## Making plots

iz = 0
plt.plot(oneway[:, iz], label='1w')
plt.plot(twoway[:, iz], label='2w')

iz = iok
ntime, nz = oneway.shape
time = np.arange(ntime)
plt.figure(figsize=(12, 12))
plt.subplot(311)
cs = plt.pcolor(time, -depth, oneway[:, iz].T, shading='auto')
plt.colorbar(cs)
plt.title('1W')
plt.subplot(312)
plt.title('2W')
cs = plt.pcolor(time, -depth, twoway[:, iz].T, shading='auto')
plt.colorbar(cs)
plt.subplot(313)
diff = twoway[:, iz].T - oneway[:, iz].T
perc = np.percentile(np.ravel(np.abs(diff)), 99)
cs = plt.pcolor(time, -depth, diff, shading='auto')
cs.set_clim(-perc, perc)
plt.colorbar(cs)
plt.title('2W - 1W')


