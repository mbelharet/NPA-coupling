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

# +
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = False

C_E_convert = 474600
varname = 'PAR'
factornemo = 1#C_E_convert * 1e-3
factorape = 1  
zmax = 200
# -

# ## Reading mesh

mesh = xr.open_dataset('one-way/mesh_mask.nc').isel(x=0, y=0)
mesh

depth = mesh['gdept_1d'].values[0]
depth.shape

iok = np.nonzero(depth <= zmax)[0]
iok

# ## NEMO outputs

data = xr.open_mfdataset('one-way/DYFAMED_[0-9]d*nc').isel(y_grid_T=1, x_grid_T=1, y=1, x=1)
data

nemo = data[varname] * factornemo
nemo

# ## One way coupling

data = xr.open_mfdataset('one-way/output_apecosm/*nc').isel(x=1, y=1)
data

oneway = data['forcing']
oneway

# ## Comparison with old forcings

# +
data = xr.open_mfdataset('/home/barrier/Work/codes/nemo/release-4.0.6/cfgs/C1D_PISCES/EXP00/DYFAMED*1d*nc').isel(y_grid_T=1, x_grid_T=1, y=1, x=1)
data

test = data[varname] * factornemo
test
# -

# ## Making plots

# +
plt.figure()
iii = 20
#ax1 = plt.subplot(311)
ntime = oneway.shape[0]
time = np.arange(ntime)
toplot = oneway
toplot = np.ma.masked_where(toplot==0, toplot)
# cs1 = plt.pcolor(time, -depth[iok], toplot.T[iok, :], shading='auto')
# plt.colorbar(cs1)
plt.plot(toplot[:, iii])
plt.title('Apecosm')

#ax2 = plt.subplot(312)
ntime = nemo.shape[0]
time = np.arange(ntime)
toplot = nemo.values
toplot = np.ma.masked_where(toplot==0, toplot)
#cs2 = plt.pcolor(time, -depth[iok], toplot.T[iok, :], shading='auto')
# plt.colorbar(cs2)
plt.plot(-toplot[:, iii])
plt.title('NEMO')
#cs1.set_clim(cs2.get_clim())
#ax2.set_xlim(ax1.get_xlim())
#ax2.set_ylim(ax1.get_ylim())

#ax3 = plt.subplot(313)
ntime = nemo.shape[0]
time = np.arange(ntime)
toplot = test.values
toplot = np.ma.masked_where(toplot==0, toplot)
# cs2 = plt.pcolor(time, -depth[iok], toplot.T[iok, :], shading='auto')
# plt.colorbar(cs2)
plt.plot(toplot[:, iii], color='blue')
plt.title('TEST')
#cs1.set_clim(cs2.get_clim())
#ax3.set_xlim(ax1.get_xlim())
#ax3.set_ylim(ax1.get_ylim())
# -
toplot[0, iii]


oneway

nemo

idepth = slice(0, 20)
nemo.isel(deptht=idepth).plot()

oneway.isel(depth=idepth).plot()


