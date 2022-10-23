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

# ## Forced outputs

data = xr.open_mfdataset('apecosm-forced/output_apecosm/*000').isel(x=1, y=1)
data

forced = data['OOPE']
forced

forced_forage = data['FORAGE']
forced_forage

# ## One way coupling

data = xr.open_mfdataset('one-way/output_apecosm/*nc').isel(x=1, y=1)
data

oneway = data['OOPE']
oneway

oneway_forage = data['FORAGE']
oneway_forage

# ## Two way coupling

data = xr.open_mfdataset('two-ways/output_apecosm/*nc').isel(x=1, y=1)
data

twoway = data['OOPE']
twoway

twoway_forage = data['FORAGE']
twoway_forage

# ## Making plots

l = 20; c = 2
forced.isel(w=l, c=c).plot(label='forced')
oneway.isel(w=l, c=c).plot(label='one way')
twoway.isel(w=l, c=c).plot(label='two way')
plt.legend()

dn = 1
c = 1
s = 0
z = slice(0, 20)
plt.figure()
ax1 = plt.subplot(211)
cs1 = forced_forage.isel(dn=dn, c=c, size_group=s, depth=z).T.plot(robust=True)
plt.title('forced')
ax2 = plt.subplot(212, sharex=ax1)
plt.title('coupled')
cs2 = oneway_forage.isel(dn=dn, c=c, size_group=s, depth=z).T.plot(robust=True)
cs2.set_clim(cs1.get_clim())


