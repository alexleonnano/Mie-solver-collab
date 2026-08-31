#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday August 31 2026

Used `calculate_spectra_table.py` to generate the Wavelength vs Radii spectra table 
for silicon from `materails/silicon.csv` (Amorphous Silicon - Pierce 1972) 
    -> generated spectra table `test_sca_amor_silicon_Pierce_1972.csv`

The reference spectra table is `silicon_σsca - [400 - 900 nm] - Air - amorphous - Pierce 1972.csv`. 
It was calculated using the same refractive index file (interpolated in a similar manner to 
`src.lambda_interpolate`) in a air medium ($n=1.0003$) uisng a `wolfram` code (Mathematica).

> Please note the reference file is the scatting spectra for the silicon particles.

The producing plot in the script show how the `miepython` package is cabable of repoducing spectra as expected.

@author: william.mcm.p
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# miepython estimation
mie_python_path = 'test_sca_amor_silicon_Pierce_1972.csv'        
miepython_df = pd.read_csv(mie_python_path, index_col=0)
miepython_df.columns = miepython_df.columns.astype(float) *1e3  # in nm
miepython_df.index = miepython_df.index.astype(float) # in nm

# William Scat. reference spectra
w_sca_path = 'silicon_σsca - [400 - 900 nm] - Air - amorphous - Pierce 1972.csv'    
w_sca_df = pd.read_csv(w_sca_path, index_col=0)
w_sca_df.columns = w_sca_df.columns.astype(float) # in nm
w_sca_df.index = w_sca_df.index.astype(float)   # in nm

#  --------------------------------------
#  -- Plot 1: Specific radii compare --
#  --------------------------------------

radius = 100 # nm
t_cross_section = (np.pi * (radius*1e-16)) # could be the wrong scale

fig, ax = plt.subplots(1)
ax.plot(w_sca_df.index, w_sca_df[float(radius)]/t_cross_section, label='William')
ax.plot(miepython_df.index, miepython_df[float(radius)], label='miepython')

ax.legend()
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Scatting Coefficient")
ax.set_xlim(400,900)
ax.set_title(f'Scatting coefficient of amorphous silicon (r={radius} nm)')
fig.tight_layout()
fig.show()


#  ----------------------------------------------
#  -- Plot 2: RMS between miepython and w_scat --
#  ----------------------------------------------

radius_values = w_sca_df.columns.astype(float)
scaling_factors = 1.0 / (np.pi * (radius_values * 1e-16))
df1_scaled = w_sca_df.multiply(scaling_factors, axis="columns")

# Aligning for common wl and radii
common_wl = df1_scaled.index.intersection(miepython_df.index)
common_radii = df1_scaled.columns.intersection(miepython_df.columns)

wl_min, wl_max = common_wl.min(), common_wl.max()
rad_min, rad_max = common_radii.min(), common_radii.max()

# Crop overlapping boundaries
high_res_df1 = (
    df1_scaled.loc[wl_min:wl_max, rad_min:rad_max].sort_index(axis=0).sort_index(axis=1)
)
high_res_df2 = miepython_df.loc[wl_min:wl_max, rad_min:rad_max].sort_index(axis=0).sort_index(axis=1)

# Get matching grid points for residuals
sub_df1 = df1_scaled.loc[common_wl, common_radii].sort_index(axis=0).sort_index(axis=1)
sub_df2 = miepython_df.loc[common_wl, common_radii].sort_index(axis=0).sort_index(axis=1)

residuals = sub_df1 - sub_df2
abs_max = float(residuals.abs().max().max())


fig, axes = plt.subplots(1, 3, figsize=(18, 5))
extent = [rad_min, rad_max, wl_max, wl_min]


im1 = axes[0].imshow(high_res_df1.values, extent=extent, aspect="auto", cmap="viridis")
axes[0].set_title(f"William Scat. {high_res_df1.shape}")
fig.colorbar(im1, ax=axes[0], label="Scatting Coefficient")

# Spectrum 2 (Original)
im2 = axes[1].imshow(high_res_df2.values, extent=extent, aspect="auto", cmap="viridis")
axes[1].set_title(f"miepython package Scat. {high_res_df2.shape}")
fig.colorbar(im2, ax=axes[1], label="Scatting Coefficient")

# Residuals imshow
im3 = axes[2].imshow(
    residuals.values, extent=extent, aspect="auto", cmap="bwr", vmin=-abs_max, vmax=abs_max
)
axes[2].set_title(f"Residuals (Intersect Grid: {residuals.shape})\nWilliam - miepython")
fig.colorbar(im3, ax=axes[2], label="Difference of Scatting Coefficients")

# Labels
for ax in axes:
    ax.set_xlabel("Particle Radius (nm)")
    ax.set_ylabel("Wavelength (nm)")
    ax.set_xlim(rad_min, rad_max)
    ax.set_ylim(wl_max, wl_min)

plt.tight_layout()
plt.show()


