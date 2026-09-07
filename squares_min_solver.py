#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday September 02 2026

Intened for a first pass of using minimisation of residual 
squares to find the best fitting spectra to experimentally measured data

@author: william.mcm.p
"""

from src import load_material as tools
from src import lambda_interpolate as interp
import numpy as np
import matplotlib.pyplot as plt
import miepython as mie
import pandas as pd

# - Helper functions
def normalised_data(data):
    """Noramlise the data between [0,1] - min-max method"""
    return (data - min(data))/(max(data)- min(data))

# Experimental Data Import
exp_data_path = 'experimental_data/SiNP_crystaline_scat.csv'
exp_data_df = pd.read_csv(exp_data_path, header=0,names=['wl', 'scat'])

# Material Refrative Index Data Import
material = "materials/Silicon - crystalline - Green 2008.txt"
refrative_index_df = pd.read_csv(material,  sep=r'\s+', header = 0, names = ["wl", "n", "k"], dtype = float)

# Setting up global wl range for interpolation
wl_step = 1 # nm step size between each wl sample
wavelengths = np.arange(min(exp_data_df['wl']), max(exp_data_df['wl']), wl_step) *1e-3 # wavelength range to me tested with

# Experimental Data Interpolation
exp_scat = np.interp(wavelengths, exp_data_df['wl']*1e-3, exp_data_df['scat'])
norm_exp_scat = normalised_data(exp_scat) # normalising experimental data

# Refrative Index Data Interpolation
n, k = interp.wvl_interpolation(wavelengths, 
                                refrative_index_df['wl'], 
                                refrative_index_df['n'],
                                refrative_index_df['k'],
                                )

m = n - 1.0j * k # complex refrative index



# range of to be fitted particles
# TODO: Radius resolution definable
radius_range = np.linspace(5, 200, 100) * 1e-3 # in micron

# Pre allocated rms array
rms = np.zeros_like(radius_range)


# Loop over each particle size, calculate the mie theory spectra
# then compare the measured spectra to what was theory. 
# Radii that has the smalled RMS to the measured data is the best fitted particle.
for i, radii in enumerate(radius_range):
    # NOTE: Could be made faster with parallel compute

    # Define the size paramter for the particle
    x = 2 * np.pi * radii / wavelengths


    # Calculate the Mie spectra for a given particle size
    qext, qsca, qback, g = mie.efficiencies_mx(m, x)

    # Assumes measured data is scattering spectra.
    norm_qscat = normalised_data(qsca)

    residuals = (norm_qscat- norm_exp_scat) ** 2
    rms_calc = np.sqrt(np.mean(residuals))

    rms[i] = rms_calc   

fitted_radii = radius_range[np.argmin(rms)] # in microns
print(f'Wavelength where RMS is minimised = {fitted_radii*1e3:.0f} nm')

# -------------------
#  - Profiling Plot -
# -------------------

fig, axes = plt.subplots(3,1, figsize=(5,7))

# Subplot 1. Measured and fitted theory spectra
ax=axes[0]

ax.plot(wavelengths*1e3, norm_exp_scat, label='Measured Scat.')

x = 2 * np.pi * fitted_radii / wavelengths
qext, qsca, qback, g = mie.efficiencies_mx(m, x)
norm_qscat = normalised_data(qsca)

ax.plot(wavelengths*1e3, norm_qscat, color="firebrick", label=f"Crystalline - r={fitted_radii*1e3} nm")

ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("Intensity (a.u.)")
ax.legend()

# Subplot 2. Residuals Profile for fitted particle radius
ax2=axes[1]

residuals = (norm_qscat- norm_exp_scat) ** 2
ax2.plot(residuals)
ax2.set_title('Residuals')

# Subplot 3. RMS score for all tested particle sizes
axes[2].plot(radius_range*1e3, rms)
axes[2].set_title('RMS')
axes[2].set_xlabel('Particle Radius (nm)')
axes[2].set_ylabel('RMS')
axes[1].set_ylabel('Residual')

fig.tight_layout()
plt.show()
