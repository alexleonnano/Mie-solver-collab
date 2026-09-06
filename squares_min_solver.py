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

exp_data_path = 'experimental_data/SiNP_crystaline_scat.csv'
exp_data_df = pd.read_csv(exp_data_path, header=0,names=['wl', 'scat'])
norm_exp = normalised_data(exp_data_df['scat'])



material = "materials/Silicon - crystalline - Green 2008.txt"
refrative_index_df = pd.read_csv(material,  sep=r'\s+', header = 0, names = ["wl", "n", "k"], dtype = float)

# NOTE: using the experimental data wl sample count is ok-ish. 
#       -> Should decouple them
n, k = interp.wvl_interpolation(exp_data_df['wl'],  
                                # Access pd.Series diretly from the df
                                refrative_index_df['wl'], 
                                refrative_index_df['n'],
                                refrative_index_df['k'],
                                )

m = n - 1.0j * k # complex refrative index



# range of to be fitted particles
# NOTE: this maybe a large array (1340)
radius_range = np.linspace(5, 200, len(norm_exp)) * 1e-3 # in micron

# Pre allocated rms array
rms = np.zeros_like(radius_range)


# Loop over each particle size, calculate the mie theory spectra
# then compare the measured spectra to what was theory. 
# Radii that has the smalled RMS to the measured data is the best fitted particle.
for i, radii in enumerate(radius_range):
    # NOTE: High N count of wl samples (1340) it takes ~5-10 s to complete all itterations. 
    #       Could be made faster with parallel compute

    # Define the size paramter for the particle
    x = 2 * np.pi * radii / exp_data_df['wl']


    # Calculate the Mie spectra for a given particle size
    qext, qsca, qback, g = mie.efficiencies_mx(m, x, 0)

    # Assumes measured data is scattering spectra.
    norm_qscat = normalised_data(qsca)

    residuals = (norm_qscat- norm_exp) ** 2
    rms_calc = np.sqrt(np.mean(residuals))

    rms[i] = rms_calc   

print(f'Wavelength where RMS is minimised = {radius_range[np.argmin(rms)]*1e3:.0f} nm')

# -------------------
#  - Profiling Plot -
# -------------------

fig, axes = plt.subplots(3,1, figsize=(5,7))

ax=axes[0]

# d = pre_df[(pre_df['wavelength (nm)'] > cryst_scat.index.min()) & (pre_df['wavelength (nm)'] < 670) ]
# cryst_filtered = cryst_scat[(cryst_scat.index >= d['wavelength (nm)'].min()) & (cryst_scat.index <= 670)]

ax.plot(exp_data_df['wl'], normalised_data(exp_data_df['scat']), label='Measured Scat.')
# ax.plot(, normalize_data(cryst_filtered[nearest_radius]), color="firebrick", label=f"Crystalline - r={radii} nm")

ax.set_xlabel("wavelength (nm)")
ax.set_ylabel("Intensity (a.u.)")
ax.legend()

ax2=axes[1]

# interp_post_sig = np.interp(
#             d['wavelength (nm)'],                # Target X-grid (pre_df wavelengths)
#             scat.index,         # Source X-grid (shifted wavelengths)
#             scat[radii]      # Source Y-data
#         )

# norm_pre = normalize_data(d['processed_sig'].values)
# norm_post = normalize_data(interp_post_sig)

# residuals = (norm_pre - norm_post) ** 2

# ax2.plot(residuals)
# ax2.set_title('Residuals')


axes[2].plot(exp_data_df['wl'], rms)
axes[2].set_title('RMS')
axes[2].set_xlabel('Offset (nm)')
axes[2].set_ylabel('RMS')
axes[1].set_ylabel('Residual')

fig.tight_layout()
plt.show()
