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

# Calaculating the intersecting wl and particle radii


