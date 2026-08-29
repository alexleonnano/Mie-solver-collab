#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Saturday August 29 2026 AEST

Computes `Wavelength vs Radii` spectra table from `miepython` package.
Generated file intended for comparision of scatting spectra calculated from `miepython` and reference spectra tables (William's)

@author: william.mcm.p
"""

# TODO: make the src import less verbose
import sys
from pathlib import Path

# Add project root (Mie-solver-collab) to python path
root_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(root_dir)) 

from src import lambda_interpolate as interp
import numpy as np
import matplotlib.pyplot as plt
import miepython as mie
import pandas as pd



material = "materials/silicon.csv" # Amorphous Silicon - Pierce 1972


refrative_index_df = pd.read_csv(material, delimiter = ",", header = 0, names = ["wl", "n", "k"], dtype = float)

# interpolating refrative index before computing Mie spectra
wl_step = 5 # nm step size between each wl sample
wl = np.arange(400, refrative_index_df['wl'].max()*1e3, wl_step)*1e-3 # in microns

# interpolate as required for high resolution
n, k = interp.wvl_interpolation(wl, 
                                # Access pd.Series diretly from the df
                                refrative_index_df['wl'], 
                                refrative_index_df['n'],
                                refrative_index_df['k'],
                                )

radius_range = np.arange(0,150, 5) * 1e-3 # in microns
m = n - 1.0j * k # complex refrative index

aDF = pd.DataFrame({"Wavelength / Radius (nm)": wl * 1e3})

for radii in radius_range:

    x = 2 * np.pi * radii / wl
    qext, qsca, qback, g = mie.efficiencies_mx(m, x)

    aDF[radii] = qsca

aDF.to_csv("tests/test_sca_amor_silicon_Pierce_1972.csv", index=False)