#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday September 02 2026

Example code for calculating the individual Mie modes

Ploted are the MD, ED, MQ and EQ modes

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

material = "../materials/silicon.csv" # Amorphous Silicon - Pierce 1972

refrative_index_df = pd.read_csv(material, delimiter = ",", header = 0, names = ["wl", "n", "k"], dtype = float)

# interpolating refrative index before computing Mie spectra
wl_step = 5 # nm step size between each wl sample
wl = np.arange(200, 1200, wl_step)*1e-3 # in microns

n, k = interp.wvl_interpolation(wl, 
                                # Access pd.Series diretly from the df
                                refrative_index_df['wl'], 
                                refrative_index_df['n'],
                                refrative_index_df['k'],
                                )

m = n - 1.0j * k # complex refrative index

# SOURCE: https://miepython.readthedocs.io/en/latest/13_resonance.html#Non-absorbing-spheres

#  NOTE: I lke plotting the MD mode first as it is always the strongest mode (findamental)

# Size of silicon particle
radii = 0.125 # in microns 
x = 2 * np.pi * radii / wl

fig, ax = plt.subplots()

qext, qsca, qback, g = mie.efficiencies_mx(m, x, 0)
qabs = qext-qsca
ax.plot(wl * 1e3, qext, color="black", label="Ext.")
ax.plot(wl * 1e3, qabs, color="orange", label="Abs.")
ax.plot(wl * 1e3, qsca, color="orangered", label="Sca.")


# Magnetic Dipole (b1)
qMDext, qMDsca, qMDback, g = mie.efficiencies_mx(m, x, 1, e_field=False)
ax.plot(wl * 1e3, qMDsca, color='blue', linestyle='--', label="MD")

# Electric Dipole (a1)
qEDext, qEDsca, qEDback, g = mie.efficiencies_mx(m, x, 1, e_field=True)
ax.plot(wl * 1e3, qEDsca, color='red', linestyle='--', label="ED")

# Magnetic Quadrupoles (b2)
qMQext, qMQsca, qMQback, g = mie.efficiencies_mx(m, x, 2, e_field=False)
ax.plot(wl * 1e3, qMQsca, color='blue', linestyle='-.', label="MQ")

# Electric Quadrupoles (a2)
qEQext, qEQsca, qEQback, g = mie.efficiencies_mx(m, x, 2, e_field=True)
ax.plot(wl * 1e3, qEQsca, color='red', linestyle='-.', label="EQ")

ax.legend()
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel(r'$\sigma$ Efficiency')
fig.suptitle(f'Mie Scatteing Modes of SiNP ($r={radii*1e3:.0f}$ nm)')

fig.tight_layout()
plt.show()