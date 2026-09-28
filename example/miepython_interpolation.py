#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday September 22 2026

@author: william.mcm.p (copied from `test_miepython_interpolation.ipynb`)
"""

# TODO: make the src import less verbose
import sys
from pathlib import Path

# Add project root (Mie-solver-collab) to python path
root_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(root_dir)) 


# Load libraries
from src import lambda_interpolate as interp
import numpy as np
import matplotlib.pyplot as plt
import miepython as mie
import pandas as pd

# Load material data and preallocate arrays
material = "materials/silicon.csv"
ref_lam = []
ref_n = []
ref_k = []


# Create daframes and double check sizes for correct calculations
df = pd.read_csv(material, delimiter = ",", header = 0, names = ["wl", "n", "k"], dtype = float)
ref_lam = df["wl"]
ref_n = df["n"]
ref_k = df["k"]

wvl_range = np.linspace(min(ref_lam), max(ref_lam), 1000)  # Wavelength range for interpolation in microns

print(len(ref_lam), len(ref_n), len(ref_k))

# Interpolate
n, k = interp.wvl_interpolation(wvl_range, ref_lam, ref_n, ref_k)

# Cross section using interpolation
radius = 0.2  # in microns
m_i = n - 1.0j * k
x_i = 2 * np.pi * radius / wvl_range
cross_section_area = np.pi * radius**2
mu_i = 4 * np.pi * k / wvl_range  # nm
qext_i, qsca_i, qback_i, g_i = mie.efficiencies_mx(m_i, x_i)
sca_cross_section_i = qsca_i * cross_section_area
abs_cross_section_i = (qext_i - qsca_i) * cross_section_area

len(wvl_range), len(mu_i), len(sca_cross_section_i), len(abs_cross_section_i), len(k)


# TEST Cross section using reference data -- checking it has the right length -- WORKING CORRECTLY
m = ref_n - 1.0j * ref_k
x = 2 * np.pi * radius / ref_lam
cross_section_area = np.pi * radius**2
mu_a = 4 * np.pi * ref_k / ref_lam  # nm
qext, qsca, qback, g = mie.efficiencies_mx(m, x)

sca_cross_section = qsca * cross_section_area
abs_cross_section = (qext - qsca) * cross_section_area

len(ref_lam), len(mu_a), len(sca_cross_section), len(abs_cross_section), len(ref_k)

# TEST cross section for different radii in Silicon -- William's collaboration -- RUNS CORRECTLY
r_range = [0.005, 0.05, 0.1, 0.15,  0.2]  # in microns
m = n - 1.0j * k

for i in r_range:
    x = 2 * np.pi * i / wvl_range
    cross_section_area = np.pi * i**2
    mu = 4 * np.pi * k / wvl_range  # nm
    qext, qsca, qback, g = mie.efficiencies_mx(m, x)
    sca_cross_section = qsca * cross_section_area
    abs_cross_section = (qext - qsca) * cross_section_area

    plt.figure(1)
    plt.plot(wvl_range*1000, abs_cross_section, label=f"{i*1000} nm")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Cross Section (μm²)")
    plt.title(r'Silicon $\sigma_{abs}$ ')
    plt.legend()
    plt.grid()

    plt.figure(2)
    plt.plot(wvl_range*1000, sca_cross_section, label=f"{i*1000} nm")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Cross Section (μm²)")
    plt.title(r'Silicon $\sigma_{sca}$ ')
    plt.legend()
    plt.grid()

    plt.figure(3)
    plt.plot(wvl_range*1000, qext, label=f"{i*1000} nm")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Extinction Efficiency")
    plt.title(r'Silicon $\sigma_{ext}$ ')
    plt.legend()
    plt.grid()

plt.show()
