#! /usr/bin/env python3

# Load libraries and functions
from src import load_material as tools
from src import lambda_interpolate as itl
import matplotlib.pyplot as plt
import miepython as mie
import numpy as np

# Load data from the refractive index file for gold (Au)
# wavelength in microns
# Load material data and interpolate to a specific wavelength range
material = "materials/Au_reformatted.txt"
df = tools.load(material)

# Radius of nanoparticle and wavelength range for calculations
radius = 0.040  # radius of the nanoparticle in microns
wvl_range = np.linspace(0.400, 0.800, 1000)  # wavelength range in microns

# Pre allocate arrays for wavelength, real part of refractive index, and imaginary part of refractive index
ref_lam = []
ref_n = []
ref_k = []

# Load arrays and turn them into floats
ref_lam = df["wl"].values;
ref_n = df["n"].values;
ref_k = df["k"].values;

# Wavelength range in microns with number of points specified
n, k = itl.wvl_interpolation(wvl_range, ref_lam, ref_n, ref_k)

# Calculate cross sections
m = n - 1.0j * k                                       # Complex refractive index of the material, miepython deals with the complex part of the refractive index as negative
x = 2 * np.pi * radius / wvl_range                     # Size parameter
area = np.pi * radius**2                               # Surface area of the nanoparticle
qext, qsca, qback, g = mie.efficiencies_mx(m, x) 
qabs = qext - qsca                                     # Absorption efficiency      


# Plots
plt.figure(1)
plt.plot(wvl_range * 1000, qsca, "--r")
plt.ylabel("Scattering efficiency")
plt.xlabel("Wavelength (nm)")
plt.legend(["Scattering"])
plt.title(f"Au nanoparticle {radius*1000:.0f} nm radius")

plt.figure(2)
plt.plot(wvl_range * 1000, qabs, '-g')
plt.ylabel("Absorption efficiency")
plt.xlabel("Wavelength (nm)")
plt.legend(["Absorption spectra"])
plt.title(f"Au nanoparticle {radius*1000:.0f} nm radius")

plt.figure(3)
plt.plot(wvl_range * 1000, qext, '-b')
plt.ylabel("Extinction efficiency")
plt.xlabel("Wavelength (nm)")
plt.legend(["Extinction spectra"])
plt.title(f"Au nanoparticle {radius*1000:.0f} nm radius")
plt.show()