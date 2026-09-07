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

def get_rsm_min(rms: np.array, 
                radius_range : np.array = np.linspace(5,200,100)*1e-3 ) -> float:
    """
    Finds the particle radius corresponding to the minimum RMS error. 
    
    The radius with the lowest RMS value is taken as the best fit between the theoretical and experimental spectra.
    Args: 
        rms (np.array): RMS error values calculated for each radius. 
        radius_range (np.array): Particle radii corresponding to the RMS values. 
        
    Returns: 
        float: Particle radius with the lowest RMS error.
    """
    return radius_range[np.argmin(rms)] # in microns


def calc_rms(theory: np.array, experimental: np.array, normalise:bool=True) -> float:
    """ 
    Calculates the root-mean-square (RMS) difference between two datasets. 
    
    The RMS value provides a measure of how closely the theoretical data matches 
    the experimental data. The datasets can optionally be normalised before calculating the RMS. 
    
    Args: 
        theory (np.array): Theoretical data to compare with the experiment. 
        experimental (np.array): Experimental data used for comparison. 
        normalise (bool): If True, normalises both datasets before calculating the RMS. Useful when comparing datasets with different amplitudes. 
    
    Returns: 
        float: RMS difference between the theoretical and experimental data. 
    """

    assert len(theory) == len(experimental), "Error: Input arrays must be the same length!"

    # Useful if amplitude are different between datasets
    if normalise:
        theory = normalised_data(theory)
        experimental = normalised_data(experimental)

    residuals = (theory- experimental) ** 2
    return np.sqrt(np.mean(residuals))
    
    

def solve_squares_min(exp_data: np.array,
                      m: np.array, # Complex refrative index of the material
                      radius_range : np.array = np.linspace(5,200,100)*1e-3,
                    ):
    """
    Finds the particle radius that best matches experimental spectral data. 
    
    Mie theory is calculated for each radius in the specified range and compared with the experimental scattering spectrum using the RMS error. The radius producing the lowest RMS error is returned as the best-fit particle radius. 
    
    Currently, only the scattering spectrum is used for the comparison. 
    
    Args: 
        exp_data (np.array): Experimental scattering data. 
        m (np.array): Complex refractive index of the particle material across the measured wavelength range. 
        radius_range (np.array): Range of particle radii to test. 
    
    Returns: 
        tuple: 
            float: Particle radius giving the lowest RMS error. 
            np.array: RMS error calculated for each tested radius.
    """

    # Checks is all data is same length
    assert len(exp_data) == len(m), "Error: Input arrays must be the same length!"

    # Pre allocated rms array
    rms_values= np.zeros_like(radius_range)

    # Loop over each particle size, calculate the mie theory spectra
    # then compare the measured spectra to what was theory. 
    # Radii that has the smalled RMS to the measured data is the best fitted particle.
    for i, radii in enumerate(radius_range):
        # NOTE: Could be made faster with parallel compute

        # Define the size paramter for the particle
        x = 2 * np.pi * radii / wavelengths

        # Calculate the Mie spectra for a given particle size
        qext, qsca, qback, g = mie.efficiencies_mx(m, x)

        # Calculate the RMS between theory and experimental data
        rms = calc_rms(qsca, exp_data)

        # Store the RMS values
        rms_values[i] = rms

    best_fit_radius = get_rsm_min(rms_values, radius_range)

    return best_fit_radius, rms_values


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

fitted_radii, rms = solve_squares_min(exp_scat, m, radius_range)

# -------------------
#  - Profiling Plot -
# -------------------

fig, axes = plt.subplots(3,1, figsize=(5,7))

# Subplot 1. Measured and fitted theory spectra
ax=axes[0]

ax.plot(wavelengths*1e3, norm_exp_scat, label='Measured Scat.')

# calc theory line
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
