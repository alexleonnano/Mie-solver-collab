#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Least-squares matching of measured spectra against Mie theory.

A grid of radii is scored by the RMS difference between the normalised theory
and measured spectra, then the best grid point is refined with a bounded 1D
minimisation.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from src import mie_model

NORMALISATIONS = ("minmax", "affine")


def normalised_data(data):
    """Normalise the data between [0,1] - min-max method. Works along the last axis."""
    data = np.asarray(data, dtype=float)
    lo = data.min(axis=-1, keepdims=True)
    hi = data.max(axis=-1, keepdims=True)
    return (data - lo) / (hi - lo)


def calc_rms(theory, experimental, normalisation="minmax"):
    """
    Calculates the root-mean-square (RMS) difference between theory and experiment.

    Args:
        theory (np.array): Theoretical spectrum, or a 2D array of spectra (one per row).
        experimental (np.array): Measured spectrum on the same wavelengths.
        normalisation (str):
            "minmax" - both spectra are min-max normalised to [0,1] before comparing.
            "affine" - the theory is least-squares fitted to the experiment as a*theory + b
                       (absorbs unknown concentration and baseline offset). The RMS is
                       divided by the experimental range so it is on the same scale as "minmax".

    Returns:
        float or np.array: RMS value(s), one per theory spectrum.
    """
    theory = np.asarray(theory, dtype=float)
    experimental = np.asarray(experimental, dtype=float)
    assert theory.shape[-1] == experimental.shape[-1], "Error: Input arrays must be the same length!"

    if normalisation == "minmax":
        residuals = normalised_data(theory) - normalised_data(experimental)
    elif normalisation == "affine":
        t_mean = theory.mean(axis=-1, keepdims=True)
        e_mean = experimental.mean()
        t_c = theory - t_mean
        a = (t_c * (experimental - e_mean)).sum(axis=-1, keepdims=True) / (t_c ** 2).sum(axis=-1, keepdims=True)
        residuals = (a * t_c + e_mean - experimental) / np.ptp(experimental)
    else:
        raise ValueError(f"Unknown normalisation '{normalisation}', expected one of {NORMALISATIONS}.")

    return np.sqrt(np.mean(residuals ** 2, axis=-1))


def local_minima(values):
    """Indices of the interior local minima of a 1D array, best (lowest) first."""
    v = np.asarray(values)
    idx = np.where((v[1:-1] < v[:-2]) & (v[1:-1] < v[2:]))[0] + 1
    return idx[np.argsort(v[idx])]


@dataclass
class FitResult:
    radius: float              # best-fit radius (microns), after refinement
    rms: float                 # RMS at the best-fit radius
    radius_grid: np.ndarray    # radii tested on the grid (microns)
    rms_grid: np.ndarray       # RMS for each grid radius
    theory: np.ndarray         # best-fit theoretical spectrum (un-normalised efficiency)
    at_grid_edge: bool         # True if the best grid point is the first/last radius tested


def fit_radius(wavelengths, signal, m_particle, radius_grid, n_medium=1.0,
               kind="sca", normalisation="minmax", refine=True):
    """
    Finds the particle radius whose Mie spectrum best matches the measured spectrum.

    Args:
        wavelengths (np.array): Vacuum wavelengths (microns).
        signal (np.array): Measured spectrum on `wavelengths`.
        m_particle (np.array): Complex refractive index of the particle on `wavelengths`.
        radius_grid (np.array): Radii to test (microns).
        n_medium (float): Refractive index of the surrounding medium.
        kind (str): Spectrum compared against the measurement: "sca", "ext" or "abs".
        normalisation (str): See `calc_rms`.
        refine (bool): Refine the best grid radius with a bounded minimisation between its neighbours.

    Returns:
        FitResult
    """
    assert len(wavelengths) == len(signal) == len(m_particle), "Error: Input arrays must be the same length!"
    radius_grid = np.asarray(radius_grid, dtype=float)

    grid = mie_model.spectra_grid(m_particle, wavelengths, radius_grid, n_medium, kind)
    rms_grid = calc_rms(grid, signal, normalisation)

    i = int(np.argmin(rms_grid))
    best_radius, best_rms = radius_grid[i], rms_grid[i]

    if refine and len(radius_grid) > 1:
        lo = radius_grid[max(i - 1, 0)]
        hi = radius_grid[min(i + 1, len(radius_grid) - 1)]

        def objective(r):
            return calc_rms(mie_model.spectrum(m_particle, wavelengths, r, n_medium, kind), signal, normalisation)

        res = minimize_scalar(objective, bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
        if res.success and res.fun < best_rms:
            best_radius, best_rms = float(res.x), float(res.fun)

    theory = mie_model.spectrum(m_particle, wavelengths, best_radius, n_medium, kind)
    return FitResult(best_radius, float(best_rms), radius_grid, rms_grid, theory,
                     at_grid_edge=i in (0, len(radius_grid) - 1))
