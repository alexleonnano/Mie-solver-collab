#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mie theory spectra for spheres embedded in a non-absorbing medium.

Wavelengths and radii are in microns. `m_particle` follows the miepython
convention of a negative imaginary part (n - 1j*k).
"""

import numpy as np
import miepython as mie

SPECTRUM_TYPES = ("sca", "ext", "abs")


def efficiencies(m_particle, wavelengths, radius, n_medium=1.0):
    """
    Calculates the extinction, scattering and absorption efficiencies of a sphere.

    The relative refractive index is m_particle / n_medium and the size parameter
    uses the wavelength inside the medium, 2*pi*r*n_medium / wavelength.

    Args:
        m_particle (np.array): Complex refractive index of the particle at each wavelength.
        wavelengths (np.array): Vacuum wavelengths (microns).
        radius (float or np.array): Particle radius (microns), broadcastable against wavelengths.
        n_medium (float): Refractive index of the surrounding medium.

    Returns:
        dict: {"ext": qext, "sca": qsca, "abs": qabs}, each shaped like the broadcast inputs.
    """
    m_rel, x = np.broadcast_arrays(np.asarray(m_particle) / n_medium,
                                   2 * np.pi * np.asarray(radius) * n_medium / np.asarray(wavelengths))
    # miepython only accepts 1D arrays, so flatten and restore the shape afterwards
    qext, qsca, _, _ = mie.efficiencies_mx(m_rel.ravel(), x.ravel())
    qext, qsca = qext.reshape(x.shape), qsca.reshape(x.shape)
    return {"ext": qext, "sca": qsca, "abs": qext - qsca}


def spectra_grid(m_particle, wavelengths, radii, n_medium=1.0, kind="sca"):
    """
    Calculates a spectrum for every radius in one call.

    Returns:
        np.array: Array of shape (len(radii), len(wavelengths)).
    """
    if kind not in SPECTRUM_TYPES:
        raise ValueError(f"Unknown spectrum type '{kind}', expected one of {SPECTRUM_TYPES}.")
    return efficiencies(m_particle, wavelengths[np.newaxis, :], np.asarray(radii)[:, np.newaxis], n_medium)[kind]


def spectrum(m_particle, wavelengths, radius, n_medium=1.0, kind="sca"):
    """Calculates the spectrum for a single radius."""
    return spectra_grid(m_particle, wavelengths, [radius], n_medium, kind)[0]
