#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnostic plot for a radius fit: spectra, residuals and RMS curve."""

import matplotlib.pyplot as plt

SPECTRUM_LABELS = {"sca": "Scat.", "ext": "Ext.", "abs": "Abs."}


def plot_fit(wavelengths, exp_norm, theory_norm, result, kind="sca", title=None):
    """
    Args:
        wavelengths (np.array): Wavelengths (microns).
        exp_norm (np.array): Normalised measured spectrum.
        theory_norm (np.array): Normalised best-fit theory spectrum.
        result (FitResult): Output of `fitting.fit_radius`.
        kind (str): Spectrum type, used for labels.
        title (str, optional): Figure title.

    Returns:
        matplotlib.figure.Figure
    """
    wl_nm = wavelengths * 1e3
    label = SPECTRUM_LABELS.get(kind, kind)
    fig, axes = plt.subplots(3, 1, figsize=(6, 8))

    # Subplot 1. Measured and fitted theory spectra
    ax = axes[0]
    ax.plot(wl_nm, exp_norm, label=f"Measured {label}")
    ax.plot(wl_nm, theory_norm, color="firebrick", label=f"Mie {label} - r={result.radius*1e3:.2f} nm")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.legend()

    # Subplot 2. Residuals profile for fitted particle radius
    ax = axes[1]
    ax.plot(wl_nm, theory_norm - exp_norm)
    ax.axhline(0, color="red", linestyle="--", alpha=0.5)
    ax.set_title("Residuals")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(r"$theory - experiment$")

    # Subplot 3. RMS score for all tested particle sizes
    ax = axes[2]
    ax.plot(result.radius_grid * 1e3, result.rms_grid)
    ax.axvline(result.radius * 1e3, color="firebrick", linestyle="--", alpha=0.7)
    ax.set_title("RMS")
    ax.set_xlabel("Particle Radius (nm)")
    ax.set_ylabel("RMS")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_batch_radii(names, radii_nm):
    """Overview of a batch: fitted radius of every sample."""
    fig, ax = plt.subplots(figsize=(max(5, 0.5 * len(names) + 2), 4))
    ax.plot(range(len(names)), radii_nm, "o", color="firebrick")
    ax.set_xticks(range(len(names)), names, rotation=45, ha="right")
    ax.set_ylabel("Fitted radius (nm)")
    ax.set_title("Batch overview")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    return fig
