#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Robust loaders for experimental spectra and refractive index files.

All wavelengths are returned in microns, which is the unit used throughout
the Mie calculations.
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd

WAVELENGTH_HINTS = ("wav", "wl", "lambda", "nm")


def _read_table(path) -> pd.DataFrame:
    """Read a delimited text table (comma, tab, semicolon or whitespace)."""
    return pd.read_csv(path, sep=r"[,;\t ]+|\s+", engine="python", comment="#")


def _to_microns(wl: np.ndarray) -> np.ndarray:
    """Wavelengths above 10 are assumed to be in nm, otherwise already in microns."""
    return wl * 1e-3 if np.nanmax(wl) > 10 else wl


def _is_row_index(col: pd.Series) -> bool:
    """True for a pandas-style 0, 1, 2, ... index column left over from `to_csv`."""
    values = col.to_numpy()
    return np.array_equal(values, np.arange(len(values)))


def load_spectrum(path, wl_column=None, signal_column=None):
    """
    Loads a measured UV-Vis spectrum.

    Columns are detected automatically: a leftover row-index column is ignored,
    the wavelength column is the one whose name looks like a wavelength
    (falling back to the first numeric column) and the signal is the next one.

    Args:
        path: Path to the spectrum file (.csv, .txt, .tsv, ...).
        wl_column (str, optional): Name of the wavelength column, overrides auto-detection.
        signal_column (str, optional): Name of the signal column, overrides auto-detection.

    Returns:
        tuple:
            np.array: Wavelengths in microns, sorted ascending.
            np.array: Measured signal.
    """
    df = _read_table(path)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.loc[:, df.notna().any()]  # drop fully non-numeric columns
    df.columns = [str(c).strip() for c in df.columns]

    candidates = [c for c in df.columns if not _is_row_index(df[c])]

    if wl_column is None:
        named = [c for c in candidates if any(h in c.lower() for h in WAVELENGTH_HINTS)]
        wl_column = named[0] if named else candidates[0]
    if signal_column is None:
        signal_column = next(c for c in candidates if c != wl_column)

    data = df[[wl_column, signal_column]].dropna().sort_values(wl_column)
    if len(data) < len(df):
        warnings.warn(f"{path}: dropped {len(df) - len(data)} rows with missing values.")

    return _to_microns(data[wl_column].to_numpy()), data[signal_column].to_numpy()


def load_refractive_index(path):
    """
    Loads a refractive index file with three columns: wavelength, n, k.

    Args:
        path: Path to the material file (wavelength in microns or nm).

    Returns:
        tuple:
            np.array: Wavelengths in microns, sorted ascending.
            np.array: Real part of the refractive index (n).
            np.array: Imaginary part of the refractive index (k).
    """
    df = _read_table(path).apply(pd.to_numeric, errors="coerce")
    df = df.iloc[:, :3].dropna()
    df.columns = ["wl", "n", "k"]
    df = df.sort_values("wl")
    return _to_microns(df["wl"].to_numpy()), df["n"].to_numpy(), df["k"].to_numpy()


def material_path(material, materials_dir="materials") -> Path:
    """Resolves a material given as a path or as a file name (with or without extension) in `materials_dir`."""
    path = Path(material)
    if path.exists():
        return path
    matches = sorted(Path(materials_dir).glob(f"{material}*"))
    if not matches:
        raise FileNotFoundError(f"Material '{material}' not found (looked in '{materials_dir}/').")
    return matches[0]
