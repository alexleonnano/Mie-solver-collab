#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the automated fitting pipeline.

Run from the project root with either:
    python tests/test_pipeline.py
    pytest tests/test_pipeline.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import fitting, io_utils, mie_model
from src import lambda_interpolate as interp
import run_fit

SI_FILE = ROOT / "materials" / "Silicon - crystalline - Green 2008.txt"
EXP_FILE = ROOT / "experimental_data" / "SiNP_crystaline_scat.csv"


def _si_index(wavelengths):
    wl, n, k = io_utils.load_refractive_index(SI_FILE)
    n, k = interp.wvl_interpolation(wavelengths, wl, n, k)
    return n - 1.0j * k


def test_synthetic_round_trip():
    """A noisy theoretical spectrum at a known radius is fitted back to that radius."""
    rng = np.random.default_rng(0)
    wavelengths = np.arange(410, 690, 1) * 1e-3
    m = _si_index(wavelengths)
    radius_grid = np.arange(5, 201, 5) * 1e-3  # deliberately coarse, refinement must recover the rest

    for true_r in (0.0453, 0.0712, 0.0987):
        clean = mie_model.spectrum(m, wavelengths, true_r, n_medium=1.333)
        noisy = 3.0 * clean + 0.2 + rng.normal(0, 0.02 * clean.max(), clean.size)
        for norm in fitting.NORMALISATIONS:
            result = fitting.fit_radius(wavelengths, noisy, m, radius_grid, n_medium=1.333, normalisation=norm)
            assert abs(result.radius - true_r) < 1e-3, (norm, true_r, result.radius)


def test_matches_original_solver():
    """In vacuum with the original grid, the refactor reproduces squares_min_solver.py (70 nm, RMS 0.1685)."""
    exp_wl, exp_signal = io_utils.load_spectrum(EXP_FILE)
    wavelengths = np.arange(exp_wl.min() * 1e3, exp_wl.max() * 1e3, 1) * 1e-3
    signal = np.interp(wavelengths, exp_wl, exp_signal)
    result = fitting.fit_radius(wavelengths, signal, _si_index(wavelengths),
                                np.linspace(5, 200, 100) * 1e-3, n_medium=1.0, refine=False)
    assert abs(result.radius * 1e3 - 70.0) < 1e-9
    assert abs(result.rms - 0.16854) < 1e-4


def test_medium_scaling():
    """A sphere in a medium equals a sphere in vacuum with m/n_med at wavelength lambda/n_med."""
    wavelengths = np.linspace(0.4, 0.7, 50)
    m = _si_index(wavelengths)
    in_water = mie_model.spectrum(m, wavelengths, 0.07, n_medium=1.333)
    equivalent = mie_model.spectrum(m / 1.333, wavelengths / 1.333, 0.07, n_medium=1.0)
    assert np.allclose(in_water, equivalent)


def test_loader_formats():
    """Spectra load the same from index+CSV, tab-separated nm and whitespace micron files."""
    wl_nm = np.linspace(400, 700, 31)
    sig = np.sin(wl_nm / 50)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "a.csv").write_text(",wavelength (nm),signal\n" +
                                   "".join(f"{i},{w},{s}\n" for i, (w, s) in enumerate(zip(wl_nm, sig))))
        (tmp / "b.txt").write_text("wl\tabs\n" + "".join(f"{w}\t{s}\n" for w, s in zip(wl_nm[::-1], sig[::-1])))
        (tmp / "c.dat").write_text("# comment line\nlambda  I\n" + "".join(f"{w/1e3}  {s}\n" for w, s in zip(wl_nm, sig)))
        for name in ("a.csv", "b.txt", "c.dat"):
            wl, s = io_utils.load_spectrum(tmp / name)
            assert np.allclose(wl, wl_nm * 1e-3), name
            assert np.allclose(s, sig), name


def test_cli_writes_txt_results():
    with tempfile.TemporaryDirectory() as tmp:
        code = run_fit.main([str(EXP_FILE), "--config", str(ROOT / "config.toml"),
                             "--results-dir", tmp, "--no-plot", "--radius-range", "50", "100", "2"])
        assert code == 0
        (out_dir,) = Path(tmp).iterdir()
        report = (out_dir / "fit_report.txt").read_text()
        assert "Best-fit radius:    70." in report
        spectrum = np.loadtxt(out_dir / "fit_spectrum.txt")
        assert spectrum.shape[1] == 4
        assert np.loadtxt(out_dir / "rms_curve.txt").shape == (26, 2)


if __name__ == "__main__":
    tests = [(name, f) for name, f in globals().items() if name.startswith("test_")]
    failed = 0
    for name, f in tests:
        try:
            f()
            print(f"PASS  {name}")
        except Exception as e:
            failed += 1
            print(f"FAIL  {name}: {e!r}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
