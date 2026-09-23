#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated nanoparticle sizing from a measured UV-Vis spectrum.

Loads a measured spectrum, calculates Mie theory spectra over a range of radii
and reports the radius with the lowest RMS difference. Settings come from
`config.toml`; command line options override them.

Usage:
    python run_fit.py                               # spectrum from config.toml
    python run_fit.py experimental_data/sample.csv  # one spectrum
    python run_fit.py experimental_data/            # every spectrum in a folder
    python run_fit.py sample.csv --medium-index 1.0 --show

Results are written to results/<spectrum name>_<timestamp>/:
    fit_report.txt    - fitted radius, RMS and the settings used
    fit_spectrum.txt  - wavelength, measured, theory and residual (normalised)
    rms_curve.txt     - RMS for every radius tested
    fit_plot.png      - diagnostic plot
A batch run also writes results/batch_summary_<timestamp>.txt.
"""

import argparse
from datetime import datetime
from pathlib import Path
import sys
import tomllib
import warnings

import numpy as np

from src import fitting, io_utils
from src import lambda_interpolate as interp

SPECTRUM_SUFFIXES = {".csv", ".txt", ".tsv", ".dat"}
MIN_MEASURED_POINTS = 20


def load_config(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Fit nanoparticle radius from a measured spectrum using Mie theory.")
    parser.add_argument("spectra", nargs="*", help="Spectrum file(s) or folder(s). Defaults to [input] spectrum in the config.")
    parser.add_argument("-c", "--config", default="config.toml", help="Config file (default: config.toml).")
    parser.add_argument("--material", help="Refractive index file, or its name inside materials/.")
    parser.add_argument("--medium-index", type=float, help="Refractive index of the surrounding medium.")
    parser.add_argument("--type", choices=["sca", "ext", "abs"], help="Spectrum type to fit.")
    parser.add_argument("--normalisation", choices=list(fitting.NORMALISATIONS))
    parser.add_argument("--radius-range", nargs=3, type=float, metavar=("MIN", "MAX", "STEP"), help="Radii to test (nm).")
    parser.add_argument("--wl-range", nargs=2, type=float, metavar=("MIN", "MAX"), help="Wavelength window to fit (nm).")
    parser.add_argument("--results-dir", help="Output folder.")
    parser.add_argument("--show", action="store_true", help="Show the plot window.")
    parser.add_argument("--no-plot", action="store_true", help="Do not save a plot.")
    return parser.parse_args(argv)


def merge_settings(cfg, args):
    """Flattens the config and applies command line overrides."""
    fit = cfg.get("fit", {})
    out = cfg.get("output", {})
    s = {
        "wl_column": cfg.get("input", {}).get("wl_column"),
        "signal_column": cfg.get("input", {}).get("signal_column"),
        "material": cfg.get("material", {}).get("file"),
        "medium_index": cfg.get("material", {}).get("medium_index", 1.333),
        "kind": fit.get("spectrum_type", "sca"),
        "radius_range": (fit.get("radius_min_nm", 5), fit.get("radius_max_nm", 200), fit.get("radius_step_nm", 1)),
        "refine": fit.get("refine", True),
        "wl_range": (fit.get("wl_min_nm"), fit.get("wl_max_nm")),
        "wl_step": fit.get("wl_step_nm", 1),
        "normalisation": fit.get("normalisation", "minmax"),
        "results_dir": out.get("results_dir", "results"),
        "save_plot": out.get("save_plot", True),
        "show_plot": out.get("show_plot", False),
    }
    overrides = {
        "material": args.material, "medium_index": args.medium_index, "kind": args.type,
        "normalisation": args.normalisation, "radius_range": args.radius_range,
        "wl_range": args.wl_range, "results_dir": args.results_dir,
    }
    s.update({k: v for k, v in overrides.items() if v is not None})
    if args.show:
        s["show_plot"] = True
    if args.no_plot:
        s["save_plot"] = False
    if s["material"] is None:
        raise SystemExit("No material given: set [material] file in the config or use --material.")
    return s


def collect_spectra(paths):
    """Expands folders into the spectrum files they contain."""
    files = []
    for p in map(Path, paths):
        if p.is_dir():
            files += sorted(f for f in p.iterdir() if f.suffix.lower() in SPECTRUM_SUFFIXES)
        elif p.exists():
            files.append(p)
        else:
            raise SystemExit(f"Spectrum not found: {p}")
    if not files:
        raise SystemExit("No spectrum files found.")
    return files


def fit_window(exp_wl, mat_wl, wl_range, notes):
    """Wavelength window (microns) covered by both the measurement and the refractive index data."""
    lo, hi = exp_wl.min(), exp_wl.max()
    if mat_wl.min() > lo or mat_wl.max() < hi:
        notes.append(f"Refractive index data ({mat_wl.min()*1e3:.0f}-{mat_wl.max()*1e3:.0f} nm) does not cover "
                     f"the whole measurement ({lo*1e3:.0f}-{hi*1e3:.0f} nm); fit window cropped.")
        lo, hi = max(lo, mat_wl.min()), min(hi, mat_wl.max())
    user_lo, user_hi = wl_range
    if user_lo is not None:
        lo = max(lo, user_lo * 1e-3)
    if user_hi is not None:
        hi = min(hi, user_hi * 1e-3)
    if hi <= lo:
        raise ValueError("Empty wavelength window: measurement, refractive index data and wl_range do not overlap.")
    return lo, hi


def fit_spectrum(spectrum_path, s):
    """Runs the full pipeline for one spectrum. Returns (result, arrays, notes)."""
    notes = []
    exp_wl, exp_signal = io_utils.load_spectrum(spectrum_path, s["wl_column"], s["signal_column"])
    mat_file = io_utils.material_path(s["material"])
    mat_wl, mat_n, mat_k = io_utils.load_refractive_index(mat_file)

    lo, hi = fit_window(exp_wl, mat_wl, s["wl_range"], notes)
    n_measured = np.count_nonzero((exp_wl >= lo) & (exp_wl <= hi))
    if n_measured < MIN_MEASURED_POINTS:
        raise ValueError(f"Only {n_measured} measured points inside the fit window "
                         f"({lo*1e3:.0f}-{hi*1e3:.0f} nm), need at least {MIN_MEASURED_POINTS}.")
    wavelengths = np.arange(lo * 1e3, hi * 1e3, s["wl_step"]) * 1e-3

    signal = np.interp(wavelengths, exp_wl, exp_signal)
    n, k = interp.wvl_interpolation(wavelengths, mat_wl, mat_n, mat_k)
    m = n - 1.0j * k  # miepython uses a negative imaginary part

    r_min, r_max, r_step = s["radius_range"]
    radius_grid = np.arange(r_min, r_max + r_step / 2, r_step) * 1e-3

    result = fitting.fit_radius(wavelengths, signal, m, radius_grid, s["medium_index"],
                                s["kind"], s["normalisation"], s["refine"])
    if result.at_grid_edge:
        notes.append("Best radius is at the edge of the radius range - widen the range and re-run.")

    exp_norm = fitting.normalised_data(signal)
    if s["normalisation"] == "affine":
        a, b = np.polyfit(result.theory, signal, 1)
        theory_norm = (a * result.theory + b - signal.min()) / np.ptp(signal)
    else:
        theory_norm = fitting.normalised_data(result.theory)

    arrays = {"wavelengths": wavelengths, "exp_norm": exp_norm, "theory_norm": theory_norm}
    return result, arrays, mat_file, notes


def write_outputs(out_dir, spectrum_path, mat_file, s, result, arrays, notes):
    out_dir.mkdir(parents=True, exist_ok=True)
    wl = arrays["wavelengths"]

    minima = fitting.local_minima(result.rms_grid)
    others = [i for i in minima if abs(result.radius_grid[i] - result.radius) > 2 * (result.radius_grid[1] - result.radius_grid[0])][:3]

    lines = [
        "Mie theory nanoparticle radius fit",
        "=" * 40,
        f"Date:               {datetime.now():%Y-%m-%d %H:%M:%S}",
        f"Spectrum file:      {spectrum_path}",
        f"Material file:      {mat_file}",
        f"Medium index:       {s['medium_index']}",
        f"Spectrum type:      {s['kind']}",
        f"Normalisation:      {s['normalisation']}",
        f"Wavelength window:  {wl[0]*1e3:.1f} - {wl[-1]*1e3:.1f} nm (step {s['wl_step']} nm, {len(wl)} points)",
        f"Radius grid:        {result.radius_grid[0]*1e3:.1f} - {result.radius_grid[-1]*1e3:.1f} nm "
        f"(step {s['radius_range'][2]} nm, refine={s['refine']})",
        "",
        "Result",
        "-" * 40,
        f"Best-fit radius:    {result.radius*1e3:.2f} nm",
        f"Best-fit diameter:  {2*result.radius*1e3:.2f} nm",
        f"RMS:                {result.rms:.4g}",
    ]
    if others:
        lines.append("Other local minima (radius nm, RMS):")
        lines += [f"    {result.radius_grid[i]*1e3:8.2f}   {result.rms_grid[i]:.4g}" for i in others]
    if notes:
        lines += ["", "Warnings"] + [f"  - {n}" for n in notes]
    (out_dir / "fit_report.txt").write_text("\n".join(lines) + "\n")

    np.savetxt(out_dir / "fit_spectrum.txt",
               np.column_stack([wl * 1e3, arrays["exp_norm"], arrays["theory_norm"],
                                arrays["theory_norm"] - arrays["exp_norm"]]),
               fmt="%.6g", delimiter="\t",
               header="wavelength_nm\tmeasured_norm\ttheory_norm\tresidual")
    np.savetxt(out_dir / "rms_curve.txt", np.column_stack([result.radius_grid * 1e3, result.rms_grid]),
               fmt="%.6g", delimiter="\t", header="radius_nm\trms")
    return "\n".join(lines)


def main(argv=None):
    args = parse_args(argv)
    cfg = load_config(args.config)
    s = merge_settings(cfg, args)

    spectra = collect_spectra(args.spectra or [cfg.get("input", {}).get("spectrum")])
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    results_dir = Path(s["results_dir"])
    summary = []

    if s["save_plot"] or s["show_plot"]:
        import matplotlib
        if not s["show_plot"]:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from src import plotting

    for path in spectra:
        print(f"Fitting {path} ...")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                result, arrays, mat_file, notes = fit_spectrum(path, s)
            except Exception as e:
                print(f"  FAILED: {e}", file=sys.stderr)
                summary.append((path, None, None, f"FAILED: {e}"))
                continue
        notes += [str(w.message) for w in caught]

        out_dir = results_dir / f"{path.stem}_{stamp}"
        report = write_outputs(out_dir, path, mat_file, s, result, arrays, notes)
        if s["save_plot"] or s["show_plot"]:
            fig = plotting.plot_fit(arrays["wavelengths"], arrays["exp_norm"], arrays["theory_norm"],
                                    result, s["kind"], title=path.name)
            if s["save_plot"]:
                fig.savefig(out_dir / "fit_plot.png", dpi=150)
        print(report)
        print(f"Results saved to {out_dir}/\n")
        summary.append((path, result.radius, result.rms, "; ".join(notes)))

    if len(spectra) > 1:
        results_dir.mkdir(parents=True, exist_ok=True)
        summary_path = results_dir / f"batch_summary_{stamp}.txt"
        rows = ["spectrum\tradius_nm\tdiameter_nm\trms\tnotes"]
        for path, r, rms, note in summary:
            rows.append(f"{path}\t{r*1e3:.2f}\t{2*r*1e3:.2f}\t{rms:.4g}\t{note}" if r is not None
                        else f"{path}\t\t\t\t{note}")
        summary_path.write_text("\n".join(rows) + "\n")
        print(f"Batch summary saved to {summary_path}")

    if s["show_plot"]:
        plt.show()
    return 0 if all(r is not None for _, r, _, _ in summary) else 1


if __name__ == "__main__":
    sys.exit(main())
