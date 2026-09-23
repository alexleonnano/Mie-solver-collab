# Mie solver implementation in python
Nanoparticle sizing from UV-Vis spectra of colloids. A measured spectrum is compared against Mie theory spectra (calculated with the `miepython` package and the refractive index data in `materials/`) over a range of particle radii, and the radius with the lowest RMS difference is reported as the particle size.

## Installation
To run locally just download the github repo, open the terminal in the project folder and depending on what you use for package management, create an environment with the dependencies and run the files.

### Using Pip
Open the terminal at the project folder location and create a new environment:

`python -m venv <environment_name>`

Then activate the environment and install the packages running the following command:

`pip install -r requirements.txt`

### Using conda
Open the conda terminal in the project folder and build the environment using the yml file:

`conda env create -f environment.yml`

After that you can activate the environment (it is called "scat"):

`conda activate scat`

## Usage
1. Measure the UV-Vis spectrum and save it as a text/CSV file (wavelength column + signal column; nm or µm, comma/tab/space separated, header optional).
2. Check the settings in `config.toml` (material, medium index, spectrum type, radius range, ...).
3. Run:

```
python run_fit.py path/to/spectrum.csv
```

Other examples:

```
python run_fit.py                              # uses [input] spectrum from config.toml
python run_fit.py experimental_data/           # batch: fits every spectrum in the folder
python run_fit.py sample.csv --medium-index 1.0 --radius-range 20 150 0.5 --show
python run_fit.py -h                           # all options
```

Command line options override `config.toml`.

Mie calculations use miepython's numba backend, which is about 50x faster with identical results. To disable it, set the environment variable `MIEPYTHON_USE_JIT=0`.

### Output
Each spectrum gets a folder `results/<spectrum name>_<timestamp>/` with:
- `fit_report.txt` - best-fit radius and diameter, RMS, other local minima of the RMS curve, warnings and all settings used.
- `fit_spectrum.txt` - wavelength, normalised measured and theory spectra and residuals (tab separated).
- `rms_curve.txt` - RMS for every radius tested.
- `fit_plot.png` - measured vs fitted spectrum, residuals and RMS curve.

### Batch processing
To size several samples at once, put their spectra in one folder, naming each file after its sample (e.g. `SiNP_A.csv`, `SiNP_B.csv`). Then run:

```
python run_fit.py path/to/folder/
python run_fit.py path/to/folder/ --pattern "SiNP_*.csv"   # only the matching files
```

All samples are fitted with the same settings and material. Files ending in `.csv`, `.txt`, `.tsv` or `.dat` are included and hidden files are skipped. If a sample fails (e.g. an unreadable file), it is recorded and the batch carries on. Results go to one folder per batch:

```
results/batch_<folder>_<timestamp>/
    batch_summary.txt     # settings, then one row per sample: sample, radius_nm, diameter_nm, rms, status, spectrum_file, notes
    batch_radii.png       # fitted radius of every sample
    SiNP_A/               # fit_report.txt, fit_spectrum.txt, rms_curve.txt, fit_plot.png
    SiNP_B/
```

`status` is `ok`, `warning` (see the notes column and that sample's `fit_report.txt`) or `FAILED`. If two files have the same name with different extensions, the extension is added to the sample name (`S2_csv`, `S2_txt`). Giving several files on the command line also runs them as a batch.

### Settings (`config.toml`)
| Setting | Meaning |
| --- | --- |
| `material.file` | Refractive index file (wavelength, n, k) |
| `material.medium_index` | Refractive index of the surrounding medium (water = 1.333) |
| `fit.spectrum_type` | `sca`, `ext` or `abs`. Si particles around 70 nm are scattering dominated, so `sca` is used |
| `fit.radius_*_nm` | Radius grid; the best grid point is refined between its neighbours when `refine = true` |
| `fit.wl_min_nm`, `fit.wl_max_nm` | Optional wavelength window, e.g. to cut noisy edges |
| `fit.normalisation` | `minmax` (both spectra scaled to [0,1]) or `affine` (theory fitted as a·theory + b, tolerant of baseline offsets) |

## Project Structure
- `run_fit.py` - command line entry point of the automated pipeline.
- `config.toml` - default settings.
- `src/io_utils.py` - loaders for measured spectra and refractive index files.
- `src/mie_model.py` - Mie spectra of spheres in a medium, for one or many radii.
- `src/fitting.py` - RMS scoring and radius fitting.
- `src/plotting.py` - diagnostic plot.
- `legacy/squares_min_solver.py` - original single-script version of the fit (vacuum, fixed grid), kept for reference. Run it from the project root: `python legacy/squares_min_solver.py`.
- `main.py` - simple example of the miepython package with gold nanoparticles.
- `tests/test_pipeline.py` - tests for the pipeline, run with `python tests/test_pipeline.py` (or `pytest`).
- `tests/` also holds the comparison between miepython and the Mathematica reference spectra, and the notebooks in the project root hold earlier interpolation/validation work.
