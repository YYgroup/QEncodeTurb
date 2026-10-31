# Reproducibility package for the initial-velocity MPS study

Extract `mps.zip` into the directory containing `GeoEncodeTurb_main.ipynb`, preserving the following structure.

```text
project/
  GeoEncodeTurb_main.ipynb
  geoencode_turb/
  mps/
    run_all.py
    MPS_reproduce.ipynb
    requirements.txt
    src/
    vendor/
    original/
    reference/
```

Running the main notebook first is unnecessary. The package does not require a copy of the original workspace's `data/` directory.
All runtime paths follow the package location or the `--output` argument.

## Recompute the complete study

From the directory containing the main notebook, run these commands in your scientific Python environment.

```bash
python -m pip install -r mps/requirements.txt
python mps/run_all.py
```

Alternatively, open `mps/MPS_reproduce.ipynb`, select a Python kernel with the required dependencies, and run all cells.
The notebook defaults to the complete numerical study, including data generation, rather than only plotting existing CSV files.

The requirements target Python 3.11-3.13; the original study used Python 3.13.5.
The first complete run performs the following steps.

1. Verify the supplied source-code and reference-data checksums.
2. Compare the Python adapter with saved outputs from the original MATLAB routines.
3. Generate all three velocity components for 15, 18, 21, 24, and 27 spatial qubits.
4. Independently compute the Schmidt spectrum at every scale cut.
5. Construct MPS representations of each component at three error tolerances, save the actual cores, reconstruct the complete fields, and verify their errors.
6. Compile and transpile the 30-qubit circuits without simulating statevectors or velocity fields.
7. Recheck every spectrum using independent physical-coordinate unfoldings and direct SVD, and independently reconstruct every saved MPS.
8. Regenerate the 15- and 27-qubit velocity fields to verify generation consistency.
9. Export summaries, error lower bounds, resource tables, and all plotting inputs, and compare them with the supplied reference results.
10. Export the latest six-panel PDF/PNG, plotting notebook, caption, and two separate LaTeX tables.

All outputs default to `mps/generated/`. Completed velocity-field and MPS calculations with matching provenance are reused, supporting interrupted runs.
To start a completely fresh calculation, specify a new output directory; deleting previous results is unnecessary.

```bash
python mps/run_all.py --output mps/generated_fresh
```

## Check the environment or redraw the figure

Recompute the complete numerical workflow for 15 and 18 qubits and compare it with the reference results.

```bash
python mps/run_all.py --check --output mps/check_run
```

This command generates the small velocity fields and 18 MPS representations, performs independent checks, and compiles the 30-qubit circuits.
It does not draw the final six-panel figure, which requires all five grids.

Verify only the code and supplied reference data, without large calculations.

```bash
python mps/run_all.py --verify-only
```

Redraw the latest six-panel figure using only the supplied reference data.

```bash
python mps/run_all.py --plot-only
```

`--plot-only` does not generate new velocity fields or MPS cores; use the default command without this argument for complete recomputation.
To study a subset of grids, use, for example, `--n 15 18 21`. The complete six-panel figure requires all five grids.

## Computing resources and plotting fonts

- The complete reference velocity arrays occupy approximately 3.43 GiB, and the final 45 MPS representations occupy approximately 13.56 GiB in core arrays.
- For complete recomputation, allow at least 25 GiB of free disk space and 64 GiB of RAM. In the original environment, the peak process RSS during 27-qubit velocity export was approximately 26 GiB.
- MPS components and grids run sequentially to avoid concurrent large constructions.
- `--threads 8` controls the MPS BLAS thread count; the generation module retains its original setting of 16 Aer threads.
- The default `--tex auto` uses the manuscript's sansmathfonts style when LaTeX/dvipng are available, and Matplotlib math fonts otherwise.
- Use `--tex required` for the original fonts. This requires LaTeX, dvipng, bm, amsmath, and sansmathfonts.
- Use `--tex off` without a TeX environment. Font selection does not affect numerical data.

Timing and memory measurements depend on the machine and need not match the reference run. Dependency versions are pinned to keep circuit-transpilation statistics and numerical results comparable.

## Data and directory layout

| Path | Contents |
|---|---|
| `src/` | Portable field generation, MPS construction, independent auditing, aggregation, and plotting-data preparation |
| `vendor/geoencode_turb/` | Necessary generating-code snapshot from the reference study; neighboring project modules, which might have changed, are not imported by default |
| `original/` | Unmodified MATLAB/C sources and author notes from the public Gourianov 2022 archive |
| `matlab/` | Optional original MATLAB/Octave validation entry points and NPY reader |
| `templates/` | Latest six-panel plotting script and notebook template |
| `reference/` | All spectra, records, summaries, metadata, and validation data from the completed study, plus complete 15-qubit velocity arrays and MPS cores as small raw-data examples |
| `generated/fields/n*/` | Recomputed complete velocity arrays and generation metadata |
| `generated/mps_uniform/n*/u*/` | Per-component spectra, actual MPS cores, errors, and rank records |
| `generated/recheck/` | Independent direct-SVD, physical-coordinate reconstruction, and original joint-error-convention results |
| `generated/results/` | Summary, 45 component/tolerance records, scale-cut table, circuit resources, and reference comparisons |
| `generated/figure/` | Six-panel figure, plotting data, notebook, and two separate LaTeX tables |
| `generated/logs/` | Stage-by-stage runtime logs |

This is a recomputable reproducibility package. To avoid placing approximately 17 GiB of large arrays in the ZIP, complete 18-27-qubit velocity arrays and MPS cores are generated at runtime; their complete parameters, spectra, and validated results are supplied.
Early exploratory calculations with equal local error budgets are not part of the final study results. The default workflow reproduces the final fixed-cap calculation with global error control.

## Scientific settings

Each grid uses one realization, with phase seeds 2025 and 2026 for the two spinor components.
The two seeds belong to the same flow realization, not two separate realizations.

- Spatial registers contain 15, 18, 21, 24, and 27 qubits, corresponding to 32³, 64³, 128³, 256³, and 512³ grids.
- `k_cutoff=15`, amplitude locality radius 3, and sparse probability-tail tolerance 1e-12.
- The first-component populations are 0.5, 0.5, 0.08, 0.02, and 0.006, respectively.
- The three real velocity components are compressed separately, with componentwise relative L2 error tolerances of 1e-2, 1e-3, and 1e-4.
- The encoding follows the coarse-to-fine octal scale order of Gourianov 2022, with physical dimension 8 at each site.
- Necessary Schmidt rank bounds and verified construction ranks are distinguished; a single-cut error lower bound is not a measured full-MPS reconstruction error.
- Cores use float64. Stored core scalars and parameter counts after removing gauge freedom are reported separately.
- The 30-qubit calculation compiles circuit parameters and transpiled resources only; it does not generate a 1024³ velocity field or MPS.

Full definitions appear in `study_config.json` and `definitions_for_SI.tex`.
The package studies static representations of initial velocity fields, without fluid time evolution.

## Optional: rerun the original MATLAB/Octave routines

The main Python workflow does not require MATLAB, Octave, or the non-open-source TNT library.
Its cross-check uses supplied inputs and outputs previously obtained by executing the original MATLAB sources.
To regenerate those original-routine validation outputs, run either of the following commands.

```bash
python mps/src/run_matlab_validation.py --octave octave-cli
```

or

```bash
python mps/src/run_matlab_validation.py --matlab matlab
```

If the executable is absent from PATH, replace the last argument with its actual path.
Alternatively, add `mps/matlab/` to the MATLAB/Octave path and execute `validate_original` and `validate_joint_original`.
`run_original_velocity_mps(15,'ux',64)` directly reads an exported array and calls the original constructor.
The authors' pointwise indexing loop remains unchanged; use the validated Python adapter for large grids.

## Sources and validation scope

Paper: [Gourianov et al., Nature Computational Science 2, 30-37 (2022)](https://doi.org/10.1038/s43588-021-00181-1).
Original code archive: [doi:10.25592/uhhfdm.14234](https://doi.org/10.25592/uhhfdm.14234).
Author notes and source checksums are retained in `original/`, `original_provenance.json`, and `source_manifest.json`.
The package includes neither the original archive's large DNS dataset nor the TNT runtime; this study uses velocity fields generated by the present project.

`PORTABILITY_TESTS.md` records the checks actually performed during packaging, distinguishing small-grid recomputation from previous large-grid calculations.
