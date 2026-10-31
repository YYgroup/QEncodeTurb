# Validation performed during packaging

## Executed checks

1. **Fresh numerical calculations.** Generated the 15- and 18-qubit velocity fields in an empty output directory, constructed 18 MPS representations, independently verified all spectra and saved-core reconstructions, and regenerated the 15-qubit fields. All six complete velocity arrays have the same SHA-256 hashes as the reference arrays.
2. **Resource comparison.** Fresh small-grid ranks, storage, gauge-reduced parameter counts, and measured errors match the reference results. Recompiling the 30-qubit instance gives 982 description parameters, 1868 CX gates, 1070 U gates, and maximum depth 194.
3. **Full-grid output stages.** Using the previously completed and independently audited 15-27-qubit records, executed this package's aggregation, plotting-data preparation, and reference comparison. Output checks passed for all 45 combinations, 90 spectra, two tables, and the six-panel figure.
4. **Original-routine cross-checks.** Reran the included original MATLAB construction, analysis, and reconstruction routines using GNU Octave 11.3.0. The Python adapter agrees with the original outputs. The maximum cut ranks for the 64³ joint-component error criterion are 475, 504, and 511.
5. **Path checks.** Runtime sources contain no hard-coded paths to the original workspace or personal directories. Historical paths in archived metadata are not runtime dependencies.
6. **Extracted-ZIP execution.** Extracted the ZIP into a new example directory with Chinese characters and spaces, beside the main notebook. Executed the entry notebook in `check` mode, redrew the six-panel figure using the required LaTeX fonts, and verified file integrity. Saving entry-notebook execution results does not invalidate the checks.

Packaging did not repeat all large 21-27-qubit MPS constructions. Those calculations and their independent validation were completed during the original study; their full records remain in `reference/`. The default complete command regenerates the large arrays and MPS cores.

The extracted-package execution results appear in `tests/relocation_test.json`. Numerical results for the small-grid and original-routine checks are stored in the JSON files under `tests/`.
