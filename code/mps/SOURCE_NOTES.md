# Source provenance

The numerical octal encoding and fixed-cap TT-SVD routines in `src/mps_ncs.py`
and `src/mps_uniform.py` are byte-for-byte copies of the audited implementation.
The packaging changes affect paths, orchestration, reporting, checkpoint checks,
and optional environment selection. They do not change the final MPS constructor.

The necessary generating modules in `vendor/geoencode_turb/` are frozen from the
implementation used to export the published reference velocity fields. The
configuration, quantum-encoding and flow-reconstruction hashes match the recorded
source hashes. The neighboring project's current Python modules are not imported.

The authors' files in `original/` are preserved unchanged from the public archive.
Their own README and the archive DOI remain included. No TNT library or original
DNS simulation dataset is redistributed in this package.

`reference/` contains the archived outcomes of the completed study. Timings,
process memory observations, source paths in old metadata, and historical hashes
describe that run. Those old absolute paths are provenance data, not runtime
dependencies. Newly generated results use the current output directory.

Default execution generates fresh numerical data in `generated/`; it does not
silently substitute supplied reference spectra for the computed spectra. The
explicit `--plot-only` mode uses the archived plotting data and says so in its
terminal output. Saved original MATLAB outputs are used for the Python/original
routine comparison; the optional MATLAB/Octave runner regenerates them.
