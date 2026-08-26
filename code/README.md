# Geometric encoding of turbulence for end-to-end quantum simulation

## Files

- `GeoEncodeTurb_main.ipynb` is the main notebook and contains concise workflow calls.
- `geoencode_turb/config.py` stores the shared workflow parameters.
- `geoencode_turb/quantum_encoding.py` contains target-distribution generation, Gray-code mapping, conditional-angle fitting, circuit construction, and statevector simulation.
- `geoencode_turb/flow_reconstruction.py` reconstructs the Pauli-spinor and primary fluid fields.
- `geoencode_turb/flow_statistics.py` contains spectra, PDFs, vorticity, velocity-gradient invariants, structure functions, and anisotropy.
- `geoencode_turb/visualization.py` contains all Matplotlib and PyVista figures.
- `geoencode_turb/affine_fit_analysis.py` contains the Gray-versus-binary conditional-angle analysis.
- `geoencode_turb/temporal_correlation.py` contains time correlation analysis.

## Use

1. Open `GeoEncodeTurb_main.ipynb` from this directory.
2. Run the setup cell.
3. Run the state-preparation and flow-reconstruction cells in order.
4. Run only the post-processing and visualization cells needed for the current analysis. The structure-function and PyVista calls remain separate because they require substantial time or memory.

All output paths remain relative to this directory. The workflow creates the original `data/3D` and `figures` directories when needed.

## Hardware requirements

The 30-qubit case corresponds to a \(1024^3=2^{30}\)-point grid. In double precision, a 30-qubit complex statevector occupies approximately 16 GiB, while each real-valued three-dimensional field occupies approximately 8 GiB. Since the main encoding, reconstruction, and post-processing stages are executed sequentially and most temporary arrays are released after use, the complete workflow runs on a workstation with 128 GB of RAM. We therefore regard 128 GB as a practical memory requirement for the current implementation. A system with 192--256 GB of RAM is recommended when additional diagnostics are enabled or multiple full-resolution intermediate fields are retained simultaneously. Approximately 16--32 CPU cores and 200 GB or more of free storage are also recommended.
