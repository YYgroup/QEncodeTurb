# Plotting data

All files describe one independently audited realization, with phase seeds 2025 and 2026.

| File | Contents |
|---|---|
| velocity_slices.npz | Middle slices along the third array axis for five grids and three velocity components; keys include n15_ux. Values retain their original scale without normalization or clipping. |
| velocity_metadata.csv | Grid size, slice location, full-field RMS, and full three-dimensional field checksums |
| schmidt_spectra.npz | Direct-SVD spectra for all 90 scale cuts; keys include n27_ux_cut5 |
| schmidt_metadata.csv | Original-field norm, entropy, cut-dimension bound, and necessary ranks at three error tolerances |
| mps_summary.csv | Resources summed over the three components; chi columns report component maxima, while storage and NVPS are summed |
| mps_components.csv | Complete records for 45 grid/component/error-tolerance combinations |
| unavoidable_error_curves.csv | Lower bound on the maximum componentwise error at each bond cap and grid; not measured full-MPS error curves |
| upper_bounds_and_error_conventions.csv | Rank bounds for octal, binary, and three-site encodings, and two error conventions |
| circuit_resources.csv | Pre-transform circuit resources for both spinor components; shared amplitude parameters are counted once |
| field_settings.json | Generation settings for each grid |
| audit_summary.json | Numerical summary of independent validation |
| provenance.json | Sources, DOIs, and checksums |

Conventions: n denotes spatial-register qubits, N=2^(n/3), and L=n/3; each MPS site has physical dimension 8.
epsilon bounds the relative L2 error of each component; 1e-2 corresponds to a squared-error budget of 1e-4.
Entropy is measured in bits from normalized velocity-tensor Schmidt coefficients, not the physical entanglement entropy of a spinor initial state.
P_store=sum(8*r_left*r_right); NVPS=P_store-sum(r_internal^2); MiB=2^20 bytes.
The previous velocity-map display clipped colors at +/-3 RMS; the stored NPZ slices remain unclipped.
The circuit and velocity-MPS tables have different output endpoints and do not define an end-to-end quantum speedup.

The added 30-qubit row includes parameter compilation and transpilation only, without statevector simulation, velocity reconstruction, or MPS calculation. Full settings, parameters, per-circuit counts, and software versions appear in circuit_n30_compilation.json. Simulation times and velocity-array sizes remain empty rather than zero.
