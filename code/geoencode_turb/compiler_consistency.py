"""Dense-sparse consistency benchmark for Table S5.

Public API: write_dense_sparse_report.
Only the report TXT persists. Large validation arrays live in an owned
temporary directory and are removed after each case.
"""

from contextlib import redirect_stdout
from dataclasses import asdict, fields
import ctypes
import gc
import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy.fft import ifftn, rfftn

__all__ = ["write_dense_sparse_report"]


def _config_from_dict(values):
    from .config import WorkflowConfig
    values = dict(values)
    for name in ("data_dir", "figure_dir"):
        if name in values:
            values[name] = Path(values[name])
    declared = {field.name: field for field in fields(WorkflowConfig)}
    unknown = set(values) - set(declared)
    if unknown:
        raise ValueError(f"Unknown configuration fields: {sorted(unknown)}")
    config = WorkflowConfig(**{
        name: value for name, value in values.items()
        if declared[name].init
    })
    # asdict includes init=False fields. Restore these only AFTER __init__
    # and __post_init__, preserving the parent's exact effective settings.
    for name, value in values.items():
        if not declared[name].init:
            setattr(config, name, value)
    return config


def _config_for_n(config, n, component_ratios):
    """Rebuild a case, allowing derived fields to follow the new grid."""
    values = {
        field.name: getattr(config, field.name)
        for field in fields(config) if field.init
    }
    values.update(nx=n // 3, ny=n // 3, nz=n // 3)
    current = type(config)(**values)
    if component_ratios is not None and n in component_ratios:
        current.component_ratio = float(component_ratios[n])
    elif n == config.n_total:
        # Preserve any deliberate population override for the current grid.
        current.component_ratio = float(config.component_ratio)
    # For another grid, init=False populations come from __post_init__;
    # ordinary init=True populations retain the user's configured value.
    return current


def _native_memory():
    """Return current RSS and native process lifetime peak, in bytes."""
    if os.name == "nt":
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = (
            wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD
        )
        psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        ok = psapi.GetProcessMemoryInfo(
            kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb
        )
        if not ok:
            raise ctypes.WinError(ctypes.get_last_error())
        return (
            int(counters.WorkingSetSize),
            int(counters.PeakWorkingSetSize),
            "Windows GetProcessMemoryInfo",
        )

    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak = int(peak if sys.platform == "darwin" else peak * 1024)
    if sys.platform.startswith("linux"):
        resident = int(Path("/proc/self/statm").read_text().split()[1])
        current = resident * os.sysconf("SC_PAGE_SIZE")
    else:
        current = None
    return current, peak, "resource.getrusage"


def _all_occupied_prefixes(bits):
    """Count all trie nodes, including the root and occupied leaves."""
    n_modes, n = bits.shape
    if n > 62:
        # A benchmark involving full grids is not practical at this size,
        # but count correctly rather than overflowing uint64.
        prefixes = [0] * n_modes
        total = 1
        for j in range(n):
            prefixes = [
                (prefix << 1) | int(bit)
                for prefix, bit in zip(prefixes, bits[:, j])
            ]
            total += len(set(prefixes))
        return total
    prefixes = np.zeros(n_modes, dtype=np.uint64)
    total = 1
    for j in range(n):
        prefixes = (prefixes << np.uint64(1)) | bits[:, j]
        total += int(np.unique(prefixes).size)
    return total


def _compile_worker(request):
    from . import quantum_encoding as qe

    cfg = _config_from_dict(request["config"])
    gc.collect()
    baseline, _, source = _native_memory()
    compiler = request["compiler"]
    start = time.perf_counter()
    if compiler == "dense":
        result = qe.compile_dense_amplitude_parameters(
            cfg, locality_radius=request["locality_radius"]
        )
    else:
        result = qe.compile_sparse_amplitude_parameters(
            cfg,
            tail_tolerance=request["tail_tolerance"],
            locality_radius=request["locality_radius"],
        )
    elapsed = time.perf_counter() - start
    current, peak, _ = _native_memory()
    diag = result.diagnostics
    if peak <= 0:
        raise RuntimeError("Native peak memory measurement failed.")

    record = {
        "compiler": compiler,
        "compilation_time_s": float(diag.compile_time_s),
        "outer_compilation_wall_time_s": elapsed,
        "peak_memory_bytes": peak,
        "peak_memory_GB": peak / 1e9,
        "peak_memory_GiB": peak / (1024**3),
        "baseline_RSS_bytes": baseline,
        "final_RSS_bytes": current,
        "peak_above_starting_RSS_bytes": (
            max(0, peak - baseline) if baseline is not None else None
        ),
        "memory_measurement": source,
        "memory_scope": (
            "native process peak in an isolated compilation-only worker; "
            "includes interpreter/import overhead; no statevector or field"
        ),
        "compiler_reported_occupied_prefixes": int(
            diag.occupied_prefix_nodes
        ),
        "fitted_prefixes": int(diag.fitted_prefix_nodes),
        "retained_array_bytes": int(diag.retained_data_bytes),
    }
    if compiler == "sparse":
        record.update({
            "retained_modes": len(result.target.probabilities),
            "occupied_prefixes": _all_occupied_prefixes(
                result.target.gray_bits
            ),
            "estimated_discarded_probability": float(
                result.target.tail_probability
            ),
            "retained_cutoff_wavenumber": float(
                result.target.cutoff_wavenumber
            ),
            "reference_tail_fraction": float(
                result.target.reference_tail_fraction
            ),
        })
    else:
        record.update({
            "dense_grid_modes": 1 << cfg.n_total,
        })
    signature = inspect.signature(qe.process_and_fit_sparse_layers)
    record["settings"] = {
        "tail_tolerance": request["tail_tolerance"],
        "locality_radius": request["locality_radius"],
        "ridge_alpha": signature.parameters["alpha"].default,
        "prefix_probability_cutoff": (
            signature.parameters["probability_cutoff"].default
        ),
        "regression_weight": (
            "sqrt(prefix probability), divided by the layer maximum"
        ),
        "bias_regularization": False,
        "controlled_Ry_pruning": "retain abs(weight) > 1e-3",
        "target_P0": 0.1,
        "target_amplitude_exponent": 0,
        "target_envelope_exponent": 10,
        "quantum_encoding_sha256": hashlib.sha256(
            Path(qe.__file__).read_bytes()
        ).hexdigest(),
    }
    # Avoid recording a stale pruning value after a future source edit.
    builder_source = inspect.getsource(qe.build_3d_coupled_circuit)
    if "np.abs(weight) > 1e-3" not in builder_source:
        raise RuntimeError(
            "The gate-pruning implementation changed. Update its recorded "
            "threshold before using this benchmark."
        )
    return {
        "record": record,
        "params": [
            [float(bias), np.asarray(weights).tolist()]
            for bias, weights in result.amp_params
        ],
    }


def _wave_worker(request):
    """Simulate ONE component and save its physical-space wavefunction."""
    from qiskit_aer import AerSimulator
    from . import quantum_encoding as qe

    cfg = _config_from_dict(request["config"])
    params = [
        (float(bias), np.asarray(weights, dtype=float))
        for bias, weights in request["params"]
    ]
    component = int(request["component"])
    seed = (
        cfg.seed_spin_up if component == 0 else cfg.seed_spin_down
    )
    circuit = qe.build_3d_coupled_circuit(
        cfg.nx, cfg.ny, cfg.nz, params, seed=seed
    )
    simulator = AerSimulator(
        method="statevector",
        max_parallel_threads=cfg.max_parallel_threads,
    )
    state = qe.simulate_fast(circuit, simulator)
    del circuit, simulator
    gc.collect()
    restored = qe.restore_amps_from_gray(
        state, cfg.nx, cfg.ny, cfg.nz
    )
    del state
    gc.collect()

    N = cfg.grid_size
    spectrum = restored.reshape((N, N, N)).T
    population = (
        cfg.component_ratio if component == 0
        else 1.0 - cfg.component_ratio
    )
    spectrum *= np.sqrt(2.0 * population)
    wave = ifftn(
        spectrum, overwrite_x=True, workers=request["fft_workers"]
    )
    np.save(request["output"], wave, allow_pickle=False)
    return {"saved": request["output"]}


def _scalar_spectrum(field, workers):
    """Same 1/2 shell sum as compute_Es; include k=0 internally."""
    N = field.shape[0]
    modes = np.fft.fftfreq(N) * N
    rmodes = np.fft.rfftfreq(N) * N
    multiplicity = np.full(rmodes.size, 2.0)
    multiplicity[0] = 1.0
    if N % 2 == 0:
        multiplicity[-1] = 1.0
    hat = rfftn(field, norm="forward", workers=workers)
    other_squared = modes[:, None]**2 + rmodes[None, :]**2
    result = np.zeros(N + 1)
    for i in range(N):
        shell = np.rint(
            np.sqrt(modes[i]**2 + other_squared)
        ).astype(np.intp)
        power = (
            hat[i].real**2 + hat[i].imag**2
        ) * (0.5 * multiplicity)
        result += np.bincount(
            shell.ravel(), weights=power.ravel(), minlength=N + 1
        )
    return result


def _close_maps(arrays):
    for array in arrays:
        mapping = getattr(array, "_mmap", None)
        if mapping is not None:
            mapping.close()


def _decode_worker(request):
    """Decode FD2 velocities to temporary files, then return Eu(k)."""
    cfg = _config_from_dict(request["config"])
    N = cfg.grid_size
    h = 2.0 * np.pi / N
    directory = Path(request["directory"])
    label = request["compiler"]
    waves = [
        np.load(directory / f"{label}_psi{a}.npy", mmap_mode="r")
        for a in (1, 2)
    ]
    velocities = [
        np.lib.format.open_memmap(
            directory / f"{label}_u{a}.npy",
            mode="w+", dtype=np.float64, shape=(N, N, N),
        )
        for a in (0, 1, 2)
    ]
    bad_velocity = 0
    bad_spin = 0
    minimum_density = np.inf

    try:
        for i in range(N):
            a, b = waves[0][i], waves[1][i]
            rho = np.abs(a)**2 + np.abs(b)**2
            minimum_density = min(minimum_density, float(rho.min()))
            bu = np.zeros((N, N), dtype=bool)
            with np.errstate(divide="ignore", invalid="ignore"):
                for axis in range(3):
                    if axis == 0:
                        da = waves[0][(i + 1) % N] - waves[0][(i - 1) % N]
                        db = waves[1][(i + 1) % N] - waves[1][(i - 1) % N]
                    else:
                        da = (
                            np.roll(a, -1, axis=axis - 1) -
                            np.roll(a, 1, axis=axis - 1)
                        )
                        db = (
                            np.roll(b, -1, axis=axis - 1) -
                            np.roll(b, 1, axis=axis - 1)
                        )
                    current = (
                        np.imag(np.conj(a) * da) +
                        np.imag(np.conj(b) * db)
                    )
                    velocities[axis][i] = current / (2 * h) / rho
                    bu |= ~np.isfinite(velocities[axis][i])
                spin = (
                    2 * (a.real * b.real + a.imag * b.imag) / rho,
                    2 * (a.real * b.imag - a.imag * b.real) / rho,
                    (np.abs(a)**2 - np.abs(b)**2) / rho,
                )
                bs = np.zeros((N, N), dtype=bool)
                for value in spin:
                    bs |= ~np.isfinite(value)
            bad_velocity += int(np.count_nonzero(bu))
            bad_spin += int(np.count_nonzero(bs))
        for velocity in velocities:
            velocity.flush()
        if bad_velocity == 0:
            energy = np.zeros(N + 1)
            for velocity in velocities:
                energy += _scalar_spectrum(
                    velocity, request["fft_workers"]
                )
        else:
            energy = np.full(N + 1, np.nan)
        return {
            "energy": energy.tolist(),
            "nonfinite_velocity_points": bad_velocity,
            "nonfinite_spin_points": bad_spin,
            "minimum_density": minimum_density,
        }
    finally:
        _close_maps(velocities)
        _close_maps(waves)


def _run_worker(request):
    """A fresh process prevents validation from contaminating compile peaks."""
    project_root = Path(__file__).resolve().parents[1]
    environment = dict(os.environ)
    environment["PYTHONPATH"] = (
        str(project_root) + os.pathsep + environment.get("PYTHONPATH", "")
    )
    process = subprocess.run(
        [sys.executable, "-m", "geoencode_turb.compiler_consistency", "--worker"],
        input=json.dumps(request, default=str),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        cwd=str(project_root),
        env=environment,
    )
    if process.returncode:
        raise RuntimeError(
            f"Benchmark worker {request['task']} failed.\n"
            + process.stderr[-12000:]
        )
    try:
        return json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            "Unexpected worker output:\n" + process.stdout[-4000:]
        ) from error


def _relative_from_sums(numerator, denominator):
    if denominator <= 0 or not np.isfinite(numerator + denominator):
        return float("nan")
    return float(np.sqrt(numerator / denominator))


def _spin_components(a, b):
    rho = np.abs(a)**2 + np.abs(b)**2
    with np.errstate(divide="ignore", invalid="ignore"):
        return (
            2 * (a.real * b.real + a.imag * b.imag) / rho,
            2 * (a.real * b.imag - a.imag * b.real) / rho,
            (np.abs(a)**2 - np.abs(b)**2) / rho,
        )


def _compare_waves(directory, N):
    """Return component fidelities and the vector-spin relative L2 error."""
    labels = ("dense", "sparse")
    arrays = {
        label: [
            np.load(directory / f"{label}_psi{a}.npy", mmap_mode="r")
            for a in (1, 2)
        ]
        for label in labels
    }
    overlaps = [0j, 0j]
    dense_norms = [0.0, 0.0]
    sparse_norms = [0.0, 0.0]
    spin_difference = spin_reference = 0.0
    try:
        for i in range(N):
            for a in (0, 1):
                reference = arrays["dense"][a][i]
                comparison = arrays["sparse"][a][i]
                overlaps[a] += np.vdot(reference, comparison)
                dense_norms[a] += float(np.vdot(reference, reference).real)
                sparse_norms[a] += float(np.vdot(comparison, comparison).real)
            sd = _spin_components(
                arrays["dense"][0][i], arrays["dense"][1][i]
            )
            ss = _spin_components(
                arrays["sparse"][0][i], arrays["sparse"][1][i]
            )
            for reference, comparison in zip(sd, ss):
                difference = comparison - reference
                spin_difference += np.sum(difference**2, dtype=np.float64)
                spin_reference += np.sum(reference**2, dtype=np.float64)

        fidelities = []
        for overlap, dn, sn in zip(overlaps, dense_norms, sparse_norms):
            fidelity = (
                float(abs(overlap)**2 / (dn * sn))
                if dn > 0 and sn > 0 else float("nan")
            )
            if np.isfinite(fidelity):
                if not -1e-10 <= fidelity <= 1 + 1e-10:
                    raise ArithmeticError("Invalid computed fidelity.")
                fidelity = float(np.clip(fidelity, 0, 1))
            fidelities.append(fidelity)
        combined_denominator = sum(dense_norms) * sum(sparse_norms)
        combined = (
            float(abs(sum(overlaps))**2 / combined_denominator)
            if combined_denominator > 0 else float("nan")
        )
        return {
            "spin_up_fidelity": fidelities[0],
            "spin_down_fidelity": fidelities[1],
            "minimum_component_fidelity": min(fidelities),
            "combined_spinor_fidelity": float(np.clip(combined, 0, 1)),
            "spin_relative_L2_difference": _relative_from_sums(
                spin_difference, spin_reference
            ),
        }
    finally:
        for label in labels:
            _close_maps(arrays[label])


def _vorticity_plane(velocity, i, N, h):
    def derivative(component, axis):
        if axis == 0:
            return (
                velocity[component][(i + 1) % N] -
                velocity[component][(i - 1) % N]
            ) / (2 * h)
        a = velocity[component][i]
        return (
            np.roll(a, -1, axis=axis - 1) -
            np.roll(a, 1, axis=axis - 1)
        ) / (2 * h)
    return (
        derivative(2, 1) - derivative(1, 2),
        derivative(0, 2) - derivative(2, 0),
        derivative(1, 0) - derivative(0, 1),
    )


def _compare_velocities(directory, N):
    arrays = {
        label: [
            np.load(directory / f"{label}_u{a}.npy", mmap_mode="r")
            for a in (0, 1, 2)
        ] for label in ("dense", "sparse")
    }
    u_diff = u_norm = w_diff = w_norm = 0.0
    h = 2 * np.pi / N
    try:
        for i in range(N):
            for reference, comparison in zip(
                arrays["dense"], arrays["sparse"]
            ):
                d = comparison[i] - reference[i]
                u_diff += np.sum(d**2, dtype=np.float64)
                u_norm += np.sum(reference[i]**2, dtype=np.float64)
            wd = _vorticity_plane(arrays["dense"], i, N, h)
            ws = _vorticity_plane(arrays["sparse"], i, N, h)
            for reference, comparison in zip(wd, ws):
                d = comparison - reference
                w_diff += np.sum(d**2, dtype=np.float64)
                w_norm += np.sum(reference**2, dtype=np.float64)
        return {
            "velocity_relative_L2_difference": _relative_from_sums(
                u_diff, u_norm
            ),
            "vorticity_relative_L2_difference": _relative_from_sums(
                w_diff, w_norm
            ),
        }
    finally:
        for label in arrays:
            _close_maps(arrays[label])


def _json_section(handle, name, values):
    handle.write("\n[" + name + "]\n")
    handle.write(json.dumps(
        values, ensure_ascii=False, indent=2, default=str, allow_nan=True
    ) + "\n")
    handle.flush()


def write_dense_sparse_report(
    config, txt_path=None, *,
    n_values=None, component_ratios=None,
    tail_tolerance=1e-12, locality_radius=3,
    fft_workers=1, scratch_dir=None, append=False,
):
    """
    Write ALL Table S5 metrics to one TXT; print nothing; return None.

    config
        Existing WorkflowConfig. Neither config nor existing flow arrays
        are modified. This is an independent benchmark, not a flow diagnostic.
    txt_path
        Default: config.data_dir / "dense_sparse_consistency.txt".
    n_values
        Optional iterable of total qubit counts, each divisible by 3.
        Default: only config.n_total. Different n use equal-size axes.
    component_ratios
        Optional {n: population} overrides. An init=False population is
        recomputed by __post_init__ for each new grid; an init=True
        population retains the configured value. The current grid preserves
        its effective population, including deliberate user overrides.
    tail_tolerance, locality_radius
        Sparse truncation and MATCHED dense/sparse locality setting.
    fft_workers
        FFT threads in the validation stage.
    scratch_dir
        Directory for owned temporary arrays, default report directory.
        Approximately 112*N^3 bytes of disk is required per case, plus margin.
        Temporary arrays are deleted after the case; only TXT persists.
    append
        True appends; False overwrites the report.

    Reported fidelity is the MINIMUM of the two component dense-sparse
    fidelities. Maximum downstream difference includes ALL FOUR quantities:
    velocity, vorticity, vector spin, and the full velocity energy spectrum.
    All relative differences use the dense compiled result as denominator.
    This does not measure fidelity with an exact analytic target.

    Supports the project's FDM velocity and vorticity definitions.
    Compilation resources are measured in separate fresh workers, excluding
    statevector emulation, decoding, comparisons and temporary-file I/O.
    """
    from .config import WorkflowConfig
    if not isinstance(config, WorkflowConfig):
        raise TypeError("config must be the current WorkflowConfig.")
    if config.velocity_method != "FDM" or config.derivative_method != "FDM":
        raise ValueError(
            "This benchmark matches the current FDM workflow only."
        )
    if not 0 < float(tail_tolerance) < 1:
        raise ValueError("tail_tolerance must lie between 0 and 1.")
    if locality_radius is not None and (
        not isinstance(locality_radius, int) or locality_radius < 0
    ):
        raise ValueError("locality_radius must be a nonnegative int or None.")
    if not isinstance(fft_workers, int) or fft_workers < 1:
        raise ValueError("fft_workers must be a positive integer.")

    if n_values is None:
        if not config.nx == config.ny == config.nz:
            raise ValueError("The benchmark requires equal-size axes.")
        n_values = (config.n_total,)
    else:
        n_values = tuple(n_values)
    if not n_values:
        raise ValueError("n_values is empty.")
    for n in n_values:
        if not isinstance(n, (int, np.integer)) or n <= 0 or n % 3:
            raise ValueError("Each n must be a positive multiple of 3.")

    path = Path(
        txt_path if txt_path is not None
        else Path(config.data_dir) / "dense_sparse_consistency.txt"
    ).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    scratch = Path(
        path.parent if scratch_dir is None else scratch_dir
    ).expanduser().resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    rows = []

    with path.open("a" if append else "w", encoding="utf-8") as handle:
        handle.write("\n" + "=" * 78 + "\n")
        _json_section(handle, "BENCHMARK_DEFINITION", {
            "n_values": n_values,
            "phase_realizations_per_n": 1,
            "reference": "dense compiled state and decoded dense fields",
            "comparison": "sparse compiled state and decoded sparse fields",
            "fidelity_in_table": "minimum of component squared overlaps",
            "maximum_difference": (
                "max(relative L2 differences of velocity, vorticity, spin, Eu)"
            ),
            "relative_difference_denominator": "dense reference",
            "spectrum": (
                "1/2 shell sum of |FFT(ux)/N^3|^2 + "
                "|FFT(uy)/N^3|^2 + |FFT(uz)/N^3|^2; k=1,...,N"
            ),
            "occupied_prefix_definition": (
                "all occupied nodes of the retained Gray-bit trie, "
                "including root and leaves; fitted-prefix count also reported"
            ),
            "memory_unit": "GB = 10^9 bytes; GiB also reported separately",
            "interpretation": (
                "compiler consistency only; not analytic-target preparation "
                "accuracy and not an end-to-end speedup test"
            ),
        })

        for requested_n in n_values:
            n = int(requested_n)
            cfg = _config_for_n(config, n, component_ratios)
            ratio = cfg.component_ratio
            if not np.isfinite(ratio) or not 0 < ratio < 1:
                raise ValueError("Component populations must be in (0,1).")
            N = cfg.grid_size
            request_common = {
                "config": asdict(cfg),
                "tail_tolerance": float(tail_tolerance),
                "locality_radius": locality_radius,
                "fft_workers": fft_workers,
            }
            _json_section(handle, f"CASE n={n}", {
                "grid": [N, N, N], "config": asdict(cfg)
            })
            try:
                # Only these fresh processes are used for reported compilation
                # times/memory. No statevector or full field is created there.
                compilers = {}
                for label in ("dense", "sparse"):
                    compilers[label] = _run_worker({
                        **request_common, "task": "compile", "compiler": label
                    })
                    _json_section(
                        handle, f"{label.upper()}_COMPILATION",
                        compilers[label]["record"],
                    )

                # 4 complex128 waves + 6 float64 velocity arrays.
                estimated_bytes = 112 * (N**3)
                required_free = int(1.15 * estimated_bytes + 1024**3)
                available = shutil.disk_usage(scratch).free
                _json_section(handle, "VALIDATION_STORAGE", {
                    "temporary_directory_parent": str(scratch),
                    "array_storage_bytes": estimated_bytes,
                    "required_free_bytes_with_margin": required_free,
                    "available_free_bytes": available,
                })
                if available < required_free:
                    raise OSError(
                        f"Insufficient scratch disk space for n={n}: "
                        f"need {required_free / 1e9:.2f} GB free, "
                        f"have {available / 1e9:.2f} GB. "
                        "Use scratch_dir on another drive."
                    )

                start_validation = time.perf_counter()
                with tempfile.TemporaryDirectory(
                    prefix=f"dense_sparse_n{n}_", dir=str(scratch)
                ) as work:
                    directory = Path(work)
                    decoded = {}
                    for label in ("dense", "sparse"):
                        for component in (0, 1):
                            _run_worker({
                                **request_common,
                                "task": "wave",
                                "compiler": label,
                                "params": compilers[label]["params"],
                                "component": component,
                                "output": str(
                                    directory /
                                    f"{label}_psi{component + 1}.npy"
                                ),
                            })
                        decoded[label] = _run_worker({
                            **request_common,
                            "task": "decode",
                            "compiler": label,
                            "directory": str(directory),
                        })
                    differences = _compare_waves(directory, N)
                    differences.update(_compare_velocities(directory, N))
                    # Match current compute_Ek: omit the zero shell.
                    dense_energy = np.asarray(decoded["dense"]["energy"])[1:]
                    sparse_energy = np.asarray(decoded["sparse"]["energy"])[1:]
                    differences["energy_spectrum_relative_L2_difference"] = (
                        _relative_from_sums(
                            float(np.sum((sparse_energy - dense_energy)**2)),
                            float(np.sum(dense_energy**2)),
                        )
                    )

                four_differences = [
                    differences[name] for name in (
                        "velocity_relative_L2_difference",
                        "vorticity_relative_L2_difference",
                        "spin_relative_L2_difference",
                        "energy_spectrum_relative_L2_difference",
                    )
                ]
                maximum = (
                    float(max(four_differences))
                    if np.all(np.isfinite(four_differences))
                    else float("nan")
                )
                differences["maximum_downstream_relative_difference"] = maximum
                differences["validation_wall_time_s_excluded_from_compilation"] = (
                    time.perf_counter() - start_validation
                )
                for label in ("dense", "sparse"):
                    differences[f"{label}_nonfinite_velocity_points"] = (
                        decoded[label]["nonfinite_velocity_points"]
                    )
                    differences[f"{label}_nonfinite_spin_points"] = (
                        decoded[label]["nonfinite_spin_points"]
                    )
                    differences[f"{label}_minimum_density"] = (
                        decoded[label]["minimum_density"]
                    )
                _json_section(handle, "DENSE_SPARSE_DIFFERENCES", differences)
                sparse_record = compilers["sparse"]["record"]
                row = {
                    "n": n,
                    "retained_modes": sparse_record["retained_modes"],
                    "occupied_prefixes": sparse_record["occupied_prefixes"],
                    "sparse_compilation_time_s": (
                        sparse_record["compilation_time_s"]
                    ),
                    "sparse_peak_memory_GB": sparse_record["peak_memory_GB"],
                    "dense_sparse_fidelity": (
                        differences["minimum_component_fidelity"]
                    ),
                    "maximum_downstream_relative_difference": maximum,
                }
                rows.append(row)
                _json_section(handle, "TABLE_S5_ROW", row)
                handle.write(
                    "\n[VELOCITY_ENERGY_SPECTRA]\n"
                    "k\tE_dense\tE_sparse\n"
                )
                np.savetxt(
                    handle,
                    np.column_stack((
                        np.arange(1, N + 1), dense_energy, sparse_energy,
                    )),
                    delimiter="\t", fmt="%.12e",
                )
                handle.flush()
                gc.collect()
            except Exception as error:
                _json_section(handle, "CASE_FAILED", {
                    "n": n,
                    "error": str(error),
                    "note": (
                        "No missing result is replaced by zero or by "
                        "dense-sparse agreement from another case."
                    ),
                })
                raise

        handle.write(
            "\n[TABLE_S5_SUMMARY]\n"
            "n\tRetained_modes\tOccupied_prefixes\t"
            "Sparse_compilation_time_s\tSparse_peak_memory_GB\t"
            "Dense_sparse_fidelity\tMaximum_downstream_relative_difference\n"
        )
        for row in rows:
            handle.write(
                f"{row['n']}\t{row['retained_modes']}\t"
                f"{row['occupied_prefixes']}\t"
                f"{row['sparse_compilation_time_s']:.12e}\t"
                f"{row['sparse_peak_memory_GB']:.12e}\t"
                f"{row['dense_sparse_fidelity']:.12e}\t"
                f"{row['maximum_downstream_relative_difference']:.12e}\n"
            )
    return None


def _worker_entry():
    request = json.loads(sys.stdin.read())
    task = request["task"]
    with redirect_stdout(sys.stderr):
        if task == "compile":
            result = _compile_worker(request)
        elif task == "wave":
            result = _wave_worker(request)
        elif task == "decode":
            result = _decode_worker(request)
        else:
            raise ValueError(f"Unknown task: {task}")
    sys.stdout.write(json.dumps(result, default=str, allow_nan=True))


if __name__ == "__main__":
    if sys.argv[1:] != ["--worker"]:
        raise SystemExit("Import write_dense_sparse_report from your notebook.")
    _worker_entry()

