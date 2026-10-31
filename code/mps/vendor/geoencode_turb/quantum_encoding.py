"""Gray-code mapping, dense/sparse compilation, and circuit simulation."""

from dataclasses import dataclass, replace
import gc
import os
from pathlib import Path
import threading
import time
from typing import Optional

import numpy as np
import pandas as pd
from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit_aer import AerSimulator

from .config import WorkflowConfig


@dataclass
class EncodingResult:
    """Outputs of the quantum state-preparation stage."""

    target_probs_phys: Optional[np.ndarray]
    target_probs_gray: Optional[np.ndarray]
    KZ: Optional[np.ndarray]
    KY: Optional[np.ndarray]
    KX: Optional[np.ndarray]
    amp_params: list
    sta_vec_1_data: Optional[np.ndarray]
    sta_vec_2_data: Optional[np.ndarray]
    weights_matrix: np.ndarray
    amp_sim_phys: Optional[np.ndarray]
    k_plot: np.ndarray
    amp_plot: np.ndarray
    c_sim: np.ndarray
    v_sim: np.ndarray
    compiler_diagnostics: Optional["CompilerDiagnostics"] = None
    sparse_target: Optional["SparseTargetDistribution"] = None


@dataclass
class CompilerDiagnostics:
    """Classical resources and support information for circuit compilation."""

    compiler: str
    n_total: int
    dense_probability_values: int
    significant_modes: int
    occupied_prefix_nodes: int
    fitted_prefix_nodes: int
    tail_probability: float
    cutoff_wavenumber: float
    locality_radius: Optional[int]
    compile_time_s: float
    baseline_rss_bytes: int
    peak_rss_bytes: int
    incremental_peak_memory_bytes: int
    retained_data_bytes: int


@dataclass
class SparseTargetDistribution:
    """Compact normalized target retained by the sparse compiler.

    The rows of ``wavenumbers`` and ``gray_bits`` refer to the same retained
    modes. ``probabilities`` sum to one after removing the discarded tail and
    renormalizing the retained state.
    """

    wavenumbers: np.ndarray
    physical_indices: np.ndarray
    gray_bits: np.ndarray
    probabilities: np.ndarray
    full_probabilities: np.ndarray
    tail_probability: float
    cutoff_wavenumber: float
    reference_tail_fraction: float


@dataclass
class DenseCompilationResult:
    """Intermediate outputs of the dense reference compiler."""

    target_probs_phys: np.ndarray
    target_probs_gray: np.ndarray
    KZ: np.ndarray
    KY: np.ndarray
    KX: np.ndarray
    amp_params: list
    diagnostics: CompilerDiagnostics


@dataclass
class SparseCompilationResult:
    """Intermediate outputs of the sparse compiler."""

    target: SparseTargetDistribution
    amp_params: list
    diagnostics: CompilerDiagnostics


@dataclass
class EncodingBenchmarkResult:
    """Dense-versus-sparse comparison used for manuscript benchmarks."""

    n_total: int
    dense_compile_time_s: float
    sparse_compile_time_s: float
    compilation_speedup: float
    dense_incremental_peak_memory_bytes: int
    sparse_incremental_peak_memory_bytes: int
    peak_memory_reduction_fraction: float
    dense_probability_values: int
    sparse_significant_modes: int
    sparse_occupied_prefix_nodes: int
    sparse_fitted_prefix_nodes: int
    sparse_tail_probability: float
    spin_up_state_fidelity: float
    spin_down_state_fidelity: float
    minimum_state_fidelity: float
    energy_spectrum_relative_l2_error: float
    velocity_relative_l2_error: float
    vorticity_relative_l2_error: float
    vortex_surface_relative_l2_error: float
    maximum_physical_relative_l2_error: float

    def as_dict(self):
        """Return a flat record suitable for a CSV row."""
        return dict(self.__dict__)


class _RSSMonitor:
    """Poll resident memory using only the Python standard library."""

    def __init__(self, interval=0.01):
        self.interval = float(interval)
        self.baseline = 0
        self.peak = 0
        self._stop_event = threading.Event()
        self._thread = None

    @staticmethod
    def _current_rss_bytes():
        try:
            with open("/proc/self/statm", "r", encoding="utf-8") as handle:
                resident_pages = int(handle.readline().split()[1])
            return resident_pages * os.sysconf("SC_PAGE_SIZE")
        except (OSError, IndexError, ValueError):
            return 0

    def _poll(self):
        while not self._stop_event.wait(self.interval):
            self.peak = max(self.peak, self._current_rss_bytes())

    def __enter__(self):
        gc.collect()
        self.baseline = self._current_rss_bytes()
        self.peak = self.baseline
        self._stop_event.clear()
        self._thread = threading.Thread(target=self._poll, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.peak = max(self.peak, self._current_rss_bytes())
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join()

    @property
    def incremental_peak(self):
        return max(0, self.peak - self.baseline)


def get_gray_indices(n):
    """Generate the length-2^n binary-to-Gray index map."""
    idxs = np.arange(2**n, dtype=np.int32)
    return np.bitwise_xor(idxs, np.right_shift(idxs, 1))


def permute_probs_to_gray(probs_physical, nx, ny, nz):
    """Reorder the three physical axes directly with fancy indexing."""
    gz_indices = get_gray_indices(nz)
    gy_indices = get_gray_indices(ny)
    gx_indices = get_gray_indices(nx)

    # The mapping is probs_gray[gray_index] = probs_physical[physical_index].
    out = np.zeros_like(probs_physical)
    out[gz_indices, :, :] = probs_physical

    res = np.zeros_like(out)
    res[:, gy_indices, :] = out

    final = np.zeros_like(res)
    final[:, :, gx_indices] = res
    return final


def restore_amps_from_gray(amps_gray_flat, nx, ny, nz):
    """Restore Gray-ordered amplitudes to the physical axis ordering."""
    Nz, Ny, Nx = 2**nz, 2**ny, 2**nx
    amps_gray_3d = amps_gray_flat.reshape((Nz, Ny, Nx))

    gz_indices = get_gray_indices(nz)
    gy_indices = get_gray_indices(ny)
    gx_indices = get_gray_indices(nx)

    # The inverse operation is physical[z] = gray[binary_to_gray(z)].
    amps_phys = amps_gray_3d[gz_indices, :, :]
    amps_phys = amps_phys[:, gy_indices, :]
    amps_phys = amps_phys[:, :, gx_indices]
    return amps_phys.flatten()


def generate_3d_target_distribution(nx, ny, nz, k_cutoff):
    """Generate the three-dimensional target probability distribution."""
    Nx, Ny, Nz = 2**nx, 2**ny, 2**nz
    kx = np.fft.fftfreq(Nx, d=1 / Nx)
    ky = np.fft.fftfreq(Ny, d=1 / Ny)
    kz = np.fft.fftfreq(Nz, d=1 / Nz)

    KZ, KY, KX = np.meshgrid(kz, ky, kx, indexing="ij")
    k_sq = KX**2 + KY**2 + KZ**2
    k_abs = np.sqrt(k_sq)

    amplitudes = np.zeros_like(k_sq)
    nonzero_mask = k_sq > 0

    amplitudes[nonzero_mask] = k_sq[nonzero_mask] ** 0

    # A super-Gaussian window suppresses hard-cutoff oscillations.
    smooth_mask = np.exp(-((k_abs / k_cutoff) ** 10))
    amplitudes *= smooth_mask

    probs = amplitudes**2
    target_zero_prob = 0.1

    tail_sum = np.sum(probs)
    if tail_sum > 0:
        probs = probs * (1.0 - target_zero_prob) / tail_sum

    probs[0, 0, 0] = target_zero_prob
    return probs, (KZ, KY, KX)


def solve_ridge_numpy(X, y, weights, alpha=1e-6):
    """Solve the weighted ridge-regression problem in closed form."""
    N, M = X.shape
    X_aug = np.empty((N, M + 1), dtype=X.dtype)
    X_aug[:, :-1] = X
    X_aug[:, -1] = 1.0

    # The weighted normal equations retain the original implementation.
    sqrt_w = np.sqrt(weights)[:, None]
    X_w = X_aug * sqrt_w
    y_w = y * weights
    _ = X_w, y_w

    X_T_W = X_aug.T * weights[None, :]
    A = X_T_W @ X_aug
    b = X_T_W @ y

    # The bias is not regularized.
    I = np.eye(M + 1)
    I[-1, -1] = 0.0
    A += alpha * I

    try:
        coefs = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        coefs = np.linalg.lstsq(A, b, rcond=None)[0]

    return coefs[-1], coefs[:-1]


def _validate_locality_radius(locality_radius, n_total):
    """Validate and normalize the backward coupling radius ``r``."""
    if locality_radius is None:
        return None
    if not isinstance(locality_radius, (int, np.integer)):
        raise TypeError("locality_radius must be an integer or None.")
    if locality_radius < 0:
        raise ValueError("locality_radius must be non-negative.")
    return min(int(locality_radius), max(0, n_total - 1))


def _active_feature_start(layer_idx, locality_radius):
    """Return the first precursor included in the r-local affine ansatz."""
    if locality_radius is None:
        return 0
    return max(0, layer_idx - locality_radius)


def process_and_fit_all_layers(
    probs_gray_3d,
    nx,
    ny,
    nz,
    locality_radius=None,
):
    """Fit dense conditional angles with an optional r-local ansatz.

    For circuit layer ``j``, ``locality_radius=r`` retains only precursor
    qubits ``m`` satisfying ``max(0, j-r) <= m <= j-1``. Passing ``None``
    preserves the original fully connected affine fit.
    """
    # Flattening follows the Z-to-Y-to-X qubit order.
    probs_flat = probs_gray_3d.flatten()
    n_total = nx + ny + nz
    locality_radius = _validate_locality_radius(locality_radius, n_total)
    cutoff = 1e-25
    params = []
    current_probs = probs_flat

    # The layers are processed from the most significant to the least
    # significant bit. Vectorized marginalization replaces recursion.
    for layer_idx in range(n_total):
        n_precursors = layer_idx
        feature_start = _active_feature_start(layer_idx, locality_radius)
        n_features = layer_idx - feature_start
        n_blocks = 2**n_precursors
        remainder = current_probs.shape[0] // (n_blocks * 2)
        reshaped = current_probs.reshape(n_blocks, 2, remainder)
        block_probs = reshaped.sum(axis=2)

        p0 = block_probs[:, 0]
        p1 = block_probs[:, 1]
        w = p0 + p1

        y_angles = np.zeros(n_blocks)
        valid_mask = w > cutoff
        if np.any(valid_mask):
            ratio = np.clip(p0[valid_mask] / w[valid_mask], 0, 1)
            y_angles[valid_mask] = 2 * np.arccos(np.sqrt(ratio))

        mask = valid_mask
        if np.sum(mask) == 0:
            params.append((0.0, []))
            continue

        y_train = y_angles[mask]
        w_train = w[mask]
        indices = np.arange(n_blocks)[mask]

        if n_precursors == 0:
            X_train = np.zeros((len(y_train), 0))
        else:
            all_bits = (
                (
                    indices[:, None]
                    & (1 << np.arange(n_precursors)[::-1])
                )
                > 0
            ).astype(float)
            X_train = all_bits[:, feature_start:]

        amp_weights = w_train**0.5
        max_w = np.max(amp_weights)
        if max_w > 0:
            amp_weights = amp_weights / max_w

        if n_features == 0:
            bias = np.average(y_train, weights=amp_weights)
            local_weights = np.empty(0, dtype=float)
        elif len(y_train) <= 1:
            bias = y_train[0]
            local_weights = np.zeros(n_features)
        else:
            bias, local_weights = solve_ridge_numpy(
                X_train,
                y_train,
                amp_weights,
                alpha=1e-9,
            )

        layer_weights = np.zeros(n_precursors, dtype=float)
        if n_features > 0:
            layer_weights[feature_start:] = local_weights
        params.append((bias, layer_weights))

    return params


def _raw_nonzero_mode_probabilities(k_squared, k_cutoff):
    """Evaluate the unnormalized nonzero probabilities used by the model."""
    k_squared = np.asarray(k_squared, dtype=np.float64)
    k_abs = np.sqrt(k_squared)
    amplitude = np.power(k_squared, 0)
    amplitude *= np.exp(-np.power(k_abs / k_cutoff, 10))
    return amplitude**2


def _axis_wavenumbers_within_radius(n_bits, radius):
    """Return valid FFT wavenumbers without allocating a length-2^n axis."""
    grid_size = 1 << int(n_bits)
    lower = max(-grid_size // 2, -int(radius))
    upper = min(grid_size // 2 - 1, int(radius))
    return np.arange(lower, upper + 1, dtype=np.int64)


def _initial_reference_radius(k_cutoff, tail_tolerance):
    """Choose a conservative n-independent search radius for the cutoff."""
    numerical_floor = min(max(tail_tolerance * 1.0e-4, 1.0e-300), 1.0e-16)
    radius_factor = (0.5 * np.log(1.0 / numerical_floor)) ** 0.1
    return max(2, int(np.ceil(k_cutoff * radius_factor)) + 2)


def _enumerate_reference_modes(
    nx,
    ny,
    nz,
    k_cutoff,
    tail_tolerance,
):
    """Enumerate only the n-independent effective spectral neighborhood."""
    radius = _initial_reference_radius(k_cutoff, tail_tolerance)
    target_boundary_fraction = max(tail_tolerance * 1.0e-4, 1.0e-18)

    for _ in range(12):
        kx = _axis_wavenumbers_within_radius(nx, radius)
        ky = _axis_wavenumbers_within_radius(ny, radius)
        kz = _axis_wavenumbers_within_radius(nz, radius)
        KZ, KY, KX = np.meshgrid(kz, ky, kx, indexing="ij")
        k_squared = KX**2 + KY**2 + KZ**2
        nonzero = k_squared > 0

        wavevectors = np.column_stack(
            (KZ[nonzero], KY[nonzero], KX[nonzero])
        ).astype(np.int64, copy=False)
        k_squared_nonzero = k_squared[nonzero].astype(np.float64, copy=False)
        raw_probabilities = _raw_nonzero_mode_probabilities(
            k_squared_nonzero,
            k_cutoff,
        )
        total_raw = float(np.sum(raw_probabilities, dtype=np.float64))
        if total_raw <= 0.0 or not np.isfinite(total_raw):
            raise FloatingPointError(
                "The sparse target has a non-positive or non-finite weight."
            )

        fully_enumerated = (
            len(kx) == (1 << nx)
            and len(ky) == (1 << ny)
            and len(kz) == (1 << nz)
        )
        if fully_enumerated:
            boundary_fraction = 0.0
            break

        boundary = np.max(np.abs(wavevectors), axis=1) >= max(0, radius - 1)
        boundary_fraction = float(
            np.sum(raw_probabilities[boundary], dtype=np.float64) / total_raw
        )
        if boundary_fraction <= target_boundary_fraction:
            break
        radius = int(np.ceil(1.25 * radius)) + 1
    else:
        raise RuntimeError(
            "The sparse reference support did not converge; increase the "
            "search-radius limit or relax tail_tolerance."
        )

    return wavevectors, k_squared_nonzero, raw_probabilities, boundary_fraction


def _fft_index_from_wavenumber(wavenumber, n_bits):
    """Map a signed FFT wavenumber to its nonnegative array index."""
    wavenumber = int(wavenumber)
    return wavenumber if wavenumber >= 0 else (1 << int(n_bits)) + wavenumber


def _gray_bits_from_wavevectors(wavenumbers, nx, ny, nz):
    """Convert retained (kz, ky, kx) modes to Gray-ordered bit strings."""
    widths = (nz, ny, nx)
    bit_offsets = np.cumsum((0,) + widths)
    n_total = nx + ny + nz
    n_modes = len(wavenumbers)
    gray_bits = np.empty((n_modes, n_total), dtype=np.uint8)
    physical_indices = np.empty((n_modes, 3), dtype=object)

    for mode_index, mode in enumerate(wavenumbers):
        for axis_index, (wavenumber, width) in enumerate(zip(mode, widths)):
            physical_index = _fft_index_from_wavenumber(wavenumber, width)
            gray_index = physical_index ^ (physical_index >> 1)
            physical_indices[mode_index, axis_index] = physical_index
            start = int(bit_offsets[axis_index])
            for local_bit in range(width):
                shift = width - 1 - local_bit
                gray_bits[mode_index, start + local_bit] = (
                    gray_index >> shift
                ) & 1

    if max(widths) <= 62:
        physical_indices = physical_indices.astype(np.int64)
    return physical_indices, gray_bits


def generate_sparse_target_distribution(
    nx,
    ny,
    nz,
    k_cutoff,
    tail_tolerance=1.0e-12,
):
    """Construct a normalized sparse target without a length-2^n array.

    The retained support is the smallest radial set, including complete
    degenerate shells, whose estimated discarded probability does not exceed
    ``tail_tolerance``. The zero mode is always retained.
    """
    if not 0.0 < tail_tolerance < 1.0:
        raise ValueError("tail_tolerance must lie strictly between zero and one.")
    if k_cutoff <= 0.0:
        raise ValueError("k_cutoff must be positive.")

    (
        nonzero_wavevectors,
        k_squared,
        raw_probabilities,
        reference_tail_fraction,
    ) = _enumerate_reference_modes(
        nx,
        ny,
        nz,
        k_cutoff,
        tail_tolerance,
    )

    order = np.argsort(k_squared, kind="stable")
    sorted_k_squared = k_squared[order]
    sorted_raw = raw_probabilities[order]
    shell_squared, shell_starts = np.unique(
        sorted_k_squared,
        return_index=True,
    )
    shell_raw = np.add.reduceat(sorted_raw, shell_starts)
    total_raw = float(np.sum(shell_raw, dtype=np.float64))
    tail_after_shell = np.zeros_like(shell_raw)
    if len(shell_raw) > 1:
        tail_after_shell[:-1] = np.cumsum(
            shell_raw[:0:-1],
            dtype=np.float64,
        )[::-1]

    reference_tail_probability = 0.9 * reference_tail_fraction
    internal_tail_budget = max(
        0.0,
        tail_tolerance - reference_tail_probability,
    )
    admissible_shells = np.flatnonzero(
        0.9 * tail_after_shell / total_raw <= internal_tail_budget
    )
    if len(admissible_shells) == 0:
        cutoff_shell = len(shell_squared) - 1
    else:
        cutoff_shell = int(admissible_shells[0])
    cutoff_squared = shell_squared[cutoff_shell]
    retained_nonzero = k_squared <= cutoff_squared
    retained_raw = raw_probabilities[retained_nonzero]
    retained_raw_sum = float(np.sum(retained_raw, dtype=np.float64))
    discarded_raw_sum = float(
        np.sum(raw_probabilities[~retained_nonzero], dtype=np.float64)
    )

    internal_tail_probability = 0.9 * discarded_raw_sum / total_raw
    tail_probability = min(
        1.0,
        internal_tail_probability + reference_tail_probability,
    )

    retained_wavevectors = nonzero_wavevectors[retained_nonzero]
    full_nonzero_probabilities = 0.9 * retained_raw / total_raw
    zero_mode = np.zeros((1, 3), dtype=np.int64)
    wavenumbers = np.vstack((zero_mode, retained_wavevectors))
    full_probabilities = np.concatenate(
        (np.array([0.1]), full_nonzero_probabilities)
    )
    retained_mass = float(np.sum(full_probabilities, dtype=np.float64))
    probabilities = full_probabilities / retained_mass

    physical_indices, gray_bits = _gray_bits_from_wavevectors(
        wavenumbers,
        nx,
        ny,
        nz,
    )
    return SparseTargetDistribution(
        wavenumbers=wavenumbers,
        physical_indices=physical_indices,
        gray_bits=gray_bits,
        probabilities=probabilities,
        full_probabilities=full_probabilities,
        tail_probability=tail_probability,
        cutoff_wavenumber=float(np.sqrt(cutoff_squared)),
        reference_tail_fraction=float(reference_tail_fraction),
    )


def _solve_ridge_from_sufficient_statistics(matrix, vector, n_features, alpha):
    """Solve a weighted affine fit from accumulated normal equations."""
    regularizer = np.eye(n_features + 1, dtype=np.float64)
    regularizer[-1, -1] = 0.0
    matrix = matrix + alpha * regularizer
    try:
        coefficients = np.linalg.solve(matrix, vector)
    except np.linalg.LinAlgError:
        coefficients = np.linalg.lstsq(matrix, vector, rcond=None)[0]
    return float(coefficients[-1]), coefficients[:-1]


def process_and_fit_sparse_layers(
    sparse_target,
    locality_radius=None,
    alpha=1.0e-9,
    probability_cutoff=1.0e-25,
):
    """Fit conditional rotations from occupied prefixes only.

    No dense probability vector is formed. For layer ``j``, the r-local
    affine ansatz is

    ``theta_j = b_j + sum_{m=max(0,j-r)}^{j-1} w_{j,m} q_m``.

    ``locality_radius=None`` includes all precursor bits and is useful for a
    direct comparison with the original dense compiler.
    """
    if not isinstance(sparse_target, SparseTargetDistribution):
        raise TypeError("sparse_target must be a SparseTargetDistribution.")

    gray_bits = sparse_target.gray_bits
    probabilities = sparse_target.probabilities
    n_modes, n_total = gray_bits.shape
    locality_radius = _validate_locality_radius(locality_radius, n_total)

    prefix_codes = [0] * n_modes
    params = []
    fitted_prefix_nodes = 0

    for layer_idx in range(n_total):
        groups = {}
        layer_bits = gray_bits[:, layer_idx]
        for mode_index in range(n_modes):
            prefix = prefix_codes[mode_index]
            branch = int(layer_bits[mode_index])
            if prefix not in groups:
                groups[prefix] = [0.0, 0.0]
            groups[prefix][branch] += float(probabilities[mode_index])

        feature_start = _active_feature_start(layer_idx, locality_radius)
        n_features = layer_idx - feature_start

        samples = []
        maximum_weight = 0.0
        for prefix, (p0, p1) in groups.items():
            total_probability = p0 + p1
            if total_probability <= probability_cutoff:
                continue
            ratio = np.clip(p0 / total_probability, 0.0, 1.0)
            angle = 2.0 * np.arccos(np.sqrt(ratio))
            amplitude_weight = np.sqrt(total_probability)
            maximum_weight = max(maximum_weight, amplitude_weight)
            samples.append((prefix, angle, amplitude_weight))

        fitted_prefix_nodes += len(samples)

        if not samples:
            params.append((0.0, np.zeros(layer_idx, dtype=float)))
        elif n_features == 0:
            normalized_weights = np.array(
                [sample[2] for sample in samples],
                dtype=np.float64,
            )
            normalized_weights /= np.max(normalized_weights)
            angles = np.array(
                [sample[1] for sample in samples],
                dtype=np.float64,
            )
            bias = float(np.average(angles, weights=normalized_weights))
            params.append((bias, np.zeros(layer_idx, dtype=float)))
        elif len(samples) == 1:
            params.append(
                (
                    float(samples[0][1]),
                    np.zeros(layer_idx, dtype=float),
                )
            )
        else:
            normal_matrix = np.zeros(
                (n_features + 1, n_features + 1),
                dtype=np.float64,
            )
            normal_vector = np.zeros(n_features + 1, dtype=np.float64)

            for prefix, angle, amplitude_weight in samples:
                feature_vector = np.ones(n_features + 1, dtype=np.float64)
                for local_feature in range(n_features):
                    shift = n_features - 1 - local_feature
                    feature_vector[local_feature] = (prefix >> shift) & 1
                sample_weight = amplitude_weight / maximum_weight
                normal_matrix += sample_weight * np.outer(
                    feature_vector,
                    feature_vector,
                )
                normal_vector += sample_weight * feature_vector * angle

            bias, local_weights = _solve_ridge_from_sufficient_statistics(
                normal_matrix,
                normal_vector,
                n_features,
                alpha,
            )
            layer_weights = np.zeros(layer_idx, dtype=float)
            layer_weights[feature_start:] = local_weights
            params.append((bias, layer_weights))

        prefix_codes = [
            (prefix << 1) | int(bit)
            for prefix, bit in zip(prefix_codes, layer_bits)
        ]

    occupied_prefix_nodes = fitted_prefix_nodes + len(set(prefix_codes))
    return params, occupied_prefix_nodes, fitted_prefix_nodes


def _retained_array_bytes(*arrays):
    """Return the total storage of the supplied NumPy arrays."""
    return int(
        sum(array.nbytes for array in arrays if isinstance(array, np.ndarray))
    )


def compile_dense_amplitude_parameters(config=None, locality_radius=None):
    """Run the original dense compiler and report its classical resources."""
    config = WorkflowConfig() if config is None else config
    n_total = config.n_total
    locality_radius = _validate_locality_radius(locality_radius, n_total)
    start = time.perf_counter()
    with _RSSMonitor() as memory:
        target_probs_phys, (KZ, KY, KX) = generate_3d_target_distribution(
            config.nx,
            config.ny,
            config.nz,
            config.k_cutoff,
        )
        target_probs_gray = permute_probs_to_gray(
            target_probs_phys,
            config.nx,
            config.ny,
            config.nz,
        )
        amp_params = process_and_fit_all_layers(
            target_probs_gray,
            config.nx,
            config.ny,
            config.nz,
            locality_radius=locality_radius,
        )
    elapsed = time.perf_counter() - start

    retained_bytes = _retained_array_bytes(
        target_probs_phys,
        target_probs_gray,
        KZ,
        KY,
        KX,
    )
    diagnostics = CompilerDiagnostics(
        compiler="dense",
        n_total=n_total,
        dense_probability_values=1 << n_total,
        significant_modes=1 << n_total,
        occupied_prefix_nodes=(1 << (n_total + 1)) - 1,
        fitted_prefix_nodes=(1 << n_total) - 1,
        tail_probability=0.0,
        cutoff_wavenumber=float("nan"),
        locality_radius=locality_radius,
        compile_time_s=elapsed,
        baseline_rss_bytes=memory.baseline,
        peak_rss_bytes=memory.peak,
        incremental_peak_memory_bytes=memory.incremental_peak,
        retained_data_bytes=retained_bytes,
    )
    return DenseCompilationResult(
        target_probs_phys=target_probs_phys,
        target_probs_gray=target_probs_gray,
        KZ=KZ,
        KY=KY,
        KX=KX,
        amp_params=amp_params,
        diagnostics=diagnostics,
    )


def compile_sparse_amplitude_parameters(
    config=None,
    tail_tolerance=1.0e-12,
    locality_radius=None,
):
    """Compile amplitude parameters from a sparse prefix tree."""
    config = WorkflowConfig() if config is None else config
    n_total = config.n_total
    locality_radius = _validate_locality_radius(locality_radius, n_total)
    start = time.perf_counter()
    with _RSSMonitor() as memory:
        sparse_target = generate_sparse_target_distribution(
            config.nx,
            config.ny,
            config.nz,
            config.k_cutoff,
            tail_tolerance=tail_tolerance,
        )
        (
            amp_params,
            occupied_prefix_nodes,
            fitted_prefix_nodes,
        ) = process_and_fit_sparse_layers(
            sparse_target,
            locality_radius=locality_radius,
        )
    elapsed = time.perf_counter() - start

    retained_bytes = _retained_array_bytes(
        sparse_target.wavenumbers,
        sparse_target.physical_indices,
        sparse_target.gray_bits,
        sparse_target.probabilities,
        sparse_target.full_probabilities,
    )
    diagnostics = CompilerDiagnostics(
        compiler="sparse",
        n_total=n_total,
        dense_probability_values=1 << n_total,
        significant_modes=len(sparse_target.probabilities),
        occupied_prefix_nodes=occupied_prefix_nodes,
        fitted_prefix_nodes=fitted_prefix_nodes,
        tail_probability=sparse_target.tail_probability,
        cutoff_wavenumber=sparse_target.cutoff_wavenumber,
        locality_radius=locality_radius,
        compile_time_s=elapsed,
        baseline_rss_bytes=memory.baseline,
        peak_rss_bytes=memory.peak,
        incremental_peak_memory_bytes=memory.incremental_peak,
        retained_data_bytes=retained_bytes,
    )
    return SparseCompilationResult(
        target=sparse_target,
        amp_params=amp_params,
        diagnostics=diagnostics,
    )


def ZZ(qc, control, target, gamma):
    """Apply the original CX-Rz-CX implementation of a ZZ interaction."""
    qc.cx(control, target)
    qc.rz(gamma, target)
    qc.cx(control, target)


def build_3d_coupled_circuit(nx, ny, nz, params, seed=None):
    """Construct the three-dimensional coupled state-preparation circuit."""
    rng = np.random.default_rng(seed)
    n_total = nx + ny + nz
    qr = QuantumRegister(n_total, "q")
    qc = QuantumCircuit(qr)

    z_qubits = [qr[i] for i in range(nx + ny + nz - 1, nx + ny - 1, -1)]
    y_qubits = [qr[i] for i in range(nx + ny - 1, nx - 1, -1)]
    x_qubits = [qr[i] for i in range(nx - 1, -1, -1)]
    all_ordered_qubits = z_qubits + y_qubits + x_qubits

    for layer_idx, (bias, weights) in enumerate(params):
        if layer_idx >= len(all_ordered_qubits):
            break
        target_qubit = all_ordered_qubits[layer_idx]
        qc.ry(bias, target_qubit)
        for i, weight in enumerate(weights):
            if np.abs(weight) > 1e-3:
                control_qubit = all_ordered_qubits[i]
                qc.cry(weight, control_qubit, target_qubit)

    phi_rand = rng.uniform(0, 2 * np.pi, n_total)
    for i in range(n_total):
        qc.rz(phi_rand[i], qr[i])

    # Retained alternative nearest-neighbor phase-scrambling construction.
    # gamma_rand = rng.uniform(0, 2*np.pi, n_total-1)
    # for i in range(n_total-1):
    #     qc.cx(qr[i], qr[i+1])
    #     qc.rz(gamma_rand[i], qr[i+1])
    #     qc.cx(qr[i], qr[i+1])

    for i in range(n_total):
        for j in range(i + 1, n_total):
            gamma_ij = rng.uniform(0, 2 * np.pi)
            qc.cx(qr[i], qr[j])
            qc.rz(gamma_ij, qr[j])
            qc.cx(qr[i], qr[j])

    return qc


def get_radial_avg_fast(k_vals, amp_vals, bins=25):
    """Compute the radial root-mean-square amplitude with vectorized bins."""
    min_k, max_k = np.min(k_vals), np.max(k_vals)
    if min_k == max_k:
        return np.array([min_k]), np.array([np.mean(amp_vals**2)])

    k_bins = np.logspace(np.log10(min_k), np.log10(max_k), bins + 1)
    inds = np.digitize(k_vals, k_bins)
    n_bins = len(k_bins)

    sq_amp = amp_vals**2
    bin_sums = np.bincount(inds, weights=sq_amp, minlength=n_bins + 1)
    bin_counts = np.bincount(inds, minlength=n_bins + 1)

    centers = []
    avgs = []
    for i in range(1, len(k_bins)):
        if bin_counts[i] > 0:
            low, high = k_bins[i - 1], k_bins[i]
            centers.append(np.sqrt(low * high))
            avgs.append(np.sqrt(bin_sums[i] / bin_counts[i]))

    return np.array(centers), np.array(avgs)


def simulate_fast(qc, simulator):
    """Compile and simulate a circuit with the statevector backend."""
    qc_sim = qc.copy()
    qc_sim.save_statevector()
    qc_compiled = transpile(qc_sim, simulator, optimization_level=3)
    result = simulator.run(qc_compiled).result()
    return result.get_statevector().data


def _build_weights_matrix(amp_params):
    """Convert layerwise affine parameters to a square weight matrix."""
    num_qubits = len(amp_params)
    weights_matrix = np.zeros((num_qubits, num_qubits))
    for target_idx, (_, weights) in enumerate(amp_params):
        if len(weights) > 0:
            weights_matrix[target_idx, : len(weights)] = weights
    return weights_matrix


def _empty_plot_arrays():
    """Return empty placeholders when full-grid plotting is disabled."""
    return tuple(np.empty(0, dtype=np.float64) for _ in range(4))


def _prepare_full_plot_data(statevector, nx, ny, nz):
    """Prepare the original full-grid amplitude arrays for plotting."""
    amp_sim_gray = np.abs(statevector)
    amp_sim_phys = restore_amps_from_gray(amp_sim_gray, nx, ny, nz)

    kx = np.fft.fftfreq(1 << nx, d=1 / (1 << nx))
    ky = np.fft.fftfreq(1 << ny, d=1 / (1 << ny))
    kz = np.fft.fftfreq(1 << nz, d=1 / (1 << nz))
    k_squared = (
        kz[:, None, None] ** 2
        + ky[None, :, None] ** 2
        + kx[None, None, :] ** 2
    )
    k_flat = np.sqrt(k_squared, out=k_squared).reshape(-1)
    mask = k_flat > 0
    k_plot = k_flat[mask]
    amp_plot = amp_sim_phys[mask]
    c_sim, v_sim = get_radial_avg_fast(k_plot, amp_plot)
    return amp_sim_phys, k_plot, amp_plot, c_sim, v_sim


def _print_compiler_diagnostics(diagnostics):
    """Print a compact summary suitable for copying into benchmark notes."""
    gib = 1024.0**3
    print(f"Compiler: {diagnostics.compiler}")
    print(f"Compilation time: {diagnostics.compile_time_s:.6f} s")
    if diagnostics.compiler == "sparse":
        print(f"Significant modes: {diagnostics.significant_modes:,}")
        print(
            "Occupied prefix-tree nodes: "
            f"{diagnostics.occupied_prefix_nodes:,}"
        )
        print(f"Estimated tail probability: {diagnostics.tail_probability:.6e}")


def run_quantum_encoding(
    config=None,
    compiler="dense",
    tail_tolerance=1.0e-12,
    locality_radius=None,
    simulate_statevectors=True,
    prepare_plot_data=True,
    report_circuit_depth=True,
):
    """Compile and optionally simulate the geometric-encoding circuit.

    Parameters
    ----------
    compiler : {"dense", "sparse"}
        ``dense`` preserves the original reference implementation.
        ``sparse`` evaluates conditional probabilities only on the retained
        spectral support and its occupied prefix tree.
    tail_tolerance : float
        Upper target for the probability discarded by sparse truncation.
    locality_radius : int or None
        Maximum backward bit distance ``r`` in the affine ansatz. ``None``
        includes all precursor qubits. Use the same value for dense and sparse
        runs when isolating the effect of sparse compilation.
    simulate_statevectors : bool
        If false, stop after classical amplitude-parameter compilation.
    prepare_plot_data : bool
        If false, avoid additional full-grid arrays used only by amplitude
        plots. This is recommended for compiler benchmarks.
    report_circuit_depth : bool
        Transpile one circuit and print its depth when statevectors are run.
    """
    config = WorkflowConfig() if config is None else config
    config.ensure_output_directories()

    nx, ny, nz = config.nx, config.ny, config.nz
    compiler = str(compiler).lower()
    if compiler == "dense":
        compilation = compile_dense_amplitude_parameters(
            config,
            locality_radius=locality_radius,
        )
        target_probs_phys = compilation.target_probs_phys
        target_probs_gray = compilation.target_probs_gray
        KZ, KY, KX = compilation.KZ, compilation.KY, compilation.KX
        sparse_target = None
    elif compiler == "sparse":
        compilation = compile_sparse_amplitude_parameters(
            config,
            tail_tolerance=tail_tolerance,
            locality_radius=locality_radius,
        )
        target_probs_phys = None
        target_probs_gray = None
        KZ = KY = KX = None
        sparse_target = compilation.target
    else:
        raise ValueError("compiler must be either 'dense' or 'sparse'.")

    amp_params = compilation.amp_params
    diagnostics = compilation.diagnostics
    _print_compiler_diagnostics(diagnostics)

    weights_matrix = _build_weights_matrix(amp_params)
    np.save(
        config.data_dir / f"weights_matrix_{compiler}.npy",
        weights_matrix,
    )
    if compiler == "dense":
        # Preserve the filename used by the original notebook.
        np.save(config.data_dir / "weights_matrix.npy", weights_matrix)

    if not simulate_statevectors:
        amp_sim_phys = None
        k_plot, amp_plot, c_sim, v_sim = _empty_plot_arrays()
        return EncodingResult(
            target_probs_phys=target_probs_phys,
            target_probs_gray=target_probs_gray,
            KZ=KZ,
            KY=KY,
            KX=KX,
            amp_params=amp_params,
            sta_vec_1_data=None,
            sta_vec_2_data=None,
            weights_matrix=weights_matrix,
            amp_sim_phys=amp_sim_phys,
            k_plot=k_plot,
            amp_plot=amp_plot,
            c_sim=c_sim,
            v_sim=v_sim,
            compiler_diagnostics=diagnostics,
            sparse_target=sparse_target,
        )

    print("Simulating circuits...")
    t_sim_start = time.time()
    qc1 = build_3d_coupled_circuit(
        nx,
        ny,
        nz,
        amp_params,
        seed=config.seed_spin_up,
    )
    qc2 = build_3d_coupled_circuit(
        nx,
        ny,
        nz,
        amp_params,
        seed=config.seed_spin_down,
    )

    simulator = AerSimulator(
        method="statevector",
        max_parallel_threads=config.max_parallel_threads,
    )
    sta_vec_1_data = simulate_fast(qc1, simulator)
    sta_vec_2_data = simulate_fast(qc2, simulator)
    print(f"   Simulation done in {time.time() - t_sim_start:.4f} s")

    if report_circuit_depth:
        print("Transpiling...")
        t_trans = time.time()
        qc_transpiled = transpile(
            qc1,
            basis_gates=["u", "cz"],
            optimization_level=3,
        )
        print(f"Transpile done in {time.time() - t_trans:.4f} s")
        print(f"Transpiled circuit depth: {qc_transpiled.depth()}")
        del qc_transpiled

    del qc1, qc2
    gc.collect()

    if prepare_plot_data:
        (
            amp_sim_phys,
            k_plot,
            amp_plot,
            c_sim,
            v_sim,
        ) = _prepare_full_plot_data(sta_vec_1_data, nx, ny, nz)
    else:
        amp_sim_phys = None
        k_plot, amp_plot, c_sim, v_sim = _empty_plot_arrays()

    return EncodingResult(
        target_probs_phys=target_probs_phys,
        target_probs_gray=target_probs_gray,
        KZ=KZ,
        KY=KY,
        KX=KX,
        amp_params=amp_params,
        sta_vec_1_data=sta_vec_1_data,
        sta_vec_2_data=sta_vec_2_data,
        weights_matrix=weights_matrix,
        amp_sim_phys=amp_sim_phys,
        k_plot=k_plot,
        amp_plot=amp_plot,
        c_sim=c_sim,
        v_sim=v_sim,
        compiler_diagnostics=diagnostics,
        sparse_target=sparse_target,
    )


def run_dense_quantum_encoding(config=None, **kwargs):
    """Run the original dense reference compiler and shared circuit backend."""
    return run_quantum_encoding(config, compiler="dense", **kwargs)


def run_sparse_quantum_encoding(config=None, **kwargs):
    """Run the sparse compiler and shared circuit backend."""
    return run_quantum_encoding(config, compiler="sparse", **kwargs)


def _normalized_state_fidelity(state_a, state_b):
    """Return the pure-state fidelity without assuming exact normalization."""
    state_a = np.asarray(state_a)
    state_b = np.asarray(state_b)
    norm_product = np.vdot(state_a, state_a).real * np.vdot(
        state_b,
        state_b,
    ).real
    if norm_product <= 0.0:
        return float("nan")
    overlap = np.vdot(state_a, state_b)
    return float(np.abs(overlap) ** 2 / norm_product)


def _relative_l2(reference_components, comparison_components):
    """Return a joint relative L2 difference for one or more fields."""
    numerator = 0.0
    denominator = 0.0
    for reference, comparison in zip(reference_components, comparison_components):
        difference = np.asarray(comparison) - np.asarray(reference)
        numerator += float(np.vdot(difference, difference).real)
        denominator += float(np.vdot(reference, reference).real)
    if denominator <= 0.0:
        return float("nan")
    return float(np.sqrt(numerator / denominator))


def _release_compiler_only_arrays(encoding):
    """Release arrays not required for statevector and flow comparisons."""
    encoding.target_probs_phys = None
    encoding.target_probs_gray = None
    encoding.KZ = None
    encoding.KY = None
    encoding.KX = None
    encoding.amp_sim_phys = None
    gc.collect()


def _compute_physical_comparison(dense_encoding, sparse_encoding, config):
    """Compare spectra, velocity, vorticity, and vortex-surface fields."""
    from .flow_reconstruction import (
        compute_Ek,
        compute_spin_fields,
        reconstruct_flow_field,
    )
    from .flow_statistics import compute_vorticity

    dense_flow = reconstruct_flow_field(dense_encoding, config)
    sparse_flow = reconstruct_flow_field(sparse_encoding, config)

    _, dense_energy = compute_Ek(
        dense_flow.ux,
        dense_flow.uy,
        dense_flow.uz,
        dense_flow.ik2,
        dense_flow.N,
    )
    _, sparse_energy = compute_Ek(
        sparse_flow.ux,
        sparse_flow.uy,
        sparse_flow.uz,
        sparse_flow.ik2,
        sparse_flow.N,
    )
    energy_error = _relative_l2((dense_energy,), (sparse_energy,))

    velocity_error = _relative_l2(
        (dense_flow.ux, dense_flow.uy, dense_flow.uz),
        (sparse_flow.ux, sparse_flow.uy, sparse_flow.uz),
    )

    dense_vorticity = compute_vorticity(
        dense_flow.ux,
        dense_flow.uy,
        dense_flow.uz,
        dense_flow.N,
        config.derivative_method,
        dense_flow.KX,
        dense_flow.KY,
        dense_flow.KZ,
    )
    sparse_vorticity = compute_vorticity(
        sparse_flow.ux,
        sparse_flow.uy,
        sparse_flow.uz,
        sparse_flow.N,
        config.derivative_method,
        sparse_flow.KX,
        sparse_flow.KY,
        sparse_flow.KZ,
    )
    vorticity_error = _relative_l2(dense_vorticity, sparse_vorticity)

    dense_spin = compute_spin_fields(
        dense_flow.psi1,
        dense_flow.psi2,
        dense_flow.rho,
    )
    sparse_spin = compute_spin_fields(
        sparse_flow.psi1,
        sparse_flow.psi2,
        sparse_flow.rho,
    )
    vortex_surface_error = _relative_l2(dense_spin, sparse_spin)

    maximum_physical_error = float(
        np.nanmax(
            [velocity_error, vorticity_error, vortex_surface_error]
        )
    )
    del dense_flow, sparse_flow, dense_vorticity, sparse_vorticity
    del dense_spin, sparse_spin
    gc.collect()
    return (
        energy_error,
        velocity_error,
        vorticity_error,
        vortex_surface_error,
        maximum_physical_error,
    )


def benchmark_dense_vs_sparse(
    config=None,
    tail_tolerance=1.0e-12,
    locality_radius=None,
    simulate_statevectors=True,
    compute_physical_diagnostics=True,
    output_csv=None,
):
    """Benchmark both compilers using identical circuits and phase seeds.

    The same ``locality_radius`` is passed to both implementations so that the
    comparison isolates sparse support construction rather than a change in
    the affine ansatz. Full-state simulation is optional because it retains
    the unavoidable O(2^n) classical validation cost.
    """
    config = WorkflowConfig() if config is None else config
    if compute_physical_diagnostics and not simulate_statevectors:
        raise ValueError(
            "compute_physical_diagnostics requires simulate_statevectors=True."
        )

    common_options = dict(
        locality_radius=locality_radius,
        simulate_statevectors=simulate_statevectors,
        prepare_plot_data=False,
        report_circuit_depth=False,
    )
    dense_encoding = run_dense_quantum_encoding(config, **common_options)
    dense_diagnostics = dense_encoding.compiler_diagnostics
    _release_compiler_only_arrays(dense_encoding)

    sparse_encoding = run_sparse_quantum_encoding(
        config,
        tail_tolerance=tail_tolerance,
        **common_options,
    )
    sparse_diagnostics = sparse_encoding.compiler_diagnostics

    if simulate_statevectors:
        spin_up_fidelity = _normalized_state_fidelity(
            dense_encoding.sta_vec_1_data,
            sparse_encoding.sta_vec_1_data,
        )
        spin_down_fidelity = _normalized_state_fidelity(
            dense_encoding.sta_vec_2_data,
            sparse_encoding.sta_vec_2_data,
        )
        minimum_fidelity = min(spin_up_fidelity, spin_down_fidelity)
    else:
        spin_up_fidelity = spin_down_fidelity = minimum_fidelity = float("nan")

    if compute_physical_diagnostics:
        (
            energy_error,
            velocity_error,
            vorticity_error,
            vortex_surface_error,
            maximum_physical_error,
        ) = _compute_physical_comparison(
            dense_encoding,
            sparse_encoding,
            config,
        )
    else:
        energy_error = float("nan")
        velocity_error = float("nan")
        vorticity_error = float("nan")
        vortex_surface_error = float("nan")
        maximum_physical_error = float("nan")

    dense_memory = dense_diagnostics.incremental_peak_memory_bytes
    sparse_memory = sparse_diagnostics.incremental_peak_memory_bytes
    if dense_memory > 0:
        memory_reduction = 1.0 - sparse_memory / dense_memory
    else:
        memory_reduction = float("nan")
    if sparse_diagnostics.compile_time_s > 0.0:
        speedup = (
            dense_diagnostics.compile_time_s
            / sparse_diagnostics.compile_time_s
        )
    else:
        speedup = float("inf")

    benchmark = EncodingBenchmarkResult(
        n_total=config.n_total,
        dense_compile_time_s=dense_diagnostics.compile_time_s,
        sparse_compile_time_s=sparse_diagnostics.compile_time_s,
        compilation_speedup=speedup,
        dense_incremental_peak_memory_bytes=dense_memory,
        sparse_incremental_peak_memory_bytes=sparse_memory,
        peak_memory_reduction_fraction=memory_reduction,
        dense_probability_values=dense_diagnostics.dense_probability_values,
        sparse_significant_modes=sparse_diagnostics.significant_modes,
        sparse_occupied_prefix_nodes=sparse_diagnostics.occupied_prefix_nodes,
        sparse_fitted_prefix_nodes=sparse_diagnostics.fitted_prefix_nodes,
        sparse_tail_probability=sparse_diagnostics.tail_probability,
        spin_up_state_fidelity=spin_up_fidelity,
        spin_down_state_fidelity=spin_down_fidelity,
        minimum_state_fidelity=minimum_fidelity,
        energy_spectrum_relative_l2_error=energy_error,
        velocity_relative_l2_error=velocity_error,
        vorticity_relative_l2_error=vorticity_error,
        vortex_surface_relative_l2_error=vortex_surface_error,
        maximum_physical_relative_l2_error=maximum_physical_error,
    )

    output_csv = (
        config.data_dir / f"dense_sparse_benchmark_n={config.n_total}.csv"
        if output_csv is None
        else Path(output_csv)
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([benchmark.as_dict()]).to_csv(output_csv, index=False)

    gib = 1024.0**3
    print("\nDense-versus-sparse benchmark")
    print("-" * 72)
    print(f"n = {benchmark.n_total}")
    print(
        f"Compilation time: dense={benchmark.dense_compile_time_s:.6f} s, "
        f"sparse={benchmark.sparse_compile_time_s:.6f} s, "
        f"speedup={benchmark.compilation_speedup:.3f}x"
    )
    print(
        f"Dense probability values: {benchmark.dense_probability_values:,}"
    )
    print(
        f"Sparse significant modes: {benchmark.sparse_significant_modes:,}"
    )
    print(
        "Sparse occupied prefix nodes: "
        f"{benchmark.sparse_occupied_prefix_nodes:,}"
    )
    print(f"Minimum state fidelity: {benchmark.minimum_state_fidelity:.12e}")
    print(
        "Energy-spectrum relative L2 error: "
        f"{benchmark.energy_spectrum_relative_l2_error:.12e}"
    )
    print(
        "Maximum velocity/vorticity/VSF relative L2 error: "
        f"{benchmark.maximum_physical_relative_l2_error:.12e}"
    )
    print(f"Benchmark saved successfully: {output_csv}")
    return benchmark


def benchmark_dense_sparse_range(
    bits_per_axis,
    base_config=None,
    tail_tolerance=1.0e-12,
    locality_radius=None,
    simulate_statevectors=True,
    compute_physical_diagnostics=True,
    output_csv=None,
):
    """Run the dense-sparse benchmark for several isotropic grid sizes."""
    base_config = WorkflowConfig() if base_config is None else base_config
    results = []
    for axis_bits in bits_per_axis:
        axis_bits = int(axis_bits)
        current_config = replace(
            base_config,
            nx=axis_bits,
            ny=axis_bits,
            nz=axis_bits,
        )
        result = benchmark_dense_vs_sparse(
            current_config,
            tail_tolerance=tail_tolerance,
            locality_radius=locality_radius,
            simulate_statevectors=simulate_statevectors,
            compute_physical_diagnostics=compute_physical_diagnostics,
            output_csv=(
                current_config.data_dir
                / f"dense_sparse_benchmark_n={current_config.n_total}.csv"
            ),
        )
        results.append(result)

    output_csv = (
        base_config.data_dir / "dense_sparse_benchmark_summary.csv"
        if output_csv is None
        else Path(output_csv)
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([result.as_dict() for result in results]).to_csv(
        output_csv,
        index=False,
    )
    print(f"Benchmark range summary saved successfully: {output_csv}")
    return results
