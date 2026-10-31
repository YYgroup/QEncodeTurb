"""Spectral, statistical, geometric, and turbulence post-processing."""

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import gc
import math
import warnings
from pathlib import Path
import json

from .config import WorkflowConfig
from .flow_reconstruction import (
    FlowFieldResult,
    compute_Ek,
    compute_Es,
    compute_spin_fields,
)


@dataclass
class VorticityResult:
    """Vorticity, helicity, and their scalar magnitudes."""

    vorx: np.ndarray
    vory: np.ndarray
    vorz: np.ndarray
    magnitude: np.ndarray
    helicity: np.ndarray


@dataclass
class EnergyFlux:
    """Shell transfer and cumulative spectral energy flux."""

    k: np.ndarray
    transfer: np.ndarray
    energy_flux: np.ndarray
    nonlinear_balance: float


@dataclass
class RQResult:
    """Joint probability data for the velocity-gradient invariants."""

    R: np.ndarray
    Q: np.ndarray
    sigma_A: float = np.nan
    val_R_star: Optional[np.ndarray] = None
    val_Q_star: Optional[np.ndarray] = None
    pdf_star: Optional[np.ndarray] = None


def compute_pdf(variable, Num, bin_num):
    """Compute a uniformly binned probability density function."""
    bin_size = (np.amax(variable) - np.amin(variable)) / bin_num
    counts, bin_edges = np.histogram(
        variable,
        bins=bin_num,
        range=(np.amin(variable), np.amax(variable)),
    )
    pdf = counts / (Num * bin_size)
    val = bin_edges[:-1] + bin_size / 2
    return val, pdf


def compute_vorticity(ux, uy, uz, N, switch, KX, KY, KZ):
    """Compute vorticity with the original SM or FDM scheme."""
    if switch == "SM":
        ux_spec = np.fft.fftn(ux)
        uy_spec = np.fft.fftn(uy)
        uz_spec = np.fft.fftn(uz)
        vorx = np.real(np.fft.ifftn(1j * (KY * uz_spec - KZ * uy_spec)))
        vory = np.real(np.fft.ifftn(1j * (KZ * ux_spec - KX * uz_spec)))
        vorz = np.real(np.fft.ifftn(1j * (KX * uy_spec - KY * ux_spec)))
    elif switch == "FDM":
        h = 2 * np.pi / N
        vorx = (
            np.roll(uz, -1, axis=1)
            - np.roll(uz, 1, axis=1)
            - np.roll(uy, -1, axis=2)
            + np.roll(uy, 1, axis=2)
        )
        vory = (
            np.roll(ux, -1, axis=2)
            - np.roll(ux, 1, axis=2)
            - np.roll(uz, -1, axis=0)
            + np.roll(uz, 1, axis=0)
        )
        vorz = (
            np.roll(uy, -1, axis=0)
            - np.roll(uy, 1, axis=0)
            - np.roll(ux, -1, axis=1)
            + np.roll(ux, 1, axis=1)
        )
        vorx = vorx / (2 * h)
        vory = vory / (2 * h)
        vorz = vorz / (2 * h)
    else:
        raise ValueError("switch must be 'SM' or 'FDM'.")
    return vorx, vory, vorz


def compute_velocity_gradient(ux, uy, uz, N, switch, KX, KY, KZ):
    """Compute all nine velocity-gradient components."""
    if switch == "SM":
        vgt11 = np.real(np.fft.ifftn(1j * KX * np.fft.fftn(ux)))
        vgt12 = np.real(np.fft.ifftn(1j * KY * np.fft.fftn(ux)))
        vgt13 = np.real(np.fft.ifftn(1j * KZ * np.fft.fftn(ux)))

        vgt21 = np.real(np.fft.ifftn(1j * KX * np.fft.fftn(uy)))
        vgt22 = np.real(np.fft.ifftn(1j * KY * np.fft.fftn(uy)))
        vgt23 = np.real(np.fft.ifftn(1j * KZ * np.fft.fftn(uy)))

        vgt31 = np.real(np.fft.ifftn(1j * KX * np.fft.fftn(uz)))
        vgt32 = np.real(np.fft.ifftn(1j * KY * np.fft.fftn(uz)))
        vgt33 = np.real(np.fft.ifftn(1j * KZ * np.fft.fftn(uz)))
    elif switch == "FDM":
        h = 2 * np.pi / N
        vgt11 = (np.roll(ux, -1, axis=2) - np.roll(ux, 1, axis=2)) / (2 * h)
        vgt12 = (np.roll(ux, -1, axis=1) - np.roll(ux, 1, axis=1)) / (2 * h)
        vgt13 = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) / (2 * h)

        vgt21 = (np.roll(uy, -1, axis=2) - np.roll(uy, 1, axis=2)) / (2 * h)
        vgt22 = (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)) / (2 * h)
        vgt23 = (np.roll(uy, -1, axis=0) - np.roll(uy, 1, axis=0)) / (2 * h)

        vgt31 = (np.roll(uz, -1, axis=2) - np.roll(uz, 1, axis=2)) / (2 * h)
        vgt32 = (np.roll(uz, -1, axis=1) - np.roll(uz, 1, axis=1)) / (2 * h)
        vgt33 = (np.roll(uz, -1, axis=0) - np.roll(uz, 1, axis=0)) / (2 * h)
    else:
        raise ValueError("switch must be 'SM' or 'FDM'.")

    return (
        vgt11,
        vgt12,
        vgt13,
        vgt21,
        vgt22,
        vgt23,
        vgt31,
        vgt32,
        vgt33,
    )

def compute_sigma_A(
    ux,
    uy,
    uz,
    N,
    KX,
    KY,
    KZ,
    switch="FDM",
    return_gradient=False,
):
    """Compute the rms magnitude of the traceless velocity-gradient tensor.

    The traceless velocity-gradient tensor is

        A'_ij = A_ij - (A_kk / 3) delta_ij,

    and its characteristic magnitude is

        sigma_A = <A'_ij A'_ij>^(1/2).

    Parameters
    ----------
    ux, uy, uz : numpy.ndarray
        Velocity components.
    N : int
        Number of grid points in each direction.
    KX, KY, KZ : numpy.ndarray
        Wavenumber arrays.
    switch : {"FDM", "SM"}, optional
        Differentiation method.
    return_gradient : bool, optional
        Return the nine components of A'_ij together with sigma_A.
        This option avoids recomputing the gradient in compute_RQ.

    Returns
    -------
    sigma_A : float
        Root-mean-square magnitude of A'_ij.
    gradient : tuple of numpy.ndarray, optional
        Nine components of A'_ij, returned only when
        return_gradient=True.
    """
    gradient = list(
        compute_velocity_gradient(
            ux,
            uy,
            uz,
            N,
            switch,
            KX,
            KY,
            KZ,
        )
    )

    # Remove the isotropic dilatational component from the diagonal
    trace_over_three = (gradient[0] + gradient[4] + gradient[8]) / 3.0
    gradient[0] = gradient[0] - trace_over_three
    gradient[4] = gradient[4] - trace_over_three
    gradient[8] = gradient[8] - trace_over_three
    del trace_over_three

    sigma_A_squared = 0.0
    for component in gradient:
        sigma_A_squared += np.mean(
            component * component,
            dtype=np.float64,
        )

    sigma_A = float(np.sqrt(sigma_A_squared))

    if not np.isfinite(sigma_A) or sigma_A <= 0.0:
        raise ValueError(
            "sigma_A must be finite and positive for R-Q normalization."
        )

    if return_gradient:
        return sigma_A, tuple(gradient)

    return sigma_A


def compute_RQ(
    ux,
    uy,
    uz,
    N,
    KX,
    KY,
    KZ,
    return_sigma_A=False,
    switch="FDM",
):
    """Compute invariants of the traceless velocity-gradient tensor."""
    sigma_A, gradient = compute_sigma_A(
        ux,
        uy,
        uz,
        N,
        KX,
        KY,
        KZ,
        switch=switch,
        return_gradient=True,
    )

    (
        a11,
        a12,
        a13,
        a21,
        a22,
        a23,
        a31,
        a32,
        a33,
    ) = gradient

    Q = -0.5 * (
        a11**2
        + a22**2
        + a33**2
        + 2.0 * a12 * a21
        + 2.0 * a13 * a31
        + 2.0 * a23 * a32
    )

    R = -(
        a11 * (a22 * a33 - a32 * a23)
        - a12 * (a21 * a33 - a31 * a23)
        + a13 * (a21 * a32 - a31 * a22)
    )

    if return_sigma_A:
        return R, Q, sigma_A

    return R, Q


def compute_jointpdf(data1, data2, bin_num):
    """
    Compute a two-dimensional probability density function using a robust
    percentile range to prevent rare extreme values from degrading the
    resolution of the central distribution.

    The returned PDF retains its physical normalization and is not divided
    by its maximum. Samples outside the displayed percentile range are
    excluded from the histogram but remain included in the normalization.
    """
    data1 = np.asarray(data1).ravel()
    data2 = np.asarray(data2).ravel()

    if data1.size != data2.size:
        raise ValueError("data1 and data2 must contain the same number of samples.")

    # Exclude 0.1% of samples from each tail when determining the bin range.
    tail_fraction = 1.0e-8

    max_range_samples = 5_000_000
    if data1.size > max_range_samples:
        rng = np.random.default_rng(1234)
        sample_indices = rng.integers(
            0,
            data1.size,
            size=max_range_samples,
        )
        sample1 = data1[sample_indices]
        sample2 = data2[sample_indices]
    else:
        sample1 = data1
        sample2 = data2

    finite_sample = np.isfinite(sample1) & np.isfinite(sample2)
    sample1 = sample1[finite_sample]
    sample2 = sample2[finite_sample]

    lower1, upper1 = np.quantile(
        sample1,
        [tail_fraction, 1.0 - tail_fraction],
    )
    lower2, upper2 = np.quantile(
        sample2,
        [tail_fraction, 1.0 - tail_fraction],
    )

    # R is expected to be approximately symmetric about zero.
    # A symmetric horizontal range facilitates comparison between panels.
    limit1 = max(abs(lower1), abs(upper1))
    lower1, upper1 = -limit1, limit1

    edges1 = np.linspace(lower1, upper1, bin_num + 1)
    edges2 = np.linspace(lower2, upper2, bin_num + 1)

    counts = np.zeros((bin_num, bin_num), dtype=np.float64)
    total_valid_samples = 0

    # Chunked accumulation avoids large temporary arrays for the 1024^3 case.
    chunk_size = 1_000_000_000

    for start in range(0, data1.size, chunk_size):
        stop = min(start + chunk_size, data1.size)

        chunk1 = data1[start:stop]
        chunk2 = data2[start:stop]

        finite = np.isfinite(chunk1) & np.isfinite(chunk2)
        chunk1 = chunk1[finite]
        chunk2 = chunk2[finite]

        total_valid_samples += chunk1.size

        chunk_counts, _, _ = np.histogram2d(
            chunk1,
            chunk2,
            bins=(edges1, edges2),
        )
        counts += chunk_counts

    bin_width1 = np.diff(edges1)[:, None]
    bin_width2 = np.diff(edges2)[None, :]

    # Normalize by all valid samples rather than only the retained samples.
    # The integral over the displayed range therefore equals the retained probability instead of being artificially rescaled to unity.
    pdf = counts / (
        total_valid_samples * bin_width1 * bin_width2
    )

    val1 = 0.5 * (edges1[:-1] + edges1[1:])
    val2 = 0.5 * (edges2[:-1] + edges2[1:])

    retained_fraction = counts.sum() / total_valid_samples
    print(
        f"Joint PDF range retains "
        f"{100.0 * retained_fraction:.3f}% of the valid samples."
    )

    return val1, val2, pdf


def calculate_structure_function(
    ux,
    uy,
    uz,
    p,
    n_samples=1_000_000,
    n_bins=50,
):
    """Estimate an absolute longitudinal moment using physical distances.

    Each box side has length 2*pi. Sampling retains the original nonperiodic
    point-pair definition. Empty bins contain NaN.
    """
    if not (ux.shape == uy.shape == uz.shape):
        raise ValueError("ux, uy, and uz must have identical shapes.")

    if p <= 0:
        raise ValueError("The structure-function order p must be positive.")

    Nx, Ny, Nz = ux.shape

    box_length = 2 * np.pi
    spacing = box_length / np.array([Nx, Ny, Nz], dtype=float)
    reference_spacing = float(spacing.min())
    spacing_ratio = spacing / reference_spacing
    uniform_spacing = Nx == Ny == Nz

    print(
        f"Computing the order-{p} structure function with "
        f"{n_samples:,} sampled point pairs..."
    )

    batch_size = 200_000

    coordinate_dtype = np.min_scalar_type(max(Nx, Ny, Nz) - 1)
    coordinates = np.empty(
        (6, n_samples),
        dtype=coordinate_dtype,
    )

    # Preserve the original coordinate draw order and random-number dtype.
    for axis, size in enumerate((Nx, Ny, Nz, Nx, Ny, Nz)):
        for start in range(0, n_samples, batch_size):
            stop = min(start + batch_size, n_samples)

            coordinates[axis, start:stop] = np.random.randint(
                0,
                size,
                size=stop - start,
            )

    # Internal distances use reference-grid units.
    if uniform_spacing:
        max_r = np.sqrt(Nx**2 + Ny**2 + Nz**2)
    else:
        max_r = np.linalg.norm(
            np.array([Nx, Ny, Nz]) * spacing_ratio
        )

    r_bins = np.logspace(
        np.log10(1.0),
        np.log10(max_r),
        n_bins + 1,
    )

    # Convert returned distances to physical length.
    r_centers = (
        0.5
        * (r_bins[:-1] + r_bins[1:])
        * reference_spacing
    )

    bin_sums = np.zeros(n_bins, dtype=np.float64)
    bin_counts = np.zeros(n_bins, dtype=np.int64)

    max_squared_distance = (
        (Nx - 1)**2
        + (Ny - 1)**2
        + (Nz - 1)**2
    )
    use_lookup = (
        uniform_spacing
        and max_squared_distance <= 8_000_000
        and n_samples > max_squared_distance
    )

    if use_lookup:
        radius_lookup = np.sqrt(
            np.arange(
                max_squared_distance + 1,
                dtype=np.float64,
            )
        )

        bin_lookup = np.searchsorted(
            r_bins,
            radius_lookup,
            side="right",
        ) - 1

        bin_lookup = bin_lookup.astype(
            np.min_scalar_type(-n_bins - 1)
        )

    velocity_dtype = np.result_type(
        ux.dtype,
        uy.dtype,
        uz.dtype,
    )
    flat_velocities = None

    if all(u.flags.c_contiguous for u in (ux, uy, uz)):
        flat_velocities = (
            ux.reshape(-1),
            uy.reshape(-1),
            uz.reshape(-1),
        )

    # Preserve the original per-bin floating-point addition order.
    sum_indices = np.empty(
        n_bins + batch_size,
        dtype=np.intp,
    )
    sum_weights = np.empty(
        n_bins + batch_size,
        dtype=np.float64,
    )
    sum_indices[:n_bins] = np.arange(n_bins)

    for start in range(0, n_samples, batch_size):
        stop = min(start + batch_size, n_samples)

        i1, j1, k1, i2, j2, k2 = coordinates[:, start:stop]

        if use_lookup:
            dr_x = i2.astype(np.int32) - i1
            dr_y = j2.astype(np.int32) - j1
            dr_z = k2.astype(np.int32) - k1

            squared_index = dr_x**2 + dr_y**2 + dr_z**2

            r = radius_lookup[squared_index]
            bin_indices = bin_lookup[squared_index]
        else:
            dr_x = (
                (i2.astype(float) - i1.astype(float))
                * spacing_ratio[0]
            )
            dr_y = (
                (j2.astype(float) - j1.astype(float))
                * spacing_ratio[1]
            )
            dr_z = (
                (k2.astype(float) - k1.astype(float))
                * spacing_ratio[2]
            )

            r = np.sqrt(dr_x**2 + dr_y**2 + dr_z**2)

            bin_indices = np.searchsorted(
                r_bins,
                r,
                side="right",
            ) - 1

        valid = (
            (r > 0)
            & (bin_indices >= 0)
            & (bin_indices < n_bins)
        )

        if not np.any(valid):
            continue

        if not np.all(valid):
            i1, j1, k1 = i1[valid], j1[valid], k1[valid]
            i2, j2, k2 = i2[valid], j2[valid], k2[valid]

            dr_x = dr_x[valid]
            dr_y = dr_y[valid]
            dr_z = dr_z[valid]

            r = r[valid]
            bin_indices = bin_indices[valid]

        if flat_velocities is not None:
            first = i1.astype(np.intp)
            first *= Ny
            first += j1
            first *= Nz
            first += k1

            second = i2.astype(np.intp)
            second *= Ny
            second += j2
            second *= Nz
            second += k2

            du = np.empty(
                (3, len(r)),
                dtype=velocity_dtype,
            )

            for axis, velocity in enumerate(flat_velocities):
                np.subtract(
                    velocity[second],
                    velocity[first],
                    out=du[axis],
                    dtype=velocity_dtype,
                )
        else:
            u1 = np.array(
                [
                    ux[i1, j1, k1],
                    uy[i1, j1, k1],
                    uz[i1, j1, k1],
                ]
            )
            u2 = np.array(
                [
                    ux[i2, j2, k2],
                    uy[i2, j2, k2],
                    uz[i2, j2, k2],
                ]
            )

            du = u2 - u1

        # A common length factor cancels in the longitudinal projection.
        du_long = np.einsum(
            "ij,ij->j",
            du,
            np.array([dr_x, dr_y, dr_z], dtype=float),
        ) / r

        sp_values = np.abs(du_long) ** p
        number = len(sp_values)

        sum_weights[:n_bins] = bin_sums
        sum_weights[n_bins:n_bins + number] = sp_values
        sum_indices[n_bins:n_bins + number] = bin_indices

        bin_sums = np.bincount(
            sum_indices[:n_bins + number],
            weights=sum_weights[:n_bins + number],
            minlength=n_bins,
        )

        bin_counts += np.bincount(
            bin_indices.astype(np.intp, copy=False),
            minlength=n_bins,
        )

    # Missing statistics must not appear as a zero structure function.
    Sp = np.full_like(r_centers, np.nan)
    populated = bin_counts > 0

    Sp[populated] = (
        bin_sums[populated]
        / bin_counts[populated]
    )

    return r_centers.reshape(-1), Sp.reshape(-1)


def run_structure_functions(flow, config=None):
    """Compute and export order-2 through order-5 structure functions."""
    config = WorkflowConfig() if config is None else config

    r, s2 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        2,
        n_samples=100_000_000,
        n_bins=100,
    )
    r, s3 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        3,
        n_samples=100_000_000,
        n_bins=100,
    )
    r, s4 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        4,
        n_samples=100_000_000,
        n_bins=100,
    )
    r, s5 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        5,
        n_samples=100_000_000,
        n_bins=100,
    )

    header_names = ["r", "s2", "s3", "s4", "s5"]
    data_dict = dict(zip(header_names, [r, s2, s3, s4, s5]))

    df = pd.DataFrame(data_dict)

    df = df.dropna(
        subset=["r", "s2", "s3", "s4", "s5"],
    ).reset_index(drop=True)

    from scipy.ndimage import median_filter
    from scipy.interpolate import make_smoothing_spline

    smoothing_r = np.geomspace(0.01, 10, 87)
    smoothing_lambda = 0.02
    source_r = df["r"].to_numpy()
    log_r = np.log(smoothing_r)
    smoothed_df = pd.DataFrame({"r": smoothing_r})
    for column in ["s2", "s3", "s4", "s5"]:
        log_sp = np.interp(log_r, np.log(source_r), np.log(df[column].to_numpy()))
        log_sp = median_filter(log_sp, size=15, mode="nearest")
        spline = make_smoothing_spline(log_r, log_sp, lam=smoothing_lambda)
        smoothed_df[column] = np.exp(spline(log_r))
    df = smoothed_df

    r = df["r"].to_numpy()
    s2 = df["s2"].to_numpy()
    s3 = df["s3"].to_numpy()
    s4 = df["s4"].to_numpy()
    s5 = df["s5"].to_numpy()

    df.to_csv(
        config.data_dir
        / f"structure_functions_n={config.nx + config.ny + config.nz}.csv",
        index=False,
        lineterminator="\n",
    )

    return r, s2, s3, s4, s5


def calculate_anisotropy_tensor(ux, uy, uz):
    """Return the six independent Reynolds-stress anisotropy components."""
    ux_prime = ux - np.mean(ux)
    uy_prime = uy - np.mean(uy)
    uz_prime = uz - np.mean(uz)

    two_K = np.mean(ux_prime**2 + uy_prime**2 + uz_prime**2)
    if two_K == 0:
        return {
            "b11": 0.0,
            "b22": 0.0,
            "b33": 0.0,
            "b12": 0.0,
            "b13": 0.0,
            "b23": 0.0,
        }

    R11 = np.mean(ux_prime * ux_prime)
    R22 = np.mean(uy_prime * uy_prime)
    R33 = np.mean(uz_prime * uz_prime)
    R12 = np.mean(ux_prime * uy_prime)
    R13 = np.mean(ux_prime * uz_prime)
    R23 = np.mean(uy_prime * uz_prime)

    return {
        "b11": R11 / two_K - 1.0 / 3.0,
        "b22": R22 / two_K - 1.0 / 3.0,
        "b33": R33 / two_K - 1.0 / 3.0,
        "b12": R12 / two_K,
        "b13": R13 / two_K,
        "b23": R23 / two_K,
    }


def export_energy_spectrum(flow, config=None):
    """Compute and export the velocity energy spectrum."""
    config = WorkflowConfig() if config is None else config
    k, Ek = compute_Ek(flow.ux, flow.uy, flow.uz, flow.ik2, flow.N)
    pd.DataFrame({"k": k, "Ek": Ek}).to_csv(
        config.data_dir / f"turb_Ek_nq={config.nx + config.ny + config.nz}.csv",
        index=False,
    )
    return k, Ek


def compute_energy_flux(flow, config=None, dealias=True):
    """Compute the scale-dependent spectral kinetic-energy flux.

    The modal transfer due to the advective nonlinearity is
    T(k) = -Re[u_i^*(k) Fourier(u_j partial_j u_i)(k)]. The forward
    energy flux through a spherical cutoff K is defined as
    Pi(K) = -sum_{|k| <= K} T(k), so positive values denote transfer from
    large to small scales. The default two-thirds truncation removes aliasing
    from the pseudospectral evaluation of the quadratic nonlinearity.

    For a non-solenoidal velocity field, this quantity represents the
    velocity-based advective flux rather than the complete compressible
    kinetic-energy budget, which additionally involves density and pressure.
    """
    config = WorkflowConfig() if config is None else config
    config.ensure_output_directories()

    velocity = (
        np.asarray(flow.ux, dtype=np.float64),
        np.asarray(flow.uy, dtype=np.float64),
        np.asarray(flow.uz, dtype=np.float64),
    )
    shape = velocity[0].shape

    if any(component.shape != shape for component in velocity):
        raise ValueError(
            "All velocity components must have identical shapes."
        )

    if len(shape) != 3 or not (
        shape[0] == shape[1] == shape[2]
    ):
        raise ValueError(
            "The energy-flux routine requires a cubic 3D grid."
        )

    Nx, Ny, Nz = shape

    kx = np.fft.fftfreq(Nx) * Nx
    ky = np.fft.fftfreq(Ny) * Ny
    kz = np.fft.fftfreq(Nz) * Nz

    wavevectors = (
        kx[:, None, None],
        ky[None, :, None],
        kz[None, None, :],
    )

    if dealias:
        cutoff = (min(shape) - 1) // 3

        dealias_mask = (
            (np.abs(wavevectors[0]) <= cutoff)
            & (np.abs(wavevectors[1]) <= cutoff)
            & (np.abs(wavevectors[2]) <= cutoff)
        )

        filtered_velocity = []

        for component in velocity:
            component_hat = np.fft.fftn(
                component,
                norm="forward",
            )
            component_hat[~dealias_mask] = 0.0

            filtered_component = np.fft.ifftn(
                component_hat,
                norm="forward",
            ).real.copy()

            filtered_velocity.append(filtered_component)
            del component_hat

        velocity = tuple(filtered_velocity)

    else:
        cutoff = min(shape) // 2
        dealias_mask = None

    k_squared = (wavevectors[0]**2 + wavevectors[1]**2 + wavevectors[2]**2)

    shell_index = np.rint(np.sqrt(k_squared)).astype(np.int32)

    valid_shell = shell_index <= cutoff
    flat_shell = shell_index[valid_shell]

    transfer = np.zeros(cutoff + 1, dtype=np.float64)
    nonlinear_balance = 0.0

    del k_squared, shell_index

    for component_index, component in enumerate(velocity):
        component_hat = np.fft.fftn(
            component,
            norm="forward",
        )

        if dealias:
            component_hat[~dealias_mask] = 0.0

        advection = np.zeros(
            shape,
            dtype=np.float64,
        )

        for advecting_velocity, wavevector in zip(
            velocity,
            wavevectors,
        ):
            derivative = np.fft.ifftn(
                1j * wavevector * component_hat,
                norm="forward",
            ).real

            advection += advecting_velocity * derivative
            del derivative

        advection_hat = np.fft.fftn(
            advection,
            norm="forward",
        )

        if dealias:
            advection_hat[~dealias_mask] = 0.0

        modal_transfer = -np.real(
            np.conj(component_hat) * advection_hat
        )

        if dealias:
            nonlinear_balance += float(
                np.sum(modal_transfer[dealias_mask])
            )
        else:
            nonlinear_balance += float(
                np.sum(modal_transfer)
            )

        transfer += np.bincount(
            flat_shell,
            weights=modal_transfer[valid_shell],
            minlength=cutoff + 1,
        )[: cutoff + 1]

        del (
            component_hat,
            advection,
            advection_hat,
            modal_transfer,
        )
        gc.collect()

        print("Energy-flux contribution completed for " f"velocity component {component_index + 1}/3.")

    k = np.arange(cutoff + 1, dtype=np.int32)

    energy_flux = -np.cumsum(transfer)

    output_path = (config.data_dir / f"energy_flux_N={Nx}.csv")

    pd.DataFrame(
        {
            "k": k[1:],
            "transfer": transfer[1:],
            "energy_flux": energy_flux[1:],
        }
    ).to_csv(
        output_path,
        index=False,
    )

    print(
        f"Energy flux saved successfully: {output_path}"
    )
    print(
        "Resolved nonlinear-transfer balance: "
        f"{nonlinear_balance:.6e}"
    )

    return EnergyFlux(
        k=k[1:],
        transfer=transfer[1:],
        energy_flux=energy_flux[1:],
        nonlinear_balance=nonlinear_balance,
    )


def export_spin_spectra(flow, config=None):
    """Compute the spin fields and export their spectra."""
    config = WorkflowConfig() if config is None else config
    s1, s2, s3 = compute_spin_fields(flow.psi1, flow.psi2, flow.rho)
    k, Es1 = compute_Es(s1, flow.ik2, flow.N)
    k, Es2 = compute_Es(s2, flow.ik2, flow.N)
    k, Es3 = compute_Es(s3, flow.ik2, flow.N)
    pd.DataFrame({"k": k, "Es1": Es1, "Es2": Es2, "Es3": Es3}).to_csv(
        config.data_dir / f"turb_Es_nq={config.nx + config.ny + config.nz}.csv",
        index=False,
    )
    return (s1, s2, s3), (k, Es1, Es2, Es3)


def export_density_pdf(flow, config=None):
    """Compute and export the density probability distribution."""
    config = WorkflowConfig() if config is None else config
    val_rho, pdf_rho = compute_pdf(flow.rho, flow.N**3, 100)
    pd.DataFrame({"rho": val_rho, "pdf": pdf_rho}).to_csv(
        config.data_dir / "turb_PDF_rho_P0=0.999.csv",
        index=False,
    )
    return val_rho, pdf_rho


def compute_vorticity_fields(flow, config=None):
    """Compute vorticity magnitude and helicity."""
    config = WorkflowConfig() if config is None else config
    vorx, vory, vorz = compute_vorticity(
        flow.ux,
        flow.uy,
        flow.uz,
        flow.N,
        config.derivative_method,
        flow.KX,
        flow.KY,
        flow.KZ,
    )
    magnitude = np.sqrt(vorx**2 + vory**2 + vorz**2)
    helicity = vorx * flow.ux + vory * flow.uy + vorz * flow.uz
    # print(np.amin(magnitude), np.amax(magnitude), np.mean(magnitude))
    # print(
    #     np.amin(np.abs(helicity)),
    #     np.amax(np.abs(helicity)),
    #     np.mean(np.abs(helicity)),
    # )
    return VorticityResult(vorx, vory, vorz, magnitude, helicity)


def export_velocity_vorticity_pdfs(flow, vorticity, config=None):
    """Compute bulk energies and export velocity and vorticity PDFs."""
    config = WorkflowConfig() if config is None else config
    u2 = flow.ux**2 + flow.uy**2 + flow.uz**2
    h = 2 * np.pi / flow.N
    tke = np.sum(u2) / 2 * h**3
    u_prime = np.sqrt(2 * tke / 3)
    print("Total kinetic energy:", tke)
    print("Mean velocity:", u_prime)

    enstrophy = np.sum(vorticity.magnitude**2) / 2 * h**3
    vor_prime = np.sqrt(2 * enstrophy / 3)
    print("Enstrophy:", enstrophy)

    vel = np.array([flow.ux, flow.uy, flow.uz])
    val_u, pdf_u = compute_pdf(vel, 3 * flow.N**3, 150)
    vor = np.array([vorticity.vorx, vorticity.vory, vorticity.vorz])
    val_vor, pdf_vor = compute_pdf(vor, 3 * flow.N**3, 150)

    pd.DataFrame({"vel": val_u / u_prime, "pdf": pdf_u}).to_csv(
        config.data_dir / f"turb_PDF_vel_n={config.nx + config.ny + config.nz}.csv",
        index=False,
    )
    pd.DataFrame({"vor": val_vor / vor_prime, "pdf": pdf_vor}).to_csv(
        config.data_dir / f"turb_PDF_vor_n={config.nx + config.ny + config.nz}.csv",
        index=False,
    )
    return {
        "tke": tke,
        "u_prime": u_prime,
        "enstrophy": enstrophy,
        "vor_prime": vor_prime,
        "val_u": val_u,
        "pdf_u": pdf_u,
        "val_vor": val_vor,
        "pdf_vor": pdf_vor,
    }


def compute_rq_statistics(flow, config=None):
    """Compute dimensional and normalized R-Q joint distributions."""
    config = WorkflowConfig() if config is None else config
    config.ensure_output_directories()

    R, Q, sigma_A = compute_RQ(
        flow.ux,
        flow.uy,
        flow.uz,
        flow.N,
        flow.KX,
        flow.KY,
        flow.KZ,
        return_sigma_A=True,
    )

    # Dimensional joint PDF retained for backward compatibility
    val_R_star, val_Q_star, pdf = compute_jointpdf(R / sigma_A**3, Q / sigma_A**2, 200)

    Nx, Ny, Nz = R.shape

    output_path = config.data_dir / f"RQ_joint_pdf_N={Nx}.npz"
    np.savez_compressed(
        output_path,
        sigma_A=sigma_A,
        val_R_star=val_R_star,
        val_Q_star=val_Q_star,
        pdf_star=pdf,
    )

    print(f"sigma_A = {sigma_A:.8e}")
    print(f"Dimensional and normalized R-Q data saved successfully: {output_path}")

    return RQResult(
        R=R,
        Q=Q,
        sigma_A=sigma_A,
        val_R_star=val_R_star,
        val_Q_star=val_Q_star,
        pdf_star=pdf,
    )


def report_anisotropy(flow, config=None):
    """Compute and print all six independent anisotropy components."""
    config = WorkflowConfig() if config is None else config
    config.ensure_output_directories()
    
    components = calculate_anisotropy_tensor(flow.ux, flow.uy, flow.uz)
    output_path = config.data_dir / f"anisotropy_n={config.nx + config.ny + config.nz}.xlsx"
    pd.DataFrame([components]).to_excel(output_path, index=False)
    return components


def compute_turbulence_scales(
    ux,
    uy,
    uz,
    *,
    L_box=2.0 * np.pi,
    k_L=1.0,
    k_eta_over_kmax=5.0,
    C_epsilon=1.0,
    nu=None,
    Re_L=None,
    component_axes=(0, 1, 2),
):
    r"""Compute the original scale-inferred turbulence diagnostics.

    The internally consistent scale-inferred definitions are

    u'                = sqrt(<u_x'^2 + u_y'^2 + u_z'^2>/3),
    eta               = 1/k_eta,
    Re_L              = C_epsilon^(-1/3) (k_eta/k_L)^(4/3),
    lambda_HIT/L      = sqrt[15/(C_epsilon Re_L)],
    Re_lambda,HIT     = sqrt(15 Re_L/C_epsilon).

    A separate kinematic Taylor scale is evaluated from the resolved velocity
    gradients as a consistency check. It is not combined with the
    scale-inferred viscosity to define the reported Re_lambda.
    """
    velocity = (np.asarray(ux), np.asarray(uy), np.asarray(uz))
    shape = velocity[0].shape
    if len(shape) != 3 or any(v.shape != shape for v in velocity):
        raise ValueError("ux, uy, and uz must be equal three-dimensional arrays.")
    if len(set(shape)) != 1:
        raise ValueError("The scale diagnostic assumes a cubic grid.")
    if sorted(component_axes) != [0, 1, 2]:
        raise ValueError("component_axes must be a permutation of (0, 1, 2).")
    if C_epsilon <= 0.0:
        raise ValueError("C_epsilon must be positive.")

    N = shape[0]
    dx = L_box / N

    # All Fourier modes are retained and no 2/3 dealiasing is applied.
    k_max = np.pi / dx
    k_eta_scale = k_eta_over_kmax * k_max
    eta = 1.0 / k_eta_scale
    kmax_eta = k_max * eta
    L_integral = 1.0 / k_L
    Re_L_scale = C_epsilon ** (-1.0 / 3.0) * (
        k_eta_scale / k_L
    ) ** (4.0 / 3.0)

    component_variances = []
    longitudinal_gradient_variances = []

    # One derivative is evaluated at a time to limit peak memory use.
    for field, axis in zip(velocity, component_axes):
        component_variances.append(float(np.var(field, dtype=np.float64)))
        derivative = (
            np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)
        ) / (2.0 * dx)
        longitudinal_gradient_variances.append(
            float(np.mean(np.asarray(derivative, dtype=np.float64) ** 2))
        )
        del derivative
        gc.collect()

    u_prime_sq = float(np.mean(component_variances))
    longitudinal_gradient_sq = float(np.mean(longitudinal_gradient_variances))
    if u_prime_sq <= 0.0:
        raise ValueError("The velocity fluctuation intensity is zero.")
    if longitudinal_gradient_sq <= 0.0:
        raise ValueError("The longitudinal velocity-gradient variance is zero.")

    u_prime = np.sqrt(u_prime_sq)
    lambda_resolved = np.sqrt(u_prime_sq / longitudinal_gradient_sq)

    if nu is not None and nu <= 0.0:
        raise ValueError("nu must be positive.")
    if Re_L is not None and Re_L <= 0.0:
        raise ValueError("Re_L must be positive.")

    lambda_taylor = L_integral * np.sqrt(15.0 / (C_epsilon * Re_L_scale))
    Re_lambda = np.sqrt(15.0 * Re_L_scale / C_epsilon)
    nu_scale = u_prime * L_integral / Re_L_scale
    epsilon_scale = C_epsilon * u_prime**3 / L_integral

    epsilon_resolved = 15.0 * nu_scale * longitudinal_gradient_sq
    resolved_dissipation_fraction = epsilon_resolved / epsilon_scale

    Re_L_from_nu = (
        u_prime * L_integral / float(nu) if nu is not None else np.nan
    )
    Re_L_input = float(Re_L) if Re_L is not None else np.nan

    result = {
        "N": N,
        "dx": dx,
        "k_max": k_max,
        "u_prime": u_prime,
        "nu_scale": nu_scale,
        "Re_L_scale": Re_L_scale,
        "C_epsilon": C_epsilon,
        "epsilon_scale": epsilon_scale,
        "eta": eta,
        "kmax_eta": kmax_eta,
        "lambda_taylor": lambda_taylor,
        "Re_lambda": Re_lambda,
        "lambda_resolved": lambda_resolved,
        "epsilon_resolved": epsilon_resolved,
        "resolved_dissipation_fraction": resolved_dissipation_fraction,
        "Re_L_from_user_nu": Re_L_from_nu,
        "Re_L_user_input": Re_L_input,
        "k_eta_scale": k_eta_scale,
        "component_variances": np.asarray(component_variances),
        "longitudinal_gradient_variances": np.asarray(
            longitudinal_gradient_variances
        ),
    }

    print("\nScale-inferred turbulence diagnostics")
    print("-" * 58)
    print(f"Grid size                         N = {N:d}")
    print(f"Maximum resolved wavenumber   k_max = {k_max:.8e}")
    print(f"k_eta = {k_eta_over_kmax:g} k_max          = {k_eta_scale:.8e}")
    print(f"Kolmogorov length       eta=1/k_eta = {eta:.8e}")
    print(f"Resolution parameter      k_max*eta = {kmax_eta:.8e}")
    print(f"Integral-scale Reynolds number Re_L = {Re_L_scale:.8e}")
    print(f"Effective viscosity        nu_scale = {nu_scale:.8e}")
    print(f"Taylor microscale       lambda_T,HIT = {lambda_taylor:.8e}")
    print(f"Taylor Reynolds number Re_lambda,HIT = {Re_lambda:.8e}")
    print("\nResolved-field consistency check")
    print("-" * 58)
    print(f"Kinematic Taylor scale lambda_resolved = {lambda_resolved:.8e}")
    print(
        "Resolved/extrapolated dissipation ratio = "
        f"{resolved_dissipation_fraction:.8e}"
    )
    print("The resolved-gradient quantities are not used to define Re_lambda.")
    if nu is not None:
        print(f"Re_L implied by user-supplied nu        = {Re_L_from_nu:.8e}")
    if Re_L is not None:
        print(f"User-supplied Re_L (comparison only)    = {Re_L_input:.8e}")

    return result


def compute_vartheta(
    ux,
    uy,
    uz,
    *,
    box_length=2 * np.pi,
    return_details=False,
):
    """
    Compute the dimensionless dilatation measure

        vartheta = <div(u)**2> / (<div(u)**2> + <|curl(u)|**2>).

    Periodic second-order central differences are used.

    Parameters
    ----------
    ux, uy, uz : real ndarray, shape (Nx, Ny, Nz)
        Velocity components on a uniform periodic grid.
    box_length : float or sequence of three floats
        Physical box lengths. A scalar specifies the same length
        in all three directions.
    return_details : bool
        If True, return a dictionary containing vartheta and
        the two mean-square quantities. Otherwise return vartheta.

    Returns
    -------
    float or dict
        The dimensionless measure, optionally with its components.
        Returns NaN when both mean-square quantities vanish.
    """
    fields = tuple(np.asarray(u) for u in (ux, uy, uz))
    ux, uy, uz = fields
    shape = ux.shape

    if (
        len(shape) != 3
        or min(shape) < 3
        or any(u.shape != shape for u in fields)
    ):
        raise ValueError(
            "ux, uy, and uz must have the same three-dimensional "
            "shape, with at least three points along each axis."
        )

    if any(
        not np.issubdtype(u.dtype, np.number)
        or np.iscomplexobj(u)
        for u in fields
    ):
        raise ValueError("Velocity arrays must contain real numbers.")

    lengths = np.broadcast_to(
        np.asarray(box_length, dtype=np.float64), (3,)
    )

    if not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("Box lengths must be finite and positive.")

    spacing = lengths / np.asarray(shape, dtype=np.float64)
    nx = shape[0]

    def derivative_slice(field, axis, i):
        """Return one x-slice of a periodic central derivative."""
        if axis == 0:
            forward = np.asarray(
                field[(i + 1) % nx], dtype=np.float64
            )
            backward = np.asarray(
                field[(i - 1) % nx], dtype=np.float64
            )
        else:
            plane = np.asarray(field[i], dtype=np.float64)
            forward = np.roll(plane, -1, axis=axis - 1)
            backward = np.roll(plane, 1, axis=axis - 1)

        return (forward - backward) / (2.0 * spacing[axis])

    sum_div2 = 0.0
    sum_curl2 = 0.0

    for i in range(nx):
        if any(not np.isfinite(u[i]).all() for u in fields):
            raise ValueError("Velocity contains NaN or infinity.")

        # div(u) = dux/dx + duy/dy + duz/dz.
        divergence = derivative_slice(ux, 0, i)
        divergence += derivative_slice(uy, 1, i)
        divergence += derivative_slice(uz, 2, i)

        sum_div2 += float(
            np.einsum("ij,ij->", divergence, divergence)
        )
        del divergence

        # curl(u): compute and release one component at a time.
        curl_terms = (
            (uz, 1, uy, 2),  # omega_x = duz/dy - duy/dz
            (ux, 2, uz, 0),  # omega_y = dux/dz - duz/dx
            (uy, 0, ux, 1),  # omega_z = duy/dx - dux/dy
        )

        for field_a, axis_a, field_b, axis_b in curl_terms:
            omega = derivative_slice(field_a, axis_a, i)
            omega -= derivative_slice(field_b, axis_b, i)

            sum_curl2 += float(
                np.einsum("ij,ij->", omega, omega)
            )
            del omega

    denominator = sum_div2 + sum_curl2

    if not np.isfinite(denominator):
        raise FloatingPointError(
            "Nonfinite derivative statistics; check velocity magnitudes."
        )

    if denominator == 0.0:
        warnings.warn(
            "Both derivative measures vanish; vartheta is undefined.",
            RuntimeWarning,
        )
        vartheta = np.nan
    else:
        vartheta = sum_div2 / denominator

    if return_details:
        n_points = ux.size
        return {
            "vartheta": float(vartheta),
            "mean_divergence_squared": sum_div2 / n_points,
            "mean_vorticity_squared": sum_curl2 / n_points,
        }

    return float(vartheta)