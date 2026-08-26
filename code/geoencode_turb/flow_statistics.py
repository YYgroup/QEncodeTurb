"""Spectral, statistical, geometric, and turbulence post-processing."""

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import gc

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


# def compute_jointpdf(data1, data2, bin_num):
#     """Compute a two-dimensional probability density function."""
#     binsize_1 = (np.max(data1) - np.min(data1)) / bin_num
#     binsize_2 = (np.max(data2) - np.min(data2)) / bin_num
#     counts, x_edges, y_edges = np.histogram2d(
#         data1.flatten(),
#         data2.flatten(),
#         bins=bin_num,
#     )
#     pdf = counts / (binsize_1 * binsize_2 * data1.size)
#     val1 = x_edges[:-1] + binsize_1 / 2
#     val2 = y_edges[:-1] + binsize_2 / 2
#     return val1, val2, pdf


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
    """Estimate a longitudinal structure function by Monte Carlo sampling.

    Parameters
    ----------
    ux, uy, uz : numpy.ndarray
        Three velocity components with shape (Nx, Ny, Nz). Unit grid spacing is
        assumed in every direction, as in the original notebook.
    p : int
        Structure-function order.
    n_samples : int, optional
        Number of randomly sampled point pairs.
    n_bins : int, optional
        Number of separation-distance bins.

    Returns
    -------
    r_centers, Sp : numpy.ndarray
        Separation-bin centers and the corresponding structure function.
    """
    if not (ux.shape == uy.shape == uz.shape):
        raise ValueError("ux, uy, and uz must have identical shapes.")
    if p <= 0:
        raise ValueError("The structure-function order p must be positive.")

    Nx, Ny, Nz = ux.shape
    print(f"Grid dimensions: {Nx}x{Ny}x{Nz}")
    print(
        f"Computing the order-{p} structure function with "
        f"{n_samples:,} sampled point pairs..."
    )

    i1 = np.random.randint(0, Nx, size=n_samples)
    j1 = np.random.randint(0, Ny, size=n_samples)
    k1 = np.random.randint(0, Nz, size=n_samples)
    i2 = np.random.randint(0, Nx, size=n_samples)
    j2 = np.random.randint(0, Ny, size=n_samples)
    k2 = np.random.randint(0, Nz, size=n_samples)

    dr_x = i2.astype(float) - i1.astype(float)
    dr_y = j2.astype(float) - j1.astype(float)
    dr_z = k2.astype(float) - k1.astype(float)
    r = np.sqrt(dr_x**2 + dr_y**2 + dr_z**2)

    # Remove zero-separation pairs.
    valid_mask = r > 0
    if not np.any(valid_mask):
        print("Warning: all sampled point pairs have zero separation.")
        return np.array([]).reshape(-1, 1), np.array([]).reshape(-1, 1)

    r = r[valid_mask]
    dr_x, dr_y, dr_z = dr_x[valid_mask], dr_y[valid_mask], dr_z[valid_mask]
    i1, j1, k1 = i1[valid_mask], j1[valid_mask], k1[valid_mask]
    i2, j2, k2 = i2[valid_mask], j2[valid_mask], k2[valid_mask]

    u1 = np.array([ux[i1, j1, k1], uy[i1, j1, k1], uz[i1, j1, k1]])
    u2 = np.array([ux[i2, j2, k2], uy[i2, j2, k2], uz[i2, j2, k2]])
    du = u2 - u1

    du_long = np.einsum(
        "ij,ij->j",
        du,
        np.array([dr_x, dr_y, dr_z]),
    ) / r
    sp_values = np.abs(du_long) ** p

    max_r = np.sqrt(Nx**2 + Ny**2 + Nz**2)
    r_bins = np.logspace(np.log10(1.0), np.log10(max_r), n_bins + 1)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    bin_indices = np.digitize(r, r_bins) - 1
    valid_bin_mask = (bin_indices >= 0) & (bin_indices < n_bins)

    bin_sums = np.bincount(
        bin_indices[valid_bin_mask],
        weights=sp_values[valid_bin_mask],
        minlength=n_bins,
    )
    bin_counts = np.bincount(
        bin_indices[valid_bin_mask],
        minlength=n_bins,
    )

    Sp = np.zeros_like(r_centers)
    non_zero_counts = bin_counts > 0
    Sp[non_zero_counts] = bin_sums[non_zero_counts] / bin_counts[non_zero_counts]
    return r_centers.reshape(-1), Sp.reshape(-1)


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
        config.data_dir / "turb_Es_nq=30.csv",
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
        config.data_dir / "turb_PDF_vel.csv",
        index=False,
    )
    pd.DataFrame({"vor": val_vor / vor_prime, "pdf": pdf_vor}).to_csv(
        config.data_dir / "turb_PDF_vor.csv",
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


def run_structure_functions(flow, config=None):
    """Compute and export the original order-2 through order-5 functions."""
    config = WorkflowConfig() if config is None else config
    r, s2 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        2,
        n_samples=1_000_000_000,
        n_bins=50,
    )
    r, s3 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        3,
        n_samples=1_000_000_000,
        n_bins=50,
    )
    r, s4 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        4,
        n_samples=1_000_000_000,
        n_bins=50,
    )
    r, s5 = calculate_structure_function(
        flow.ux,
        flow.uy,
        flow.uz,
        5,
        n_samples=1_000_000_000,
        n_bins=50,
    )

    header_names = ["r", "s2", "s3", "s4", "s5"]
    data_dict = dict(zip(header_names, [r, s2, s3, s4, s5]))
    df = pd.DataFrame(data_dict)
    df.to_csv(config.data_dir / "structure_functions.csv", index=False)
    return r, s2, s3, s4, s5


def report_anisotropy(flow):
    """Compute and print all six independent anisotropy components."""
    components = calculate_anisotropy_tensor(flow.ux, flow.uy, flow.uz)
    for component, value in components.items():
        print(f"{component}: {value:.6f}")
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