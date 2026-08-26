"""Fixed-position, two-time velocity correlations under free spinor evolution."""

from dataclasses import dataclass
import gc

import numpy as np
import pandas as pd
from scipy import fft as scipy_fft

from .config import WorkflowConfig
from .flow_reconstruction import FlowFieldResult, compute_velocity


@dataclass
class TemporalCorrelationResult:
    """Single-component Eulerian two-time velocity correlation."""

    tau: np.ndarray
    correlation: np.ndarray
    mean_initial: np.ndarray
    mean_evolved: np.ndarray
    variance_initial: np.ndarray
    variance_evolved: np.ndarray
    retained_fraction: np.ndarray


def _evolve_spinor_spectrum(
    spectrum,
    phase_x,
    phase_y,
    phase_z,
    workers,
):
    """Apply the separable free-Schr\"odinger phase and transform to space."""
    work = np.array(spectrum, copy=True, order="C")
    work *= phase_x[:, None, None]
    work *= phase_y[None, :, None]
    work *= phase_z[None, None, :]
    return scipy_fft.ifftn(
        work,
        workers=workers,
        overwrite_x=True,
    )


def _accumulate_density_and_current(
    psi,
    density,
    current,
    axis,
):
    """Accumulate one spinor component using plane-wise finite differences."""
    psi_axis = np.moveaxis(psi, axis, 0)
    density_axis = np.moveaxis(density, axis, 0)
    current_axis = np.moveaxis(current, axis, 0)
    n_axis = psi_axis.shape[0]

    for index in range(n_axis):
        center = psi_axis[index]
        forward = psi_axis[(index + 1) % n_axis]
        backward = psi_axis[(index - 1) % n_axis]

        center_real = center.real
        center_imag = center.imag

        density_axis[index] += (
            center_real * center_real
            + center_imag * center_imag
        )
        current_axis[index] += (
            center_real * (forward.imag - backward.imag)
            - center_imag * (forward.real - backward.real)
        )


def _correlation_moments(
    reference_velocity,
    density,
    current,
    grid_spacing,
    density_floor_relative,
    chunk_size,
):
    """Evaluate correlation moments without constructing the evolved velocity."""
    reference_flat = np.ravel(reference_velocity)
    density_flat = np.ravel(density)
    current_flat = np.ravel(current)

    mean_density = float(np.mean(density_flat, dtype=np.float64))
    density_floor = density_floor_relative * mean_density

    count = 0
    sum_reference = 0.0
    sum_evolved = 0.0
    sum_reference_squared = 0.0
    sum_evolved_squared = 0.0
    sum_cross = 0.0

    for start in range(0, density_flat.size, chunk_size):
        stop = min(start + chunk_size, density_flat.size)

        rho_chunk = density_flat[start:stop]
        current_chunk = current_flat[start:stop]
        reference_chunk = reference_flat[start:stop]

        valid = (
            np.isfinite(rho_chunk)
            & np.isfinite(current_chunk)
            & np.isfinite(reference_chunk)
            & (rho_chunk > density_floor)
        )
        if not np.any(valid):
            continue

        evolved_chunk = current_chunk[valid] / (
            2.0 * grid_spacing * rho_chunk[valid]
        )
        reference_chunk = reference_chunk[valid]

        count += evolved_chunk.size
        sum_reference += float(np.sum(reference_chunk, dtype=np.float64))
        sum_evolved += float(np.sum(evolved_chunk, dtype=np.float64))
        sum_reference_squared += float(
            np.sum(reference_chunk * reference_chunk, dtype=np.float64)
        )
        sum_evolved_squared += float(
            np.sum(evolved_chunk * evolved_chunk, dtype=np.float64)
        )
        sum_cross += float(
            np.sum(reference_chunk * evolved_chunk, dtype=np.float64)
        )

    if count == 0:
        raise FloatingPointError(
            "No grid points remain after applying the density threshold."
        )

    mean_reference = sum_reference / count
    mean_evolved = sum_evolved / count
    variance_reference = (
        sum_reference_squared / count - mean_reference**2
    )
    variance_evolved = (
        sum_evolved_squared / count - mean_evolved**2
    )
    covariance = sum_cross / count - mean_reference * mean_evolved

    roundoff_scale = np.finfo(np.float64).eps
    variance_reference = max(variance_reference, roundoff_scale)
    variance_evolved = max(variance_evolved, roundoff_scale)
    correlation = covariance / np.sqrt(
        variance_reference * variance_evolved
    )

    return (
        correlation,
        mean_reference,
        mean_evolved,
        variance_reference,
        variance_evolved,
        count / density_flat.size,
    )


def _initial_velocity_moments(reference_velocity, chunk_size):
    """Evaluate initial moments in chunks without a full-size mask or copy."""
    reference_flat = np.ravel(reference_velocity)
    count = 0
    total = 0.0
    total_squared = 0.0

    for start in range(0, reference_flat.size, chunk_size):
        stop = min(start + chunk_size, reference_flat.size)
        chunk = reference_flat[start:stop]
        valid = np.isfinite(chunk)
        if not np.any(valid):
            continue

        chunk = chunk[valid]
        count += chunk.size
        total += float(np.sum(chunk, dtype=np.float64))
        total_squared += float(np.sum(chunk * chunk, dtype=np.float64))

    if count == 0:
        raise FloatingPointError("The initial velocity contains no finite values.")

    mean = total / count
    variance = max(
        total_squared / count - mean**2,
        np.finfo(np.float64).eps,
    )
    return mean, variance, count / reference_flat.size


def compute_two_time_velocity_correlation(
    flow,
    tau,
    config=None,
    component="x",
    kinetic_coefficient=0.5,
    density_floor_relative=0.0,
    chunk_size=2_000_000,
    output_filename="two_time_velocity_correlation.csv",
):
    """Compute a fixed-position, two-time velocity correlation.

    The Pauli-spinor spectra evolve analytically according to

        psi_hat(k, tau) = exp[-i * kinetic_coefficient * |k|^2 * tau]
                            psi_hat(k, 0).

    For the equation i partial_t psi = -(1/2) Laplacian(psi), retain the
    default ``kinetic_coefficient=0.5``. Only the requested velocity component
    is reconstructed. The evolved three-dimensional velocity field is never
    stored; its moments and cross-correlation with the initial field are
    accumulated in chunks.
    """
    if not isinstance(flow, FlowFieldResult):
        raise TypeError("flow must be a FlowFieldResult instance.")

    config = WorkflowConfig() if config is None else config
    config.ensure_output_directories()
    if config.velocity_method != "FDM":
        raise ValueError(
            "This memory-efficient implementation uses the same FDM velocity "
            "definition as the default reconstruction. Set "
            "config.velocity_method='FDM' before reconstructing the flow."
        )

    component_to_axis = {"x": 0, "y": 1, "z": 2}
    component_to_field = {
        "x": flow.ux,
        "y": flow.uy,
        "z": flow.uz,
    }
    if component not in component_to_axis:
        raise ValueError("component must be 'x', 'y', or 'z'.")

    axis = component_to_axis[component]
    reference_velocity = component_to_field[component]
    if reference_velocity is None:
        raise ValueError(f"The initial {component}-velocity field is missing.")

    spectra = (flow.psi1_spec, flow.psi2_spec)
    if any(spectrum is None for spectrum in spectra):
        raise ValueError("Both initial Pauli-spinor spectra are required.")
    if any(spectrum.shape != spectra[0].shape for spectrum in spectra):
        raise ValueError("The two spinor spectra must have identical shapes.")

    tau = np.atleast_1d(np.asarray(tau, dtype=np.float64))
    if tau.ndim != 1 or not np.all(np.isfinite(tau)):
        raise ValueError("tau must be a finite one-dimensional array.")
    if kinetic_coefficient <= 0.0:
        raise ValueError("kinetic_coefficient must be positive.")
    if density_floor_relative < 0.0:
        raise ValueError("density_floor_relative must be nonnegative.")

    workers = max(1, int(config.max_parallel_threads))
    grid_spacing = 2.0 * np.pi / flow.N
    real_dtype = np.empty((), dtype=spectra[0].dtype).real.dtype

    kx_squared = np.asarray(flow.kx, dtype=np.float64) ** 2
    ky_squared = np.asarray(flow.ky, dtype=np.float64) ** 2
    kz_squared = np.asarray(flow.kz, dtype=np.float64) ** 2

    correlation = np.empty(tau.size, dtype=np.float64)
    mean_initial = np.empty(tau.size, dtype=np.float64)
    mean_evolved = np.empty(tau.size, dtype=np.float64)
    variance_initial = np.empty(tau.size, dtype=np.float64)
    variance_evolved = np.empty(tau.size, dtype=np.float64)
    retained_fraction = np.empty(tau.size, dtype=np.float64)

    for time_index, time_lag in enumerate(tau):
        if time_lag == 0.0:
            (
                reference_mean,
                reference_variance,
                reference_retained_fraction,
            ) = _initial_velocity_moments(
                reference_velocity,
                chunk_size,
            )

            correlation[time_index] = 1.0
            mean_initial[time_index] = reference_mean
            mean_evolved[time_index] = reference_mean
            variance_initial[time_index] = reference_variance
            variance_evolved[time_index] = reference_variance
            retained_fraction[time_index] = reference_retained_fraction
            print("Two-time velocity correlation completed: tau=0, C=1.")
            continue

        phase_x = np.exp(
            -1j * kinetic_coefficient * time_lag * kx_squared
        ).astype(spectra[0].dtype, copy=False)
        phase_y = np.exp(
            -1j * kinetic_coefficient * time_lag * ky_squared
        ).astype(spectra[0].dtype, copy=False)
        phase_z = np.exp(
            -1j * kinetic_coefficient * time_lag * kz_squared
        ).astype(spectra[0].dtype, copy=False)

        density = np.zeros(spectra[0].shape, dtype=real_dtype)
        current = np.zeros(spectra[0].shape, dtype=real_dtype)

        for spectrum in spectra:
            psi = _evolve_spinor_spectrum(
                spectrum,
                phase_x,
                phase_y,
                phase_z,
                workers,
            )
            _accumulate_density_and_current(
                psi,
                density,
                current,
                axis,
            )
            del psi
            gc.collect()

        (
            correlation[time_index],
            mean_initial[time_index],
            mean_evolved[time_index],
            variance_initial[time_index],
            variance_evolved[time_index],
            retained_fraction[time_index],
        ) = _correlation_moments(
            reference_velocity,
            density,
            current,
            grid_spacing,
            density_floor_relative,
            chunk_size,
        )

        del density, current, phase_x, phase_y, phase_z
        gc.collect()
        print(
            "Two-time velocity correlation completed: "
            f"tau={time_lag:.8e}, "
            f"C={correlation[time_index]:.8e}, "
            f"retained={retained_fraction[time_index]:.6%}."
        )

    result = TemporalCorrelationResult(
        tau=tau,
        correlation=correlation,
        mean_initial=mean_initial,
        mean_evolved=mean_evolved,
        variance_initial=variance_initial,
        variance_evolved=variance_evolved,
        retained_fraction=retained_fraction,
    )

    output_path = config.data_dir / output_filename
    pd.DataFrame(
        {
            "tau": result.tau,
            "correlation": result.correlation,
            "mean_initial": result.mean_initial,
            "mean_evolved": result.mean_evolved,
            "variance_initial": result.variance_initial,
            "variance_evolved": result.variance_evolved,
            "retained_fraction": result.retained_fraction,
        }
    ).to_csv(output_path, index=False)
    print(f"Two-time velocity correlation saved successfully: {output_path}")

    return result


def estimate_spinor_dephasing_time(flow, kinetic_coefficient=0.5):
    """Estimate the spectral dephasing time without forming a 3D K2 array."""
    kx2 = np.asarray(flow.kx, dtype=np.float64) ** 2
    ky2 = np.asarray(flow.ky, dtype=np.float64) ** 2
    kz2 = np.asarray(flow.kz, dtype=np.float64) ** 2

    k2_yz = ky2[:, None] + kz2[None, :]

    weight_sum = 0.0
    omega_sum = 0.0
    omega_squared_sum = 0.0

    for ix, kx_squared in enumerate(kx2):
        weight = np.abs(flow.psi1_spec[ix]) ** 2
        weight += np.abs(flow.psi2_spec[ix]) ** 2

        omega = kinetic_coefficient * (kx_squared + k2_yz)

        weight_sum += np.sum(weight, dtype=np.float64)
        omega_sum += np.sum(weight * omega, dtype=np.float64)
        omega_squared_sum += np.sum(
            weight * omega**2,
            dtype=np.float64,
        )

    omega_mean = omega_sum / weight_sum
    omega_variance = max(
        omega_squared_sum / weight_sum - omega_mean**2,
        0.0,
    )
    delta_omega = np.sqrt(omega_variance)
    tau_phi = 1.0 / delta_omega

    return omega_mean, delta_omega, tau_phi


def evolve_flow_field(
    flow,
    tau,
    config=None,
    kinetic_coefficient=0.5,
):
    """Reconstruct the flow at time tau under free spinor evolution.

    The two Pauli-spinor components evolve exactly according to

        psi_hat_s(k, tau)
        = exp[-i * kinetic_coefficient * |k|^2 * tau]
          * psi_hat_s(k, 0).

    The default kinetic_coefficient=0.5 corresponds to

        i * partial_t psi_s = -(1/2) * Laplacian(psi_s).

    The evolved spinors are transformed to physical space and subsequently
    mapped to the density and velocity fields using the generalized
    Madelung transform. Evolved spectral copies are not retained in the
    returned object in order to reduce memory usage.
    """
    if not isinstance(flow, FlowFieldResult):
        raise TypeError("flow must be a FlowFieldResult instance.")

    config = WorkflowConfig() if config is None else config

    if config.velocity_method != "FDM":
        raise ValueError(
            "The memory-efficient evolution currently uses the FDM velocity "
            "definition. Set config.velocity_method='FDM'."
        )

    tau = float(tau)

    if not np.isfinite(tau):
        raise ValueError("tau must be finite.")

    if kinetic_coefficient <= 0.0:
        raise ValueError("kinetic_coefficient must be positive.")

    if flow.psi1_spec is None or flow.psi2_spec is None:
        raise ValueError(
            "Both initial Pauli-spinor spectra are required."
        )

    physical_fields = (
        flow.psi1,
        flow.psi2,
        flow.rho,
        flow.ux,
        flow.uy,
        flow.uz,
    )

    # Return the existing initial field without recomputation.
    if tau == 0.0 and all(
        field is not None for field in physical_fields
    ):
        return flow

    workers = max(1, int(config.max_parallel_threads))

    # Since |k|^2 = k_x^2 + k_y^2 + k_z^2, the evolution phase separates as
    # exp(-i a |k|^2 tau)
    # = exp(-i a k_x^2 tau)
    #   exp(-i a k_y^2 tau)
    #   exp(-i a k_z^2 tau).
    # This avoids constructing a full three-dimensional complex phase array.
    phase_x = np.exp(
        -1j
        * kinetic_coefficient
        * tau
        * np.asarray(flow.kx, dtype=np.float64) ** 2
    ).astype(flow.psi1_spec.dtype, copy=False)

    phase_y = np.exp(
        -1j
        * kinetic_coefficient
        * tau
        * np.asarray(flow.ky, dtype=np.float64) ** 2
    ).astype(flow.psi1_spec.dtype, copy=False)

    phase_z = np.exp(
        -1j
        * kinetic_coefficient
        * tau
        * np.asarray(flow.kz, dtype=np.float64) ** 2
    ).astype(flow.psi1_spec.dtype, copy=False)

    psi1 = _evolve_spinor_spectrum(
        flow.psi1_spec,
        phase_x,
        phase_y,
        phase_z,
        workers,
    )

    psi2 = _evolve_spinor_spectrum(
        flow.psi2_spec,
        phase_x,
        phase_y,
        phase_z,
        workers,
    )

    del phase_x, phase_y, phase_z
    gc.collect()

    # Generalized Madelung reconstruction of density and velocity.
    rho, ux, uy, uz = compute_velocity(
        psi1,
        psi2,
        None,
        None,
        "FDM",
        flow.KX,
        flow.KY,
        flow.KZ,
        flow.N,
    )

    print(f"Flow field reconstructed successfully at t={tau:.8e}.")

    return FlowFieldResult(
        N=flow.N,
        kx=flow.kx,
        ky=flow.ky,
        kz=flow.kz,
        KX=flow.KX,
        KY=flow.KY,
        KZ=flow.KZ,
        K2=flow.K2,
        ik2=flow.ik2,
        psi1_spec=None,
        psi2_spec=None,
        psi1=psi1,
        psi2=psi2,
        rho=rho,
        ux=ux,
        uy=uy,
        uz=uz,
    )