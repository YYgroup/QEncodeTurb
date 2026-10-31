"""Single-realization density, dilatation and decoding stability diagnostics."""

import gc
import json
from pathlib import Path
import numpy as np
from scipy.fft import rfftn


def _decoding_fd_spectra(velocity, rho, mean_rho, epsilon, box_length, workers):
    """Full-vector spectra and FD2 divergence/curl moments, via Parseval."""
    N = rho.shape[0]
    h = box_length / N
    modes = np.fft.fftfreq(N) * N
    rmodes = np.fft.rfftfreq(N) * N
    kap = np.sin(2 * np.pi * modes / N) / h
    rkap = np.sin(2 * np.pi * rmodes / N) / h
    if N % 2 == 0:
        kap[N // 2] = 0.0
        rkap[-1] = 0.0
    multiplicity = np.full(rmodes.size, 2.0)
    multiplicity[0] = 1.0
    if N % 2 == 0:
        multiplicity[-1] = 1.0
    div_hat = np.zeros((N, N, rmodes.size), dtype=np.complex128)
    Eu = np.zeros(N + 1)
    Egrad = np.zeros(N + 1)
    m12 = modes[:, None]**2 + rmodes[None, :]**2
    kap12 = kap[:, None]**2 + rkap[None, :]**2

    for c, component in enumerate(velocity):
        if epsilon == 0.0:
            field = np.asarray(component, dtype=np.float64)
        else:
            field = np.empty(rho.shape, dtype=np.float64)
            for i in range(N):
                q = rho[i] / mean_rho
                field[i] = component[i] * (q / (q + epsilon))
        hat = rfftn(field, norm="forward", workers=workers)
        del field
        for i in range(N):
            shell = np.rint(np.sqrt(modes[i]**2 + m12)).astype(np.intp)
            power = (hat[i].real**2 + hat[i].imag**2) * (0.5 * multiplicity)
            Eu += np.bincount(shell.ravel(), weights=power.ravel(),
                              minlength=N + 1)
            Egrad += np.bincount(
                shell.ravel(),
                weights=(power * (kap[i]**2 + kap12)).ravel(),
                minlength=N + 1,
            )
            symbol = (kap[i], kap[:, None], rkap[None, :])[c]
            div_hat[i] += symbol * hat[i]
        del hat
        gc.collect()

    Ediv = np.zeros(N + 1)
    for i in range(N):
        shell = np.rint(np.sqrt(modes[i]**2 + m12)).astype(np.intp)
        power = (div_hat[i].real**2 + div_hat[i].imag**2) * (0.5 * multiplicity)
        Ediv += np.bincount(shell.ravel(), weights=power.ravel(),
                           minlength=N + 1)
    del div_hat
    Ew = Egrad - Ediv
    scale = max(float(Egrad.max()), np.finfo(float).tiny)
    if np.min(Ew) < -1e-10 * scale:
        raise ArithmeticError("Unexpected negative curl energy.")
    np.maximum(Ew, 0.0, out=Ew)  # remove roundoff-level negative values
    div2 = 2.0 * Ediv.sum()
    curl2 = 2.0 * Ew.sum()
    theta = div2 / (div2 + curl2) if div2 + curl2 > 0 else np.nan
    return {
        "k": np.arange(N + 1) * (2 * np.pi / box_length),
        "Eu": Eu, "Ew": Ew,
        "mean_u_squared": float(2 * Eu.sum()),
        "mean_divergence_squared": float(div2),
        "mean_vorticity_squared": float(curl2),
        "vartheta": float(theta),
    }


def _decoding_regularization_changes(velocity, rho, mean_rho, epsilon, h):
    """Squared differences, evaluated plane by plane with periodic FD2."""
    N = rho.shape[0]
    u_difference_sum = 0.0
    w_difference_sum = 0.0
    low_count = 0

    def delta_plane(i):
        i %= N
        q = rho[i] / mean_rho
        factor = q / (q + epsilon)
        return tuple(
            np.asarray(v[i], dtype=np.float64) * factor - v[i]
            for v in velocity
        )

    def dd(a, axis):
        return (np.roll(a, -1, axis=axis) -
                np.roll(a, 1, axis=axis)) / (2 * h)

    for i in range(N):
        d = delta_plane(i)
        previous = delta_plane(i - 1)
        following = delta_plane(i + 1)
        low_count += np.count_nonzero(rho[i] / mean_rho < epsilon)
        for a in d:
            u_difference_sum += np.sum(a * a, dtype=np.float64)

        # Original 3-D axes 1 and 2 become plane axes 0 and 1.
        dwx = dd(d[2], 0) - dd(d[1], 1)
        dwy = dd(d[0], 1) - (following[2] - previous[2]) / (2 * h)
        dwz = (following[1] - previous[1]) / (2 * h) - dd(d[0], 0)
        for a in (dwx, dwy, dwz):
            w_difference_sum += np.sum(a * a, dtype=np.float64)
        del d, previous, following, dwx, dwy, dwz

    return (
        float(u_difference_sum / rho.size),
        float(w_difference_sum / rho.size),
        int(low_count),
    )


def _decoding_slope(k, spectrum, fit_range):
    lo, hi = fit_range
    good = ((k >= lo) & (k <= hi) & (k > 0) &
            np.isfinite(spectrum) & (spectrum > 0))
    kfit = k[good]
    result = {
        "slope": np.nan, "r_squared": np.nan,
        "number_of_shells": int(kfit.size),
        "actual_k_min": float(kfit.min()) if kfit.size else np.nan,
        "actual_k_max": float(kfit.max()) if kfit.size else np.nan,
    }
    if kfit.size >= 3:
        x, y = np.log(kfit), np.log(spectrum[good])
        xc, yc = x - x.mean(), y - y.mean()
        slope = np.dot(xc, yc) / np.dot(xc, xc)
        residual = yc - slope * xc
        sst = np.dot(yc, yc)
        result["slope"] = float(slope)
        result["r_squared"] = float(
            1 - np.dot(residual, residual) / sst
        ) if sst > 0 else np.nan
    return result


def _decoding_density_and_nonfinite(flow):
    rho = np.asarray(flow.rho)
    U = tuple(np.asarray(getattr(flow, name))
              for name in ("ux", "uy", "uz"))
    P = (np.asarray(flow.psi1), np.asarray(flow.psi2))
    if rho.ndim != 3 or len(set(rho.shape)) != 1:
        raise ValueError("A cubic 3-D grid is required.")
    if any(a.shape != rho.shape for a in U + P):
        raise ValueError("rho, velocity and spinor shapes must match.")
    if np.iscomplexobj(rho) or any(np.iscomplexobj(a) for a in U):
        raise ValueError("Density and velocity must be real arrays.")
    bad_u = bad_s = bad_union = zero_rho = 0

    for i in range(rho.shape[0]):
        r = rho[i]
        if not np.all(np.isfinite(r)) or np.any(r < 0):
            raise ValueError("Density must be finite and nonnegative.")
        zero_rho += np.count_nonzero(r == 0)
        bu = np.zeros(r.shape, dtype=bool)
        bs = np.zeros(r.shape, dtype=bool)
        for a in U:
            bu |= ~np.isfinite(a[i])
        a, b = P[0][i], P[1][i]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            numerators = (
                2 * (a.real * b.real + a.imag * b.imag),
                2 * (a.real * b.imag - a.imag * b.real),
                a.real**2 + a.imag**2 - b.real**2 - b.imag**2,
            )
            for numerator in numerators:
                bs |= ~np.isfinite(numerator / r)
        bad_u += np.count_nonzero(bu)
        bad_s += np.count_nonzero(bs)
        bad_union += np.count_nonzero(bu | bs)

    mean_rho = float(np.mean(rho, dtype=np.float64))
    if not np.isfinite(mean_rho) or mean_rho <= 0:
        raise ValueError("Mean density must be positive and finite.")
    minimum = float(np.min(rho)) / mean_rho
    # Exact quantiles: one temporary density copy, no full q array.
    buffer = np.array(rho, dtype=np.float64, order="C", copy=True).ravel()
    quantiles = np.quantile(
        buffer, (0.0001, 0.001, 0.01),
        method="linear", overwrite_input=True,
    ) / mean_rho
    del buffer
    gc.collect()
    return rho, U, {
        "mean_rho": mean_rho,
        "minimum_q": minimum,
        "q_quantile_0.01_percent": float(quantiles[0]),
        "q_quantile_0.1_percent": float(quantiles[1]),
        "q_quantile_1_percent": float(quantiles[2]),
        "zero_density_points": int(zero_rho),
        "nonfinite_velocity_points": int(bad_u),
        "nonfinite_spin_points": int(bad_s),
        "nonfinite_velocity_or_spin_points": int(bad_union),
    }



__all__ = ["write_flow_decoding_report"]


def _scalar_shell_spectrum(field, box_length, workers):
    """Return 0.5 * shell sum of squared normalized scalar FFT."""
    N = field.shape[0]
    modes = np.fft.fftfreq(N) * N
    rmodes = np.fft.rfftfreq(N) * N
    mult = np.full(rmodes.size, 2.0)
    mult[0] = 1.0
    if N % 2 == 0:
        mult[-1] = 1.0
    hat = rfftn(field, norm="forward", workers=workers)
    spectrum = np.zeros(N + 1)
    transverse_m2 = modes[:, None]**2 + rmodes[None, :]**2
    for i in range(N):
        shell = np.rint(
            np.sqrt(modes[i]**2 + transverse_m2)
        ).astype(np.intp)
        power = (
            hat[i].real**2 + hat[i].imag**2
        ) * (0.5 * mult)
        spectrum += np.bincount(
            shell.ravel(), weights=power.ravel(), minlength=N + 1
        )
    del hat
    return spectrum


def _spin_spectra(flow, rho, mean_rho, epsilon, box_length, workers):
    """
    Scalar spectra of the three spin-defined surface fields, separately.

    Reference: s_a = (Psi^dagger sigma_a Psi) / rho.
    Diagnostic: s_a,epsilon = numerator / (rho + epsilon * mean_rho).
    No subsequent unit-length normalization is applied.
    """
    N = rho.shape[0]
    spectra = {}
    difference_squared_sum = 0.0
    for component in range(3):
        field = np.empty(rho.shape, dtype=np.float64)
        for i in range(N):
            a, b = flow.psi1[i], flow.psi2[i]
            if component == 0:
                numerator = 2 * (a.real * b.real + a.imag * b.imag)
            elif component == 1:
                numerator = 2 * (a.real * b.imag - a.imag * b.real)
            else:
                numerator = (
                    a.real**2 + a.imag**2 - b.real**2 - b.imag**2
                )
            reference = numerator / rho[i]
            if epsilon == 0:
                field[i] = reference
            else:
                field[i] = numerator / (rho[i] + epsilon * mean_rho)
                difference = field[i] - reference
                difference_squared_sum += np.sum(
                    difference * difference, dtype=np.float64
                )
        spectra[f"Es{component + 1}"] = _scalar_shell_spectrum(
            field, box_length, workers
        )
        del field
        gc.collect()
    spectra["mean_s_squared"] = float(
        2 * sum(spectra[f"Es{a}"].sum() for a in (1, 2, 3))
    )
    spectra["mean_spin_difference_squared"] = float(
        difference_squared_sum / rho.size
    )
    return spectra


def write_flow_decoding_report(
    flow, txt_path=None, *, config=None,
    epsilons=(1e-6, 1e-5, 1e-4, 1e-3),
    phase_seeds=None, fit_range=(15.0, None),
    spin_fit_range=None, box_length=2 * np.pi,
    fft_workers=1, append=False, parameters=None,
):
    """
    Write all diagnostics for ONE realization to ONE TXT; return None.

    Parameters
    ----------
    flow
        Existing FlowFieldResult. Its rho, ux, uy, uz, psi1 and psi2
        must still exist. The function never modifies these arrays.
    txt_path : str or Path, optional
        Defaults to config.data_dir / "decoding_diagnostics.txt",
        or ./decoding_diagnostics.txt without config.
    config : WorkflowConfig, optional
        Used only for seeds, provenance and default output directory.
        derivative_method must be FDM if config is supplied.
    epsilons : iterable of positive float
        All diagnostic regularizations evaluated by this one call.
    phase_seeds : optional
        Actual phase seeds, recorded only; never generate a new field.
    fit_range : (float, float or None)
        Velocity-spectrum fit interval. None means N/2 for box 2*pi.
    spin_fit_range : (float, float or None), optional
        Shared interval for Es1, Es2, Es3; defaults to fit_range.
    box_length : float
        Periodic cubic box length; default 2*pi.
    fft_workers : int
        Number of SciPy FFT workers; default 1.
    append : bool
        False overwrites the TXT; True appends a realization.
    parameters : dict, optional
        Extra provenance, e.g. kd or component population.

    Definitions
    -----------
    Eu(k) = 1/2 sum_shell (|FFT(ux)/N^3|^2
                        + |FFT(uy)/N^3|^2 + |FFT(uz)/N^3|^2).
    Esa(k) = 1/2 sum_shell |FFT(sa)/N^3|^2, a=1,2,3.
    This is NOT FFT(sqrt(ux^2+uy^2+uz^2)), NOT the spectrum of
    |s|, and NOT a compensated or shell-averaged spectrum.

    Divergence and curl use periodic centered second-order differences
    with ux along array axis 0, uy along axis 1, uz along axis 2.
    All epsilon cases keep rho and the decoding numerators fixed.
    No projection, density floor or cross-seed statistics are applied.

    Notes
    -----
    Exact quantiles use one temporary density copy. FFTs are evaluated
    sequentially; full fields for all epsilons are never retained.
    """
    epsilons = tuple(dict.fromkeys(float(e) for e in epsilons))
    if not epsilons or any(
        not np.isfinite(e) or e <= 0 for e in epsilons
    ):
        raise ValueError("epsilons must contain positive finite values.")
    if not np.isfinite(box_length) or box_length <= 0:
        raise ValueError("box_length must be positive and finite.")
    if not isinstance(fft_workers, int) or fft_workers < 1:
        raise ValueError("fft_workers must be a positive integer.")
    if config is not None and getattr(
        config, "derivative_method", "FDM"
    ) != "FDM":
        raise ValueError(
            "This diagnostic uses FDM (periodic centered FD2). "
            "It must not be labeled as an SM derivative test."
        )

    # Fail clearly if cleanup has already released required arrays.
    for name in ("rho", "ux", "uy", "uz", "psi1", "psi2"):
        if getattr(flow, name, None) is None:
            raise ValueError(
                f"flow.{name} is unavailable. Run before MEMORY_RELEASE N03."
            )

    rho, velocity, density = _decoding_density_and_nonfinite(flow)
    N = rho.shape[0]
    mean_rho = density["mean_rho"]

    def interval(bounds):
        lo = float(bounds[0])
        hi = (
            np.pi * N / box_length
            if bounds[1] is None else float(bounds[1])
        )
        if not np.isfinite(lo + hi) or lo <= 0 or hi <= lo:
            raise ValueError("Invalid fitting interval.")
        return (lo, hi)

    velocity_interval = interval(fit_range)
    spin_interval = interval(
        fit_range if spin_fit_range is None else spin_fit_range
    )
    if txt_path is None:
        directory = (
            Path(config.data_dir) if config is not None else Path.cwd()
        )
        txt_path = directory / "decoding_diagnostics.txt"
    path = Path(txt_path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)

    provenance = {}
    if config is not None:
        for key in (
            "k_cutoff", "component_ratio", "velocity_method",
            "derivative_method",
        ):
            if hasattr(config, key):
                provenance[key] = getattr(config, key)
        if phase_seeds is None:
            phase_seeds = (
                getattr(config, "seed_spin_up", None),
                getattr(config, "seed_spin_down", None),
            )
    if parameters is not None:
        provenance.update(parameters)

    k = np.arange(N + 1) * (2 * np.pi / box_length)
    nan_spectrum = np.full(N + 1, np.nan)

    def dump(handle, title, values):
        handle.write("\n[" + title + "]\n")
        handle.write(json.dumps(
            values, indent=2, ensure_ascii=False, default=str,
            allow_nan=True,
        ) + "\n")

    def fits(u_data, s_data):
        answer = {"velocity": _decoding_slope(
            k, nan_spectrum if u_data is None else u_data["Eu"],
            velocity_interval,
        )}
        for a in (1, 2, 3):
            answer[f"s{a}"] = _decoding_slope(
                k, nan_spectrum if s_data is None else s_data[f"Es{a}"],
                spin_interval,
            )
        return answer

    def save_spectra(handle, title, u_data, s_data):
        arrays = [
            k,
            nan_spectrum if u_data is None else u_data["Eu"],
        ]
        arrays.extend(
            nan_spectrum if s_data is None else s_data[f"Es{a}"]
            for a in (1, 2, 3)
        )
        handle.write("\n[" + title + "]\n")
        handle.write("k\tE_velocity\tE_s1\tE_s2\tE_s3\n")
        np.savetxt(
            handle, np.column_stack(arrays), delimiter="\t", fmt="%.12e"
        )

    def relative(change_squared, reference_squared):
        return (
            float(np.sqrt(change_squared / reference_squared))
            if reference_squared > 0 else np.nan
        )

    u_valid = density["nonfinite_velocity_points"] == 0
    s_valid = density["nonfinite_spin_points"] == 0
    with path.open("a" if append else "w", encoding="utf-8") as handle:
        handle.write("\n" + "=" * 72 + "\n")
        dump(handle, "REALIZATION", {
            "grid": f"{N} x {N} x {N}",
            "n_qubits": (
                3 * (N.bit_length() - 1) if N & (N - 1) == 0 else None
            ),
            "number_of_phase_realizations": 1,
            "phase_seeds": phase_seeds,
            "prescribed_P0": 0.1,
            "parameters": provenance,
            "box_length": box_length,
            "derivative": "periodic centered second-order difference",
            "component_axes": "ux:axis0, uy:axis1, uz:axis2",
            "velocity_fit_interval": velocity_interval,
            "spin_fit_interval": spin_interval,
            "epsilon_values": epsilons,
            "velocity_spectrum": (
                "1/2 shell sum of |ux_hat|^2+|uy_hat|^2+|uz_hat|^2; "
                "hat=FFT/N^3"
            ),
            "surface_field_spectra": (
                "Es1, Es2, Es3 separately; 1/2 shell sum of "
                "|FFT(sa)/N^3|^2"
            ),
            "shell_definition": "round(|integer Fourier mode|)",
            "regularization": (
                "u_epsilon=J/(rho+epsilon*mean_rho); "
                "sa_epsilon=spin_numerator_a/(rho+epsilon*mean_rho); "
                "no subsequent spin normalization"
            ),
            "nonfinite_point_rule": (
                "count each grid point once even if several components "
                "are nonfinite"
            ),
            "NaN": "undefined or not evaluated; never assumed zero",
        })
        dump(handle, "TABLE_S3_DENSITY_AND_NONFINITE_COUNTS", density)

        base_u = (
            _decoding_fd_spectra(
                velocity, rho, mean_rho, 0.0, box_length, fft_workers
            ) if u_valid else None
        )
        base_s = (
            _spin_spectra(
                flow, rho, mean_rho, 0.0, box_length, fft_workers
            ) if s_valid else None
        )
        dump(handle, "DILATATION_AND_REFERENCE_NORMS", {
            name: (
                base_u[name] if base_u is not None else np.nan
            ) for name in (
                "vartheta", "mean_u_squared",
                "mean_divergence_squared", "mean_vorticity_squared",
            )
        })
        dump(handle, "REFERENCE_SPIN_NORM", {
            "mean_s_squared": (
                base_s["mean_s_squared"] if base_s is not None else np.nan
            )
        })
        base_fits = fits(base_u, base_s)
        dump(handle, "REFERENCE_FITS", base_fits)
        save_spectra(handle, "REFERENCE_SPECTRA", base_u, base_s)
        handle.flush()

        for epsilon in epsilons:
            if u_valid:
                du2, dw2, low_count = _decoding_regularization_changes(
                    velocity, rho, mean_rho, epsilon, box_length / N
                )
                reg_u = _decoding_fd_spectra(
                    velocity, rho, mean_rho, epsilon,
                    box_length, fft_workers,
                )
                delta_u = relative(du2, base_u["mean_u_squared"])
                delta_w = relative(
                    dw2, base_u["mean_vorticity_squared"]
                )
            else:
                low_count = int(sum(
                    np.count_nonzero(r / mean_rho < epsilon) for r in rho
                ))
                reg_u = None
                delta_u = delta_w = np.nan
            reg_s = (
                _spin_spectra(
                    flow, rho, mean_rho, epsilon,
                    box_length, fft_workers,
                ) if s_valid else None
            )
            delta_s = (
                relative(
                    reg_s["mean_spin_difference_squared"],
                    base_s["mean_s_squared"],
                ) if s_valid else np.nan
            )
            reg_fits = fits(reg_u, reg_s)
            row = {
                "epsilon": epsilon,
                "low_density_points": int(low_count),
                "f_epsilon": low_count / rho.size,
                "Delta_u": delta_u,
                "Delta_omega": delta_w,
                "Delta_s": delta_s,
                "velocity_spectrum_slope_change": (
                    reg_fits["velocity"]["slope"] -
                    base_fits["velocity"]["slope"]
                ),
            }
            for a in (1, 2, 3):
                row[f"s{a}_spectrum_slope_change"] = (
                    reg_fits[f"s{a}"]["slope"] -
                    base_fits[f"s{a}"]["slope"]
                )
            row["regularized_vartheta"] = (
                reg_u["vartheta"] if reg_u is not None else np.nan
            )
            dump(handle, f"TABLE_S4 epsilon={epsilon:.12e}", row)
            dump(handle, "REGULARIZED_FITS", reg_fits)
            save_spectra(
                handle, f"REGULARIZED_SPECTRA epsilon={epsilon:.12e}",
                reg_u, reg_s,
            )
            handle.flush()
            del reg_u, reg_s
            gc.collect()
    return None

