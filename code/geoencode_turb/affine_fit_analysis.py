"""Conditional-angle affine-fit diagnostics for Gray and binary orderings."""

from dataclasses import dataclass
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import WorkflowConfig
from .quantum_encoding import EncodingResult


CUTOFF = 1.0e-25
RIDGE_ALPHA = 1.0e-9
CHUNK_SIZE = 200_000
WEIGHT_THRESHOLD = 1.0e-3
DO_BINARY_REFERENCE = True
REPRESENTATIVE_LAYER = None
MAX_SCATTER_POINTS = 20_000

OUT_DIR = Path("data/conditional_angle_affine_fit")
FIG_DIR = Path("figures")


@dataclass
class AffineFitResult:
    """Complete Gray-versus-binary conditional-angle diagnostics."""

    gray_params: list
    gray_records: list
    binary_params: list
    binary_records: list
    gray_summary: dict
    binary_summary: dict
    diagnostic_table: pd.DataFrame


def layer_samples(probs_flat, layer):
    """Yield precursor indices, target angles, and omega=sqrt[P(c_j)]."""
    probs_flat = np.asarray(probs_flat, dtype=np.float64).reshape(-1)
    n_blocks = 1 << layer
    remainder = probs_flat.size // (2 * n_blocks)
    view = probs_flat.reshape(n_blocks, 2, remainder)

    for start in range(0, n_blocks, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, n_blocks)
        p0 = view[start:stop, 0].sum(axis=1, dtype=np.float64)
        p1 = view[start:stop, 1].sum(axis=1, dtype=np.float64)
        marginal = p0 + p1
        valid = marginal > CUTOFF
        if not np.any(valid):
            continue

        idx = np.arange(start, stop, dtype=np.uint64)[valid]
        marginal = marginal[valid]
        theta = 2.0 * np.arccos(
            np.sqrt(np.clip(p0[valid] / marginal, 0.0, 1.0))
        )
        omega = np.sqrt(marginal)
        yield idx, theta, omega


def design_matrix(indices, layer):
    """Return [q_0,...,q_{j-1},1] in the original bit convention."""
    if layer == 0:
        return np.ones((indices.size, 1), dtype=np.float64)
    shifts = np.arange(layer - 1, -1, -1, dtype=np.uint64)
    bits = ((indices[:, None] >> shifts[None]) & 1).astype(np.float64)
    return np.column_stack((bits, np.ones(indices.size)))


def fit_parameter_set(probs_3d, name, n_total):
    """Perform a chunked affine ridge fit for every layer."""
    probs_flat = np.asarray(probs_3d, dtype=np.float64).reshape(-1)
    params = []

    for layer in range(n_total):
        tic = time.perf_counter()
        p = layer + 1
        A = np.zeros((p, p), dtype=np.float64)
        b = np.zeros(p, dtype=np.float64)
        max_omega = 0.0

        for idx, theta, omega in layer_samples(probs_flat, layer):
            X = design_matrix(idx, layer)
            A += X.T @ (omega[:, None] * X)
            b += X.T @ (omega * theta)
            max_omega = max(max_omega, float(omega.max()))

        # This is equivalent to the original omega/max(omega) weighting.
        penalty = np.eye(p)
        penalty[-1, -1] = 0.0
        system = A + RIDGE_ALPHA * max_omega * penalty
        try:
            beta = np.linalg.solve(system, b)
        except np.linalg.LinAlgError:
            beta = np.linalg.lstsq(system, b, rcond=None)[0]

        params.append((float(beta[-1]), beta[:-1].copy()))
        print(
            f"{name:>6s} fit layer {layer:2d}: "
            f"time={time.perf_counter() - tic:.1f} s"
        )
    return params


def evaluate_parameter_set(probs_3d, params, name):
    """Evaluate K_j and all supporting diagnostics."""
    probs_flat = np.asarray(probs_3d, dtype=np.float64).reshape(-1)
    records = []

    for layer, (bias, weights) in enumerate(params):
        tic = time.perf_counter()
        beta = np.r_[np.asarray(weights, dtype=np.float64), float(bias)]

        # First pass: stable weighted mean and effective sample size.
        sw = sw2 = swy = 0.0
        n_valid = 0
        target_min, target_max = np.inf, -np.inf
        for _, theta, omega in layer_samples(probs_flat, layer):
            sw += float(omega.sum())
            sw2 += float(omega @ omega)
            swy += float(omega @ theta)
            n_valid += theta.size
            target_min = min(target_min, float(theta.min()))
            target_max = max(target_max, float(theta.max()))

        theta_mean = swy / sw
        n_eff = sw**2 / sw2

        # Second pass: direct TSS, SSE, and predicted-angle range.
        sse = tss = outside_weight = 0.0
        outside_count = 0
        pred_min, pred_max = np.inf, -np.inf
        for idx, theta, omega in layer_samples(probs_flat, layer):
            pred = design_matrix(idx, layer) @ beta
            residual = theta - pred
            centered = theta - theta_mean
            sse += float(omega @ (residual * residual))
            tss += float(omega @ (centered * centered))
            pred_min = min(pred_min, float(pred.min()))
            pred_max = max(pred_max, float(pred.max()))
            outside = (pred < 0.0) | (pred > np.pi)
            outside_count += int(outside.sum())
            outside_weight += float(omega[outside].sum())

        weighted_mse = sse / sw
        weighted_rmse = np.sqrt(weighted_mse)
        K = weighted_mse / np.pi**2
        sqrt_K = weighted_rmse / np.pi

        tss_tol = 100.0 * np.finfo(float).eps * max(sw * np.pi**2, 1.0)
        R2 = np.nan if tss <= tss_tol else 1.0 - sse / tss
        weights = np.asarray(weights, dtype=np.float64)

        record = {
            "layer": layer,
            "n_parameters": layer + 1,
            "n_valid": n_valid,
            "sample_to_parameter_ratio": n_valid / (layer + 1),
            "effective_sample_size": n_eff,
            "effective_sample_to_parameter_ratio": n_eff / (layer + 1),
            "sum_weight": sw,
            "target_angle_mean": theta_mean,
            "target_angle_std": np.sqrt(tss / sw),
            "target_angle_min": target_min,
            "target_angle_max": target_max,
            "prediction_angle_min": pred_min,
            "prediction_angle_max": pred_max,
            "outside_angle_count": outside_count,
            "outside_angle_fraction": outside_count / n_valid,
            "outside_angle_weight_fraction": outside_weight / sw,
            "n_affine_weights": weights.size,
            "n_significant_weights": int(
                (np.abs(weights) >= WEIGHT_THRESHOLD).sum()
            ),
            "SSE": sse,
            "TSS": tss,
            "weighted_MSE": weighted_mse,
            "weighted_RMSE": weighted_rmse,
            "K": K,
            "sqrt_K": sqrt_K,
            "R2_weighted": R2,
        }
        records.append(record)

        print(
            f"{name:>6s} layer {layer:2d}: "
            f"K={K:.6e}, sqrt(K)={sqrt_K:.6e}, "
            f"N_valid={n_valid:,}, N_eff={n_eff:.2f}, "
            f"outside={outside_weight / sw:.3e}, "
            f"time={time.perf_counter() - tic:.1f} s"
        )
    return records


def summarize(records):
    """Summarize layer-resolved and global affine-fit diagnostics."""
    nontrivial = [record for record in records if record["layer"] > 0]
    K = np.asarray([record["K"] for record in nontrivial])
    layers = np.asarray([record["layer"] for record in nontrivial])
    imax = int(np.argmax(K))
    total_sse = sum(record["SSE"] for record in records)
    total_weight = sum(record["sum_weight"] for record in records)
    global_K = total_sse / (np.pi**2 * total_weight)

    informative = [
        record
        for record in nontrivial
        if record["target_angle_std"] > 1.0e-12
        and record["n_valid"] > record["n_parameters"]
    ]
    informative_median = (
        float(np.median([record["K"] for record in informative]))
        if informative
        else np.nan
    )
    return {
        "K_min": float(K.min()),
        "K_max": float(K.max()),
        "K_median": float(np.median(K)),
        "K_max_layer": int(layers[imax]),
        "global_K": float(global_K),
        "global_sqrt_K": float(np.sqrt(global_K)),
        "informative_K_median": informative_median,
        "max_outside_weight_fraction": float(
            max(record["outside_angle_weight_fraction"] for record in records)
        ),
    }


def choose_representative_layer(records):
    """Select a nontrivial layer near the informative-layer median K."""
    candidates = [
        record
        for record in records
        if record["layer"] > 0
        and record["target_angle_std"] > 1.0e-12
        and record["n_valid"] > record["n_parameters"]
    ]
    if not candidates:
        candidates = [record for record in records if record["layer"] > 0]
    median_K = max(float(np.median([r["K"] for r in candidates])), 1.0e-30)
    return min(
        candidates,
        key=lambda record: abs(
            np.log10(max(record["K"], 1.0e-30) / median_K)
        ),
    )["layer"]


def scatter_data(probs_3d, params, layer):
    """Collect a bounded set of target and fitted conditional angles."""
    probs_flat = np.asarray(probs_3d, dtype=np.float64).reshape(-1)
    bias, weights = params[layer]
    beta = np.r_[np.asarray(weights), bias]
    target_list, pred_list, omega_list = [], [], []

    for idx, theta, omega in layer_samples(probs_flat, layer):
        target_list.append(theta)
        pred_list.append(design_matrix(idx, layer) @ beta)
        omega_list.append(omega)

    target = np.concatenate(target_list)
    pred = np.concatenate(pred_list)
    omega = np.concatenate(omega_list)
    if target.size > MAX_SCATTER_POINTS:
        take = np.linspace(0, target.size - 1, MAX_SCATTER_POINTS, dtype=int)
        target, pred, omega = target[take], pred[take], omega[take]
    return target, pred, omega


def _format_summary(gray_summary, binary_summary):
    summary_lines = [
        "K_j = SSE_j/[pi^2 sum_c omega_c]; sqrt(K_j)=weighted RMSE/pi.",
        "",
        (
            "Gray ordering: "
            f"K range=[{gray_summary['K_min']:.6e}, "
            f"{gray_summary['K_max']:.6e}], "
            f"median={gray_summary['K_median']:.6e}, "
            f"maximum at layer {gray_summary['K_max_layer']}, "
            f"global K={gray_summary['global_K']:.6e}, "
            f"global sqrt(K)={gray_summary['global_sqrt_K']:.6e}, "
            "informative-layer median="
            f"{gray_summary['informative_K_median']:.6e}, "
            "maximum out-of-range weight fraction="
            f"{gray_summary['max_outside_weight_fraction']:.6e}."
        ),
    ]

    if binary_summary:
        ratio = (
            binary_summary["global_K"] / gray_summary["global_K"]
            if gray_summary["global_K"] > 0.0
            else np.inf
        )
        summary_lines += [
            (
                "Binary ordering: "
                f"K range=[{binary_summary['K_min']:.6e}, "
                f"{binary_summary['K_max']:.6e}], "
                f"median={binary_summary['K_median']:.6e}, "
                f"maximum at layer {binary_summary['K_max_layer']}, "
                f"global K={binary_summary['global_K']:.6e}, "
                f"global sqrt(K)={binary_summary['global_sqrt_K']:.6e}, "
                "informative-layer median="
                f"{binary_summary['informative_K_median']:.6e}, "
                "maximum out-of-range weight fraction="
                f"{binary_summary['max_outside_weight_fraction']:.6e}."
            ),
            f"Global-K ratio (binary/Gray)={ratio:.6e}.",
        ]
    return "\n".join(summary_lines)


def run_conditional_angle_analysis(encoding, config=None):
    """Run the complete Gray-versus-binary K_j analysis and export the CSV."""
    if not isinstance(encoding, EncodingResult):
        raise TypeError("encoding must be an EncodingResult instance.")
    config = WorkflowConfig() if config is None else config
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    gray_params = encoding.amp_params
    gray_records = evaluate_parameter_set(
        encoding.target_probs_gray,
        gray_params,
        "Gray",
    )

    if DO_BINARY_REFERENCE:
        binary_params = fit_parameter_set(
            encoding.target_probs_phys,
            "binary",
            config.n_total,
        )
        binary_records = evaluate_parameter_set(
            encoding.target_probs_phys,
            binary_params,
            "binary",
        )
    else:
        binary_params = binary_records = None

    rows = []
    for layer, gray in enumerate(gray_records):
        row = {"layer": layer}
        row.update({f"{key}_gray": value for key, value in gray.items() if key != "layer"})
        if binary_records is not None:
            row.update(
                {
                    f"{key}_binary": value
                    for key, value in binary_records[layer].items()
                    if key != "layer"
                }
            )
        rows.append(row)

    diagnostic_table = pd.DataFrame(rows)
    diagnostic_table.to_csv(
        OUT_DIR / "conditional_angle_affine_fit_error.csv",
        index=False,
    )

    gray_summary = summarize(gray_records)
    binary_summary = summarize(binary_records) if binary_records else None
    summary_text = _format_summary(gray_summary, binary_summary)
    (OUT_DIR / "conditional_angle_affine_fit_error_summary.txt").write_text(
        summary_text + "\n",
        encoding="utf-8",
    )
    print("\n" + summary_text)

    return AffineFitResult(
        gray_params=gray_params,
        gray_records=gray_records,
        binary_params=binary_params,
        binary_records=binary_records,
        gray_summary=gray_summary,
        binary_summary=binary_summary,
        diagnostic_table=diagnostic_table,
    )


def plot_affine_fit_error_comparison(result, pdf_path=None):
    """Plot the three normalized Gray-versus-binary error metrics."""
    gray_records = result.gray_records
    binary_records = result.binary_records

    gray_layers = np.array([record["layer"] for record in gray_records], dtype=int)
    binary_layers = np.array(
        [record["layer"] for record in binary_records],
        dtype=int,
    )
    gray_K = np.array([record["K"] for record in gray_records], dtype=float)
    binary_K = np.array([record["K"] for record in binary_records], dtype=float)

    if not np.array_equal(gray_layers, binary_layers):
        raise ValueError("Gray and binary records use different circuit layers.")

    nontrivial = gray_layers > 0
    gray_median_K = np.median(gray_K[nontrivial])
    binary_median_K = np.median(binary_K[nontrivial])

    gray_total_SSE = sum(record["SSE"] for record in gray_records)
    gray_total_weight = sum(record["sum_weight"] for record in gray_records)
    binary_total_SSE = sum(record["SSE"] for record in binary_records)
    binary_total_weight = sum(record["sum_weight"] for record in binary_records)

    gray_global_K = gray_total_SSE / (np.pi**2 * gray_total_weight)
    binary_global_K = binary_total_SSE / (np.pi**2 * binary_total_weight)
    gray_global_sqrt_K = np.sqrt(gray_global_K)
    binary_global_sqrt_K = np.sqrt(binary_global_K)

    gray_zero_layers = gray_layers[nontrivial & (gray_K == 0.0)]
    binary_zero_layers = binary_layers[nontrivial & (binary_K == 0.0)]
    print("Gray exact-zero layers:", gray_zero_layers.tolist())
    print("Binary exact-zero layers:", binary_zero_layers.tolist())

    median_reduction = (1.0 - gray_median_K / binary_median_K) * 100.0
    global_reduction = (1.0 - gray_global_K / binary_global_K) * 100.0
    global_rmse_reduction = (
        1.0 - gray_global_sqrt_K / binary_global_sqrt_K
    ) * 100.0

    print(f"Gray median K       = {gray_median_K:.6e}")
    print(f"Binary median K     = {binary_median_K:.6e}")
    print(f"Median reduction    = {median_reduction:.2f}%")
    print(f"Gray global K       = {gray_global_K:.6e}")
    print(f"Binary global K     = {binary_global_K:.6e}")
    print(f"Global-K reduction  = {global_reduction:.2f}%")
    print(f"Gray global sqrt(K) = {gray_global_sqrt_K:.6e}")
    print(f"Binary global sqrt(K) = {binary_global_sqrt_K:.6e}")
    print(f"Global-RMSE reduction = {global_rmse_reduction:.2f}%")

    metric_labels = [
        "Median\nlayerwise $\\mathcal{K}_j$",
        "Global\n$\\mathcal{K}$",
        "Global\n$\\sqrt{\\mathcal{K}}$",
    ]
    gray_absolute = np.array(
        [gray_median_K, gray_global_K, gray_global_sqrt_K]
    )
    binary_absolute = np.array(
        [binary_median_K, binary_global_K, binary_global_sqrt_K]
    )
    gray_relative = gray_absolute / binary_absolute
    binary_relative = np.ones_like(binary_absolute)
    reductions = (1.0 - gray_relative) * 100.0

    fig_width = 18 / 2.54
    fig_height = 18 / 2.54
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax_width = 10 / 2.54 / fig_width
    ax_height = 5 / 2.54 / fig_height
    ax = fig.add_axes([0, 0, ax_width, ax_height])

    x = np.arange(len(metric_labels))
    bar_width = 0.34
    gray_bars = ax.bar(
        x - bar_width / 2,
        gray_relative,
        width=bar_width,
        color="#2878B5",
        edgecolor="black",
        linewidth=0.5,
        label="Gray code",
    )
    binary_bars = ax.bar(
        x + bar_width / 2,
        binary_relative,
        width=bar_width,
        color="#F28E2B",
        edgecolor="black",
        linewidth=0.5,
        label="Standard binary",
    )
    ax.axhline(1.0, color="0.35", linestyle="--", linewidth=0.9)
    ax.set_xticks(x, metric_labels)
    ax.set_ylabel("Error relative to standard binary")
    ax.set_ylim(0.0, 1.22)
    ax.legend(
        ncol=1,
        frameon=True,
        labelspacing=0.2,
        handlelength=2,
        handletextpad=0.5,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        fontsize=8,
    )

    for index, bar in enumerate(gray_bars):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + 0.01,
            rf"$-{reductions[index]:.1f}\%$",
            ha="center",
            va="bottom",
            color="#2878B5",
            fontsize=8,
        )
    for bar in binary_bars:
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            "Reference",
            ha="center",
            va="bottom",
            color="#F28E2B",
            fontsize=8,
        )
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(
        axis="x",
        which="major",
        direction="out",
        top=False,
        right=False,
        length=3,
        width=0.5,
        pad=1.5,
    )
    ax.tick_params(
        axis="y",
        which="major",
        direction="out",
        top=False,
        right=False,
        length=3,
        width=0.5,
        pad=1,
    )
    ax.tick_params(
        which="minor",
        direction="out",
        top=False,
        right=False,
        length=1.5,
        width=0.5,
    )

    pdf_path = (
        Path("figures/gray_binary_error_histogram.pdf")
        if pdf_path is None
        else Path(pdf_path)
    )
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.show()
    print(f"Saved PDF: {pdf_path}")
    return fig
