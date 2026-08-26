from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, hex2color
import numpy as np
from scipy.interpolate import interp1d

from .config import WorkflowConfig


def create_smooth_cmap(colors, name="smooth_cmap", N=256):
    """Interpolate a smooth color map through the supplied hex colors."""
    rgb_colors = np.array([hex2color(color) for color in colors])
    original_nodes = np.linspace(0, 1, len(rgb_colors))
    interpolator = interp1d(original_nodes, rgb_colors, kind="cubic", axis=0)
    new_nodes = np.linspace(0, 1, N)
    new_colors = interpolator(new_nodes)
    new_colors = np.clip(new_colors, 0, 1)
    return LinearSegmentedColormap.from_list(name, new_colors)


cmap_BuRd = create_smooth_cmap(
    ["#0f5fb5", "#86c5e5", "#ffffff", "#f69d74", "#b80b1f"],
    "cmap_BuRd",
)


def plot_encoded_amplitudes(encoding, config=None, output_path="figures/encoded_amplitudes.pdf"):
    """Plot modewise and radially averaged encoded amplitudes."""
    config = WorkflowConfig() if config is None else config
    k_plot = encoding.k_plot
    amp_plot = encoding.amp_plot
    c_sim = encoding.c_sim
    v_sim = encoding.v_sim

    fig_width = 18 / 2.54
    fig_height = 18 / 2.54
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax_width = 4.5 / 2.54 / fig_width
    ax_height = 3 / 2.54 / fig_height
    ax = fig.add_axes([0, 0, ax_width, ax_height])
    ax.grid(color="k", linestyle="-", linewidth=0.5, alpha=0.1)

    # Downsample the scatter points to reduce rendering time.
    step = 2 ** (config.nx - 5) if len(k_plot) > 100000 else 1
    is_valid_mode = k_plot <= config.k_cutoff
    plt.scatter(
        k_plot[is_valid_mode][::step],
        amp_plot[is_valid_mode][::step],
        marker=".",
        alpha=1,
        s=3,
        color="lightblue",
        rasterized=True,
    )
    plt.scatter(
        k_plot[~is_valid_mode][::step],
        amp_plot[~is_valid_mode][::step],
        marker=".",
        alpha=1,
        s=3,
        color="lightgray",
        rasterized=True,
    )
    plt.plot(
        c_sim,
        v_sim,
        linestyle="-",
        color="k",
        linewidth=1,
        label=r"Averaged amplitude",
    )

    if len(c_sim) > 0:
        valid_idxs = np.where(c_sim < config.k_cutoff)[0]
        if len(valid_idxs) > 0:
            ref_idx = valid_idxs[0]
            ref_k = c_sim[ref_idx]
            ref_v = v_sim[ref_idx]
            k_theory = np.logspace(
                np.log10(np.min(c_sim)),
                np.log10(np.max(config.k_cutoff)),
                100,
            )
            v_theory = ref_v * (k_theory / ref_k) ** 0
            plt.plot(
                k_theory,
                0.1 * v_theory,
                linestyle="-",
                color="gray",
                linewidth=1,
                label=r"Theory",
            )
            plt.axvline(
                x=config.k_cutoff,
                color="gray",
                linestyle="--",
                dashes=(2, 2),
                linewidth=1,
            )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1, 1e2)
    plt.xlabel(
        r"Wavevector magnitude, $|\boldsymbol{k}|$",
        labelpad=0,
        size=config.fontsize,
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
        length=0,
        width=0.5,
    )
    plt.savefig(
        output_path,
        transparent=True,
        orientation="portrait",
        format="pdf",
        bbox_inches="tight",
        dpi=600,
    )
    return fig


def plot_probability_histogram(encoding, config=None):
    """Compare target and prepared probability distributions."""
    config = WorkflowConfig() if config is None else config
    target_values = encoding.target_probs_phys.flatten()
    simulated_values = encoding.amp_sim_phys.flatten() ** 2
    return probability_histogram_from_arrays(
        target_values,
        simulated_values,
        config.fontsize,
    )


def probability_histogram_from_arrays(t_vals, s_vals, fontsize=8):
    """Plot the original probability-density comparison from full arrays."""
    mask = t_vals > 1e-20
    t_vals = t_vals[mask]
    s_vals = s_vals[mask]

    fig_width = 18 / 2.54
    fig_height = 18 / 2.54
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax_width = 6 / 2.54 / fig_width
    ax_height = 4 / 2.54 / fig_height
    ax = fig.add_axes([0, 0, ax_width, ax_height])
    ax.grid(color="k", linestyle="-", linewidth=0.5, alpha=0.1)
    bins = np.logspace(np.log10(min(t_vals)), np.log10(max(t_vals)), 50)

    plt.hist(
        t_vals,
        bins=bins,
        alpha=0.5,
        color="gray",
        label="Target",
        density=True,
        log=True,
    )
    plt.hist(
        s_vals,
        bins=bins,
        alpha=0.5,
        color="dodgerblue",
        histtype="step",
        linewidth=1,
        label="Linear ansatz",
        density=True,
        log=True,
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
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Probability amplitude, $P(k)$", labelpad=0, fontsize=fontsize)
    plt.ylabel("Probability density function", labelpad=1, fontsize=fontsize)
    plt.legend()
    return fig

def _add_background_plane(plotter, position, normal, size, color="black"):
    import pyvista as pv

    plane = pv.Plane(
        center=position,
        direction=normal,
        i_size=size[0],
        j_size=size[1],
    )
    plotter.add_mesh(plane, color=color, opacity=1)


def _add_standard_light(plotter):
    import pyvista as pv

    light = pv.Light()
    light.intensity = 0.5
    light.positional = True
    light.position = (6, 9, 9)
    plotter.add_light(light)


# def render_vortex_isosurface(flow, vorticity):
#     """Render the vorticity isosurface colored by helicity."""
#     import pyvista as pv

#     N = flow.N
#     spacing = tuple(2 * np.pi / (N - 1) for _ in range(3))
#     grid = pv.ImageData(dimensions=(N, N, N), spacing=spacing, origin=(0, 0, 0))
#     grid.point_data["vor"] = vorticity.magnitude.ravel(order="F")
#     grid.point_data["helicity"] = vorticity.helicity.ravel(order="F")

#     plotter = pv.Plotter(notebook=False, off_screen=True)
#     contours = grid.contour(isosurfaces=[150], scalars="vor")
#     if hasattr(contours, "smooth"):
#         contours = contours.smooth(n_iter=100, relaxation_factor=0.1)
#     bounds = contours.bounds
#     plotter.add_mesh(
#         contours,
#         scalars="helicity",
#         cmap=cmap_BuRd,
#         clim=[-1e3, 1e3],
#         opacity=0.7,
#     )
#     plotter.add_mesh(
#         pv.Box(bounds=bounds),
#         color="white",
#         style="wireframe",
#         line_width=10,
#     )
#     plotter.remove_scalar_bar()

#     planes = [
#         ((2 * np.pi, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2), (1, 0, 0)),
#         (((bounds[0] + bounds[1]) / 2, 0, (bounds[4] + bounds[5]) / 2), (0, -1, 0)),
#         (((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, 0), (0, 0, -1)),
#     ]
#     size = (
#         bounds[1] - bounds[0],
#         bounds[3] - bounds[2],
#         bounds[5] - bounds[4],
#     )
#     for position, normal in planes:
#         _add_background_plane(plotter, position, normal, size)

#     plotter.enable_anti_aliasing()
#     _add_standard_light(plotter)
#     plotter.camera.azimuth = 70
#     plotter.camera.elevation = -15
#     plotter.disable_parallel_projection()
#     plotter.show(
#         screenshot="figures/iso-vor_try.png",
#         window_size=[3000, 3000],
#     )

def render_vortex_isosurface(
    flow,
    vorticity,
    iso_value=22.0,
    helicity_limit=70.0,
    output_path="figures/iso-vor_nq=15.png",
):
    """Render the vorticity isosurface colored by helicity."""
    import pyvista as pv

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    N = flow.N
    spacing = tuple(2 * np.pi / (N - 1) for _ in range(3))
    grid = pv.ImageData(dimensions=(N, N, N), spacing=spacing, origin=(0, 0, 0))
    grid.point_data["vor"] = vorticity.magnitude.ravel(order="F")
    grid.point_data["helicity"] = vorticity.helicity.ravel(order="F")

    plotter = pv.Plotter(notebook=False, off_screen=True)
    contours = grid.contour(isosurfaces=[iso_value], scalars="vor")
    if contours.n_points == 0:
        raise ValueError(
            f"No vorticity isosurface exists at iso_value={iso_value}."
        )
    if hasattr(contours, "smooth"):
        contours = contours.smooth(n_iter=100, relaxation_factor=0.1)
    bounds = contours.bounds
    plotter.add_mesh(
        contours,
        scalars="helicity",
        cmap=cmap_BuRd,
        clim=[-helicity_limit, helicity_limit],
        opacity=0.7,
    )
    plotter.add_mesh(
        pv.Box(bounds=bounds),
        color="white",
        style="wireframe",
        line_width=10,
    )
    plotter.remove_scalar_bar()

    planes = [
        ((2 * np.pi, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2), (1, 0, 0)),
        (((bounds[0] + bounds[1]) / 2, 0, (bounds[4] + bounds[5]) / 2), (0, -1, 0)),
        (((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, 0), (0, 0, -1)),
    ]
    size = (
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
    )
    for position, normal in planes:
        _add_background_plane(plotter, position, normal, size)

    plotter.enable_anti_aliasing()
    _add_standard_light(plotter)
    plotter.camera.azimuth = 70
    plotter.camera.elevation = -15
    plotter.disable_parallel_projection()
    plotter.show(
        screenshot=str(output_path),
        window_size=[3000, 3000],
    )
    print(f"Vorticity isosurface saved successfully: {output_path}")
    return output_path


def render_spin_isosurface(flow, spin_fields, vorticity):
    """Render the s1 isosurface colored by vorticity magnitude."""
    import pyvista as pv

    N = flow.N
    spacing = tuple(2 * np.pi / (N - 1) for _ in range(3))
    grid = pv.ImageData(dimensions=(N, N, N), spacing=spacing, origin=(0, 0, 0))
    grid.point_data["s"] = spin_fields[0].ravel(order="F")
    grid.point_data["vor"] = vorticity.magnitude.ravel(order="F")

    plotter = pv.Plotter(notebook=False, off_screen=True)
    contours = grid.contour(isosurfaces=[-0.8], scalars="s")
    if hasattr(contours, "smooth"):
        contours = contours.smooth(n_iter=100, relaxation_factor=0.1)
    bounds = contours.bounds
    plotter.add_mesh(
        contours,
        scalars="vor",
        cmap="cividis",
        clim=[0, 350],
        opacity=0.7,
    )
    plotter.add_mesh(
        pv.Box(bounds=bounds),
        color="white",
        style="wireframe",
        line_width=10,
    )
    plotter.remove_scalar_bar()
    size = (
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
    )
    planes = [
        ((2 * np.pi, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2), (1, 0, 0)),
        (((bounds[0] + bounds[1]) / 2, 0, (bounds[4] + bounds[5]) / 2), (0, -1, 0)),
        (((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, 0), (0, 0, -1)),
    ]
    for position, normal in planes:
        _add_background_plane(plotter, position, normal, size)

    plotter.enable_anti_aliasing()
    _add_standard_light(plotter)
    plotter.camera.azimuth = 70
    plotter.camera.elevation = -15
    plotter.disable_parallel_projection()
    plotter.show(
        screenshot="figures/iso-s_nq=30.png",
        window_size=[3000, 3000],
    )


def _render_spectral_volume(data, clim, screenshot):
    import pyvista as pv

    N = data.shape[0]
    grid = pv.ImageData()
    grid.dimensions = (N, N, N)
    grid.spacing = (1, 1, 1)
    grid["amplitude"] = data.flatten(order="F")
    plotter = pv.Plotter(notebook=False, off_screen=False)
    plotter.add_volume(
        grid,
        scalar_bar_args={"title": "Amplitude"},
        cmap="Blues",
        clim=clim,
        log_scale=True,
        show_scalar_bar=False,
    )
    _add_standard_light(plotter)
    plotter.camera.azimuth = 70
    plotter.camera.elevation = -15
    plotter.disable_parallel_projection()
    plotter.show(screenshot=screenshot, window_size=[3000, 3000])


def render_spinor_spectral_volume(flow):
    """Render the normalized spin-up spectral probability."""
    data = np.abs(np.fft.fftshift(flow.psi1_spec)) ** 2
    data = data / np.amax(data)
    _render_spectral_volume(data, [1e-6, 0.8], "figures/psi1_vol.png")


def render_density_spectral_volume(flow):
    """Render the normalized density spectrum."""
    data = np.abs(np.fft.fftshift(np.fft.fftn(flow.rho))) ** 2
    data = data / np.unique(data)[-2]
    _render_spectral_volume(data, [1e-6, 0.3], "figures/psi1_vol.png")


def render_velocity_spectral_volume(flow):
    """Render the normalized x-velocity spectrum."""
    data = np.abs(np.fft.fftshift(np.fft.fftn(flow.ux)) / flow.N**3) ** 2
    data = data / np.unique(data)[-1]
    _render_spectral_volume(data, [1e-6, 0.1], "figures/psi1_vol.png")


def _render_physical_volume(data, cmap, clim, opacity, screenshot):
    import pyvista as pv

    N = data.shape[0]
    grid = pv.ImageData()
    grid.dimensions = (N, N, N)
    grid.spacing = (1, 1, 1)
    grid["amplitude"] = data.flatten(order="F")
    plotter = pv.Plotter(notebook=False, off_screen=True)
    volume = plotter.add_volume(
        grid,
        scalar_bar_args={"title": "Amplitude"},
        cmap=cmap,
        clim=clim,
        opacity=opacity,
        show_scalar_bar=False,
    )
    bounds = volume.bounds
    size = (
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
    )
    planes = [
        ((N, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2), (-1, 0, 0)),
        (((bounds[0] + bounds[1]) / 2, 0, (bounds[4] + bounds[5]) / 2), (0, -1, 0)),
        (((bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, 0), (0, 0, -1)),
    ]
    for position, normal in planes:
        _add_background_plane(plotter, position, normal, size)
    plotter.add_mesh(grid.outline(), color="grey", line_width=5)
    _add_standard_light(plotter)
    plotter.camera.azimuth = 70
    plotter.camera.elevation = -15
    plotter.disable_parallel_projection()
    plotter.show(screenshot=screenshot, window_size=[3000, 3000])


def render_vorticity_volume(vorticity):
    """Render the physical-space vorticity magnitude."""
    _render_physical_volume(
        vorticity.magnitude,
        "Blues",
        [0, 500],
        "linear",
        "figures/vor_volume_render_nq=27.png",
    )


def render_density_volume(flow):
    """Render the physical-space density field."""
    _render_physical_volume(
        flow.rho,
        "Reds",
        [0, 1e-16],
        "sigmoid",
        "figures/rho_volume_render_nq=27.png",
    )
