"""Central configuration for the geometric-encoding workflow."""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=False)
class WorkflowConfig:
    """Parameters retained from the original notebook."""

    nx: int = 10
    ny: int = 10
    nz: int = 10
    k_cutoff: float = 15.0
    seed_spin_up: int = 2025
    seed_spin_down: int = 2026
    max_parallel_threads: int = 16
    velocity_method: str = "FDM"
    derivative_method: str = "FDM"
    fontsize: float = 8.0
    data_dir: Path = Path("data/3D")
    figure_dir: Path = Path("figures")

    @property
    def n_total(self):
        return self.nx + self.ny + self.nz

    @property
    def grid_size(self):
        return 2**self.nx

    def ensure_output_directories(self):
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.figure_dir.mkdir(parents=True, exist_ok=True)
