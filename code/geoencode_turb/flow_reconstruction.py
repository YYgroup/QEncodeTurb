"""Pauli-spinor reconstruction and primary spectral observables."""

from dataclasses import dataclass
import gc

import numpy as np

from .config import WorkflowConfig
from .quantum_encoding import EncodingResult, restore_amps_from_gray


@dataclass
class FlowFieldResult:
    """Reconstructed Pauli-spinor and fluid fields."""

    N: int
    kx: np.ndarray
    ky: np.ndarray
    kz: np.ndarray
    KX: np.ndarray
    KY: np.ndarray
    KZ: np.ndarray
    K2: np.ndarray
    ik2: np.ndarray
    psi1_spec: np.ndarray
    psi2_spec: np.ndarray
    psi1: np.ndarray
    psi2: np.ndarray
    rho: np.ndarray
    ux: np.ndarray
    uy: np.ndarray
    uz: np.ndarray


def compute_velocity(
    psi1,
    psi2,
    psi1_spec,
    psi2_spec,
    switch,
    KX,
    KY,
    KZ,
    N,
):
    """Compute density and velocity with the original SM or FDM scheme."""
    rho = np.abs(psi1) ** 2 + np.abs(psi2) ** 2

    if switch == "SM":
        dpsi1_x = np.fft.fftn(1j * KX * psi1_spec) / np.sqrt(N**3)
        dpsi1_y = np.fft.fftn(1j * KY * psi1_spec) / np.sqrt(N**3)
        dpsi1_z = np.fft.fftn(1j * KZ * psi1_spec) / np.sqrt(N**3)
        dpsi2_x = np.fft.fftn(1j * KX * psi2_spec) / np.sqrt(N**3)
        dpsi2_y = np.fft.fftn(1j * KY * psi2_spec) / np.sqrt(N**3)
        dpsi2_z = np.fft.fftn(1j * KZ * psi2_spec) / np.sqrt(N**3)
        ux = np.real(
            np.real(psi1) * np.imag(dpsi1_x)
            - np.imag(psi1) * np.real(dpsi1_x)
            + np.real(psi2) * np.imag(dpsi2_x)
            - np.imag(psi2) * np.real(dpsi2_x)
        ) / rho
        uy = np.real(
            np.real(psi1) * np.imag(dpsi1_y)
            - np.imag(psi1) * np.real(dpsi1_y)
            + np.real(psi2) * np.imag(dpsi2_y)
            - np.imag(psi2) * np.real(dpsi2_y)
        ) / rho
        uz = np.real(
            np.real(psi1) * np.imag(dpsi1_z)
            - np.imag(psi1) * np.real(dpsi1_z)
            + np.real(psi2) * np.imag(dpsi2_z)
            - np.imag(psi2) * np.real(dpsi2_z)
        ) / rho
    elif switch == "FDM":
        Jx = np.imag(
            np.conj(psi1)
            * (np.roll(psi1, -1, axis=0) - np.roll(psi1, 1, axis=0))
        ) + np.imag(
            np.conj(psi2)
            * (np.roll(psi2, -1, axis=0) - np.roll(psi2, 1, axis=0))
        )
        Jy = np.imag(
            np.conj(psi1)
            * (np.roll(psi1, -1, axis=1) - np.roll(psi1, 1, axis=1))
        ) + np.imag(
            np.conj(psi2)
            * (np.roll(psi2, -1, axis=1) - np.roll(psi2, 1, axis=1))
        )
        Jz = np.imag(
            np.conj(psi1)
            * (np.roll(psi1, -1, axis=2) - np.roll(psi1, 1, axis=2))
        ) + np.imag(
            np.conj(psi2)
            * (np.roll(psi2, -1, axis=2) - np.roll(psi2, 1, axis=2))
        )
        h = 2 * np.pi / N
        ux = Jx / (2 * h) / rho
        uy = Jy / (2 * h) / rho
        uz = Jz / (2 * h) / rho
        del Jx, Jy, Jz
        gc.collect()
    else:
        raise ValueError("switch must be 'SM' or 'FDM'.")

    return rho, ux, uy, uz


def compute_Ek(ux, uy, uz, ik2, N):
    """Compute the velocity energy spectrum with the original shell sum."""
    ux_spec = np.fft.fftn(ux) / N**3
    uy_spec = np.fft.fftn(uy) / N**3
    uz_spec = np.fft.fftn(uz) / N**3
    energy_spec = (
        np.abs(ux_spec) ** 2 + np.abs(uy_spec) ** 2 + np.abs(uz_spec) ** 2
    ) / 2

    nek = N
    Ek = np.zeros(nek)
    k = np.linspace(1, nek, nek)
    for i in range(1, nek + 1):
        Ek[i - 1] = k[i - 1] * np.sum(energy_spec[ik2 == i])
    return k, Ek


def compute_spin_fields(psi1, psi2, rho):
    """Compute the three components of the normalized spin-vector field."""
    s1 = 2 * (
        np.real(psi1) * np.real(psi2) + np.imag(psi1) * np.imag(psi2)
    ) / rho
    s2 = 2 * (
        np.real(psi1) * np.imag(psi2) - np.imag(psi1) * np.real(psi2)
    ) / rho
    s3 = (np.abs(psi1) ** 2 - np.abs(psi2) ** 2) / rho
    return s1, s2, s3


def compute_Es(s, ik2, N):
    """Compute a spin-component spectrum with the original shell sum."""
    s_spec = np.fft.fftn(s) / N**3
    energy_spec = np.abs(s_spec) ** 2 / 2

    nek = N
    Es = np.zeros(nek)
    k = np.linspace(1, nek, nek)
    for i in range(1, nek + 1):
        Es[i - 1] = k[i - 1] * np.sum(energy_spec[ik2 == i])
    return k, Es


def reconstruct_flow_field(encoding, config=None):
    """Restore both spinor components and evaluate the fluid observables."""
    if not isinstance(encoding, EncodingResult):
        raise TypeError("encoding must be an EncodingResult instance.")
    config = WorkflowConfig() if config is None else config

    nx, ny, nz = config.nx, config.ny, config.nz
    N = 2**nx
    kx = np.fft.fftfreq(N) * N
    ky = np.fft.fftfreq(N) * N
    kz = np.fft.fftfreq(N) * N
    KX, KY, KZ = np.meshgrid(kx, ky, kz)
    K2 = KX**2 + KY**2 + KZ**2
    ik2 = np.round(np.sqrt(K2))

    psi1_spec = restore_amps_from_gray(
        encoding.sta_vec_1_data,
        nx,
        ny,
        nz,
    ).reshape((2**nz, 2**ny, 2**nx)).T
    psi2_spec = restore_amps_from_gray(
        encoding.sta_vec_2_data,
        nx,
        ny,
        nz,
    ).reshape((2**nz, 2**ny, 2**nx)).T

    # Release the statevectors after their physical ordering has been restored.
    encoding.sta_vec_1_data = None
    encoding.sta_vec_2_data = None
    gc.collect()

    psi1 = np.fft.ifftn(psi1_spec)
    psi2 = np.fft.ifftn(psi2_spec)
    rho, ux, uy, uz = compute_velocity(
        psi1,
        psi2,
        psi1_spec,
        psi2_spec,
        config.velocity_method,
        KX,
        KY,
        KZ,
        N,
    )
    # print(np.amin(rho), np.amax(rho), np.mean(rho))

    return FlowFieldResult(
        N=N,
        kx=kx,
        ky=ky,
        kz=kz,
        KX=KX,
        KY=KY,
        KZ=KZ,
        K2=K2,
        ik2=ik2,
        psi1_spec=psi1_spec,
        psi2_spec=psi2_spec,
        psi1=psi1,
        psi2=psi2,
        rho=rho,
        ux=ux,
        uy=uy,
        uz=uz,
    )
