# Copyright 2026 Ivan Lobato / NeuralSoftX
# SPDX-License-Identifier: Apache-2.0
# Email: ivan.lobato@neuralsoftx.com
"""Evaluate the released scattering factors, densities and potentials.

Author: Ivan Lobato
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
DATA_ROOT = HERE / "data"
DEFAULT_COEFFICIENTS = DATA_ROOT / "coefficients.h5"
DEFAULT_REFERENCE = DATA_ROOT / "reference" / "reference_dirac_fock.h5"
DEFAULT_CERTIFICATE = DATA_ROOT / "global_positivity_certificate.json"
A_0 = 0.52917721077817892
INVERSE_KAPPA = 47.877642131544313031


@dataclass(frozen=True)
class Coefficients:
    """Immutable coefficient set for one element and one released model."""

    model: str
    atomic_number: int
    symbol: str
    a: np.ndarray
    b: np.ndarray
    is_dirac_pade: np.ndarray

    @property
    def n_terms(self) -> int:
        """Return the number of basis terms."""
        return int(self.a.size)


def available_models(path: str | Path = DEFAULT_COEFFICIENTS) -> tuple[str, ...]:
    """Return the model names present in a coefficient archive."""
    archive = _archive_path(path)
    with h5py.File(archive, "r") as handle:
        return tuple(handle.keys())


def available_elements(
    model: str = "adaptive_dirac_pade",
    path: str | Path = DEFAULT_COEFFICIENTS,
) -> tuple[tuple[int, str], ...]:
    """Return ``(atomic_number, symbol)`` pairs for one model."""
    archive = _archive_path(path)
    with h5py.File(archive, "r") as handle:
        group = _model_group(handle, model, archive)
        elements = [
            (int(group[key].attrs["Z"]), str(group[key].attrs["symbol"]))
            for key in group
        ]
    return tuple(sorted(elements))


def load_coefficients(
    atomic_number: int,
    model: str = "adaptive_dirac_pade",
    path: str | Path = DEFAULT_COEFFICIENTS,
) -> Coefficients:
    """Load one element's released float32 coefficients as float64 arrays."""
    if isinstance(atomic_number, bool) or not isinstance(
        atomic_number, (int, np.integer)
    ):
        raise TypeError("atomic_number must be an integer in the released range")
    z = int(atomic_number)
    archive = _archive_path(path)
    with h5py.File(archive, "r") as handle:
        group = _model_group(handle, model, archive)
        matches = [key for key in group if int(group[key].attrs["Z"]) == z]
        if len(matches) != 1:
            available = sorted(int(group[key].attrs["Z"]) for key in group)
            if not available:
                raise ValueError(f"model {model!r} contains no elements in {archive}")
            raise ValueError(
                f"atomic_number={z} is unavailable for model {model!r}; "
                f"released range is {available[0]}..{available[-1]}"
            )
        item = group[matches[0]]
        a = np.asarray(item["a"][...], dtype=np.float64)
        b = np.asarray(item["b"][...], dtype=np.float64)
        term_types = np.asarray(item["type"][...])
        symbol = str(item.attrs["symbol"])
        if "n_t" not in item.attrs:
            raise ValueError(
                f"missing n_t attribute for model={model!r}, Z={z} in {archive}"
            )
        declared_n_terms_value = item.attrs["n_t"]
        if isinstance(declared_n_terms_value, bool) or not isinstance(
            declared_n_terms_value,
            (int, np.integer),
        ):
            raise ValueError(
                f"n_t must be an integer for model={model!r}, Z={z}; "
                f"got {declared_n_terms_value!r}"
            )
        declared_n_terms = int(declared_n_terms_value)

    if a.ndim != 1 or b.shape != a.shape or term_types.shape != a.shape:
        raise ValueError(
            f"invalid coefficient schema for model={model!r}, Z={z}: "
            f"a={a.shape}, b={b.shape}, type={term_types.shape}"
        )
    if a.size == 0:
        raise ValueError(f"empty coefficient set for model={model!r}, Z={z}")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError(f"non-finite coefficients for model={model!r}, Z={z}")
    if not np.all(b > 0.0):
        raise ValueError(f"non-positive width coefficient for model={model!r}, Z={z}")
    if not np.all(np.isin(term_types, (b"NR", b"DP"))):
        raise ValueError(f"unknown basis type for model={model!r}, Z={z}")
    if declared_n_terms != a.size:
        raise ValueError(
            f"n_t={declared_n_terms} does not match {a.size} coefficients "
            f"for model={model!r}, Z={z}"
        )
    if model == "adaptive_dirac_pade":
        valid_layout = (
            term_types.size > 0
            and np.all(term_types[:-1] == b"NR")
            and term_types[-1] == b"DP"
        )
        if not valid_layout:
            raise ValueError(
                f"model={model!r}, Z={z} must contain NR terms followed by "
                "exactly one final DP term"
            )
    elif model in ("adaptive_hydrogenic", "five_term_hydrogenic") and not np.all(term_types == b"NR"):
        raise ValueError(f"model={model!r}, Z={z} must contain only NR terms")

    is_dirac_pade = np.asarray(term_types == b"DP", dtype=bool)
    for array in (a, b, is_dirac_pade):
        array.setflags(write=False)
    return Coefficients(model, z, symbol, a, b, is_dirac_pade)


def electron_scattering_factor(
    coefficients: Coefficients,
    g: float | np.ndarray,
) -> float | np.ndarray:
    """Evaluate the electron scattering factor ``f_e(g)`` in ångström."""
    values, shape, scalar = _coordinates(g, name="g", allow_zero=True)
    g2 = values[:, None] ** 2
    a, b = coefficients.a, coefficients.b
    denominator = 1.0 + b * g2
    nr = a * (2.0 + b * g2) / denominator**2
    u = b * g2
    dp = a * (3.0 + 3.0 * u + u * u) / denominator**3
    result = np.where(coefficients.is_dirac_pade, dp, nr).sum(axis=1)
    return _restore_shape(result, shape, scalar)


def xray_scattering_factor(
    coefficients: Coefficients,
    g: float | np.ndarray,
) -> float | np.ndarray:
    """Evaluate the X-ray scattering factor ``f_x(g)`` in electrons."""
    values, shape, scalar = _coordinates(g, name="g", allow_zero=True)
    g2 = values[:, None] ** 2
    a, b = coefficients.a, coefficients.b
    denominator = 1.0 + b * g2
    amplitude = 2.0 * np.pi**2 * A_0 * a / b
    nr = amplitude / denominator**2
    dp = amplitude / denominator**3
    result = np.where(coefficients.is_dirac_pade, dp, nr).sum(axis=1)
    return _restore_shape(result, shape, scalar)


def electron_density(
    coefficients: Coefficients,
    r: float | np.ndarray,
) -> float | np.ndarray:
    """Evaluate the electron density ``rho(r)`` in electrons per ångström cubed."""
    values, shape, scalar = _coordinates(r, name="r", allow_zero=True)
    a, b = coefficients.a, coefficients.b
    decay = 2.0 * np.pi / np.sqrt(b)
    amplitude = 2.0 * np.pi**4 * A_0 * a / b**2.5
    x = values[:, None] * decay
    exponential = np.exp(-x)
    nr = amplitude * exponential
    dp = (amplitude / 4.0) * (1.0 + x) * exponential
    result = np.where(coefficients.is_dirac_pade, dp, nr).sum(axis=1)
    return _restore_shape(result, shape, scalar)


def electrostatic_potential(
    coefficients: Coefficients,
    r: float | np.ndarray,
) -> float | np.ndarray:
    """Evaluate the radial electrostatic potential ``V(r)``."""
    values, shape, scalar = _coordinates(r, name="r", allow_zero=False)
    a, b = coefficients.a, coefficients.b
    decay = 2.0 * np.pi / np.sqrt(b)
    amplitude = np.pi**2 * INVERSE_KAPPA * a / b**1.5
    x = values[:, None] * decay
    exponential = np.exp(-x)
    nr = amplitude * exponential * (2.0 / x + 1.0)
    dp = (amplitude / 4.0) * exponential * (8.0 / x + 5.0 + x)
    result = np.where(coefficients.is_dirac_pade, dp, nr).sum(axis=1)
    return _restore_shape(result, shape, scalar)


def projected_potential(
    coefficients: Coefficients,
    radius: float | np.ndarray,
) -> float | np.ndarray:
    """Evaluate the projected potential ``V(R)`` for strictly positive ``R``."""
    try:
        from scipy.special import kv
    except ImportError as exc:
        raise ImportError(
            "projected_potential requires scipy; install the declared package "
            "dependencies"
        ) from exc

    values, shape, scalar = _coordinates(
        radius,
        name="radius",
        allow_zero=False,
    )
    a, b = coefficients.a, coefficients.b
    decay = 2.0 * np.pi / np.sqrt(b)
    amplitude = 2.0 * np.pi**2 * INVERSE_KAPPA * a / b**1.5
    radius_column = values[:, None]
    x = radius_column * decay
    k0 = kv(0, x)
    k1 = kv(1, x)
    nr = amplitude * (2.0 * k0 / decay + radius_column * k1)
    dp = (amplitude / 8.0) * (
        (16.0 / decay + decay * radius_column**2) * k0
        + 10.0 * radius_column * k1
        + decay * radius_column**2 * kv(2, x)
    )
    result = np.where(coefficients.is_dirac_pade, dp, nr).sum(axis=1)
    return _restore_shape(result, shape, scalar)


def _archive_path(path: str | Path) -> Path:
    archive = Path(path).expanduser().resolve()
    if not archive.is_file():
        raise FileNotFoundError(f"coefficient archive not found: {archive}")
    return archive


def _model_group(
    handle: h5py.File,
    model: str,
    archive: Path,
) -> h5py.Group:
    if not isinstance(model, str) or not model:
        raise TypeError("model must be a non-empty string")
    if model not in handle:
        raise ValueError(
            f"model {model!r} is unavailable in {archive}; "
            f"choose from {tuple(handle.keys())}"
        )
    return handle[model]


def _coordinates(
    value: float | np.ndarray,
    *,
    name: str,
    allow_zero: bool,
) -> tuple[np.ndarray, tuple[int, ...], bool]:
    array = np.asarray(value, dtype=np.float64)
    scalar = array.ndim == 0
    shape = array.shape
    flat = array.reshape(-1)
    if not np.all(np.isfinite(flat)):
        raise ValueError(f"{name} must contain only finite values")
    if allow_zero:
        valid = flat >= 0.0
        requirement = "non-negative"
    else:
        valid = flat > 0.0
        requirement = "strictly positive"
    if not np.all(valid):
        raise ValueError(f"{name} must contain only {requirement} values")
    return flat, shape, scalar


def _restore_shape(
    value: np.ndarray,
    shape: tuple[int, ...],
    scalar: bool,
) -> float | np.ndarray:
    output = np.asarray(value, dtype=np.float64).reshape(shape)
    return float(output) if scalar else output
