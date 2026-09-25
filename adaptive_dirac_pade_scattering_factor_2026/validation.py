#!/usr/bin/env python3
# Copyright 2026 Ivan Lobato / NeuralSoftX
# SPDX-License-Identifier: Apache-2.0
# Email: ivan.lobato@neuralsoftx.com
"""Certify global positivity of every released analytic parametrisation.

Author: Ivan Lobato
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp
from numpy.polynomial import polynomial
from scipy.optimize import brentq
from scipy.special import logsumexp

from .scattering_factors import (
    A_0,
    DEFAULT_CERTIFICATE,
    DEFAULT_COEFFICIENTS,
    available_elements,
    load_coefficients,
)


MODEL_NAMES = ('adaptive_dirac_pade', 'adaptive_hydrogenic', 'five_term_hydrogenic')
QUANTITIES_RECIPROCAL = ('fe', 'fx')
QUANTITIES_REAL = ('rho', 'w')
CERTIFICATE_NAME = 'global_positivity_certificate.json'
CERTIFICATE_SCHEMA = 1
TAIL_RELATIVE_TARGET = 0.24
TAIL_FINITE_OVERLAP_FACTOR = 1.05
FINITE_RELATIVE_TOLERANCE = 1.0e-12
MAX_INTERVAL_DEPTH = 60
RECIPROCAL_GRID_POINTS = 40001
RECIPROCAL_X_MIN = 1.0e-20
RECIPROCAL_X_MAX = 1.0e24
EXPECTED_ATOMIC_NUMBERS = tuple(range(1, 119))


@dataclass(frozen=True)
class component_group:
    """One exact-width exponential-polynomial component group."""

    name: str
    width: float
    decay: float
    rho_polynomial: np.ndarray
    w_polynomial: np.ndarray


def _sha256(path: Path) -> str:
    """Return the lowercase SHA-256 digest of one file."""
    if not path.is_file():
        raise FileNotFoundError(f'file not found for SHA-256: {path}')
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate_release_inventory(coefficient_path: Path) -> None:
    """Require every released model to contain each element exactly once."""
    expected = set(EXPECTED_ATOMIC_NUMBERS)
    for model in MODEL_NAMES:
        elements = available_elements(model, path=coefficient_path)
        atomic_numbers = tuple(atomic_number for atomic_number, _ in elements)
        if atomic_numbers == EXPECTED_ATOMIC_NUMBERS:
            continue
        actual = set(atomic_numbers)
        duplicates = sorted(
            atomic_number
            for atomic_number in actual
            if atomic_numbers.count(atomic_number) > 1
        )
        raise ValueError(
            f'invalid atomic-number inventory for model={model!r}: '
            f'missing={sorted(expected - actual)}, '
            f'unexpected={sorted(actual - expected)}, duplicates={duplicates}'
        )


def _float32_rational(value: float) -> sp.Rational:
    """Return the exact dyadic rational represented by one float32 value."""
    return sp.Rational(float(np.float32(value)))


def _reciprocal_numerator(
    coefficient_path: Path,
    model: str,
    atomic_number: int,
    quantity: str,
) -> sp.Poly:
    """Return the exact common-denominator numerator in ``x = g**2``."""
    if quantity not in QUANTITIES_RECIPROCAL:
        raise ValueError(f'unknown reciprocal quantity: {quantity!r}')
    x = sp.symbols('x')
    coefficients = load_coefficients(
        atomic_number,
        model=model,
        path=coefficient_path,
    )
    denominators: list[sp.Poly] = []
    numerators: list[sp.Poly] = []
    for amplitude_value, width_value, is_dirac_pade in zip(
        coefficients.a,
        coefficients.b,
        coefficients.is_dirac_pade,
        strict=True,
    ):
        amplitude = _float32_rational(amplitude_value)
        width = _float32_rational(width_value)
        denominator = sp.Poly(1 + width * x, x, domain=sp.QQ)
        if quantity == 'fx':
            power = 3 if is_dirac_pade else 2
            denominators.append(denominator**power)
            numerators.append(sp.Poly(amplitude / width, x, domain=sp.QQ))
        elif is_dirac_pade:
            denominators.append(denominator**3)
            numerators.append(
                sp.Poly(
                    amplitude
                    * (3 + 3 * width * x + width**2 * x**2),
                    x,
                    domain=sp.QQ,
                )
            )
        else:
            denominators.append(denominator**2)
            numerators.append(
                sp.Poly(amplitude * (2 + width * x), x, domain=sp.QQ)
            )

    common_denominator = sp.Poly(1, x, domain=sp.QQ)
    for denominator in denominators:
        common_denominator *= denominator
    common_numerator = sp.Poly(0, x, domain=sp.QQ)
    for numerator, denominator in zip(numerators, denominators, strict=True):
        common_numerator += numerator * common_denominator.exquo(denominator)
    return common_numerator


def _reciprocal_root_task(
    task: tuple[str, str, int, str],
) -> dict[str, int | str | bool]:
    """Isolate every exact real root of one reciprocal-space numerator."""
    coefficient_path_text, model, atomic_number, quantity = task
    numerator = _reciprocal_numerator(
        Path(coefficient_path_text),
        model,
        atomic_number,
        quantity,
    )
    positive_roots = 0
    zero_roots = 0
    for interval, multiplicity in sp.polys.polytools.intervals(numerator):
        left, right = interval
        if left == 0 and right == 0:
            zero_roots += int(multiplicity)
        elif right <= 0:
            continue
        elif left >= 0:
            positive_roots += int(multiplicity)
        else:
            raise RuntimeError(
                f'{model} Z={atomic_number} {quantity} has a root interval '
                f'that straddles x=0: {interval}'
            )
    return {
        'model': model,
        'atomic_number': atomic_number,
        'quantity': quantity,
        'degree': int(numerator.degree()),
        'positive_at_origin': bool(numerator.eval(0) > 0),
        'positive_leading_coefficient': bool(numerator.LC() > 0),
        'positive_root_count': positive_roots,
        'zero_root_count': zero_roots,
    }


def _reciprocal_components(
    coefficient_path: Path,
    model: str,
    atomic_number: int,
    quantity: str,
    x_values: np.ndarray,
) -> np.ndarray:
    """Return reciprocal-space components on an ``x = g**2`` grid."""
    coefficients = load_coefficients(
        atomic_number,
        model=model,
        path=coefficient_path,
    )
    x_column = np.asarray(x_values, dtype=np.float64).reshape(-1, 1)
    amplitude = coefficients.a.reshape(1, -1)
    width = coefficients.b.reshape(1, -1)
    u = x_column * width
    denominator = 1.0 + u
    if quantity == 'fe':
        nr = amplitude * (2.0 + u) / denominator**2
        dp = amplitude * (3.0 + 3.0 * u + u**2) / denominator**3
        return np.where(coefficients.is_dirac_pade.reshape(1, -1), dp, nr)
    if quantity == 'fx':
        scaled_amplitude = amplitude / width
        power = np.where(
            coefficients.is_dirac_pade.reshape(1, -1),
            3.0,
            2.0,
        )
        return scaled_amplitude / denominator**power
    raise ValueError(f'unknown reciprocal quantity: {quantity!r}')


def _reciprocal_tail_ratio(
    coefficient_path: Path,
    model: str,
    atomic_number: int,
    quantity: str,
) -> float:
    """Return the component-normalised positive asymptotic coefficient."""
    coefficients = load_coefficients(
        atomic_number,
        model=model,
        path=coefficient_path,
    )
    if quantity == 'fe':
        terms = coefficients.a / coefficients.b
    elif quantity == 'fx':
        use = ~coefficients.is_dirac_pade
        terms = coefficients.a[use] / coefficients.b[use] ** 3
    else:
        raise ValueError(f'unknown reciprocal quantity: {quantity!r}')
    return float(np.sum(terms) / np.sum(np.abs(terms)))


def _reciprocal_diagnostics(
    coefficient_path: Path,
    model: str,
    quantity: str,
) -> dict[str, float | int | str]:
    """Return the smallest sampled relative component-sum margin."""
    x_values = np.concatenate(
        (
            np.asarray([0.0]),
            np.geomspace(
                RECIPROCAL_X_MIN,
                RECIPROCAL_X_MAX,
                RECIPROCAL_GRID_POINTS,
            ),
        )
    )
    worst = (float('inf'), -1, float('nan'), 'grid')
    for atomic_number, _ in available_elements(model, path=coefficient_path):
        components = _reciprocal_components(
            coefficient_path,
            model,
            atomic_number,
            quantity,
            x_values,
        )
        ratio = np.sum(components, axis=1) / np.sum(np.abs(components), axis=1)
        index = int(np.argmin(ratio))
        candidate = (
            float(ratio[index]),
            atomic_number,
            float(np.sqrt(x_values[index])),
            'grid',
        )
        tail_ratio = _reciprocal_tail_ratio(
            coefficient_path,
            model,
            atomic_number,
            quantity,
        )
        if tail_ratio < candidate[0]:
            candidate = (tail_ratio, atomic_number, float('inf'), 'analytic_tail')
        if candidate < worst:
            worst = candidate
    return {
        'relative_margin': worst[0],
        'atomic_number': worst[1],
        'g_inverse_angstrom': worst[2],
        'location': worst[3],
    }


def _component_groups(
    coefficient_path: Path,
    model: str,
    atomic_number: int,
) -> list[component_group]:
    """Return exact-width NR groups and the independent Dirac--Padé group."""
    coefficients = load_coefficients(
        atomic_number,
        model=model,
        path=coefficient_path,
    )
    amplitude = coefficients.a
    width = coefficients.b
    decay = 2.0 * np.pi / np.sqrt(width)
    n_nr = amplitude.size - 1 if model == 'adaptive_dirac_pade' else amplitude.size
    groups: list[component_group] = []
    for group_width in np.sort(np.unique(width[:n_nr]))[::-1]:
        indices = np.flatnonzero(width[:n_nr] == group_width)
        effective_amplitude = float(np.sum(amplitude[indices], dtype=np.float64))
        if effective_amplitude == 0.0:
            continue
        index = int(indices[0])
        rho_amplitude = (
            2.0 * np.pi**4 * A_0 * effective_amplitude / width[index] ** 2.5
        )
        w_amplitude = effective_amplitude / width[index] ** 1.5
        groups.append(
            component_group(
                name=f'adaptive_hydrogenic:b={group_width:.9g}',
                width=float(group_width),
                decay=float(decay[index]),
                rho_polynomial=np.asarray([rho_amplitude]),
                w_polynomial=np.asarray(
                    [
                        2.0 * w_amplitude / decay[index],
                        w_amplitude,
                    ]
                ),
            )
        )
    if model == 'adaptive_dirac_pade' and amplitude[-1] != 0.0:
        index = amplitude.size - 1
        rho_amplitude = 2.0 * np.pi**4 * A_0 * amplitude[index] / width[index] ** 2.5
        w_amplitude = amplitude[index] / width[index] ** 1.5
        groups.append(
            component_group(
                name=f'adaptive_dirac_pade:b={width[index]:.9g}',
                width=float(width[index]),
                decay=float(decay[index]),
                rho_polynomial=np.asarray(
                    [rho_amplitude / 4.0, rho_amplitude * decay[index] / 4.0]
                ),
                w_polynomial=np.asarray(
                    [
                        2.0 * w_amplitude / decay[index],
                        5.0 * w_amplitude / 4.0,
                        w_amplitude * decay[index] / 4.0,
                    ]
                ),
            )
        )
    if not groups:
        raise ValueError(f'{model} Z={atomic_number} has no non-zero component groups')
    return groups


def _group_polynomial(group: component_group, quantity: str) -> np.ndarray:
    """Return one group's polynomial for ``rho`` or ``w``."""
    if quantity == 'rho':
        return group.rho_polynomial
    if quantity == 'w':
        return group.w_polynomial
    raise ValueError(f'unknown real-space quantity: {quantity!r}')


def _leading_group(
    groups: list[component_group],
    quantity: str,
) -> tuple[component_group, list[component_group]]:
    """Select the slowest group, using polynomial degree for an exact tie."""
    ordered = sorted(
        groups,
        key=lambda group: (
            group.width,
            _group_polynomial(group, quantity).size - 1,
        ),
        reverse=True,
    )
    leading = ordered[0]
    if float(_group_polynomial(leading, quantity)[0]) <= 0.0:
        raise ValueError(
            f'leading {quantity} group is not positive: {leading.name}'
        )
    return leading, ordered[1:]


def _relative_monotonic_onset(
    leading: component_group,
    other: component_group,
    quantity: str,
) -> float:
    """Return a radius after which ``abs(other) / leading`` cannot increase."""
    leading_polynomial = np.abs(_group_polynomial(leading, quantity))
    other_polynomial = np.abs(_group_polynomial(other, quantity))
    delta_decay = other.decay - leading.decay
    if delta_decay < 0.0:
        raise ValueError('a nominally sub-leading group decays more slowly')
    numerator = polynomial.polysub(
        polynomial.polymul(polynomial.polyder(other_polynomial), leading_polynomial),
        polynomial.polymul(polynomial.polyder(leading_polynomial), other_polynomial),
    )
    numerator = polynomial.polysub(
        numerator,
        delta_decay * polynomial.polymul(other_polynomial, leading_polynomial),
    )
    scale = float(np.max(np.abs(numerator)))
    if scale == 0.0:
        return 0.0
    tolerance = 64.0 * np.finfo(np.float64).eps * scale
    numerator[np.abs(numerator) < tolerance] = 0.0
    numerator = polynomial.polytrim(numerator)
    roots = polynomial.polyroots(numerator)
    real_roots = [
        float(root.real)
        for root in roots
        if root.real >= 0.0
        and abs(root.imag) <= 1.0e-7 * max(1.0, abs(root.real))
    ]
    onset = max(real_roots, default=0.0)
    probe = onset + max(1.0e-9, 1.0e-7 * max(1.0, onset))
    if polynomial.polyval(probe, numerator) > tolerance:
        raise ValueError(
            f'relative contribution remains increasing after r={onset:.9g}'
        )
    return onset


def _tail_certificate(
    coefficient_path: Path,
    model: str,
    atomic_number: int,
    quantity: str,
) -> dict[str, float | str]:
    """Return a positive-tail dominance radius and analytic margin."""
    leading, others = _leading_group(
        _component_groups(coefficient_path, model, atomic_number),
        quantity,
    )
    monotonic_onset = max(
        (
            _relative_monotonic_onset(leading, other, quantity)
            for other in others
        ),
        default=0.0,
    )

    def log_ratio_sum(radius: float) -> float:
        leading_value = abs(
            polynomial.polyval(radius, _group_polynomial(leading, quantity))
        )
        log_ratios = [
            np.log(
                abs(polynomial.polyval(radius, _group_polynomial(other, quantity)))
            )
            - np.log(leading_value)
            - (other.decay - leading.decay) * radius
            for other in others
        ]
        return float(logsumexp(log_ratios)) if log_ratios else -np.inf

    target = float(np.log(TAIL_RELATIVE_TARGET))
    if log_ratio_sum(monotonic_onset) <= target:
        radius = monotonic_onset
    else:
        upper = max(1.0, 2.0 * monotonic_onset)
        while log_ratio_sum(upper) > target:
            upper *= 2.0
            if upper > 1.0e8:
                raise RuntimeError(
                    f'failed to bracket the {quantity} dominance radius for '
                    f'{model} Z={atomic_number}'
                )
        radius = brentq(
            lambda value: log_ratio_sum(value) - target,
            monotonic_onset,
            upper,
            xtol=1.0e-13,
            rtol=1.0e-13,
        )
    ratio_sum = float(np.exp(log_ratio_sum(radius)))
    return {
        'radius_angstrom': float(radius),
        'monotonic_onset_angstrom': float(monotonic_onset),
        'leading_group': leading.name,
        'subleading_over_leading': ratio_sum,
        'leading_relative_margin': 1.0 - ratio_sum,
    }


def _real_components(
    coefficients: Any,
    quantity: str,
    radius: float,
) -> np.ndarray:
    """Return signed real-space components at one radius."""
    amplitude = coefficients.a
    width = coefficients.b
    decay = 2.0 * np.pi / np.sqrt(width)
    exponential = np.exp(-decay * radius)
    values: np.ndarray
    if quantity == 'rho':
        scaled_amplitude = 2.0 * np.pi**4 * A_0 * amplitude / width**2.5
        values = scaled_amplitude * exponential
        use = coefficients.is_dirac_pade
        values[use] = (
            scaled_amplitude[use]
            / 4.0
            * (1.0 + decay[use] * radius)
            * exponential[use]
        )
        return values
    if quantity == 'w':
        scaled_amplitude = amplitude / width**1.5
        values = scaled_amplitude * (2.0 / decay + radius) * exponential
        use = coefficients.is_dirac_pade
        values[use] = (
            scaled_amplitude[use]
            / 4.0
            * (
                8.0 / decay[use]
                + 5.0 * radius
                + decay[use] * radius**2
            )
            * exponential[use]
        )
        return values
    raise ValueError(f'unknown real-space quantity: {quantity!r}')


def _finite_interval_certificate(
    coefficient_path: Path,
    model: str,
    atomic_number: int,
    quantity: str,
    radius_end: float,
) -> dict[str, float | int]:
    """Certify positivity by monotone component enclosures on ``[0, R]``."""
    coefficients = load_coefficients(
        atomic_number,
        model=model,
        path=coefficient_path,
    )
    positive = coefficients.a >= 0.0
    stack = [(0.0, radius_end, 0)]
    accepted_intervals = 0
    maximum_depth = 0
    smallest_relative_lower_bound = float('inf')
    while stack:
        left, right, depth = stack.pop()
        left_values = _real_components(
            coefficients,
            quantity,
            left,
        )
        right_values = _real_components(
            coefficients,
            quantity,
            right,
        )
        lower_bound = float(
            np.sum(np.where(positive, right_values, left_values), dtype=np.float64)
        )
        component_scale = float(
            np.sum(
                np.maximum(np.abs(left_values), np.abs(right_values)),
                dtype=np.float64,
            )
        )
        if not np.isfinite(component_scale) or component_scale <= 0.0:
            raise FloatingPointError(
                f'invalid {quantity} component scale for {model} '
                f'Z={atomic_number} on [{left}, {right}]'
            )
        relative_lower_bound = lower_bound / component_scale
        if relative_lower_bound > FINITE_RELATIVE_TOLERANCE:
            accepted_intervals += 1
            maximum_depth = max(maximum_depth, depth)
            smallest_relative_lower_bound = min(
                smallest_relative_lower_bound,
                relative_lower_bound,
            )
            continue
        if depth >= MAX_INTERVAL_DEPTH:
            raise ValueError(
                f'{model} Z={atomic_number} {quantity} cannot be certified on '
                f'[{left:.17g}, {right:.17g}]: relative lower bound '
                f'{relative_lower_bound:.6e}'
            )
        midpoint = 0.5 * (left + right)
        stack.append((midpoint, right, depth + 1))
        stack.append((left, midpoint, depth + 1))
    return {
        'radius_end_angstrom': radius_end,
        'accepted_intervals': accepted_intervals,
        'maximum_depth': maximum_depth,
        'smallest_relative_lower_bound': smallest_relative_lower_bound,
        'required_relative_tolerance': FINITE_RELATIVE_TOLERANCE,
    }


def _reciprocal_certificate(
    coefficient_path: Path,
    workers: int,
) -> dict[str, Any]:
    """Return exact root-isolation and sampled-margin results."""
    tasks = [
        (str(coefficient_path), model, atomic_number, quantity)
        for model in MODEL_NAMES
        for atomic_number, _ in available_elements(model, path=coefficient_path)
        for quantity in QUANTITIES_RECIPROCAL
    ]
    if workers == 1:
        root_results = [_reciprocal_root_task(task) for task in tasks]
    else:
        with mp.Pool(processes=workers) as pool:
            root_results = list(
                pool.imap_unordered(_reciprocal_root_task, tasks, chunksize=1)
            )
    root_results.sort(
        key=lambda item: (
            str(item['model']),
            int(item['atomic_number']),
            str(item['quantity']),
        )
    )
    failures = [
        item
        for item in root_results
        if not item['positive_at_origin']
        or not item['positive_leading_coefficient']
        or item['positive_root_count']
        or item['zero_root_count']
    ]
    if failures:
        raise ValueError(f'reciprocal-space positivity failures: {failures}')

    result: dict[str, Any] = {}
    for model in MODEL_NAMES:
        model_roots = [item for item in root_results if item['model'] == model]
        result[model] = {
            'certified_polynomials': len(model_roots),
            'positive_root_count': sum(
                int(item['positive_root_count']) for item in model_roots
            ),
            'zero_root_count': sum(int(item['zero_root_count']) for item in model_roots),
            'maximum_numerator_degree': max(
                int(item['degree']) for item in model_roots
            ),
            'minimum_sampled_relative_margin': {
                quantity: _reciprocal_diagnostics(
                    coefficient_path,
                    model,
                    quantity,
                )
                for quantity in QUANTITIES_RECIPROCAL
            },
        }
    return result


def _real_space_certificate(coefficient_path: Path) -> dict[str, Any]:
    """Return finite-interval and positive-tail certificates."""
    result: dict[str, Any] = {}
    for model in MODEL_NAMES:
        per_quantity: dict[str, Any] = {}
        for quantity in QUANTITIES_REAL:
            records = []
            for atomic_number, _ in available_elements(model, path=coefficient_path):
                tail = _tail_certificate(
                    coefficient_path,
                    model,
                    atomic_number,
                    quantity,
                )
                finite = _finite_interval_certificate(
                    coefficient_path,
                    model,
                    atomic_number,
                    quantity,
                    max(
                        TAIL_FINITE_OVERLAP_FACTOR
                        * float(tail['radius_angstrom']),
                        float(tail['radius_angstrom']) + 1.0e-9,
                    ),
                )
                records.append(
                    {
                        'atomic_number': atomic_number,
                        'tail': tail,
                        'finite': finite,
                    }
                )
            worst_finite = min(
                records,
                key=lambda item: float(
                    item['finite']['smallest_relative_lower_bound']
                ),
            )
            largest_tail_radius = max(
                records,
                key=lambda item: float(item['tail']['radius_angstrom']),
            )
            per_quantity[quantity] = {
                'certified_records': len(records),
                'worst_finite_enclosure': {
                    'atomic_number': worst_finite['atomic_number'],
                    **worst_finite['finite'],
                },
                'largest_tail_radius': {
                    'atomic_number': largest_tail_radius['atomic_number'],
                    **largest_tail_radius['tail'],
                },
                'total_accepted_intervals': sum(
                    int(item['finite']['accepted_intervals']) for item in records
                ),
            }
        result[model] = per_quantity
    return result


def build_certificate(coefficient_path: Path, workers: int) -> dict[str, Any]:
    """Build the complete deterministic positivity certificate."""
    coefficient_path = coefficient_path.expanduser().resolve()
    if not coefficient_path.is_file():
        raise FileNotFoundError(f'coefficient archive not found: {coefficient_path}')
    if workers < 1:
        raise ValueError(f'workers must be positive, got {workers}')
    _validate_release_inventory(coefficient_path)
    return {
        'schema_version': CERTIFICATE_SCHEMA,
        'coefficient_file': coefficient_path.name,
        'coefficient_sha256': _sha256(coefficient_path),
        'delivery_precision': 'float32',
        'scope': {
            'models': list(MODEL_NAMES),
            'elements_per_model': len(EXPECTED_ATOMIC_NUMBERS),
            'reciprocal_domains': {
                'fe': 'finite g >= 0; limit 0 from above as g -> infinity',
                'fx': 'finite g >= 0; limit 0 from above as g -> infinity',
            },
            'real_domains': {
                'rho': 'finite r >= 0; limit 0 from above as r -> infinity',
                'w': 'W=rV on finite r >= 0, equivalent to V>0 for r>0; '
                'limit 0 from above as r -> infinity',
            },
        },
        'method': {
            'reciprocal': 'exact rational common-denominator polynomial root isolation',
            'real_finite': 'adaptive monotone-component interval lower bounds',
            'real_tail': 'positive leading exponential-polynomial dominance',
            'tail_subleading_over_leading_target': TAIL_RELATIVE_TARGET,
            'tail_finite_overlap_factor': TAIL_FINITE_OVERLAP_FACTOR,
            'finite_relative_tolerance': FINITE_RELATIVE_TOLERANCE,
            'reciprocal_diagnostic_grid': {
                'points': RECIPROCAL_GRID_POINTS + 1,
                'x_equals_g_squared_min': RECIPROCAL_X_MIN,
                'x_equals_g_squared_max': RECIPROCAL_X_MAX,
            },
        },
        'reciprocal': _reciprocal_certificate(coefficient_path, workers),
        'real_space': _real_space_certificate(coefficient_path),
        'projected_potential_corollary': (
            'V(R)>0 for R>0 because it is the line integral of V(r)>0.'
        ),
    }


def _compare_certificate(actual: Any, expected: Any, path: str = 'root') -> None:
    """Compare certificates with tight tolerance for platform-dependent roots."""
    if isinstance(actual, dict) and isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError(
                f'certificate keys differ at {path}: '
                f'actual={sorted(actual)}, expected={sorted(expected)}'
            )
        for key in actual:
            _compare_certificate(actual[key], expected[key], f'{path}.{key}')
        return
    if isinstance(actual, list) and isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(
                f'certificate list lengths differ at {path}: '
                f'{len(actual)} != {len(expected)}'
            )
        for index, (actual_item, expected_item) in enumerate(
            zip(actual, expected, strict=True)
        ):
            _compare_certificate(actual_item, expected_item, f'{path}[{index}]')
        return
    if isinstance(actual, float) or isinstance(expected, float):
        if np.isinf(actual) and np.isinf(expected):
            return
        if not np.isclose(actual, expected, rtol=2.0e-10, atol=1.0e-14):
            raise ValueError(
                f'certificate float differs at {path}: {actual} != {expected}'
            )
        return
    if actual != expected:
        raise ValueError(f'certificate value differs at {path}: {actual!r} != {expected!r}')


def _print_summary(certificate: dict[str, Any]) -> None:
    """Print a concise certificate summary."""
    print('global positivity certification passed')
    elements_per_model = int(certificate['scope']['elements_per_model'])
    for model in MODEL_NAMES:
        reciprocal = certificate['reciprocal'][model]
        real_space = certificate['real_space'][model]
        print(
            f'  {model}: reciprocal '
            f'{reciprocal["certified_polynomials"]}/'
            f'{len(QUANTITIES_RECIPROCAL) * elements_per_model}; '
            f'rho {real_space["rho"]["certified_records"]}/{elements_per_model}; '
            f'V {real_space["w"]["certified_records"]}/{elements_per_model}'
        )
        for quantity in QUANTITIES_RECIPROCAL:
            margin = reciprocal['minimum_sampled_relative_margin'][quantity]
            print(
                f'    worst {quantity} relative margin: '
                f'{margin["relative_margin"]:.6e} at Z='
                f'{margin["atomic_number"]}'
            )


def _write_text_atomic(path: Path, text: str) -> None:
    """Write UTF-8 text through an exclusive sibling temporary file."""
    temporary_path = path.with_name(f'.{path.name}.tmp')
    if temporary_path.exists():
        raise FileExistsError(f'temporary output already exists: {temporary_path}')
    try:
        with temporary_path.open('x', encoding='utf-8') as handle:
            handle.write(text)
        temporary_path.replace(path)
    finally:
        temporary_path.unlink(missing_ok=True)


def _require_external_certificate_path(path: Path) -> None:
    """Reject generated certificates inside the release checkout or package."""
    package_root = Path(__file__).resolve().parent
    checkout_root = package_root.parent
    release_root = (
        checkout_root
        if (checkout_root / 'pyproject.toml').is_file()
        else package_root
    )
    if path.resolve().is_relative_to(release_root):
        raise ValueError(
            f'--write-certificate path must be outside the release: {path}'
        )


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        '--coefficients',
        type=Path,
        default=DEFAULT_COEFFICIENTS,
        help='released coefficient archive',
    )
    parser.add_argument(
        '--certificate',
        type=Path,
        default=DEFAULT_CERTIFICATE,
        help='certificate to verify or replace',
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=min(6, os.cpu_count() or 1),
        help='worker processes for exact reciprocal root isolation',
    )
    parser.add_argument(
        '--write-certificate',
        action='store_true',
        help=f'write {CERTIFICATE_NAME} instead of checking it',
    )
    args = parser.parse_args()
    coefficient_path = args.coefficients.expanduser().resolve()
    certificate_path = args.certificate.expanduser().resolve()
    if args.write_certificate:
        _require_external_certificate_path(certificate_path)
    actual = build_certificate(coefficient_path, args.workers)
    if args.write_certificate:
        _write_text_atomic(
            certificate_path,
            json.dumps(actual, indent=2, sort_keys=True, allow_nan=False) + '\n',
        )
        print(f'wrote {certificate_path}')
    else:
        if not certificate_path.is_file():
            raise FileNotFoundError(
                f'positivity certificate not found: {certificate_path}'
            )
        expected = json.loads(certificate_path.read_text(encoding='utf-8'))
        _compare_certificate(actual, expected)
        print(f'certificate matches {certificate_path}')
    _print_summary(actual)


if __name__ == '__main__':
    main()
