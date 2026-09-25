#!/usr/bin/env python3
# Copyright 2026 Ivan Lobato / NeuralSoftX
# SPDX-License-Identifier: Apache-2.0
# Email: ivan.lobato@neuralsoftx.com
"""Per-atom diagnostic PDF from the released data.

Author: Ivan Lobato
"""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages

from .scattering_factors import (
    Coefficients,
    DEFAULT_COEFFICIENTS,
    DEFAULT_REFERENCE,
    available_elements,
    available_models,
    electron_density as rho,
    electron_scattering_factor as feg,
    electrostatic_potential as vr,
    load_coefficients,
    xray_scattering_factor as fxg,
)


def _release_root() -> Path:
    """Return the protected release root for a checkout or installation."""
    package_root = Path(__file__).resolve().parent
    checkout_root = package_root.parent
    if (checkout_root / 'pyproject.toml').is_file():
        return checkout_root
    return package_root

def _inset_hi_g(ax, g, ys, colours, lws):
    """High-g (g>=20) inset to show the low-amplitude tail."""
    m = g >= 20.0
    if not m.any():
        return
    axi = ax.inset_axes([0.45, 0.22, 0.50, 0.32])
    for y, c, lw in zip(ys, colours, lws):
        axi.plot(g[m], y[m], c, lw=lw)
    axi.tick_params(labelsize=6)
    axi.grid(True, alpha=0.25)
    axi.set_title('g >= 20', fontsize=7)


def _tidy_ticks(fig):
    axes = []
    for a in fig.axes:
        axes.append(a)
        axes.extend(getattr(a, 'child_axes', []))
    for a in axes:
        fmt = mticker.ScalarFormatter()
        fmt.set_powerlimits((0, 0))
        a.yaxis.set_major_formatter(fmt)
        a.yaxis.get_offset_text().set_fontsize(6)


def draw_atom(fig, coefficients: Coefficients, ref):
    """Draw the 8-panel page for one atom into ``fig``."""
    fig.clf()
    z = coefficients.atomic_number
    sym = coefficients.symbol
    model = coefficients.model
    reference_label = 'Dirac 1s' if z == 1 else 'Dirac–Fock'
    g = ref['g']
    r = ref['r']
    feg_ref, fxg_ref, pr_ref = ref['feg'], ref['fxg'], ref['pr']
    rmax = float(r[-1])
    r_log = np.logspace(-5.0, np.log10(rmax), 600)
    r_v_small = np.linspace(1e-2, 2e-1, 400)
    r_v_full = np.logspace(-2.0, np.log10(rmax), 800)

    feg_fit = feg(coefficients, g)
    fxg_fit = fxg(coefficients, g)
    pr_fit = rho(coefficients, r)
    pr_fit_log = rho(coefficients, r_log)
    vr_small = vr(coefficients, r_v_small)
    vr_full = vr(coefficients, r_v_full)

    # (1) f_e(g)
    ax = fig.add_subplot(2, 4, 1)
    ax.plot(g, feg_ref, '-r', lw=1.2, label=reference_label)
    ax.plot(g, feg_fit, '-b', lw=1.2, label=model)
    ax.set_xlabel('g (1/A)'); ax.set_ylabel('f_e(g) (A)'); ax.set_title('f_e(g)')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, alpha=0.25)
    _inset_hi_g(ax, g, [feg_ref, feg_fit], ['-r', '-b'], [1.0, 1.0])

    # (2) f_x(g)
    ax = fig.add_subplot(2, 4, 2)
    ax.plot(g, fxg_ref, '-r', lw=1.2, label=reference_label)
    ax.plot(g, fxg_fit, '-b', lw=1.2, label=model)
    ax.set_xlabel('g (1/A)'); ax.set_ylabel('f_x(g) (e)'); ax.set_title('f_x(g)')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, alpha=0.25)
    _inset_hi_g(ax, g, [fxg_ref, fxg_fit], ['-r', '-b'], [1.0, 1.0])

    # (3) rho(r) small-r + full-r (log) inset
    ax = fig.add_subplot(2, 4, 3)
    ax.plot(r, pr_ref, '-r', lw=1.2, label=reference_label)
    ax.plot(r, pr_fit, '-b', lw=1.2, label=model)
    ax.set_xlim(0.0, 0.03); ax.set_xlabel('r (A)')
    ax.set_ylabel('rho(r) (e/A^3)'); ax.set_title('rho(r) - small r')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, alpha=0.25)
    axi = ax.inset_axes([0.45, 0.25, 0.50, 0.45])
    axi.plot(r, pr_ref, '-r', lw=1.0); axi.plot(r, pr_fit, '-b', lw=1.0)
    axi.set_xscale('log'); axi.set_xlim(1e-4, rmax)
    axi.tick_params(labelsize=6); axi.grid(True, alpha=0.25)
    axi.set_title('full r (log)', fontsize=7)

    # (4) V(r) zoom + full-r (log) inset (model only)
    ax = fig.add_subplot(2, 4, 4)
    ax.plot(r_v_small, vr_small, '-b', lw=1.2, label=model)
    ax.set_xlabel('r (A)'); ax.set_ylabel('V(r)')
    ax.set_title('V(r) zoom [1e-2, 2e-1]')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, alpha=0.25)
    axi = ax.inset_axes([0.45, 0.25, 0.50, 0.45])
    axi.plot(r_v_full, vr_full, '-b', lw=1.0)
    axi.set_xscale('log'); axi.set_xlim(1e-2, rmax)
    axi.tick_params(labelsize=6); axi.grid(True, alpha=0.25)
    axi.set_title('full r (log)', fontsize=7)

    # (5) g^2 f_e(g)
    g2 = g ** 2
    ax = fig.add_subplot(2, 4, 5)
    ax.plot(g, g2 * feg_ref, '-r', lw=1.2, label=reference_label)
    ax.plot(g, g2 * feg_fit, '-b', lw=1.2, label=model)
    ax.set_xlabel('g (1/A)'); ax.set_ylabel('g^2 f_e(g)'); ax.set_title('g^2 f_e(g)')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, alpha=0.25)
    _inset_hi_g(ax, g, [g2 * feg_ref, g2 * feg_fit], ['-r', '-b'], [1.0, 1.0])

    # (6) g^2 f_x(g)
    ax = fig.add_subplot(2, 4, 6)
    ax.plot(g, g2 * fxg_ref, '-r', lw=1.2, label=reference_label)
    ax.plot(g, g2 * fxg_fit, '-b', lw=1.2, label=model)
    ax.set_xlabel('g (1/A)'); ax.set_ylabel('g^2 f_x(g)'); ax.set_title('g^2 f_x(g)')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, alpha=0.25)
    _inset_hi_g(ax, g, [g2 * fxg_ref, g2 * fxg_fit], ['-r', '-b'], [1.0, 1.0])

    # (7) r^2 rho(r), log-r
    ax = fig.add_subplot(2, 4, 7)
    ax.plot(r, r ** 2 * pr_ref, '-r', lw=1.2, label=reference_label)
    ax.plot(r_log, r_log ** 2 * pr_fit_log, '-b', lw=1.2, label=model)
    ax.set_xscale('log'); ax.set_xlim(1e-5, rmax)
    ax.set_xlabel('r (A)'); ax.set_ylabel('r^2 rho(r)'); ax.set_title('r^2 rho(r)')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, which='both', alpha=0.25)

    # (8) r^2 V(r), log-r (model only)
    ax = fig.add_subplot(2, 4, 8)
    ax.plot(r_v_full, r_v_full ** 2 * vr_full, '-b', lw=1.2, label=model)
    ax.set_xscale('log'); ax.set_xlim(1e-2, rmax)
    ax.set_xlabel('r (A)'); ax.set_ylabel('r^2 V(r)'); ax.set_title('r^2 V(r)')
    ax.legend(loc='upper right', fontsize='small'); ax.grid(True, which='both', alpha=0.25)

    _tidy_ticks(fig)
    fig.suptitle(
        f'Z = {z}  {sym}   |   model = {model}   |   '
        f'n_t = {coefficients.n_terms}',
        fontsize=11,
        fontweight='bold',
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))


def _resolve_paths(
    coefficient_value: str,
    reference_value: str,
    output_value: str,
) -> tuple[Path, Path, Path]:
    """Resolve and validate the CLI input and output paths."""
    coefficient_path = Path(coefficient_value).expanduser().resolve()
    reference_path = Path(reference_value).expanduser().resolve()
    output_path = Path(output_value).expanduser().resolve()
    if not reference_path.is_file():
        raise FileNotFoundError(f'reference archive not found: {reference_path}')
    if output_path.suffix.lower() != '.pdf':
        raise ValueError(f'output path must end in .pdf: {output_path}')
    if output_path.is_relative_to(_release_root()):
        raise ValueError(f'output path must be outside the release: {output_path}')
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f'output directory not found: {output_path.parent}')
    if output_path in (coefficient_path, reference_path):
        raise ValueError(
            f'output path must differ from both input archives: {output_path}'
        )
    if output_path.exists() and not output_path.is_file():
        raise ValueError(f'output path is not a regular file: {output_path}')
    return coefficient_path, reference_path, output_path


def _select_elements(
    elements: tuple[tuple[int, str], ...],
    requested: list[int] | None,
) -> tuple[tuple[int, str], ...]:
    """Return the requested elements after exact availability validation."""
    if not elements:
        raise ValueError('the selected model contains no elements')
    if requested is None:
        return elements
    if len(set(requested)) != len(requested):
        raise ValueError(f'--only contains duplicate atomic numbers: {requested}')
    available = {atomic_number for atomic_number, _ in elements}
    unavailable = sorted(set(requested) - available)
    if unavailable:
        raise ValueError(f'--only contains unavailable atomic numbers: {unavailable}')
    selected = set(requested)
    return tuple(item for item in elements if item[0] in selected)


def _write_pdf(
    coefficient_path: Path,
    reference_path: Path,
    output_path: Path,
    model: str,
    elements: tuple[tuple[int, str], ...],
) -> None:
    """Write the diagnostic PDF atomically."""
    temporary_path = output_path.with_name(f'.{output_path.name}.tmp')
    if temporary_path.exists():
        raise FileExistsError(f'temporary output already exists: {temporary_path}')
    try:
        with h5py.File(reference_path, 'r') as reference_handle:
            figure = plt.figure(figsize=(14, 7))
            try:
                with PdfPages(temporary_path) as pdf:
                    for atomic_number, symbol in elements:
                        coefficients = load_coefficients(
                            atomic_number,
                            model,
                            coefficient_path,
                        )
                        reference_group = reference_handle[
                            f'Z{atomic_number:03d}_{symbol}'
                        ]
                        reference = {
                            key: np.asarray(reference_group[key][...], float)
                            for key in ('g', 'feg', 'fxg', 'r', 'pr')
                        }
                        draw_atom(figure, coefficients, reference)
                        pdf.savefig(figure)
                        print(f'  Z={atomic_number:3d} {symbol}')
            finally:
                plt.close(figure)
        temporary_path.replace(output_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Per-atom diagnostic PDF from the released data.'
    )
    ap.add_argument('--coef', default=str(DEFAULT_COEFFICIENTS))
    ap.add_argument(
        '--reference',
        default=str(DEFAULT_REFERENCE),
    )
    ap.add_argument('--model', default='adaptive_dirac_pade',
                    help='fitted model: adaptive_dirac_pade, adaptive_hydrogenic or '
                         'five_term_hydrogenic')
    ap.add_argument('--only', type=int, nargs='+', default=None,
                    help='atomic numbers to render (default: all in the model)')
    ap.add_argument('--out', required=True, help='output PDF path')
    args = ap.parse_args()

    coefficient_path, reference_path, output_path = _resolve_paths(
        args.coef,
        args.reference,
        args.out,
    )
    models = available_models(coefficient_path)
    if args.model not in models:
        raise SystemExit(
            f"model {args.model!r} not in {coefficient_path}; "
            f"choose from {models}"
        )
    elements = available_elements(args.model, coefficient_path)
    elements = _select_elements(elements, args.only)
    _write_pdf(
        coefficient_path,
        reference_path,
        output_path,
        args.model,
        elements,
    )
    print(f'wrote {len(elements)} page(s) -> {output_path}')


if __name__ == '__main__':
    main()
