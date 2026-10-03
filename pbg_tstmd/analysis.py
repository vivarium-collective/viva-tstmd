"""Quantitative ECM-remodeling metrics for TST-MD runs.

These reproduce, in 2D, the measures defined in Keijzer & Merks 2026
(arXiv:2609.02375) for the 3D hybrid CPM + bead-spring model:

- ``giant_component_fraction`` — network percolation (their Fig 2B): the
  fraction of ECM beads in the largest connected component of the
  bond+crosslink graph.
- ``reorientation_q`` / ``q_by_distance`` — fiber reorientation toward the
  cell (their Fig 4A-C, eq. 7): the fraction of fibers nearly radial to the
  cell. In 2D the baseline distribution differs from 3D, so the *ratio*
  q/q0 (final vs initial) is the directional readout, not the absolute q.
- ``densification_factor`` — ECM accumulation near the cell (their Fig 4D-E,
  eqs. 8-9): ratio of bead density close to the cell vs near the domain
  boundary.

All functions are pure and operate on plain numpy arrays extracted from the
engine's pickle state dumps; none require the simulator or Docker.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

__all__ = [
    "giant_component_fraction",
    "reorientation_q",
    "q_by_distance",
    "ratio_with_nan_guard",
    "densification_factor",
]


def giant_component_fraction(n_nodes: int, edges: Iterable[Sequence[int]]) -> float:
    """Largest-connected-component size / n_nodes, via union-find.

    Nodes are ``0..n_nodes-1``; ``edges`` is an iterable of ``(i, j)`` pairs
    (ECM bonds and crosslinks). With no edges every node is its own
    component, so the fraction is ``1 / n_nodes`` (never a crash).
    """
    if n_nodes <= 0:
        return float("nan")
    parent = list(range(n_nodes))

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:  # path compression
            parent[x], x = root, parent[x]
        return root

    for i, j in edges:
        i, j = int(i), int(j)
        if i < 0 or j < 0 or i >= n_nodes or j >= n_nodes:
            continue
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    sizes = np.zeros(n_nodes, dtype=np.int64)
    for x in range(n_nodes):
        sizes[find(x)] += 1
    return float(sizes.max()) / float(n_nodes)


def _alignment(fiber_centers, fiber_dirs, cell_centroid):
    """gamma_i = |c_i . v_i| / (||c_i|| ||v_i||) in [0, 1]; NaN where degenerate."""
    c = np.asarray(fiber_centers, dtype=float) - np.asarray(cell_centroid, dtype=float)
    v = np.asarray(fiber_dirs, dtype=float)
    cn = np.linalg.norm(c, axis=1)
    vn = np.linalg.norm(v, axis=1)
    denom = cn * vn
    dots = np.abs(np.sum(c * v, axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        gamma = np.where(denom > 0, dots / denom, np.nan)
    return gamma


def reorientation_q(fiber_centers, fiber_dirs, cell_centroid, lo: float = 0.9) -> float:
    """Fraction of fibers with alignment gamma >= ``lo`` (their q; eq. 7).

    gamma = |cos(angle between the cell->fiber vector and the fiber
    direction)|; gamma near 1 means the fiber points toward/away from the
    cell. Degenerate fibers (zero-length direction or a fiber centered on the
    cell centroid) are excluded from the fraction.
    """
    gamma = _alignment(fiber_centers, fiber_dirs, cell_centroid)
    valid = ~np.isnan(gamma)
    if not valid.any():
        return float("nan")
    return float(np.mean(gamma[valid] >= lo))


def q_by_distance(fiber_centers, fiber_dirs, cell_centroid, bins, lo: float = 0.9):
    """q computed per distance-from-cell bin; NaN for bins with no fibers.

    ``bins`` are the bin edges; returns an array of length ``len(bins)-1``.
    The q(r)/q0(r) ratio of Fig 4C is ``ratio_with_nan_guard(q_final, q0)``.
    """
    c = np.asarray(fiber_centers, dtype=float) - np.asarray(cell_centroid, dtype=float)
    dist = np.linalg.norm(c, axis=1)
    gamma = _alignment(fiber_centers, fiber_dirs, cell_centroid)
    edges = np.asarray(bins, dtype=float)
    out = np.full(len(edges) - 1, np.nan)
    for b in range(len(edges) - 1):
        in_bin = (dist >= edges[b]) & (dist < edges[b + 1]) & ~np.isnan(gamma)
        if in_bin.any():
            out[b] = float(np.mean(gamma[in_bin] >= lo))
    return out


def ratio_with_nan_guard(numer, denom):
    """Elementwise numer/denom; NaN wherever denom is 0, NaN, or missing."""
    a = np.asarray(numer, dtype=float)
    b = np.asarray(denom, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((b != 0) & ~np.isnan(b) & ~np.isnan(a), a / b, np.nan)


def densification_factor(
    bead_coords,
    cell_centroid,
    domain_size: float,
    close_px: float = 30.0,
    far_px: float = 50.0,
) -> float:
    """rho_close / rho_far (their rho_dens; eqs. 8-9), 2D areas for volumes.

    ``rho_close`` is the bead density within ``close_px`` of the cell
    centroid; ``rho_far`` the density within ``far_px`` of the square
    domain's boundary. Returns NaN if the far region is empty (no beads or
    zero area), never +/-inf.
    """
    coords = np.asarray(bead_coords, dtype=float)
    if coords.size == 0:
        return float("nan")
    centroid = np.asarray(cell_centroid, dtype=float)
    L = float(domain_size)

    dist_cell = np.linalg.norm(coords - centroid, axis=1)
    n_close = int(np.count_nonzero(dist_cell <= close_px))
    area_close = float(np.pi * close_px ** 2)

    x, y = coords[:, 0], coords[:, 1]
    dist_boundary = np.minimum.reduce([x, L - x, y, L - y])
    n_far = int(np.count_nonzero(dist_boundary <= far_px))
    inner = max(0.0, L - 2.0 * far_px)
    area_far = float(L * L - inner * inner)

    if n_far == 0 or area_far <= 0 or area_close <= 0:
        return float("nan")
    rho_close = n_close / area_close
    rho_far = n_far / area_far
    if rho_far == 0:
        return float("nan")
    return float(rho_close / rho_far)
