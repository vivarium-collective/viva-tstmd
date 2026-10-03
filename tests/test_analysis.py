"""Unit tests for the ECM-remodeling metrics on synthetic networks.

These are deterministic and require no simulator or Docker.
"""
import numpy as np
import pytest

from pbg_tstmd.analysis import (
    densification_factor,
    giant_component_fraction,
    q_by_distance,
    ratio_with_nan_guard,
    reorientation_q,
)


# --- giant component / percolation (Fig 2B) ----------------------------------

def test_giant_component_two_triangles_joined():
    # nodes 0-2 triangle, 3-5 triangle, joined by edge (2,3) -> all connected
    edges = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (2, 3)]
    assert giant_component_fraction(6, edges) == pytest.approx(1.0)


def test_giant_component_fully_disconnected():
    assert giant_component_fraction(5, []) == pytest.approx(1 / 5)


def test_giant_component_zero_edges_no_crash():
    # Review Focus: degenerate (zero-crosslink) network -> 1/n, not a crash.
    assert giant_component_fraction(10, iter([])) == pytest.approx(1 / 10)


def test_giant_component_two_components():
    # one 4-clique-ish chain + an isolated pair -> largest = 4 of 6
    edges = [(0, 1), (1, 2), (2, 3), (4, 5)]
    assert giant_component_fraction(6, edges) == pytest.approx(4 / 6)


# --- reorientation q / q0 (Fig 4A-C, eq. 7) ----------------------------------

def _radial_fibers(centroid, n=36, radius=20.0):
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    dirs = np.c_[np.cos(angles), np.sin(angles)]
    centers = centroid + radius * dirs
    return centers, dirs  # each fiber points radially from the centroid


def test_reorientation_all_radial_q_is_one():
    centroid = np.array([100.0, 100.0])
    centers, dirs = _radial_fibers(centroid)
    assert reorientation_q(centers, dirs, centroid, lo=0.9) == pytest.approx(1.0)


def test_reorientation_random_lower_than_radial():
    rng = np.random.default_rng(0)
    centroid = np.array([100.0, 100.0])
    centers = centroid + rng.uniform(-50, 50, size=(2000, 2))
    thetas = rng.uniform(0, 2 * np.pi, size=2000)
    dirs = np.c_[np.cos(thetas), np.sin(thetas)]
    q = reorientation_q(centers, dirs, centroid, lo=0.9)
    assert 0.0 < q < 0.5  # far from the all-radial value of 1.0


def test_q_by_distance_empty_bin_is_nan():
    # Review Focus: a distance bin with no fibers -> NaN, not a spurious 0.
    centroid = np.array([100.0, 100.0])
    centers, dirs = _radial_fibers(centroid, radius=20.0)
    # all fibers sit at r=20; bins [0,10) and [40,50) are empty
    bins = [0, 10, 30, 40, 50]
    q = q_by_distance(centers, dirs, centroid, bins, lo=0.9)
    assert np.isnan(q[0])          # [0,10) empty
    assert q[1] == pytest.approx(1.0)  # [10,30) holds the radial fibers
    assert np.isnan(q[3])          # [40,50) empty


def test_ratio_with_nan_guard():
    out = ratio_with_nan_guard([2.0, 1.0, np.nan], [1.0, 0.0, 1.0])
    assert out[0] == pytest.approx(2.0)
    assert np.isnan(out[1])  # divide-by-zero guarded
    assert np.isnan(out[2])  # NaN numerator propagates


# --- densification factor (Fig 4D-E, eqs. 8-9) -------------------------------

def test_densification_packed_near_cell_above_one():
    rng = np.random.default_rng(1)
    centroid = np.array([100.0, 100.0])
    near = centroid + rng.uniform(-10, 10, size=(500, 2))   # dense near cell
    far = rng.uniform(0, 200, size=(50, 2))                 # sparse elsewhere
    coords = np.vstack([near, far])
    rho = densification_factor(coords, centroid, domain_size=200.0)
    assert rho > 1.0


def test_densification_uniform_near_one():
    rng = np.random.default_rng(2)
    coords = rng.uniform(0, 200, size=(20000, 2))
    rho = densification_factor(coords, [100.0, 100.0], domain_size=200.0)
    assert rho == pytest.approx(1.0, abs=0.25)


def test_densification_empty_far_is_nan():
    # Review Focus: no beads in the far (near-boundary) region -> NaN, not inf.
    centroid = np.array([100.0, 100.0])
    coords = centroid + np.random.default_rng(3).uniform(-5, 5, size=(100, 2))
    rho = densification_factor(coords, centroid, domain_size=200.0,
                               close_px=30.0, far_px=10.0)
    assert np.isnan(rho)
