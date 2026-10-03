"""Process-bigraph surface + sweep driver over the real TST-MD engine.

- ``run_once`` runs one coupled simulation and returns the per-dump
  trajectory of the paper's observables (cell area, percolation, fiber
  reorientation q/q0, densification).
- ``run_sweep`` sweeps one parameter and is resilient: a failed point is
  recorded and the sweep continues.
- ``TstmdEcmProcess`` exposes a single equilibrium run to process-bigraph /
  the dashboard.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from process_bigraph import Process

from .analysis import (
    densification_factor,
    giant_component_fraction,
    reorientation_q,
)
from .runtime import (
    TstmdSession,
    cell_area,
    cell_centroid,
    fiber_beads,
)


def metrics_for_snapshot(snap: dict, n_beads: int) -> dict:
    """The paper's observables for one state dump (2D analogs)."""
    L = float(snap["cpm"].shape[0])  # domain/lattice size in pixels (e.g. 200)
    centroid = cell_centroid(snap)
    centers, dirs = fiber_beads(snap, n_beads)
    q = reorientation_q(centers, dirs, centroid) if len(centers) else float("nan")
    n_nodes = int(snap["positions"].shape[0])
    gc = giant_component_fraction(n_nodes, snap["bond_groups"])
    # paper: 30 um close to cell, 50 um from boundary on a 200 um domain.
    dens = densification_factor(
        snap["positions"], centroid, L, close_px=30.0, far_px=40.0,
    )
    return {
        "mcs": int(snap["mcs"]),
        "cell_area": cell_area(snap),
        "reorientation_q": float(q),
        "giant_component": float(gc),
        "densification": float(dens),
    }


def run_once(settings: Optional[dict] = None, n_beads: int = 9, **session_kw) -> dict:
    """Run one coupled sim; return trajectory + derived q/q0 and equilibrium.

    Returns ``{settings, trajectory:[metrics...], q0, equilibrium:{...},
    crosslinks_created}``. ``q0`` is the first snapshot's reorientation; each
    trajectory point also carries ``q_ratio = q / q0``.
    """
    settings = dict(settings or {})
    with TstmdSession(settings=settings, **session_kw) as sess:
        snaps = sess.run()
    traj = [metrics_for_snapshot(s, n_beads) for s in snaps]
    q0 = traj[0]["reorientation_q"] if traj else float("nan")
    for m in traj:
        m["q_ratio"] = (m["reorientation_q"] / q0) if q0 else float("nan")
    # crosslinks actually created = bonds that are not plain fibers (type != 0)
    crosslinks = int((snaps[-1]["bond_types"] != 0).sum()) if snaps else 0
    eq = traj[-1] if traj else {}
    return {
        "settings": settings,
        "trajectory": traj,
        "q0": q0,
        "equilibrium": eq,
        "crosslinks_created": crosslinks,
    }


def run_sweep(param: str, values, base_settings: Optional[dict] = None,
              n_beads: int = 9, seeds=(12345678,), **session_kw) -> list[dict]:
    """Sweep ``param`` over ``values`` (x seeds); resilient to per-point failures.

    Each result is ``run_once(...)`` plus ``{param, value, seed}``, or
    ``{param, value, seed, error}`` if that point raised.
    """
    base = dict(base_settings or {})
    results = []
    for v in values:
        for seed in seeds:
            settings = {**base, param: v, "md_seed": seed}
            try:
                out = run_once(settings, n_beads=n_beads, **session_kw)
                out.update({"param": param, "value": v, "seed": seed})
            except Exception as exc:  # noqa: BLE001 - record and continue
                out = {"param": param, "value": v, "seed": seed, "error": repr(exc)}
            results.append(out)
    return results


class TstmdEcmProcess(Process):
    """One equilibrium run of the real TST-MD coupled model.

    This is a one-shot bridge: the engine runs a full ``mcs``-step trajectory
    per invocation and cannot resume, so each ``update`` runs a complete
    simulation and emits the equilibrium (final-snapshot) observables. The
    observables are absolute sensor-style readings of the finished run, hence
    ``overwrite[float]`` ports.
    """

    config_schema = {
        "n_cross": {"_type": "integer", "_default": 2000},
        "fiber_stiffness": {"_type": "float", "_default": 200.0},
        "target_area": {"_type": "integer", "_default": 50},
        "mcs": {"_type": "integer", "_default": 300},
        "strands": {"_type": "integer", "_default": 800},
        "n_beads": {"_type": "integer", "_default": 9},
        "output_interval": {"_type": "integer", "_default": 20},
        "seed": {"_type": "integer", "_default": 12345678},
    }

    def inputs(self):
        return {}

    def outputs(self):
        return {
            "cell_area": "overwrite[float]",
            "giant_component": "overwrite[float]",
            "reorientation_q": "overwrite[float]",
            "densification": "overwrite[float]",
            "mcs": "overwrite[float]",
        }

    def _settings(self) -> dict:
        c = self.config
        return {
            "mcs": c["mcs"],
            "state_output_interval": c["output_interval"],
            "make_ecm.strands": c["strands"],
            "make_ecm.num_init_crosslinks": c["n_cross"],
            "make_ecm.spring_k": c["fiber_stiffness"],
            "make_ecm.crosslink_k": c["fiber_stiffness"],
            "make_ecm.bend_k": c["fiber_stiffness"],
            "cellular_potts.target_area": c["target_area"],
            "md_seed": c["seed"],
        }

    def update(self, state, interval):
        out = run_once(self._settings(), n_beads=self.config["n_beads"])
        eq = out["equilibrium"]
        return {
            "cell_area": float(eq.get("cell_area", float("nan"))),
            "giant_component": float(eq.get("giant_component", float("nan"))),
            "reorientation_q": float(eq.get("reorientation_q", float("nan"))),
            "densification": float(eq.get("densification", float("nan"))),
            "mcs": float(eq.get("mcs", float("nan"))),
        }
