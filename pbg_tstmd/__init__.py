"""viva-tstmd: a process-bigraph bridge to the real TST-MD engine.

TST-MD is the hybrid cellular-Potts + bead-spring ECM simulator of
Tsingos et al. 2023 (Biophys. J.), the published 2D ancestor of Merks'
3D model (Keijzer & Merks 2026, arXiv:2609.02375). This package drives the
genuine engine CPU-only inside a Docker image and exposes it as a
process-bigraph Process, plus pure-Python metrics reproducing the paper's
network-remodeling measures.
"""
from .analysis import (
    densification_factor,
    giant_component_fraction,
    q_by_distance,
    ratio_with_nan_guard,
    reorientation_q,
)

__all__ = [
    "giant_component_fraction",
    "reorientation_q",
    "q_by_distance",
    "ratio_with_nan_guard",
    "densification_factor",
]
