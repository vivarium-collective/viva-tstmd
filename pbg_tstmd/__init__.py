"""viva-tstmd: a process-bigraph bridge to the real TST-MD engine + the
merks-ecm-reciprocity-2d investigation workspace.

TST-MD is the hybrid cellular-Potts + bead-spring ECM simulator of
Tsingos et al. 2023 (Biophys. J.), the published 2D ancestor of Merks'
3D model (Keijzer & Merks 2026, arXiv:2609.02375). This package drives the
genuine engine CPU-only inside Docker and exposes it as a process-bigraph
Process, with pure-Python ECM-remodeling metrics and the investigation's
visualizations.
"""
from .analysis import (
    densification_factor,
    giant_component_fraction,
    q_by_distance,
    ratio_with_nan_guard,
    reorientation_q,
)
from .processes import TstmdEcmProcess
from .core import build_core
from . import composites  # noqa: F401  (fires @composite_generator decorators)

__all__ = [
    "TstmdEcmProcess",
    "build_core",
    "giant_component_fraction",
    "reorientation_q",
    "q_by_distance",
    "ratio_with_nan_guard",
    "densification_factor",
]
