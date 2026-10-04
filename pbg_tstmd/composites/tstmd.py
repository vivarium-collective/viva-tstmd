"""Composite generator for the real TST-MD bridge (merks-ecm-reciprocity-2d).

The genuine coupled cellular-Potts + bead-spring ECM engine is driven by
``TstmdEcmProcess`` in the sibling ``viva-tstmd`` package, which runs the real
simulator CPU-only inside Docker (TST C++ + HOOMD + MUSCLE3). This generator
wires that Process into a composite by ADDRESS (``local:TstmdEcmProcess``), so
the document builds structurally without importing the Docker-backed package;
actually RUNNING it requires ``viva-tstmd`` installed and a Docker daemon.

The study figures are baked from real runs of this engine (see
``viva_cpm_studies/visualizations/_tstmd_data.py`` and
``viva-tstmd/docs/engine-io.md``); this generator documents their provenance.
"""
from __future__ import annotations

from process_bigraph.composite_generator import composite_generator


@composite_generator(
    name="tstmd_coupled",
    description="Real TST-MD coupled CPM + bead-spring ECM run (via viva-tstmd / Docker).",
    parameters={
        "n_cross": {"type": "integer", "default": 5000,
                    "description": "Requested ECM crosslinks"},
        "fiber_stiffness": {"type": "float", "default": 200.0,
                            "description": "Fiber/crosslink spring constant"},
        "target_area": {"type": "integer", "default": 10,
                        "description": "Cell target area (contractility sink)"},
        "mcs": {"type": "integer", "default": 300,
                "description": "Monte Carlo steps"},
    },
    emitters=[{"address": "local:ParquetEmitter",
               "paths": ["stores/cell_area", "stores/reorientation_q",
                         "stores/giant_component"]}],
)
def tstmd_coupled(core=None, *, n_cross=5000, fiber_stiffness=200.0,
                  target_area=10, mcs=300):
    return {
        "tstmd": {
            "_type": "process",
            "address": "local:TstmdEcmProcess",   # resolved by viva-tstmd at run time
            "config": {"n_cross": n_cross, "fiber_stiffness": fiber_stiffness,
                       "target_area": target_area, "mcs": mcs},
            "interval": float(mcs),
            "outputs": {
                "cell_area": ["stores", "cell_area"],
                "giant_component": ["stores", "giant_component"],
                "reorientation_q": ["stores", "reorientation_q"],
                "densification": ["stores", "densification"],
                "mcs": ["stores", "mcs"],
            },
        },
        "stores": {"cell_area": 0.0, "giant_component": 0.0,
                   "reorientation_q": 0.0, "densification": 0.0, "mcs": 0.0},
        "emitter": {
            "_type": "step",
            "address": "local:RAMEmitter",
            "config": {"emit": {"cell_area": "float", "giant_component": "float",
                                "reorientation_q": "float"}},
            "inputs": {"cell_area": ["stores", "cell_area"],
                       "giant_component": ["stores", "giant_component"],
                       "reorientation_q": ["stores", "reorientation_q"]},
        },
    }
