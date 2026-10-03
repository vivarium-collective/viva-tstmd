"""Tests for the sweep driver and the process-bigraph surface."""
import numpy as np
import pytest

from pbg_tstmd import processes
from pbg_tstmd.processes import TstmdEcmProcess, run_sweep
from pbg_tstmd.runtime import docker_available, image_present


def test_process_schema_and_ports():
    proc = TstmdEcmProcess(config={}, core=None) if False else None  # noqa
    # inspect class-level schema without instantiation (no engine needed)
    assert "n_cross" in TstmdEcmProcess.config_schema
    assert "fiber_stiffness" in TstmdEcmProcess.config_schema
    inst = TstmdEcmProcess.__new__(TstmdEcmProcess)
    outs = inst.outputs()
    assert set(outs) >= {"cell_area", "giant_component", "reorientation_q",
                         "densification", "mcs"}


def test_process_registers_in_core():
    from process_bigraph import allocate_core
    core = allocate_core()
    core.register_link("TstmdEcmProcess", TstmdEcmProcess)
    assert "TstmdEcmProcess" in core.list_processes()


def test_run_sweep_continues_on_failure(monkeypatch):
    # Review Focus: a failed sweep point is recorded and the sweep continues.
    def fake_run_once(settings, n_beads=9, **kw):
        if settings["make_ecm.num_init_crosslinks"] == 999:
            raise RuntimeError("boom")
        return {"trajectory": [], "equilibrium": {"cell_area": 42},
                "q0": float("nan"), "crosslinks_created": 0}
    monkeypatch.setattr(processes, "run_once", fake_run_once)

    results = run_sweep("make_ecm.num_init_crosslinks", [100, 999, 300])
    assert len(results) == 3
    errs = [r for r in results if "error" in r]
    oks = [r for r in results if "error" not in r]
    assert len(errs) == 1 and errs[0]["value"] == 999
    assert len(oks) == 2 and all(r["equilibrium"]["cell_area"] == 42 for r in oks)


_HAVE_ENGINE = docker_available() and image_present()


@pytest.mark.docker
@pytest.mark.skipif(not _HAVE_ENGINE, reason="Docker daemon + TST-MD image required")
def test_run_once_trajectory():
    out = processes.run_once(
        {"mcs": 20, "state_output_interval": 10,
         "make_ecm.strands": 300, "make_ecm.num_init_crosslinks": 300,
         "equilibrate_ecm.md_its": 150, "simulate_ecm.md_its": 30},
    )
    assert len(out["trajectory"]) >= 2
    m = out["trajectory"][0]
    assert set(m) >= {"mcs", "cell_area", "reorientation_q",
                      "giant_component", "densification", "q_ratio"}
    assert m["cell_area"] > 0
    assert 0.0 <= m["giant_component"] <= 1.0
