"""Tests for the TST-MD Docker bridge runtime.

The end-to-end test is marked ``docker`` and skips when no daemon/image is
available, so CI without Docker stays green.
"""
import numpy as np
import pytest

from pbg_tstmd import runtime
from pbg_tstmd.runtime import (
    TstmdSession,
    cell_area,
    cell_centroid_physical,
    docker_available,
    image_present,
    write_override,
)


def test_run_without_docker_raises(monkeypatch, tmp_path):
    # Review Focus: daemon down -> clear, actionable error, not a hang.
    monkeypatch.setattr(runtime, "docker_available", lambda: False)
    sess = TstmdSession(settings={"mcs": 10}, workroot=tmp_path)
    with pytest.raises(RuntimeError, match="(?i)docker|colima"):
        sess.run()


def test_write_override(tmp_path):
    p = tmp_path / "o.ymmsl"
    write_override(p, {"mcs": 50, "cellular_potts.graphics": True,
                       "make_ecm.spring_k": 25.0})
    text = p.read_text()
    assert "ymmsl_version: v0.1" in text
    assert "  mcs: 50" in text
    assert "  cellular_potts.graphics: true" in text  # bool -> lowercase
    assert "  make_ecm.spring_k: 25.0" in text


def test_missing_image_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(runtime, "docker_available", lambda: True)
    monkeypatch.setattr(runtime, "image_present", lambda tag=runtime.IMAGE: False)
    sess = TstmdSession(settings={"mcs": 10}, workroot=tmp_path)
    with pytest.raises(RuntimeError, match="(?i)pull|image"):
        sess.run()


_HAVE_ENGINE = docker_available() and image_present()


@pytest.mark.docker
@pytest.mark.skipif(not _HAVE_ENGINE, reason="Docker daemon + TST-MD image required")
def test_short_run_returns_snapshots():
    # No workroot override: must stay under $HOME so Colima mounts it.
    sess = TstmdSession(
        settings={
            "mcs": 20, "state_output_interval": 10,
            "make_ecm.strands": 300, "make_ecm.num_init_crosslinks": 300,
            "equilibrate_ecm.md_its": 150, "simulate_ecm.md_its": 30,
        },
        timeout=1800,
    )
    snaps = sess.run()
    assert len(snaps) >= 2
    s0 = snaps[0]
    assert s0["cpm"].shape == (runtime.LATTICE, runtime.LATTICE)
    assert cell_area(s0) > 0                      # a cell exists
    c = cell_centroid_physical(s0)
    assert 0 <= c[0] <= s0["Lx"] and 0 <= c[1] <= s0["Ly"]
    assert s0["positions"].shape[1] == 2


def test_float_settings_coerced():
    # Review Focus: float-typed engine settings written as int crash MUSCLE3.
    import tempfile, pathlib
    p = pathlib.Path(tempfile.mktemp())
    write_override(p, {"make_ecm.spring_k": 15, "make_ecm.bend_k": 4,
                       "make_ecm.strands": 1000, "mcs": 300})
    text = p.read_text()
    assert "make_ecm.spring_k: 15.0" in text   # coerced to float
    assert "make_ecm.bend_k: 4.0" in text
    assert "make_ecm.strands: 1000" in text     # int stays int
    assert "mcs: 300" in text
