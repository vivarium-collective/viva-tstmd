"""Drive the real TST-MD engine (Docker, CPU-only) and read its state dumps.

``TstmdSession`` injects parameters via a YMMSL override, runs the coupled
cellular-Potts + bead-spring ECM simulation headless inside the pre-built
Docker image, and returns the per-dump snapshots as plain numpy (no
libmuscle on the host). See ``docs/engine-io.md`` for the confirmed recipe
and pickle schema.
"""
from __future__ import annotations

import os
import pickle
import shutil
import subprocess
import uuid
from pathlib import Path
from typing import Optional

import numpy as np

IMAGE = "sbrshakibi/tst2-docker:latest"
PLATFORM = "linux/amd64"
TST = "/home/mduser/Tissue-Simulation-Toolkit"
LATTICE = 200  # cpm sigma grid is LATTICE x LATTICE

# Colima only mounts paths under $HOME into its VM, so run dirs must live there.
DEFAULT_WORKROOT = Path.home() / ".cache" / "viva-tstmd" / "runs"

# Default settings for a CPU-tractable single contractile cell in a fibrous ECM.
# Overridden per-run by `settings`. Sizes are reduced vs the paper's HPC sweeps
# (documented per study); the mechanisms are size-robust.
DEFAULT_SETTINGS = {
    "mcs": 300,
    "state_output_interval": 20,
    "make_ecm.strands": 800,
    "make_ecm.num_init_crosslinks": 2000,
    "equilibrate_ecm.md_its": 500,
    "simulate_ecm.md_its": 60,
    # one viable contractile cell: initial size well above target so it
    # contracts toward equilibrium instead of vanishing (sic=50 dies; 120 lives).
    "cellular_potts.n_init_cells": 1,
    "cellular_potts.size_init_cells": 120,
    "cellular_potts.divisions": 0,
    "cellular_potts.target_area": 50,
}

# In-container parser: loads muscle3 pickles, converts to plain numpy, and
# dumps one plain pickle the host can read without libmuscle. Must NOT be named
# inspect.py (shadows stdlib inspect and breaks numpy/ymmsl imports).
_PARSER = r'''
import pickle, glob, numpy as np
def arr(x): return np.asarray(x.array) if hasattr(x, "array") else np.asarray(x)
files = sorted(glob.glob(
    "/work/run_adhesions_*/instances/state_dumper/workdir/state_*.pickle"))
snaps = []
for fn in files:
    d = pickle.load(open(fn, "rb"))
    cpm = arr(d["cpm_state"]["cpm"])
    P = d["ecm_state"]["particles"]; B = d["ecm_state"]["bonds"]
    snaps.append({
        "mcs": int(d["mcs"]),
        "Lx": float(d.get("Lx", d.get("sizex", cpm.shape[0]))),
        "Ly": float(d.get("Ly", d.get("sizey", cpm.shape[1]))),
        "cpm": cpm.astype(np.int32),
        "positions": arr(P["positions"]).astype(np.float64),
        "ptypes": arr(P["types"]).astype(np.int32),
        "bond_groups": arr(B["groups"]).astype(np.int32),
        "bond_types": arr(B["types"]).astype(np.int32),
    })
pickle.dump(snaps, open("/work/parsed.pkl", "wb"))
print(f"PARSED {len(snaps)} snapshots")
'''


def _run(cmd, timeout=None):
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def docker_available() -> bool:
    """True if a Docker daemon is reachable (on this machine: Colima)."""
    if shutil.which("docker") is None:
        return False
    try:
        return _run(["docker", "info"], timeout=20).returncode == 0
    except Exception:
        return False


def image_present(tag: str = IMAGE) -> bool:
    try:
        return _run(["docker", "image", "inspect", tag], timeout=20).returncode == 0
    except Exception:
        return False


# Settings the engine declares as float — an int here raises a MUSCLE3 TypeError
# ("is of type int, where float was expected"), so coerce them to float literals.
_FLOAT_SETTINGS = frozenset({
    "contour_length", "md_kT", "md_dt", "viscosity",
    "make_ecm.spring_r0", "make_ecm.spring_k", "make_ecm.crosslink_k",
    "make_ecm.helix_angle", "make_ecm.bend_t0", "make_ecm.bend_k",
    "make_ecm.crosslink_max_r", "make_ecm.crosslink_quant_step",
    "make_ecm.crosslink_bin_size",
    "cellular_potts.lambda2", "cellular_potts.saturation",
    "cellular_potts.dt", "cellular_potts.dx",
})


def _ymmsl_value(v, as_float: bool = False) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if as_float and isinstance(v, int):
        return repr(float(v))  # 15 -> "15.0"
    return str(v)


def write_override(path: Path, settings: dict) -> None:
    lines = ["ymmsl_version: v0.1", "settings:"]
    for k, v in settings.items():
        lines.append(f"  {k}: {_ymmsl_value(v, as_float=(k in _FLOAT_SETTINGS))}")
    path.write_text("\n".join(lines) + "\n")


class TstmdSession:
    """One coupled TST-MD run with a given parameter override.

    Parameters
    ----------
    settings:
        YMMSL setting overrides (merged over ``DEFAULT_SETTINGS``), e.g.
        ``{"make_ecm.num_init_crosslinks": 4000, "make_ecm.spring_k": 50}``.
    workroot:
        Parent dir for the run (must be under $HOME for Colima). Defaults to
        ``~/.cache/viva-tstmd/runs``.
    """

    def __init__(
        self,
        settings: Optional[dict] = None,
        workroot: Optional[os.PathLike] = None,
        image: str = IMAGE,
        platform: str = PLATFORM,
        timeout: int = 3600,
        keep: bool = False,
    ):
        self.settings = {**DEFAULT_SETTINGS, **(settings or {})}
        self.image = image
        self.platform = platform
        self.timeout = timeout
        self.keep = keep
        root = Path(workroot) if workroot else DEFAULT_WORKROOT
        self.rundir = Path(root) / uuid.uuid4().hex[:12]

    def _preflight(self) -> None:
        if not docker_available():
            raise RuntimeError(
                "Docker daemon is not reachable. Start it with `colima start` "
                "(or your Docker runtime) before running a TstmdSession."
            )
        if not image_present(self.image):
            raise RuntimeError(
                f"TST-MD image '{self.image}' is not present. Pull it once with:\n"
                f"  docker pull --platform {self.platform} {self.image}"
            )

    def run(self) -> list[dict]:
        """Run the coupled sim and return the list of parsed snapshots.

        Each snapshot is a dict with plain numpy arrays:
        ``{mcs, Lx, Ly, cpm, positions, ptypes, bond_groups, bond_types}``.
        """
        self._preflight()
        home = str(Path.home())
        if not str(self.rundir.resolve()).startswith(home):
            raise RuntimeError(
                f"Run dir {self.rundir} is not under $HOME ({home}). Colima only "
                "mounts paths under $HOME into its VM, so the bind mount would be "
                "empty. Use the default workroot or one under $HOME."
            )
        self.rundir.mkdir(parents=True, exist_ok=True)
        write_override(self.rundir / "override.ymmsl", self.settings)
        (self.rundir / "parse_states.py").write_text(_PARSER)

        inner = (
            f". {TST}/venv/bin/activate && "
            f"muscle_manager --start-all "
            f"{TST}/ymmsl/adhesions.ymmsl {TST}/ymmsl/dump_state.ymmsl "
            f"/work/override.ymmsl > /work/mm.log 2>&1 ; "
            f"python /work/parse_states.py"
        )
        cmd = [
            "docker", "run", "--rm", "--platform", self.platform, "--user", "root",
            "-e", "QT_QPA_PLATFORM=offscreen",
            "-v", f"{self.rundir}:/work", "-w", "/work",
            self.image, "bash", "-lc", inner,
        ]
        proc = _run(cmd, timeout=self.timeout)
        parsed = self.rundir / "parsed.pkl"
        if not parsed.exists():
            log = (self.rundir / "mm.log").read_text()[-2000:] if (self.rundir / "mm.log").exists() else ""
            raise RuntimeError(
                "TST-MD run produced no parsed snapshots.\n"
                f"docker stdout/stderr:\n{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}\n"
                f"muscle manager log tail:\n{log}"
            )
        snaps = pickle.loads(parsed.read_bytes())
        if not self.keep:
            # keep parsed.pkl + logs, drop the bulky run dir
            for p in self.rundir.glob("run_adhesions_*"):
                shutil.rmtree(p, ignore_errors=True)
        return snaps

    def close(self) -> None:
        if not self.keep:
            shutil.rmtree(self.rundir, ignore_errors=True)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


# --- extraction helpers: engine snapshot -> analysis-ready quantities --------

def cell_area(snap: dict) -> int:
    """Number of lattice sites occupied by cells (sigma >= 1)."""
    return int((snap["cpm"] >= 1).sum())


def cell_centroid(snap: dict) -> np.ndarray:
    """Cell centroid in the shared pixel frame used by bead ``positions``.

    Empirically, bead positions live in the same [0, L] pixel frame as the
    cpm lattice (fibers spill outside it), so the mapping is the identity:
    the centroid is just the mean nonzero-sigma pixel. (x = columns,
    y = rows, matching positions[:,0]=x, positions[:,1]=y.)
    """
    ys, xs = np.nonzero(snap["cpm"] >= 1)
    L = snap["cpm"].shape[0]
    if xs.size == 0:
        return np.array([L / 2.0, L / 2.0])
    return np.array([xs.mean(), ys.mean()])


# backward-compat alias
cell_centroid_physical = cell_centroid


def fiber_beads(snap: dict, n_beads: int):
    """Return (centers, directions) for fibers, from free-particle chains.

    Strands are ``n_beads`` consecutive free particles; the center is the
    middle bead and the direction the middle bond vector. Boundary/adhesion
    particles are excluded.
    """
    pos = snap["positions"]
    free = snap["ptypes"] == 0
    idx = np.nonzero(free)[0]
    # free particles are laid out as contiguous strands of n_beads
    n_full = (idx.size // n_beads) * n_beads
    idx = idx[:n_full].reshape(-1, n_beads)
    mid = n_beads // 2
    centers = pos[idx[:, mid]]
    directions = pos[idx[:, mid]] - pos[idx[:, mid - 1]]
    return centers, directions
