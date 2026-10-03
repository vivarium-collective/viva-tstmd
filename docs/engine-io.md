# TST-MD engine I/O (empirically confirmed)

Facts the bridge (`pbg_tstmd/runtime.py`, `processes.py`, `analysis.py`)
depends on, confirmed against the real engine on 2026-10-03. This is the
M0 Task 0.2 deliverable from the implementation plan.

## The engine

- **Image:** `sbrshakibi/tst2-docker:latest` — a pre-built TST 2.0 + HOOMD +
  MUSCLE3 stack (Debian buster, Python 3.10.13). Pulled, not built; skips the
  multi-hour HOOMD compile. Built for **linux/amd64** → runs under emulation
  on Apple Silicon (slower, but CPU-only and correct).
- **Build source:** the image clones `sbr-shakibi/Tissue-Simulation-Toolkit`
  branch `TST2.0` (a fork). Its dump format differs slightly from the
  `mathbioleiden` branch (e.g. `Lx/Ly` keys, not `sizex/sizey`), so trust
  *this* document over upstream source reading.
- **Install root in image:** `/home/mduser/Tissue-Simulation-Toolkit` (`$TST`).
  Run as `--user root` (the image sets `USER mduser` but never created that
  passwd entry). venv at `$TST/venv`; ymmsl configs at `$TST/ymmsl/`.

## Headless run recipe (the one that works)

```bash
docker run --rm --platform linux/amd64 --user root \
  -e QT_QPA_PLATFORM=offscreen \
  -v <host_rundir>:/work -w /work sbrshakibi/tst2-docker:latest \
  bash -lc '. /home/mduser/Tissue-Simulation-Toolkit/venv/bin/activate && \
    muscle_manager --start-all \
      /home/mduser/Tissue-Simulation-Toolkit/ymmsl/adhesions.ymmsl \
      /home/mduser/Tissue-Simulation-Toolkit/ymmsl/dump_state.ymmsl \
      /work/override.ymmsl'
```

- **`QT_QPA_PLATFORM=offscreen` is mandatory** and graphics must stay **on**.
  The MCS loop (`for t in 0..par.mcs: TimeStep()`) lives in the Qt graphics
  backend (`src/graphics/graph.cpp`). Setting `cellular_potts.graphics: false`
  takes a different path that runs the loop ~once → `simulate_ecm` loses its
  peer and the whole run collapses with no pickles. Do NOT override graphics;
  render offscreen instead.
- **Mounts:** Colima only mounts paths under `$HOME` into its VM. Run dirs
  must live under `$HOME` (we use `viva-tstmd/.runs/…`), not `/tmp`.
- **Overrides:** an extra YMMSL file (later file wins) is the parameter-
  injection surface. The sweep knobs:
  - `make_ecm.num_init_crosslinks` — crosslinking (Figs 2, 4).
    NOTE: *requested* ≠ *created* (only pairs within `crosslink_max_r` are
    linked); the engine logs "Number of crosslinks created: N" — use the
    created count as the x-axis, read from the make_ecm instance log or the
    bond graph.
  - `make_ecm.spring_k`, `make_ecm.crosslink_k`, `make_ecm.bend_k` — fiber
    stiffness `k` (Fig 3).
  - `cellular_potts.target_area`, `cellular_potts.lambda` — contractility.
  - `cellular_potts.n_init_cells`, `size_init_cells`, `divisions` — cell(s).
    `divisions:0` + `n_init_cells:1` = one cell (the paper's single
    contractile cell); default `n_init_cells:100` = a blob. A single tiny
    cell with a small target area can vanish immediately and end the run —
    size the initial cell well above its target.
  - `mcs`, `state_output_interval`, `make_ecm.strands`,
    `equilibrate_ecm.md_its`, `simulate_ecm.md_its` — runtime/cost.
- **Outputs:** pickle dumps every `state_output_interval` MCS at
  `<rundir>/run_adhesions_<ts>/instances/state_dumper/workdir/state_<mcs:05d>.pickle`.

## Pickle schema (one snapshot)

```python
{
  'Lx': float, 'Ly': float,          # physical domain size (um); e.g. 100.0
  'mcs': int,
  'cpm_state': {
     'cpm':  int32  (200, 200),      # spin field sigma. lattice is 200x200.
                                     #   sigma == 0 -> medium; sigma >= 1 -> cell id;
                                     #   sigma == -1 seen (wall/fixed). CELL AREA = (cpm >= 1).sum()
     'pde':  float32 (1, 200, 200),  # chemical field (unused by our metrics)
     'act_field': dict,              # sparse Act activity {"i,j": int}; empty when Act off
  },
  'ecm_state': {
     'particles': {
        'positions': float64 (N, 2), # bead coords in MD/physical units (0..Lx)
        'types':     int32  (N,),    # 0=free, 1=boundary, 2=adhesion, 3=excluded
     },
     'bond_types': {'r0': float64 (T,), 'k': float64 (T,)},   # per-type rest len + stiffness
     'bonds': {
        'groups': int32 (M, 2),      # particle-id pairs (one row per bond)
        'types':  int32 (M,),        # index into bond_types; NamedBondTypes.fiber == 0
     },
     'angle_cst_types': {'t0': float64 (A,), 'k': float64 (A,)},
     'angle_csts': {'groups': int32 (P, 3), 'types': int32 (P,)},  # bending triples
  }
}
```

### Coordinate systems (important)

`cpm` is on the **200×200 lattice**; bead `positions` are in **physical units
`0..Lx`** (Lx=100 → lattice spacing 0.5, i.e. lattice = positions × (200/Lx)).
Reorientation and densification need the cell centroid and bead positions in
the *same* frame — the bridge's extraction layer converts the lattice centroid
to physical units (× Lx/200) before calling `analysis.py` (whose functions are
frame-agnostic and take already-consistent coordinates).

### `.array` wrappers → parse inside the container

Pickle leaves are MUSCLE3 `Data`/Grid objects (`.array` yields the numpy
array), and unpickling imports `libmuscle`/`ymmsl`. **Parse pickles inside the
container** (venv has libmuscle) and re-emit plain numpy (`.npz`) for the host.
Do NOT name the parser `inspect.py` — it shadows the stdlib `inspect` and
breaks numpy/ymmsl imports.

## Fiber extraction

A "fiber"/strand is a run of `make_ecm.beads` consecutive free particles
(sequential ids per strand). Fiber center = middle bead; direction = middle
bond vector. Boundary (type 1) and adhesion (type 2) particles are not part of
fiber chains. Confirm the exact bead-id layout against a real dump before
trusting it (recorded as a follow-up in `runtime.py`'s extraction).
