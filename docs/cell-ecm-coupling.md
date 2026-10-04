# Cell–ECM coupling: why it fails here, and what must change

Deep-dive into *why* the bridged TST-MD engine does not reproduce Merks'
cell-driven ECM remodeling (the cell contracts but never drags its adhesions /
the fiber network inward), and a prioritized plan to make it behave like the
paper. Companion to [`engine-io.md`](engine-io.md). Written 2026-10-04.

## 1. The symptom (measured)

Across **every** regime tried — fiber stiffness (stiff→soft), crosslink
density, contractility `λ` (up to 400), MD noise `md_kT`, adhesion-annihilation
penalty, and `gradient` vs `uniform` displacement selection — the result is the
same:

- the cell contracts (area 13000 → ~0), but
- the **adhesion ring stays frozen** (mean radius ~58 px, <2% change) and
- near-fiber radius is unchanged (~33 px); densification ~0; q/q0 only ~1.09.

The cell **detaches from its adhesions** instead of pulling them (and the
crosslinked ECM) inward. This is the opposite of Keijzer & Merks 2026 Fig 4.

## 2. How the coupling actually works in the engine

(Source: `src/cellular_potts/ca.cpp`, `src/adhesions/*.cpp`,
`src/parameters/parameters.hpp`, on the `sbrshakibi/tst2-docker` build.)

Per CPM copy-attempt `(source → target)`, `CellularPotts::DeltaH` (ca.cpp:493)
sums the usual CPM terms **plus**, when `adhesions_enabled`:

```cpp
// ca.cpp:571-581
if (par.adhesions_enabled) {
    double adh_dh = adhesion_mover.move_dh({xp,yp},{x,y}, *adh_disp);
    DH += static_cast<int>(round(adh_dh));   // <-- rounded to int!
}
```

- `adh_dh` (adhesion_index.cpp:50 `move_dh`) is the **harmonic-spring energy
  change** of displacing the adhesion particle(s) in the retracted/extended
  voxel: `Σ_bonds (k/2)(‖to−other‖−r0)² − (k/2)(‖from−other‖−r0)²` (+ angle
  terms). It is proportional to the fiber spring constant `k`.
- On a **retraction** over an adhesion, `retraction_displacements`
  (adhesion_movement.cpp:25) offers the cell-neighbor pixels of the retracted
  voxel (i.e. move the adhesion ~1 px inward); `select_displacement` picks one.
- Selection mode is `par.adhesion_displacement_selection`
  (parameters.hpp:217, **default `"uniform"`** = random) vs `"gradient"` =
  the min-ΔH choice.

So the cell pulls the ECM *only* by retracting over adhesion voxels; each such
accepted retraction drags that adhesion ~1 px inward, stretching its fibers
(cost `adh_dh`). The crosslinked network should then transmit that pull.

The **3D paper uses argmin** (gradient): eq. (5)–(6),
`v = argmin_v Σ ΔH_ECM,α,v` — a coherent, energy-minimizing inward migration,
**not** the stock `uniform` random choice.

## 3. Root causes (why it never engages)

### R1 — Contraction energy dwarfs ECM resistance (the dominant cause)

`DeltaH` is **integer**, Metropolis temperature `T=50`. The area term
(ca.cpp:537-548) for a retraction with `σ'=MEDIUM` is

```
DH += λ · (1 − 2·(Area − TargetArea))
```

For a big, fully-contractile cell (`Area≈13000, TargetArea=0, λ=50`):
`DH ≈ 50·(1 − 2·13000) ≈ −1.3×10⁶`. The ECM work `adh_dh ≈ (k/2)Δr²` is only
**O(10–100)** per adhesion. So every retraction is accepted *regardless of the
ECM* — the cell collapses to ~0 in ~75 MCS, ignoring the matrix, and the
adhesions get swept past before they can migrate. The ECM can **never** balance
contraction at these settings.

The balance the paper achieves requires `λ·(Area−Target)` to be **comparable to
`adh_dh`** near equilibrium. With `λ=50` that only happens at `Area≈1`
(collapse). Lowering `λ` to ~1–5 moves the balance point to a *positive*
equilibrium area where the ECM resistance can halt contraction — the regime in
which adhesions are dragged gradually inward and held. **Higher `λ` (what was
tried) is exactly backwards.**

### R2 — Integer rounding of `adh_dh` erases soft-fiber coupling

`round(adh_dh)` (ca.cpp:575): with soft fibers (`spring_k≈15`) a 1-px adhesion
move costs `adh_dh < 0.5` → **rounds to 0** → the ECM exerts *zero* influence on
acceptance. With stiff fibers it is a large int that (combined with R1) still
loses to the λ term until collapse. There is no `spring_k` that both (a) avoids
rounding to zero and (b) is small enough for the cell to drag — unless `λ` is
also brought down (R1). For faithful mechanics the engine should accumulate
`adh_dh` as a **float** (drop the `round`).

### R3 — Stock uses `uniform`, the paper uses `gradient`

`adhesion_displacement_selection` defaults to `"uniform"` (random inward pixel);
the paper's argmin (`"gradient"`) gives coherent inward migration that tracks
the retracting boundary. Random selection scatters adhesion moves and weakens
coherent dragging.

### R4 — No physical→internal unit scaling

The paper's parameters (Keijzer 2025 Frontiers **Table 1**: `λ=4.96×10⁷ Nm⁻³`,
`K=3.1×10⁻² Nm⁻¹`, `K_cyto=3.1×10⁻⁴ Nm⁻¹`, fiber density `0.48 µm⁻²`, …) are
**scaled** so `adh_dh` and the CPM energies are commensurate (they state
`K = Y·A/L`, `Y=10⁶ Pa`, fiber diameter `0.125 µm`). Our arbitrary `spring_k`
(15/200/600) and `λ` (50) never land in that commensurate window — so we only
ever see "ECM ignored" (R1/R2) or "ECM rigid, cell detaches".

### R5 — Fast collapse vs gradual equilibrium

Because of R1, the cell collapses in ~75 MCS instead of relaxing to equilibrium
over the paper's ~2000 MCS. Adhesions need *time* (many small boundary steps +
MD relaxation between them) to migrate inward with the boundary; the instant
collapse leaves them behind.

## 4. What needs to happen next (prioritized)

1. **Lower `λ` so the ECM can balance contraction (R1).** Sweep
   `cellular_potts.lambda ∈ {1,2,5,10}` with `target_area=0`, moderate
   `spring_k≈200`, and measure whether the cell halts at a *positive*
   equilibrium area with the adhesion-ring radius shrinking. This is the single
   highest-leverage, cheapest test. *(Result below.)*
2. **Use `adhesion_displacement_selection: "gradient"` (R3)** — match the
   paper's argmin.
3. **Run to equilibrium (R5):** `mcs ≈ 1500–2000`, `simulate_ecm.md_its ≈ 100`,
   `md_kT = 0.001`.
4. **Do the unit scaling (R4):** derive internal `spring_k`, `crosslink_k`,
   `λ`, `target_area` from Frontiers Table 1 via `K=Y·A/L` and the paper's
   CPM–MD conversion, so `adh_dh ≈ λ·(equilibrium area)`. This is the
   principled fix that should reproduce Fig 2/3 magnitudes too.
5. **Patch the integer rounding (R2)** — carry `adh_dh` as a float in
   `DeltaH` (a small engine change in the fork) so soft-fiber coupling isn't
   quantized away. Rebuild the Docker image from a patched source.
6. **Get the authors' exact configuration** — the definitive answer. Check the
   Zenodo deposit (10.5281/zenodo.7906973) for `.par`/`.ymmsl` figure configs,
   or contact the Merks lab. Every public config we have is a demo blob, not the
   paper's contraction experiment.
7. **Percolation (separate study):** reaching giant-component ≈1 (Fig 2B) needs
   near-paper-scale fiber/crosslink density (`strands`, `num_init_crosslinks`)
   beyond the CPU-tractable sizes used here — likely GPU/HPC via HOOMD.

## 5. Evidence

- Measured adhesion-ring vs cell-area trajectories (every regime: ring frozen).
- `ca.cpp:537-581` (λ term vs `round(adh_dh)`); `parameters.hpp:217`
  (`uniform` default); paper eq. (5)–(6) (argmin) and Frontiers Table 1
  (scaled parameters).
- **Low-λ + gradient test (500 MCS, measured):**
  | λ   | cell area trajectory          | adhesion-ring radius |
  |-----|-------------------------------|----------------------|
  | 2   | 13505 → 1875 → **199** (halts)| 58.0 → 58.0 (frozen) |
  | 5   | 13635 → 1802 → **91** (halts) | 57.7 → 57.7 (frozen) |
  | 10  | 13093 → … → 23 → 0            | 57.4 → 57.4 (frozen) |

  **Half-confirmed, half-refuted.** R1 is real: lowering λ *does* halt the cell
  at a positive, λ-dependent equilibrium area (the Fig-2 "contraction reaches
  equilibrium" phenomenon appears). **But the adhesion ring stays frozen to the
  decimal in every case** — the adhesions genuinely never relocate, so the halt
  is almost certainly CPM-internal (surface tension `J` balancing a weak λ),
  **not** ECM-driven, and the Fig-4 *dragging/remodeling* still does not engage.

## 6. The real blocker, refined

The adhesion-ring radius being **constant to the decimal** across *every* test
(stiffness, λ, selection mode, penalty, md_kT) means the 200 adhesions created
at MCS 0 are **never moved after creation** — the cell contracts *past* them and
leaves them as orphaned `type=2` beads. The `retraction_displacements` path
(adhesion_movement.cpp:25) that should migrate an adhesion inward when the
boundary retracts over it is, in practice, **not relocating adhesions at all**
on this build. The highest-value next step is therefore a **mechanism probe
inside the engine**, not more parameter sweeps:

- Instrument `CellularPotts::DeltaH` / the copy-attempt loop (ca.cpp ~494, the
  `AmoebaeMove`/`ConvertSpin` path that passes `adh_disp`) to log, per MCS:
  how many copy-attempts touch an adhesion voxel, how many are accepted, and
  whether `adh_disp` is actually applied to the particle positions sent back to
  the MD. Determine whether (a) retractions over adhesions are always rejected,
  (b) the chosen displacement is always `(0,0)`/no-op, or (c) the applied move
  is overwritten by the MD step.
- Cross-check against the Merks lab's **actual** figure configuration (Zenodo
  10.5281/zenodo.7906973 contents, or direct contact) — the public image is a
  demo build and may differ from the paper's.

Items 2–5 of §4 (unit scaling, `gradient`, long runs, float `adh_dh`) remain
necessary, but they are **downstream** of fixing why adhesions don't relocate —
that is the gate.

> Bridge fix landed while investigating: float-typed settings (`spring_k`,
> `md_kT`, …) were written as ints and silently crashed every varied-stiffness
> run (`viva-tstmd` c6ff070), which had masked the whole parameter space.
