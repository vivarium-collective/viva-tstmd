# Keijzer & Merks 2026 — 3D hybrid CPM + deformable fiber network

**THE REPRODUCTION TARGET** for the `merks-ecm-reciprocity-2d` investigation.
arXiv:2609.02375 (in revision, Biomech. Model. Mechanobiol.).

A 3D hybrid cellular Potts model coupled to a discrete, deformable bead-spring
ECM. A contractile cell attached to the fiber network via **static adhesion
sites** contracts; the ECM resists via fiber elasticity and crosslinking.

## Headline results (the figures we reproduce in 2D)
- **Fig 2** — ECM crosslinking limits contraction: equilibrium cell volume
  rises with crosslink density; a network *percolation* (giant-component)
  threshold emerges (~15,000 crosslinks in 3D).
- **Fig 3** — fiber stiffness `k` introduces an effective drag; equilibrium
  volume rises with `k` (and a spurious non-zero residual at very high `k`,
  an unrealistic regime the authors avoid by using `k≈10`).
- **Fig 4** — contractile cells remodel the ECM: fiber **reorientation**
  `q/q0 > 1` (stronger with crosslinking) and **biphasic densification**
  (peaks at intermediate crosslinking).

## Why no code / why we bridge the 2D ancestor
No public code or data; a custom 3D CPM + HOOMD run on an HPC/GPU cluster. We
reproduce the **2D analogs** of Figs 2-4 with the *real* published 2D
static-adhesion engine [tsingos2023beadspring] via `viva-tstmd`. Fidelity is
directional (trend/threshold/biphasic shape), not quantitative/3D-identical.

Uses **static** adhesions; names dynamic focal adhesions
[keijzer2025focaladhesions] as future work.
