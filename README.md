# viva-tstmd

A process-bigraph bridge to the **real** TST-MD engine — the hybrid
cellular-Potts + bead-spring ECM simulator of Tsingos et al. 2023
(Biophys. J.), the published 2D ancestor of Merks' 3D model
(Keijzer & Merks 2026, arXiv:2609.02375).

The engine runs **CPU-only inside a Docker image** (`sbrshakibi/tst2-docker`,
linux/amd64; no GPU). This package injects parameters via a YMMSL override,
runs the coupled simulation headless, parses its pickle state dumps, and
exposes it as a process-bigraph Process plus pure-Python metrics that
reproduce the paper's ECM-remodeling measures (percolation, fiber
reorientation q/q0, densification).

See `docs/engine-io.md` for the confirmed run recipe and pickle schema.

Built to drive the `merks-ecm-reciprocity-2d` reproduction investigation in
the viva-cpm workspace.
