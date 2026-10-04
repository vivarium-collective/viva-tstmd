# Tsingos et al. 2023 — Hybrid CPM + bead-spring ECM (2D)

**THE ENGINE** bridged by `viva-tstmd` and consumed by this investigation.
Biophys. J. 122(13):2609-2622. doi:10.1016/j.bpj.2023.05.013.
bioRxiv 2022.06.10.495667.

The published 2D hybrid cellular Potts + bead-spring ECM model that the 3D
target paper [keijzer2026cpmecm3d] extends. A contractile CPM cell is coupled
to a crosslinked fiber network (beads + linear/angular springs) integrated to
mechanical equilibrium by HOOMD-blue, via **static adhesion particles**.
Establishes the same mechanical-reciprocity phenomena (crosslinking-limited
contraction, fiber reorientation, densification) in 2D.

## Code
- Zenodo deposit: `10.5281/zenodo.7906973` (`TST-MD` + `PyTST`).
- Tissue Simulation Toolkit: github.com/mathbioleiden/Tissue-Simulation-Toolkit
  (branch `TST2.0`), built with `make with_adhesions` (TST C++ + HOOMD +
  MUSCLE3).
- We run it CPU-only via the pre-built image `sbrshakibi/tst2-docker`; see
  `viva-tstmd/docs/engine-io.md`.

*(Registered bib-only — the preprint PDF download is Cloudflare rate-limited;
drop the PDF into `papers/tsingos-2023-cpm-beadspring.pdf` + add a
`references_pdfs` entry later if wanted.)*
