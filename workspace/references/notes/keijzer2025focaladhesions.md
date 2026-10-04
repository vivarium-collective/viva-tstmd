# Keijzer et al. 2025 — Dynamic focal adhesions (2D)

**Context / background** for the `merks-ecm-reciprocity-2d` investigation —
NOT the engine we bridge. Front. Cell Dev. Biol. 12:1462277.
doi:10.3389/fcell.2024.1462277.

Extends the 2D hybrid CPM + bead-spring model [tsingos2023beadspring] with
**dynamic mechanosensitive focal adhesions**: an ODE (after Novikova & Storm)
builds/breaks adhesion clusters under mechanical tension. Shows stiffness-
dependent spreading, local ECM remodeling, and ECM-alignment-dependent cell
elongation from mechanical reciprocity.

## Why it is background, not the bridge target
The 3D target paper [keijzer2026cpmecm3d] deliberately uses **static**
adhesions and names these dynamic focal adhesions as a *future* extension. So
to faithfully reproduce the 3D paper's Figs 2-4 in 2D we bridge the
**static-adhesion** engine [tsingos2023beadspring], not this dynamic-FA model.
Kept here for the mechanistic lineage and as the natural next step if the
investigation later explores adhesion dynamics.
