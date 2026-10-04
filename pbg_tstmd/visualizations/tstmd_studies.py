"""Visualizations for the merks-ecm-reciprocity-2d investigation.

Every figure is built from REAL output of the bridged TST-MD engine (the
published 2D static-adhesion hybrid CPM + bead-spring ECM of Tsingos et al.
2023, driven CPU-only via viva-tstmd / Docker+HOOMD), baked into
``_tstmd_data.py``. These reproduce, at the mechanism level and in 2D, the
phenomena of Keijzer & Merks 2026 (arXiv:2609.02375): cell contraction
deforming a crosslinked ECM, fibers reorienting toward the cell (q/q0 > 1),
and crosslinking building network connectivity.
"""
from __future__ import annotations

import base64
import json
import zlib

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from viva_superpowers.visualization import as_visualization

from ._tstmd_data import TSTMD_DATA

# palette (readable in light + dark)
C_AREA = "#e4572e"   # contraction
C_Q = "#1b9e77"      # reorientation q/q0
C_GC = "#3773b8"     # connectivity
C_FIBER = "rgba(130,130,130,0.45)"
C_CELL = "rgba(228,87,46,0.85)"


def _snapshot_trace(snap, name):
    fib = snap["fibers"]
    cell = snap["cell"]
    fx = [p[0] for p in fib]; fy = [p[1] for p in fib]
    cx = [p[0] for p in cell]; cy = [p[1] for p in cell]
    return (
        go.Scatter(x=fx, y=fy, mode="markers", name="ECM fibers",
                   marker=dict(size=2.5, color=C_FIBER), showlegend=False,
                   hoverinfo="skip"),
        go.Scatter(x=cx, y=cy, mode="markers", name="cell",
                   marker=dict(size=3, color=C_CELL), showlegend=False,
                   hoverinfo="skip"),
    )


def _build_contraction_reorientation() -> str:
    dyn = TSTMD_DATA["dynamics"]
    traj = dyn["trajectory"]
    mcs = [m["mcs"] for m in traj]
    area = [m["cell_area"] for m in traj]
    qratio = [m["q_ratio"] for m in traj]
    start, end = dyn["start_snapshot"], dyn["end_snapshot"]
    L = start["L"]

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"colspan": 2, "secondary_y": True}, None], [{}, {}]],
        row_heights=[0.52, 0.48], vertical_spacing=0.14, horizontal_spacing=0.08,
        subplot_titles=(
            "Cell contracts while fibers reorient toward it",
            f"ECM + cell at start (MCS {start['mcs']})",
            f"ECM + cell at end (MCS {end['mcs']})",
        ),
    )

    fig.add_trace(go.Scatter(x=mcs, y=area, mode="lines+markers", name="cell area",
                             line=dict(color=C_AREA, width=3),
                             marker=dict(size=6)), row=1, col=1, secondary_y=False)
    fig.add_trace(go.Scatter(x=mcs, y=qratio, mode="lines+markers", name="reorientation q/q₀",
                             line=dict(color=C_Q, width=3, dash="dot"),
                             marker=dict(size=6)), row=1, col=1, secondary_y=True)
    fig.add_hline(y=1.0, line=dict(color=C_Q, width=1, dash="dash"),
                  opacity=0.4, row=1, col=1, secondary_y=True)

    for tr in _snapshot_trace(start, "start"):
        fig.add_trace(tr, row=2, col=1)
    for tr in _snapshot_trace(end, "end"):
        fig.add_trace(tr, row=2, col=2)

    fig.update_xaxes(title_text="Monte Carlo step (MCS)", row=1, col=1)
    fig.update_yaxes(title_text="cell area (lattice sites)", color=C_AREA,
                     row=1, col=1, secondary_y=False)
    fig.update_yaxes(title_text="q/q₀  (fibers toward cell)", color=C_Q,
                     row=1, col=1, secondary_y=True)
    for c in (1, 2):
        fig.update_xaxes(range=[-10, L + 10], row=2, col=c, scaleanchor=f"y{3 if c==1 else 4}",
                         showticklabels=False)
        fig.update_yaxes(range=[-10, L + 10], row=2, col=c, showticklabels=False)

    fig.update_layout(
        template="plotly_white", height=760,
        title=dict(text="<b>Mechanical reciprocity (2D):</b> a contractile cell "
                        "remodels a real crosslinked ECM", x=0.02, font=dict(size=17)),
        legend=dict(orientation="h", y=1.07, x=0.0),
        margin=dict(l=70, r=70, t=90, b=50),
    )
    return fig.to_html(full_html=True, include_plotlyjs="cdn",
                       config={"displayModeBar": False})


def _build_network_connectivity() -> str:
    conn = TSTMD_DATA["connectivity"]
    created = [c["created"] for c in conn]
    gc = [c["giant_component"] for c in conn]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=created, y=gc, mode="lines+markers",
        line=dict(color=C_GC, width=3), marker=dict(size=10, color=C_GC),
        hovertemplate="crosslinks: %{x}<br>giant component: %{y:.3f}<extra></extra>",
        name="giant component"))
    fig.add_annotation(x=created[0], y=gc[0], text="no crosslinks →<br>disconnected",
                       showarrow=True, arrowhead=2, ax=50, ay=-30, font=dict(size=11))
    fig.update_layout(
        template="plotly_white", height=460,
        title=dict(text="<b>Crosslinking builds ECM network connectivity</b><br>"
                        "<span style='font-size:12px;color:#666'>Giant component rises sharply "
                        "once crosslinks appear, then saturates (~0.3) at these reduced densities "
                        "— a percolation-onset signature, not full percolation.</span>",
                   x=0.02, font=dict(size=16)),
        xaxis=dict(title="crosslinks created", type="log"),
        yaxis=dict(title="giant-component fraction", range=[0, 1]),
        margin=dict(l=70, r=40, t=90, b=55),
    )
    return fig.to_html(full_html=True, include_plotlyjs="cdn",
                       config={"displayModeBar": False})


C_XLINK = "rgba(34,160,94,0.75)"   # crosslinks (paper draws these green)


def _load_frames():
    from ._tstmd_spatial import TSTMD_FRAMES_B64  # lazy: baked after capture
    return json.loads(zlib.decompress(base64.b64decode(TSTMD_FRAMES_B64)))


def _grid_from(f):
    """Decode the baked cell-occupancy grid and 2x-downsample (stays solid)."""
    g = np.frombuffer(zlib.decompress(base64.b64decode(f["grid"])),
                      np.uint8).reshape(f["gshape"])
    ny, nx = g.shape
    g = g[: ny // 2 * 2, : nx // 2 * 2].reshape(ny // 2, 2, nx // 2, 2).max(axis=(1, 3))
    return g


C_ADH_LINK = "rgba(255,140,0,0.75)"   # adhesion springs: the cell's grip on the ECM
C_ADH = "#8e44ad"                      # adhesion sites (paper draws these purple)


def _frame_traces(f):
    g = _grid_from(f)
    L = f["L"]
    ny, nx = g.shape
    return [
        go.Scatter(x=f["fx"], y=f["fy"], mode="lines", name="ECM fibers",
                   line=dict(color="rgba(150,150,150,0.30)", width=0.8), hoverinfo="skip"),
        go.Scatter(x=f["xx"], y=f["xy"], mode="lines", name="crosslinks",
                   line=dict(color="rgba(34,160,94,0.30)", width=0.8), hoverinfo="skip"),
        go.Scatter(x=f.get("lx", []), y=f.get("ly", []), mode="lines",
                   name="cell–ECM adhesions",
                   line=dict(color=C_ADH_LINK, width=1.4), hoverinfo="skip"),
        go.Heatmap(z=g, x0=L / nx / 2, dx=L / nx, y0=L / ny / 2, dy=L / ny,
                   colorscale=[[0, "rgba(255,255,255,0)"], [1, "rgba(228,87,46,0.95)"]],
                   zmin=0, zmax=1, showscale=False, hoverinfo="skip", name="cell"),
        go.Scatter(x=f.get("ax", []), y=f.get("ay", []), mode="markers",
                   name="adhesion sites",
                   marker=dict(size=5, color=C_ADH, line=dict(width=0)), hoverinfo="skip"),
    ]


def _build_contraction_video() -> str:
    frames = _load_frames()
    L = frames[0]["L"]
    lo, hi = 0, L   # zoom to the simulation domain (crop the outer fiber buffer)

    def ann(f):
        return [dict(x=0.02, y=0.98, xref="paper", yref="paper", showarrow=False,
                     align="left", bgcolor="rgba(255,255,255,0.7)",
                     text=f"<b>MCS {f['mcs']}</b>  ·  cell area {f['area']}",
                     font=dict(size=13))]

    fig = go.Figure(
        data=_frame_traces(frames[0]),
        frames=[go.Frame(data=_frame_traces(f), name=str(f["mcs"]),
                         layout=go.Layout(annotations=ann(f))) for f in frames],
    )
    slider_steps = [dict(method="animate", label=str(f["mcs"]),
                         args=[[str(f["mcs"])],
                               dict(mode="immediate", frame=dict(duration=0, redraw=True),
                                    transition=dict(duration=0))]) for f in frames]
    fig.update_layout(
        template="plotly_white", height=720, width=720, plot_bgcolor="#0e1116",
        title=dict(text="<b>The simulation unfolding:</b> a cell (red) grips the ECM at "
                        "adhesions (purple) and contracts<br>"
                        "<span style='font-size:11px;color:#888'>orange = cell–ECM adhesion "
                        "links · gray = fibers · green = crosslinks</span>",
                   x=0.02, font=dict(size=15)),
        xaxis=dict(range=[lo, hi], showticklabels=False, scaleanchor="y",
                   showgrid=False, zeroline=False, constrain="domain"),
        yaxis=dict(range=[lo, hi], showticklabels=False, showgrid=False, zeroline=False),
        legend=dict(orientation="h", y=1.10, x=0.0),
        annotations=ann(frames[0]),
        margin=dict(l=20, r=20, t=90, b=60),
        updatemenus=[dict(type="buttons", direction="left", showactive=False,
                          x=0.0, xanchor="left", y=-0.04, yanchor="top",
                          buttons=[
                              dict(label="▶ Play", method="animate",
                                   args=[None, dict(frame=dict(duration=450, redraw=True),
                                                    fromcurrent=True,
                                                    transition=dict(duration=0))]),
                              dict(label="⏸ Pause", method="animate",
                                   args=[[None], dict(mode="immediate",
                                                      frame=dict(duration=0, redraw=True))])])],
        sliders=[dict(active=0, x=0.12, len=0.88, y=-0.02, yanchor="top",
                      currentvalue=dict(prefix="MCS ", font=dict(size=13)),
                      steps=slider_steps)],
    )
    return fig.to_html(full_html=True, include_plotlyjs="cdn",
                       config={"displayModeBar": False})


@as_visualization(inputs={"mcs": "list[float]"},
                  name="TstmdContractionReorientation", demo={"mcs": [0.0]})
def update_tstmd_contraction_reorientation(state):
    """Real TST-MD: cell contracts while fibers reorient toward it (q/q0 > 1)."""
    return {"html": _build_contraction_reorientation()}


@as_visualization(inputs={"mcs": "list[float]"},
                  name="TstmdContractionVideo", demo={"mcs": [0.0]})
def update_tstmd_contraction_video(state):
    """Real TST-MD time-lapse: the cell contracting and remodeling the ECM, frame by frame."""
    return {"html": _build_contraction_video()}


def _load_conn_frames():
    from ._tstmd_spatial_conn import TSTMD_CONN_FRAMES_B64  # lazy: baked after capture
    return json.loads(zlib.decompress(base64.b64decode(TSTMD_CONN_FRAMES_B64)))


def _net_traces(f):
    return [
        go.Scatter(x=f["fx"], y=f["fy"], mode="lines", name="ECM fibers",
                   line=dict(color=C_FIBER, width=1), hoverinfo="skip"),
        go.Scatter(x=f["xx"], y=f["xy"], mode="lines", name="crosslinks",
                   line=dict(color=C_XLINK, width=1.3), hoverinfo="skip"),
    ]


def _build_connectivity_video() -> str:
    frames = _load_conn_frames()
    L = frames[0]["L"]

    def ann(f):
        return [dict(x=0.02, y=0.98, xref="paper", yref="paper", showarrow=False,
                     align="left", bgcolor="rgba(255,255,255,0.7)",
                     text=f"<b>{f['ncross_created']} crosslinks</b>  ·  "
                          f"giant component {f['gc']:.2f}", font=dict(size=13))]

    fig = go.Figure(
        data=_net_traces(frames[0]),
        frames=[go.Frame(data=_net_traces(f), name=str(f["ncross_created"]),
                         layout=go.Layout(annotations=ann(f))) for f in frames],
    )
    steps = [dict(method="animate", label=str(f["ncross_created"]),
                  args=[[str(f["ncross_created"])],
                        dict(mode="immediate", frame=dict(duration=0, redraw=True),
                             transition=dict(duration=0))]) for f in frames]
    fig.update_layout(
        template="plotly_white", height=680, width=720,
        title=dict(text="<b>Building the network:</b> crosslinking knits loose "
                        "fibers into a connected ECM", x=0.02, font=dict(size=16)),
        xaxis=dict(range=[-5, L + 5], showticklabels=False, scaleanchor="y",
                   showgrid=False, zeroline=False),
        yaxis=dict(range=[-5, L + 5], showticklabels=False, showgrid=False, zeroline=False),
        legend=dict(orientation="h", y=1.10, x=0.0),
        annotations=ann(frames[0]),
        margin=dict(l=20, r=20, t=90, b=60),
        updatemenus=[dict(type="buttons", direction="left", showactive=False,
                          x=0.0, xanchor="left", y=-0.04, yanchor="top",
                          buttons=[
                              dict(label="▶ Play", method="animate",
                                   args=[None, dict(frame=dict(duration=900, redraw=True),
                                                    fromcurrent=True,
                                                    transition=dict(duration=0))]),
                              dict(label="⏸ Pause", method="animate",
                                   args=[[None], dict(mode="immediate",
                                                      frame=dict(duration=0, redraw=True))])])],
        sliders=[dict(active=0, x=0.12, len=0.88, y=-0.02, yanchor="top",
                      currentvalue=dict(prefix="crosslinks ", font=dict(size=13)),
                      steps=steps)],
    )
    return fig.to_html(full_html=True, include_plotlyjs="cdn",
                       config={"displayModeBar": False})


@as_visualization(inputs={"mcs": "list[float]"},
                  name="TstmdNetworkVideo", demo={"mcs": [0.0]})
def update_tstmd_network_video(state):
    """Real TST-MD: the fiber network gaining connectivity as crosslinks increase."""
    return {"html": _build_connectivity_video()}


@as_visualization(inputs={"mcs": "list[float]"},
                  name="TstmdNetworkConnectivity", demo={"mcs": [0.0]})
def update_tstmd_network_connectivity(state):
    """Real TST-MD: crosslinking builds network connectivity (percolation onset)."""
    return {"html": _build_network_connectivity()}
