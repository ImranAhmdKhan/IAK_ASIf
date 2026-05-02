from __future__ import annotations

import numpy as np

try:
    import matplotlib
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - registers 3D projection

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def _draw_wavy(ax, x1, y1, x2, y2, n_waves=5, amp_frac=0.06, **kw):
    t = np.linspace(0, 1, 300)
    dx, dy = x2 - x1, y2 - y1
    length = max(np.hypot(dx, dy), 1e-9)
    px, py = -dy / length, dx / length
    amp = amp_frac * length
    wave = amp * np.sin(n_waves * 2 * np.pi * t)
    ax.plot(x1 + t * dx + wave * px, y1 + t * dy + wave * py, **kw)


def _draw_zigzag(ax, x1, y1, x2, y2, n_zigs=7, amp_frac=0.05, **kw):
    n = n_zigs * 2 + 1
    t = np.linspace(0, 1, n)
    dx, dy = x2 - x1, y2 - y1
    length = max(np.hypot(dx, dy), 1e-9)
    px, py = -dy / length, dx / length
    amp = amp_frac * length
    displace = np.array([amp * ((-1) ** i) if i % 2 == 1 else 0 for i in range(n)])
    ax.plot(x1 + t * dx + displace * px, y1 + t * dy + displace * py, **kw)


def _draw_curled(ax, x1, y1, x2, y2, rad=0.35, **kw):
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    length = max(np.hypot(dx, dy), 1e-9)
    px, py = -dy / length, dx / length
    cp = (mx + rad * length * px, my + rad * length * py)
    path = mpath.Path(
        [(x1, y1), cp, (x2, y2)],
        [mpath.Path.MOVETO, mpath.Path.CURVE3, mpath.Path.CURVE3],
    )
    lw = kw.pop("lw", kw.pop("linewidth", 1.0))
    color = kw.pop("color", "k")
    alpha = kw.pop("alpha", 1.0)
    ls = kw.pop("ls", kw.pop("linestyle", "-"))
    zorder = kw.pop("zorder", 3)
    ax.add_patch(
        mpatches.PathPatch(
            path,
            fill=False,
            edgecolor=color,
            linewidth=lw,
            alpha=alpha,
            linestyle=ls,
            zorder=zorder,
        )
    )


def _draw_fancy_arrow(
    ax,
    x1,
    y1,
    x2,
    y2,
    style="Straight",
    color="k",
    lw=1.0,
    alpha=0.75,
    zorder=10,
    head_size=0.012,
):
    kw = dict(color=color, lw=lw, alpha=alpha, zorder=zorder)
    dx, dy = x2 - x1, y2 - y1
    length = max(np.hypot(dx, dy), 1e-9)
    ux, uy = dx / length, dy / length
    gap = head_size * length
    bx2, by2 = x2 - gap * ux, y2 - gap * uy

    if style == "Wavy":
        _draw_wavy(ax, x1, y1, bx2, by2, n_waves=5, amp_frac=0.06, **kw)
    elif style == "Zigzag":
        _draw_zigzag(ax, x1, y1, bx2, by2, n_zigs=6, amp_frac=0.05, **kw)
    elif style in ("Curled (Arc)", "Curled"):
        _draw_curled(ax, x1, y1, bx2, by2, rad=0.30, **kw)
    elif style in ("Smooth (Bezier)", "Smooth"):
        xs = np.linspace(x1, bx2, 60)
        ts = (xs - x1) / max(bx2 - x1, 1e-9)
        ax.plot(xs, y1 + (by2 - y1) * (3 * ts**2 - 2 * ts**3), **kw)
    else:
        ax.plot([x1, bx2], [y1, by2], **kw)

    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(bx2, by2),
        arrowprops=dict(arrowstyle="->", color=color, lw=lw, alpha=alpha, mutation_scale=10 * lw),
        zorder=zorder + 1,
    )

