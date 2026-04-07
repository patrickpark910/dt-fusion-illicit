"""
1D blanket schematics (three rows: FLiBe, pebble bed, DCLL).
Layer thicknesses from ``Python/parameters.py``; horizontal extent is proportional
to thickness (0.2 cm W as reference). Over-thick layers render as compressed strips

"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from Python.parameters import (
    DCLL_FW_CM,
    DCLL_FWF_CM,
    DCLL_FWC_CM,
    DCLL_FWB_CM,
    DCLL_BR1_CM,
    DCLL_D1_CM,
    DCLL_BR2_CM,
    DCLL_D2_CM,
    DCLL_BR3_CM,
    DCLL_IM_CM,
    DCLL_BP_CM,
    DCLL_SS_CM,
    FLIBE_FW_CM,
    FLIBE_BE_CM,
    FLIBE_ST1_CM,
    FLIBE_BR1_CM,
    FLIBE_ST2_CM,
    FLIBE_BR2_CM,
    FLIBE_ST3_CM,
    HCPB_FW_CM,
    HCPB_ST1_CM,
    HCPB_BR1_I_CM,
    HCPB_BR1_O_CM,
    HCPB_ST2_CM,
)


GRAY = "#7F7F7F"
BE_GRAY = "#B0B0B0"

COLORS = {
    "w": "#2a2a2a",
    "structure": GRAY,
    "flibe": "#90EE90",
    "pebble": "#F08080",
    "pbli": "#87CEEB",
    "coolant": "#B3E5FC",
    "be": BE_GRAY,
    "back_plate": BE_GRAY,
}

PLASMA_CMAP = LinearSegmentedColormap.from_list(
    "plasma_sym",
    ["#FFF5FB", "#FF74D2", "#FFF5FB"],
    N=256,
)

REF_CM = FLIBE_FW_CM
W_MIN_REF = 0.2
UNIT_PER_CM = W_MIN_REF / REF_CM
W_RAW_BREAK_THRESHOLD = 5.5
BREAK_STRIP_EQUIV_CM = 2.0
HALF_WIDTH_BUDGET = 24.0


def _fmt_cm(x: float) -> str:
    return f"{x:g}"


def _raw_width(th_cm: float) -> float:
    return th_cm * UNIT_PER_CM


def _layer_width(
    th_cm: float,
    break_strip_equiv_cm: float = BREAK_STRIP_EQUIV_CM,
) -> tuple[float, bool]:
    """Pre-scale width; second value True if layer is drawn as a compressed ``//`` strip."""
    raw = _raw_width(th_cm)
    w_break_strip = _raw_width(break_strip_equiv_cm)
    if raw > W_RAW_BREAK_THRESHOLD:
        return w_break_strip, True
    return raw, False


def _stack_extent(
    layers: list,
    break_strip_equiv_cm: float = BREAK_STRIP_EQUIV_CM,
) -> float:
    return sum(_layer_width(t, break_strip_equiv_cm)[0] for t, _, _ in layers)


def _label_text(th_cm: float, name: str) -> str:
    return f"{_fmt_cm(th_cm)} cm\n{name}"


def _estimate_label_width(text: str, fs_pt: float) -> float:
    lines = text.split("\n")
    m = max(len(line) for line in lines) if lines else len(text)
    return max(m * 0.018 * fs_pt * 0.4, 0.6)


def _assign_stagger_levels(
    xs: list[float],
    texts: list[str],
    fs: float,
    min_gap: float = 0.55,
) -> list[int]:
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    levels = [0] * len(xs)
    occupied: list[list[tuple[float, float]]] = []

    for oi in order:
        cx = xs[oi]
        half = _estimate_label_width(texts[oi], fs)
        for lev in range(14):
            while len(occupied) <= lev:
                occupied.append([])
            lo, hi = cx - half, cx + half
            bad = any(
                not (hi < lo2 - min_gap or lo > hi2 + min_gap)
                for lo2, hi2 in occupied[lev]
            )
            if not bad:
                occupied[lev].append((lo, hi))
                levels[oi] = lev
                break
        else:
            levels[oi] = 13

    return levels


def plot_emma_blanket_schematic(
    out_path: str | None = None,
    *,
    bar_h: float = 1,
    plasma_w: float = 3.85,
    dpi: int = 200,
    svg_dpi: int = 400,
    label_fs: float = 8.5,
    stagger_dy: float = 0.44,
    save_svg: bool = True,
) -> str:
    """Write PNG (and optionally SVG); return path to the PNG."""
    flibe_sym = [
        (FLIBE_FW_CM, "W", "w"),
        (FLIBE_BE_CM, "Be", "be"),
        (FLIBE_ST1_CM, "V-4Cr-4Ti", "structure"),
        (FLIBE_BR1_CM, "FLiBe", "flibe"),
        (FLIBE_ST2_CM, "V-4Cr-4Ti", "structure"),
        (FLIBE_BR2_CM, "FLiBe", "flibe"),
        (FLIBE_ST3_CM, "V-4Cr-4Ti", "structure"),
    ]
    specs = {
        "FLiBe": {
            "title": "FLiBe",
            "left": flibe_sym,
            "right": list(flibe_sym),
            "symmetric_note": True,
            "break_strip_equiv_cm": 3.0,
        },
        "Pebble bed": {
            "title": "Pebble bed",
            "left": [
                (HCPB_FW_CM, "W", "w"),
                (HCPB_ST1_CM, "Eurofer", "structure"),
                (HCPB_BR1_I_CM, "pebbles", "pebble"),
                (HCPB_ST2_CM, "Eurofer", "structure"),
            ],
            "right": [
                (HCPB_FW_CM, "W", "w"),
                (HCPB_ST1_CM, "Eurofer", "structure"),
                (HCPB_BR1_O_CM, "pebbles", "pebble"),
                (HCPB_ST2_CM, "Eurofer", "structure"),
            ],
            "symmetric_note": False,
            "break_strip_equiv_cm": 4.0,
        },
        "Lead-lithium": {
            "title": "Lead-lithium",
            "left": [
                (DCLL_FW_CM, "W", "w"),
                (DCLL_FWF_CM, "F82H", "structure"),
                (DCLL_FWC_CM, "coolant", "coolant"),
                (DCLL_FWB_CM, "F82H", "structure"),
                (DCLL_BR1_CM, "Pb-Li", "pbli"),
                (DCLL_D1_CM, "divider", "structure"),
                (DCLL_BR2_CM, "Pb-Li", "pbli"),
                (DCLL_IM_CM, "manifold", "structure"),
                (DCLL_BP_CM, "back plate", "back_plate"),
                (DCLL_SS_CM, "shield", "structure"),
            ],
            "right": [
                (DCLL_FW_CM, "W", "w"),
                (DCLL_FWF_CM, "F82H", "structure"),
                (DCLL_FWC_CM, "coolant", "coolant"),
                (DCLL_FWB_CM, "F82H", "structure"),
                (DCLL_BR1_CM, "Pb-Li", "pbli"),
                (DCLL_D1_CM, "divider", "structure"),
                (DCLL_BR2_CM, "Pb-Li", "pbli"),
                (DCLL_D2_CM, "divider", "structure"),
                (DCLL_BR3_CM, "Pb-Li", "pbli"),
                (DCLL_IM_CM, "manifold", "structure"),
                (DCLL_BP_CM, "back plate", "back_plate"),
                (DCLL_SS_CM, "shield", "structure"),
            ],
            "symmetric_note": False,
        },
    }

    rows = ["FLiBe", "Pebble bed", "Lead-lithium"]

    max_half = 0.0
    for r in rows:
        pl, pr = specs[r]["left"], specs[r]["right"]
        beq = specs[r].get("break_strip_equiv_cm", BREAK_STRIP_EQUIV_CM)
        max_half = max(max_half, _stack_extent(pl, beq), _stack_extent(pr, beq))
    max_half += plasma_w / 2
    global_scale = min(1.0, HALF_WIDTH_BUDGET / max_half) if max_half > 0 else 1.0

    pw = plasma_w * global_scale

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
        }
    )

    row_gap = 3.25
    y0 = 1.0
    y_positions = [y0 + 2 * row_gap, y0 + row_gap, y0]
    label_anchor_below_bar = 0.64
    leader_drop = 0.06

    x_max = HALF_WIDTH_BUDGET + 2.2

    fig_w, fig_h = 13.5, 7.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_aspect("auto")
    ax.axis("off")

    y_top = max(y_positions)
    y_bot = min(y_positions)
    ymax = y_top + bar_h + 1.45
    ymin = y_bot - 1.65 - 14 * stagger_dy

    ax.set_xlim(-x_max, x_max)
    ax.set_ylim(ymin, ymax)

    fs_header = 17
    fs_row = 17
    fs_plasma = 15

    ax.text(
        -x_max + 8.0,
        ymax - 0.6,
        "Inboard",
        ha="left",
        va="center",
        fontsize=fs_header,
        color=GRAY,
    )
    ax.text(
        x_max - 8.0,
        ymax - 0.6,
        "Outboard",
        ha="right",
        va="center",
        fontsize=fs_header,
        color=GRAY,
    )

    plasma_h = bar_h * 1.06
    plasma_xpad = 0.06

    def draw_breaks(cx: float, cy: float, h: float) -> None:
        for dx in (-0.07, 0.07):
            ax.plot(
                [cx + dx - 0.11, cx + dx + 0.11],
                [cy - h / 2 + 0.05, cy + h / 2 - 0.05],
                color="white",
                lw=1.7,
                zorder=5,
            )

    for row_name, cy in zip(rows, y_positions):
        spec = specs[row_name]
        symmetric = spec.get("symmetric_note", False)
        break_equiv = spec.get("break_strip_equiv_cm", BREAK_STRIP_EQUIV_CM)

        def build_rects(layers: list, side: str) -> list[dict]:
            x_edge = pw / 2 if side == "right" else -pw / 2
            direction = 1.0 if side == "right" else -1.0
            out = []
            for th_cm, name, key in layers:
                w_pre, brk = _layer_width(th_cm, break_equiv)
                w = w_pre * global_scale
                x0 = x_edge
                x1 = x_edge + direction * w
                out.append(
                    {
                        "x0": x0,
                        "x1": x1,
                        "cx": (x0 + x1) / 2,
                        "th_cm": th_cm,
                        "name": name,
                        "key": key,
                        "brk": brk,
                    }
                )
                x_edge = x1
            return out

        ax.text(
            -x_max + 0.85,
            cy + 0.28,
            spec["title"],
            ha="left",
            va="center",
            fontsize=fs_row,
            color=GRAY,
        )

        left_rects = build_rects(spec["left"], "left")
        right_rects = build_rects(spec["right"], "right")

        for r in left_rects + right_rects:
            face = COLORS.get(r["key"], GRAY)
            ax.add_patch(
                mpatches.Rectangle(
                    (min(r["x0"], r["x1"]), cy - bar_h / 2),
                    abs(r["x1"] - r["x0"]),
                    bar_h,
                    facecolor=face,
                    edgecolor="none",
                    linewidth=0,
                    zorder=1,
                )
            )
            if r["brk"]:
                draw_breaks(r["cx"], cy, bar_h)

        grad = np.linspace(0.0, 1.0, 512, dtype=np.float64).reshape(1, -1)
        im = ax.imshow(
            grad,
            extent=[
                -pw / 2 - plasma_xpad,
                pw / 2 + plasma_xpad,
                cy - plasma_h / 2,
                cy + plasma_h / 2,
            ],
            aspect="auto",
            cmap=PLASMA_CMAP,
            zorder=4,
            interpolation="bilinear",
            interpolation_stage="rgba",
            filternorm=False,
        )
        im.set_clip_on(True)

        ax.text(
            0,
            cy,
            "plasma",
            ha="center",
            va="center",
            fontsize=fs_plasma,
            color="white",
            zorder=5,
        )

        base_y = cy - bar_h / 2 - label_anchor_below_bar

        def annotate_stack(rects: list[dict]) -> int:
            if not rects:
                return 0
            texts = [_label_text(r["th_cm"], r["name"]) for r in rects]
            xs = [r["cx"] for r in rects]
            levels = _assign_stagger_levels(xs, texts, label_fs)
            for r, txt, lev in zip(rects, texts, levels):
                y_t = base_y - lev * stagger_dy
                ax.annotate(
                    txt,
                    xy=(r["cx"], cy - bar_h / 2 - leader_drop),
                    xytext=(r["cx"], y_t),
                    ha="center",
                    va="top",
                    ma="center",
                    fontsize=label_fs,
                    color=GRAY,
                    linespacing=1.08,
                    arrowprops=dict(
                        arrowstyle="-",
                        color=GRAY,
                        lw=0.35,
                        shrinkA=0,
                        shrinkB=2,
                    ),
                )
            return max(levels) if levels else 0

        pair_w = (
            len(left_rects) > 0
            and len(right_rects) > 0
            and left_rects[0]["key"] == "w"
            and right_rects[0]["key"] == "w"
        )
        left_ann = left_rects[1:] if pair_w else left_rects
        right_ann = right_rects[1:] if pair_w else right_rects

        max_lev_l = annotate_stack(left_ann)
        if not symmetric:
            annotate_stack(right_ann)

        if pair_w:
            wl = left_rects[0]
            txt = _label_text(wl["th_cm"], wl["name"])
            y_anchor = cy - bar_h / 2 - leader_drop
            y_t = base_y
            ax.plot(
                [-pw / 2, 0],
                [y_anchor, y_t],
                color=GRAY,
                lw=0.3,
                zorder=6,
                clip_on=False,
                solid_capstyle="round",
            )
            ax.plot(
                [pw / 2, 0],
                [y_anchor, y_t],
                color=GRAY,
                lw=0.3,
                zorder=6,
                clip_on=False,
                solid_capstyle="round",
            )
            ax.text(
                0,
                y_t,
                txt,
                ha="center",
                va="top",
                ma="center",
                fontsize=label_fs,
                color=GRAY,
                linespacing=1.08,
                zorder=7,
            )

        if symmetric:
            sym_note_x = 4.2
            sym_note_y = base_y - (max_lev_l + 1.35) * stagger_dy + stagger_dy
            ax.text(
                sym_note_x,
                sym_note_y,
                "same as inboard",
                ha="left",
                va="top",
                fontsize=10,
                style="italic",
                color=GRAY,
            )

    fig.tight_layout()

    save_kw = dict(bbox_inches="tight", facecolor="white", edgecolor="none", pad_inches=0.02)

    if out_path is not None:
        fig.savefig(out_path, dpi=dpi, **save_kw)
        plt.close(fig)
        return out_path

    os.makedirs("./Figures/PNG", exist_ok=True)
    png_path = os.path.join(".", "Figures", "PNG", "fig_emmaBlanket.png")
    fig.savefig(png_path, dpi=dpi, format="png", **save_kw)

    if save_svg:
        os.makedirs("./Figures/SVG", exist_ok=True)
        svg_path = os.path.join(".", "Figures", "SVG", "fig_emmaBlanket.svg")
        with plt.rc_context(
            {
                "svg.fonttype": "none",
                "svg.hashsalt": None,
            }
        ):
            fig.savefig(
                svg_path,
                format="svg",
                dpi=svg_dpi,
                **save_kw,
            )

    plt.close(fig)
    return png_path


if __name__ == "__main__":
    png = plot_emma_blanket_schematic()
    print(f"Saved {png}")
    print(f"Saved {os.path.join('.', 'Figures', 'SVG', 'fig_emmaBlanket.svg')}")
