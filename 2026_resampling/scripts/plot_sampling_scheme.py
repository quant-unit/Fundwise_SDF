#!/usr/bin/env python3
"""Create the within-vintage resampling scheme figure.

This script mirrors the TikZ figure in ../figure/sampling_scheme.tex, but is
kept as a reusable plotting template for later empirical versions with actual
fund counts and realized random partition assignments.

Usage:
    python plot_sampling_scheme.py

Optional:
    python plot_sampling_scheme.py --out-dir ../figure
"""

from __future__ import annotations

import argparse
from pathlib import Path


VINTAGE_COUNTS = {
    2005: 8,
    2006: 10,
    2007: 14,
    2008: 13,
    2009: 11,
    2010: 9,
    2011: 7,
    2012: 6,
}

SEASONED_VINTAGES = set(range(2005, 2011))
UNSEASONED_VINTAGES = {2011, 2012}
APPROACH_BOUNDARY_YEAR = max(SEASONED_VINTAGES)
APPROACH_BOUNDARY_LABEL = "approach boundary after 2010:\ncash-flow above, normalized Euler below"

FOLD_SEQUENCE = [
    ("estimation_fold_1", "Estimation fold 1", "#1B7837"),
    ("estimation_fold_2", "Estimation fold 2", "#4DAF4A"),
    ("estimation_fold_3", "Estimation fold 3", "#A6D96A"),
    ("validation_fold", "Validation fold", "#2C7FB8"),
]


def build_points(vintage_counts: dict[int, int]) -> list[dict[str, object]]:
    """Return one plotting record per fund dot."""
    points: list[dict[str, object]] = []
    for vintage, count in vintage_counts.items():
        for fund_index in range(1, count + 1):
            fold_key, fold_label, color = FOLD_SEQUENCE[(fund_index - 1) % len(FOLD_SEQUENCE)]
            points.append(
                {
                    "vintage": vintage,
                    "fund_index": fund_index,
                    "fold_key": fold_key,
                    "fold_label": fold_label,
                    "color": color,
                }
            )
    return points


def add_approach_boundary(ax, y_pos: dict[int, int]) -> None:
    ax.axhline(
        y_pos[APPROACH_BOUNDARY_YEAR] - 0.5,
        color="#303030",
        linewidth=0.8,
        linestyle=(0, (7, 5)),
    )
    ax.text(
        9.55,
        y_pos[APPROACH_BOUNDARY_YEAR] + 0.05,
        APPROACH_BOUNDARY_LABEL,
        fontsize=8.7,
        color="#303030",
        va="top",
        bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )


def plot_sampling_scheme(out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for this reusable Python plot. "
            "Install it in the active environment, or compile the TikZ version "
            "in ../figure/sampling_scheme.tex."
        ) from exc

    out_dir.mkdir(parents=True, exist_ok=True)
    points = build_points(VINTAGE_COUNTS)
    vintages = sorted(VINTAGE_COUNTS)
    y_pos = {vintage: len(vintages) - 1 - i for i, vintage in enumerate(vintages)}

    fig, ax = plt.subplots(figsize=(11.8, 5.2), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.axhspan(y_pos[2012] - 0.45, y_pos[2011] + 0.5, color="#F5F5F5", zorder=0)

    for point in points:
        ax.scatter(
            point["fund_index"],
            y_pos[point["vintage"]],
            s=56,
            color=point["color"],
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )

    add_approach_boundary(ax, y_pos)

    ax.set_title(
        "Example split: three estimation folds plus one validation fold",
        loc="left",
        fontsize=15,
        fontweight="bold",
        color="#303030",
        pad=18,
    )
    ax.text(
        0,
        1.02,
        "Each dot is one fund; folds are assigned within vintage year across both pricing approaches.",
        transform=ax.transAxes,
        fontsize=11,
        color="#303030",
        va="bottom",
    )

    ax.set_yticks([y_pos[vintage] for vintage in vintages])
    ax.set_yticklabels([str(vintage) for vintage in vintages], fontsize=11)
    ax.set_ylabel("Vintage year", fontsize=11, labelpad=34)
    ax.set_xlabel("Fund index within vintage year", fontsize=11, labelpad=10)

    max_count = max(VINTAGE_COUNTS.values())
    ax.set_xlim(0.5, max_count + 0.5)
    ax.set_xticks(range(1, max_count + 1))
    ax.tick_params(axis="x", labelsize=9, length=0)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="both", color="#ECECEC", linewidth=0.7)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_visible(False)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=color, markeredgecolor="white")
        for _, _, color in FOLD_SEQUENCE
    ]
    labels = [label for _, label, _ in FOLD_SEQUENCE]
    ax.legend(
        handles,
        labels,
        title="Fold assignment",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.03, 0.58),
        borderaxespad=0,
        labelspacing=1.15,
        fontsize=11,
        title_fontsize=12,
    )

    ax.text(
        1.03,
        0.23,
        "Schematic counts are illustrative. In the empirical figure,\n"
        "the x-axis can use actual fund counts and the colors can\n"
        "use the realized random partition.",
        transform=ax.transAxes,
        fontsize=9,
        color="#303030",
        va="top",
    )

    for suffix in ("pdf", "png"):
        fig.savefig(out_dir / f"sampling_scheme_from_python.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_sampling_scheme_nav_pricing(out_dir: Path) -> None:
    """Hybrid-estimator variant: unseasoned vintages (2011, 2012) are included
    and receive fold assignments like seasoned vintages rather than being ignored."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for this reusable Python plot. "
            "Install it in the active environment, or compile the TikZ version "
            "in ../figure/sampling_scheme.tex."
        ) from exc

    out_dir.mkdir(parents=True, exist_ok=True)
    points = build_points(VINTAGE_COUNTS)
    vintages = sorted(VINTAGE_COUNTS)
    y_pos = {vintage: len(vintages) - 1 - i for i, vintage in enumerate(vintages)}

    fig, ax = plt.subplots(figsize=(11.8, 5.2), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.axhspan(y_pos[2012] - 0.45, y_pos[2011] + 0.5, color="#F5F5F5", zorder=0)

    for point in points:
        ax.scatter(
            point["fund_index"],
            y_pos[point["vintage"]],
            s=56,
            color=point["color"],
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )

    add_approach_boundary(ax, y_pos)

    ax.set_title(
        "Example split under the hybrid estimator",
        loc="left",
        fontsize=15,
        fontweight="bold",
        color="#303030",
        pad=30,
    )
    ax.text(
        0,
        1.02,
        "Each dot is one fund; folds are assigned within vintage year across both pricing approaches.",
        transform=ax.transAxes,
        fontsize=11,
        color="#303030",
        va="bottom",
    )

    ax.set_yticks([y_pos[vintage] for vintage in vintages])
    ax.set_yticklabels([str(vintage) for vintage in vintages], fontsize=11)
    ax.set_ylabel("Vintage year", fontsize=11, labelpad=34)
    ax.set_xlabel("Fund index within vintage year", fontsize=11, labelpad=10)

    max_count = max(VINTAGE_COUNTS.values())
    ax.set_xlim(0.5, max_count + 0.5)
    ax.set_xticks(range(1, max_count + 1))
    ax.tick_params(axis="x", labelsize=9, length=0)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="both", color="#ECECEC", linewidth=0.7)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_visible(False)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=color, markeredgecolor="white")
        for _, _, color in FOLD_SEQUENCE
    ]
    labels = [label for _, label, _ in FOLD_SEQUENCE]
    ax.legend(
        handles,
        labels,
        title="Fold assignment",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.03, 0.58),
        borderaxespad=0,
        labelspacing=1.15,
        fontsize=11,
        title_fontsize=12,
    )

    ax.text(
        1.03,
        0.23,
        "Schematic counts are illustrative. In the empirical figure,\n"
        "the x-axis can use actual fund counts and the colors can\n"
        "use the realized random partition.",
        transform=ax.transAxes,
        fontsize=9,
        color="#303030",
        va="top",
    )

    for suffix in ("pdf", "png"):
        fig.savefig(out_dir / f"sampling_scheme_nav_pricing.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "figure",
        help="Directory where sampling_scheme_from_python.pdf/png should be written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plot_sampling_scheme(args.out_dir)
    plot_sampling_scheme_nav_pricing(args.out_dir)


if __name__ == "__main__":
    main()
