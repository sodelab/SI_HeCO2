#!/usr/bin/env python3
"""Create a publication-ready He-CO2 energy-level diagram."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


# Energies in cm^-1 relative to the variational ground state.
D0 = 16.106

BOUND_LEVELS = [
    (7.700, r"$\Delta E(0,1)$"),
    (8.637, r"$\Delta E(0,2)$"),
    (12.339, r"$\Delta E(0,3)$"),
    (15.001, r"$\Delta E(0,4)$"),
]

# Selected local-mode excitation energies from the VCI calculation.
VCI_LEVELS = [
    (12.44, r"$\nu_b$"),
    (18.79, r"$\nu_s$"),
    (26.95, r"$2\nu_b$"),
    (32.05, r"$\nu_b+\nu_s$"),
    (41.19, r"$2\nu_s$"),
]


def draw_level(ax, x, width, energy, color, linestyle, linewidth=2.0):
    ax.hlines(
        energy,
        x - width / 2,
        x + width / 2,
        colors=color,
        linestyles=linestyle,
        linewidth=linewidth,
        zorder=3,
    )


def main():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
            "font.size": 9,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig, ax = plt.subplots(figsize=(6.7, 5.6))

    bound_x = 0.0
    vci_x = 1.75
    level_width = 0.58
    ymax = 44.0
    xmin = -0.72
    xmax = 2.42

    bound_color = "#1F4E79"
    vci_color = "#9A4F00"
    # Lightly distinguish the region above dissociation.
    ax.axhspan(D0, ymax, color="#D9D9D9", alpha=0.24, zorder=0)
    # Variational bound-state ladder.
    for energy, label in BOUND_LEVELS:
        draw_level(ax, bound_x, level_width, energy, bound_color, "solid")
        # Separate the two closely spaced levels so their labels remain legible.
        label_y = {7.700: 7.30, 8.637: 9.05}.get(energy, energy)
        for x_text, text, align in [
            (bound_x - level_width / 2 - 0.07, f"{energy:.3f}", "right"),
            (bound_x + level_width / 2 + 0.07, label, "left"),
        ]:
            ax.annotate(
                text,
                xy=(bound_x - level_width / 2 if align == "right" else bound_x + level_width / 2, energy),
                xytext=(x_text, label_y),
                textcoords="data",
                ha=align,
                va="center",
                color="#222222",
                arrowprops=(
                    {"arrowstyle": "-", "color": "#777777", "lw": 0.6}
                    if label_y != energy
                    else None
                ),
            )

    # Selected VCI energies. Dashed lines reinforce that these are nominal
    # local-mode excitations rather than variationally bound levels.
    for energy, label in VCI_LEVELS:
        draw_level(ax, vci_x, level_width, energy, vci_color, (0, (4, 2)))
        ax.text(
            vci_x - level_width / 2 - 0.07,
            energy,
            f"{energy:.2f}",
            ha="right",
            va="center",
            color="#222222",
        )
        ax.text(
            vci_x + level_width / 2 + 0.07,
            energy,
            label,
            ha="left",
            va="center",
            color="#222222",
        )

    # Both columns contain excitation energies, so use one common zero-energy
    # reference rather than plotting a ground state in only one column.
    ax.hlines(
        0.0,
        xmin,
        xmax,
        colors="#4A4A4A",
        linestyles="solid",
        linewidth=1.4,
        zorder=2,
    )
    ax.text(
        0.87,
        0.55,
        "Ground-state energy reference",
        ha="center",
        va="bottom",
        color="#4A4A4A",
        fontsize=8.5,
    )

    # Dissociation threshold across both ladders.
    ax.axhline(
        D0,
        color=bound_color,
        linestyle=(0, (7, 3)),
        linewidth=1.8,
        zorder=2,
    )
    ax.text(
        0.16,
        D0 + 0.55,
        rf"Dissociation threshold  $D_0={D0:.3f}\ \mathrm{{cm}}^{{-1}}$",
        ha="center",
        va="bottom",
        color=bound_color,
        fontsize=9,
    )

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-1.4, ymax)
    ax.set_ylabel(r"Energy relative to the ground state (cm$^{-1}$)")
    ax.set_xticks([bound_x, vci_x])
    ax.set_xticklabels(
        [r"Variational $J=0$ bound-state excitations", "Selected VCI local-mode excitations"]
    )
    ax.tick_params(axis="x", length=0, pad=7)
    ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35, 40])
    ax.tick_params(axis="y", direction="out", length=3, width=0.8)

    # Quiet publication-style frame.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.grid(False)

    fig.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.14)

    stem = Path("heco2_energy_levels")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.05)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.05)
    fig.savefig(
        stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.05
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
