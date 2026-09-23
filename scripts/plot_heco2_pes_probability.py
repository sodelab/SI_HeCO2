#!/usr/bin/env python3
"""This script reads a bound-state wavefunction and a 2D intermolecular PES, then
combines them to visualize the ground-state probability distribution in the
(R, theta) plane.

Overview
--------
The script is built for quantum-bound-state analysis of the He--CO2 complex. It takes
as input:

1. A formatted BOUND wavefunction file produced by a bound-state quantum chemistry or
   scattering calculation.
2. A PES file containing the effective intermolecular potential energy surface as a
   function of intermolecular distance R and angle theta.

Its purpose is to compare the spatial probability density of the bound state with the
underlying potential energy landscape. The output is a publication-style two-panel
figure showing the PES contour map and the corresponding ground-state probability
distribution.

Why the script exists
---------------------
For a weakly bound van der Waals complex such as He--CO2, the wavefunction describes
where the complex is most likely to be found. The PES tells us the energetic landscape
that traps the complex. By combining the two, the script helps answer questions like:

- Where is the molecule most likely to be found in the angular coordinate theta?
- What is the average radial separation <R> of the complex?
- Where is the PES minimum located in geometry space?
- How does the most probable quantum distribution compare to the equilibrium geometry
  of the interaction potential?

The script is therefore a practical analysis and plotting utility for molecular
spectroscopy and intermolecular dynamics work.

Wavefunction parsing
--------------------
The function parse_bound_state() reads a text file containing a bound-state solution.
Quantum codes often print structured wavefunction data with a header that includes the
state number, energy, basis size, and number of radial points.

The parser looks for a pattern like:

    # WAVEFUNCTION FOR STATE 1 AT ENERGY ...

This identifies the state and extracts its energy. It then finds the section of the file
that reports the total number of radial points, as well as metadata telling it how many
basis functions are used in the angular expansion.

The output is a NumPy array of shape (n_points, n_channels + 1):

    [R, channel_1, channel_2, ..., channel_n]

where R is the radial grid and the remaining values are amplitudes for each angular
channel. These amplitudes are the fundamental representation of the bound-state wavefunction.

PES loading
-----------
The load_pes() function reads the PES file with np.loadtxt(). It assumes each row has
three values:

    R, theta, V

The code extracts the unique grid values in R and theta, validates that the grid is
complete, and reshapes the interaction energy values into a 2D array of shape
(len(R), len(theta)).

The energy values are converted from Hartree to wavenumbers using:

    HARTREE_TO_WAVENUMBER = 219474.6313705

This conversion matches the unit conventions of spectroscopic and quantum-chemistry
analysis where energies are often reported in cm^-1.

Reconstructing the probability density
--------------------------------------
The function reconstruct_probability() converts the radial amplitudes into a probability
in the continuum of R and theta. The angular dependence is expanded in Legendre
polynomials:

    P_j(cos(theta))

with a normalization factor:

    sqrt((2j + 1) / 2)

This creates an angular basis. The radial amplitudes are projected onto that basis via a
matrix multiplication, producing a reduced wavefunction:

    reduced_wavefunction = radial_amplitudes @ angular_basis

The 2D probability density per element dR dtheta is then:

    probability = reduced_wavefunction^2 * sin(theta)

The factor sin(theta) is the Jacobian associated with spherical/angular coordinates. This
ensures that the density is properly weighted when integrating over theta.

Normalization and observable quantities
--------------------------------------
Inside main(), the script calculates the radial probability density by summing the squared
amplitudes across all angular channels:

    radial_density = sum(amplitudes**2, axis=1)

Using a trapezoidal integration over the radial grid, it normalizes the 1D radial density
and computes the average radius:

    <R> = trapz(R * radial_density, R) / trapz(radial_density, R)

It also computes the radial variance and standard deviation:

    sigma_R = sqrt(<R^2> - <R>^2)

These values describe how spread out the ground-state density is along the radial coordinate.

The script also evaluates the full 2D probability density on a finely sampled theta grid,
normalizes it over both coordinates, and finds the location (R, theta) of maximum probability.
This identifies the most likely geometry of the bound state in the intermolecular coordinate system.

PES minimum and geometry analysis
--------------------------------
The script interpolates the PES using RectBivariateSpline(), a smooth 2D cubic spline,
and then minimizes the energy numerically with scipy.optimize.minimize.

This yields the equilibrium geometry of the PES:

    (R_e, theta_e, V_min)

which corresponds to the minimum energy geometry of the interaction potential. The location
is marked on the PES plot as a black star and annotated with the equilibrium distance and
angle.

Plot generation
---------------
The script creates a two-panel Matplotlib figure.

Left panel (a): PES contour plot
- R is plotted on the x-axis.
- theta is plotted on the y-axis, converted from radians to degrees.
- The PES is shown as filled contours using a red-blue diverging colormap.
- The zero-energy contour is drawn as a dashed black line.
- The PES minimum is marked with a star and labeled with R_e and theta_e.

Right panel (b): J = 0 probability density
- The x-axis is the radial coordinate R.
- The y-axis is theta in degrees.
- The reconstructed probability density is drawn as filled contours.
- Selected contour levels (0.1, 0.5, 0.9 relative probability) are overlaid.
- The mean radius <R> is shown as a red dashed vertical line.
- The absolute probability maximum is marked with a circle.

The plotting style is tuned for a clean publication-quality figure, with carefully chosen
axis limits, fonts, tick marks, and colorbars.

Output files
------------
The script saves the figure to three formats:

- PDF
- SVG
- PNG

It writes the outputs with a user-specified stem via the --output flag. These files are
typically used directly in manuscripts, presentations, or analysis reports.

Command-line usage
------------------
The script is invoked as:

    python scripts/plot_heco2_pes_probability.py <wavefunction_file> <pes_file> \
        --state 1 --output heco2_pes_probability

where:
- wavefunction_file is the formatted BOUND wavefunction output,
- pes_file is the R-theta effective PES data,
- --state selects the bound state to analyze,
- --output defines the output file prefix.

At the end, the script prints summary information to the terminal, including:
- bound-state energy,
- wavefunction norm,
- average radius <R>,
- radial standard deviation sigma_R,
- maximum probability location,
- PES minimum geometry and well depth.

Overall significance
--------------------
This script is a bridge between quantum structure calculations and geometric visualization.
It turns a bound-state wavefunction into a physically interpretable probability map and then
compares that map with the effective interaction potential. The result is a compact, visual
representation of where the He--CO2 complex is most likely to reside and how that distribution
relates to the energetically preferred geometry.

In short, the program reads quantum-chemistry output, reconstructs the probability density,
locates the important geometric features of both the wavefunction and the PES, and exports a
publication-ready figure summarizing the physics of the bound complex.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize
from scipy.special import eval_legendre


HARTREE_TO_WAVENUMBER = 219474.6313705


def parse_bound_state(path: Path, state: int = 1) -> tuple[float, np.ndarray, np.ndarray]:
    """Return energy, radial grid, and channel amplitudes for one BOUND state."""
    lines = path.read_text().splitlines()
    marker = re.compile(
        rf"^# WAVEFUNCTION FOR STATE\s+{state}\s+AT ENERGY\s+([-+0-9.EeDd]+)"
    )
    start = None
    energy = None
    for i, line in enumerate(lines):
        match = marker.match(line)
        if match:
            start = i
            energy = float(match.group(1).replace("D", "E"))
            break
    if start is None or energy is None:
        raise ValueError(f"State {state} was not found in {path}")

    npoints_line = next(
        i for i in range(start, min(start + 12, len(lines)))
        if "TOTAL NUMBER OF POINTS" in lines[i]
    )
    npoints = int(lines[npoints_line].split()[-1])

    header = "\n".join(lines[:start])
    basis_match = re.search(r"JTOT, IBLOCK, NBASIS:\s+\d+\s+\d+\s+(\d+)", header)
    row_match = re.search(
        r"COMPONENTS OF WAVEFUNCTION VECTOR OCCUPY\s+(\d+)\s+LINES", header
    )
    if basis_match is None or row_match is None:
        raise ValueError("Could not determine the BOUND basis or record layout")
    nbasis = int(basis_match.group(1))
    lines_per_point = int(row_match.group(1))

    rows = []
    cursor = npoints_line + 1
    for _ in range(npoints):
        values: list[float] = []
        for _ in range(lines_per_point):
            values.extend(
                float(token.replace("D", "E")) for token in lines[cursor].split()
            )
            cursor += 1
        if len(values) != nbasis + 1:
            raise ValueError(
                f"Expected {nbasis + 1} values per radial point; found {len(values)}"
            )
        rows.append(values)

    data = np.asarray(rows)
    return energy, data[:, 0], data[:, 1:]


def load_pes(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return R, theta, and the effective PES in wavenumbers."""
    raw = np.loadtxt(path)
    r_grid = np.unique(raw[:, 0])
    theta_grid = np.unique(raw[:, 1])
    expected = r_grid.size * theta_grid.size
    if raw.shape[0] != expected:
        raise ValueError(f"PES grid is incomplete: expected {expected} rows")
    values = raw[:, 2].reshape(r_grid.size, theta_grid.size)
    return r_grid, theta_grid, values * HARTREE_TO_WAVENUMBER


def reconstruct_probability(
    radial_amplitudes: np.ndarray,
    theta: np.ndarray,
) -> np.ndarray:
    """Return probability per dR dtheta from the J=0 channel amplitudes."""
    x = np.cos(theta)
    angular_basis = np.asarray(
        [
            np.sqrt((2 * j + 1) / 2) * eval_legendre(j, x)
            for j in range(radial_amplitudes.shape[1])
        ]
    )
    reduced_wavefunction = radial_amplitudes @ angular_basis
    return reduced_wavefunction**2 * np.sin(theta)[None, :]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("wavefunction", type=Path, help="Formatted BOUND wavefunction file")
    parser.add_argument("pes", type=Path, help="R, theta, V effective-PES grid")
    parser.add_argument("--state", type=int, default=1)
    parser.add_argument("--output", type=Path, default=Path("heco2_pes_probability"))
    args = parser.parse_args()

    energy, r_wave, amplitudes = parse_bound_state(args.wavefunction, args.state)
    r_pes, theta_pes, potential = load_pes(args.pes)

    radial_density = np.sum(amplitudes**2, axis=1)
    norm = np.trapezoid(radial_density, r_wave)
    mean_r = np.trapezoid(r_wave * radial_density, r_wave) / norm
    mean_r2 = np.trapezoid(r_wave**2 * radial_density, r_wave) / norm
    sigma_r = np.sqrt(mean_r2 - mean_r**2)

    theta_plot = np.linspace(0.0, np.pi, 721)
    probability = reconstruct_probability(amplitudes, theta_plot)
    probability /= np.trapezoid(
        np.trapezoid(probability, theta_plot, axis=1), r_wave
    )
    probability_relative = probability / probability.max()
    max_index = np.unravel_index(np.argmax(probability), probability.shape)
    r_probability_max = r_wave[max_index[0]]
    theta_probability_max = np.degrees(theta_plot[max_index[1]])

    spline = RectBivariateSpline(r_pes, theta_pes, potential, kx=3, ky=3)
    minimum = minimize(
        lambda x: float(spline(x[0], x[1])[0, 0]),
        x0=(3.1, np.pi / 2),
        bounds=((r_pes.min(), 6.0), (theta_pes.min(), theta_pes.max())),
    )
    r_e, theta_e = minimum.x
    v_min = minimum.fun

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.35), sharey=True)
    r_limits = (2.45, 6.25)
    theta_degrees_pes = np.degrees(theta_pes)
    theta_degrees_probability = np.degrees(theta_plot)

    pes_levels = np.asarray(
        [-50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 20, 40]
    )
    pes_norm = BoundaryNorm(pes_levels, plt.get_cmap("RdYlBu_r").N)
    pes_contours = axes[0].contourf(
        r_pes,
        theta_degrees_pes,
        potential.T,
        levels=pes_levels,
        norm=pes_norm,
        cmap="RdYlBu_r",
        extend="both",
    )
    axes[0].contour(
        r_pes,
        theta_degrees_pes,
        potential.T,
        levels=[0.0],
        colors="black",
        linewidths=1.1,
        linestyles="--",
    )
    axes[0].plot(
        r_e,
        np.degrees(theta_e),
        marker="*",
        markersize=9,
        color="black",
        markeredgecolor="white",
        markeredgewidth=0.7,
        zorder=5,
    )
    axes[0].annotate(
        rf"$R_e={r_e:.3f}$ Å" + "\n" + rf"$\theta_e={np.degrees(theta_e):.0f}^\circ$",
        xy=(r_e, np.degrees(theta_e)),
        xytext=(3.72, 61),
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "-", "color": "black", "lw": 0.7},
    )
    axes[0].set_title(r"(a) Effective intermolecular PES ($v_3=0$)")
    axes[0].set_xlabel(r"$R$ (Å)")
    axes[0].set_ylabel(r"$\theta$ (degrees)")
    cbar0 = fig.colorbar(pes_contours, ax=axes[0], orientation="horizontal", pad=0.17)
    cbar0.set_label(r"Interaction energy (cm$^{-1}$)")
    cbar0.set_ticks([-40, -20, 0, 20, 40])

    probability_levels = np.linspace(0.0, 1.0, 21)
    probability_contours = axes[1].contourf(
        r_wave,
        theta_degrees_probability,
        probability_relative.T,
        levels=probability_levels,
        cmap="Blues",
    )
    axes[1].contour(
        r_wave,
        theta_degrees_probability,
        probability_relative.T,
        levels=[0.1, 0.5, 0.9],
        colors="#333333",
        linewidths=0.55,
    )
    axes[1].axvline(mean_r, color="#7A1F1F", linestyle=(0, (5, 2)), linewidth=1.3)
    axes[1].plot(
        r_probability_max,
        theta_probability_max,
        marker="o",
        markersize=5,
        markerfacecolor="white",
        markeredgecolor="black",
        markeredgewidth=0.8,
        zorder=5,
    )
    axes[1].text(
        mean_r + 0.07,
        166,
        rf"$\langle R\rangle={mean_r:.3f}$ Å",
        ha="left",
        va="top",
        color="#7A1F1F",
    )
    axes[1].set_title(r"(b) $J=0$ ground-state probability")
    axes[1].set_xlabel(r"$R$ (Å)")
    cbar1 = fig.colorbar(
        probability_contours, ax=axes[1], orientation="horizontal", pad=0.17
    )
    cbar1.set_label("Relative probability density")
    cbar1.set_ticks([0.0, 0.5, 1.0])

    for ax in axes:
        ax.set_xlim(*r_limits)
        ax.set_ylim(0, 180)
        ax.set_xticks([2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
        ax.set_yticks([0, 45, 90, 135, 180])
        ax.tick_params(direction="out", length=3, width=0.8)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)

    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.18, top=0.92, wspace=0.13)

    stem = args.output
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.04)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.04)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    print(f"State {args.state} energy: {energy:.11f} cm^-1")
    print(f"Wavefunction norm: {norm:.10f}")
    print(f"<R>: {mean_r:.6f} A")
    print(f"sigma_R: {sigma_r:.6f} A")
    print(
        f"Probability maximum: R={r_probability_max:.6f} A, "
        f"theta={theta_probability_max:.2f} deg"
    )
    print(f"PES minimum: R={r_e:.6f} A, theta={np.degrees(theta_e):.3f} deg")
    print(f"PES well depth: {-v_min:.6f} cm^-1")


if __name__ == "__main__":
    main()
