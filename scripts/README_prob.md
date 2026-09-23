# He--CO2 effective PES and probability-density figure

`plot_heco2_pes_probability.py` creates a two-panel figure showing:

1. the effective He--CO2 intermolecular potential energy surface (PES) for
   CO2 in its `v3 = 0` state; and
2. the probability density reconstructed from a formatted BOUND
   wavefunction file for a selected `J = 0` bound state.

The script writes publication-ready PDF and SVG files and a 600 dpi PNG
preview.

## Requirements

- Python 3.9 or newer
- NumPy 2.0 or newer
- SciPy
- Matplotlib

Install the required Python packages, if necessary, with:

```bash
python3 -m pip install "numpy>=2.0" scipy matplotlib
```

## Input files

### Formatted BOUND wavefunction

The first positional argument is a formatted BOUND wavefunction file, such as
`fort.10`. The calculation used for the manuscript requested formatted
wavefunctions on unit 10 with settings equivalent to:

```text
IWAVE  = 10
IWAVEF = .TRUE.
IWVSTP = 1
```

The parser reads the BOUND header to determine the basis size and record
layout. It then extracts the requested state number, radial grid, energy, and
channel amplitudes.

### Effective PES grid

The second positional argument is a three-column text file such as
`Veff0_grid.dat`:

```text
R (angstrom)    theta (radian)    V (hartree)
```

The grid must be rectangular and ordered in blocks of constant `R`, with
`theta` varying within each block. Potential energies are converted internally
from hartree to inverse centimeters using
`1 hartree = 219474.6313705 cm^-1`.

## Usage

From the directory containing the script and input files, run:

```bash
python3 plot_heco2_pes_probability.py fort.10 Veff0_grid.dat
```

To choose a BOUND state and a different output stem:

```bash
python3 plot_heco2_pes_probability.py fort.10 Veff0_grid.dat \
    --state 1 --output figures/heco2_pes_probability
```

Command-line options:

- `--state N`: sequential state number in the formatted BOUND file; default
  is state 1.
- `--output PATH`: output path without an extension; default is
  `heco2_pes_probability`.

The manuscript figure uses state 1. The panel title identifies the plotted
state as the ground state, so the title should be edited if the script is used
to display an excited state.

## Outputs

For the default output stem, the script creates:

- `heco2_pes_probability.pdf`
- `heco2_pes_probability.svg`
- `heco2_pes_probability.png`

It also prints the selected state energy, wavefunction normalization,
`<R>`, the standard deviation of `R`, the probability-density maximum, the
interpolated PES minimum, and the PES well depth.

## Probability-density reconstruction

For `J = 0`, the angular channel functions are represented with normalized
Legendre polynomials,

```text
phi_j(theta) = sqrt[(2j + 1)/2] P_j(cos(theta)).
```

The plotted density is the probability per `dR dtheta`,

```text
rho(R,theta) = |sum_j F_j(R) phi_j(theta)|^2 sin(theta),
```

where `F_j(R)` are the radial channel amplitudes written by BOUND. The
`sin(theta)` factor is the angular-coordinate measure. The density is first
normalized by numerical integration and is then divided by its maximum for
the plotted color scale. The color bar therefore reports relative probability
density from 0 to 1.

## Reference results for the supplied He--CO2 data

Using state 1 from the supplied `fort.10` file and `Veff0_grid.dat` gives:

| Quantity | Value |
|---|---:|
| BOUND state energy | -16.09868560855 cm^-1 |
| Wavefunction norm | 0.9999999970 |
| `<R>` | 3.619363 angstrom |
| Standard deviation of `R` | 0.463356 angstrom |
| Probability-density maximum | `R = 3.353384` angstrom, `theta = 90.00` degrees |
| Effective-PES minimum | `R_e = 3.072879` angstrom, `theta_e = 90.000` degrees |
| Effective-PES well depth | 47.586141 cm^-1 |

The potential minimum and probability-density maximum are both T-shaped.
Because the two ends of homonuclear CO2 are equivalent, the density is
symmetric about 90 degrees and the average over the full 0--180 degree range
is necessarily 90 degrees. If the angle is folded onto the unique 0--90
degree range using
`theta_fold = min(theta, 180 degrees - theta)`, the state-1 distribution gives
`<theta_fold> = 78.26 degrees`. This folded expectation value is distinct from
an effective angle obtained by fitting experimental rotational constants.

## Notes

- The PES minimum is obtained from a bicubic spline interpolation followed by
  numerical minimization.
- The dashed contour in the PES panel is the zero-interaction-energy contour.
- The open circle in the probability panel marks the global density maximum.
- The vertical dashed line marks `<R>`, not the radial-density mode and not an
  experimentally fitted `R_0` value.
- The state energy printed by this wavefunction-generating run can differ
  slightly from an energy obtained with another BOUND propagator. The figure
  depends on the wavefunction stored in the supplied file.

