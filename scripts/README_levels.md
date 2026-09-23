# He--CO2 energy-level diagram

`plot_heco2_energy_levels.py` creates a publication-ready comparison of the
variational `J = 0` intermolecular bound-state excitations and selected VCI
local-mode excitation energies of He--CO2. It also displays the dissociation
threshold relative to the ground-state energy.

## Requirements

- Python 3.9 or newer
- Matplotlib

Install Matplotlib, if necessary, with:

```bash
python3 -m pip install matplotlib
```

## Usage

The numerical values are defined near the beginning of the script; no external
data file or command-line argument is required. Run:

```bash
python3 plot_heco2_energy_levels.py
```

The script creates:

- `heco2_energy_levels.pdf`
- `heco2_energy_levels.svg`
- `heco2_energy_levels.png`

The PNG file is written at 600 dpi. The PDF and SVG files are preferable for
manuscript preparation because they retain vector text and lines.

## Data represented in the figure

All plotted values are excitation energies in inverse centimeters relative to
the corresponding ground-state reference at zero.

### Variational bound-state excitations

| Assignment | Energy (cm^-1) |
|---|---:|
| `Delta E(0,1)` | 7.700 |
| `Delta E(0,2)` | 8.637 |
| `Delta E(0,3)` | 12.339 |
| `Delta E(0,4)` | 15.001 |

These levels are stored in `BOUND_LEVELS` and are drawn as solid blue lines.

### Selected VCI local-mode excitations

| Assignment | Energy (cm^-1) |
|---|---:|
| `nu_b` | 12.44 |
| `nu_s` | 18.79 |
| `2nu_b` | 26.95 |
| `nu_b + nu_s` | 32.05 |
| `2nu_s` | 41.19 |

These values are stored in `VCI_LEVELS` and are drawn as dashed brown lines.
The dashed style emphasizes that they are nominal local-mode excitations, not
additional variationally bound levels.

## Dissociation threshold

The dissociation energy is defined in the script as:

```python
D0 = 16.106
```

Because the vertical axis is energy relative to the bound ground state, the
dissociation limit appears at positive `D0 = 16.106 cm^-1`. Equivalently, the
ground-state energy is `-16.106 cm^-1` relative to separated He + CO2.

The gray region above the threshold identifies excitations above the bound
region. Thus, the VCI stretching fundamental and the displayed VCI overtones
and combination band lie above dissociation, whereas all four variational
levels shown remain below the threshold.

## Modifying the figure

Edit these objects near the top of the script:

- `D0` changes the dissociation threshold.
- `BOUND_LEVELS` changes the variational levels and labels.
- `VCI_LEVELS` changes the selected VCI levels and labels.

Each level is represented by a tuple containing its energy and a Matplotlib
math-text label, for example:

```python
(12.44, r"$\nu_b$")
```

The output stem is set near the end of the script:

```python
stem = Path("heco2_energy_levels")
```

The common horizontal line at zero is an excitation-energy reference for both
columns. It is not an additional VCI level. The dissociation threshold is
derived from the variational ground-state binding energy and is therefore
drawn using the same color as the bound-state ladder.

## Interpretation

The figure compares two different descriptions of the intermolecular motion.
The variational calculation produces a finite ladder of levels below
dissociation. In contrast, several oscillator-like VCI excitation energies lie
above the dissociation threshold and should not be interpreted as bound van
der Waals vibrational states. The diagram is intended to make this distinction
visually explicit; it does not imply a one-to-one assignment between every
variational level and every VCI local-mode excitation.

