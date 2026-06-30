# Supporting Information for He--CO2 Vibrational Structure Study

This repository contains supporting data for the manuscript:

**Anharmonic Vibrational Structure and Bound States of the He–CO$_2$ van der Waals Complex**

Larry Rodriguez, Michael Natalizio, John Oliver Lobo, Angel Vasquez, and Olaseni Sode

Department of Chemistry and Biochemistry, California State University, Los Angeles

## Overview

This repository provides the numerical data and supporting files associated with the development and analysis of the He--CO2 potential energy surface and the vibrational structure calculations discussed in the manuscript. The data include optimized geometries, potential-energy-surface information, harmonic/VSCF/VCI vibrational results, VCI wavefunction coefficients, effective intermolecular potentials, and bound-state calculation files.

The accompanying Supporting Information document summarizes the dominant VCI wavefunction contributions and the J = 0 bound-state energies most directly relevant to the assignments discussed in the manuscript. This repository contains the more complete numerical data used to generate those summaries.

## Directory Structure and Repository Contents

```text
SI_HeCO2/
│
├── README.md
├── LICENSE
├── CITATION.cff
│
├── supporting_information/
├── geometries/
├── PES/
├── vibrational_calculations/
├── bound_state_calculations/
└── scripts/
```

### `supporting_information/`

Contains the Supporting Information PDF and LaTeX source associated with the manuscript.

### `geometries/`

Contains optimized Cartesian coordinates provided in XYZ format and internal-coordinate representations of the He--CO<sub>2</sub> complex. The equilibrium structure corresponds to the T-shaped He--CO<sub>2</sub> minimum discussed in the manuscript.

### `PES/`

Contains the fitted He--CO<sub>2</sub> potential-energy-surface coefficients, fitting statistics, and associated training/test-set data. 

### `vibrational_calculations/`

Contains harmonic, VSCF, and VCI calculations for all isotopologues considered in the manuscript and the VCI wavefunction expansion coefficients. VCI configurations are reported in the mode order:

```text
(nu_b, nu_s, nu_2^o, nu_2^i, nu_1, nu_3)
```

where `nu_b` and `nu_s` are the intermolecular bend and stretch, `nu_2^o` and `nu_2^i` are the out-of-plane and in-plane CO2 bends, `nu_1` is the symmetric stretch, and `nu_3` is the asymmetric stretch.

```text
vibrational_calculations/
│
├── He-16O12C16O/
│   ├── harmonic/
│   ├── VSCF/
│   └── VCI/
│
├── He-16O13C16O/
├── He-18O12C18O/
├── He-18O13C18O/
├── He-16O12C18O/
└── He-16O13C18O/
```

Each isotopologue directory follows the same internal organization as `He-16O12C16O`.

### `bound_state_calculations/`

Contains effective intermolecular potentials, BOUND input files, output files, and processed bound-state energies. 

```text
bound_state_calculations/
│
├── He-16O12C16O/
│   ├── Veff_v0/
│   ├── Veff_v1/
│   ├── BOUND_inputs/
│   └── outputs/
│
├── He-16O13C16O/
├── He-18O12C18O/
├── He-18O13C18O/
├── He-16O12C18O/
└── He-16O13C18O/
```

Each isotopologue directory follows the same internal organization as `He-16O12C16O`.

### `scripts/`

Contains utility scripts used for data analysis and generation of Supporting Information tables, including extraction of dominant VCI wavefunction contributions. 


## Notes on Reproducibility

The files in this repository are intended to allow readers to inspect the numerical data supporting the vibrational assignments, potential-energy-surface fitting, and bound-state calculations reported in the manuscript. Some calculations require external software packages, including NITROGEN for the VSCF/VCI calculations and BOUND for the intermolecular bound-state calculations.

Where possible, processed summary files are provided alongside raw output files to make the data easier to inspect without rerunning the full calculations.

## Citation

If you use this data, please cite the associated manuscript:

Rodriguez, L.; Natalizio, M.; Lobo, J. O.; Vasquez, A.; Sode, O. *Anharmonic Vibrational Structure and Bound States of the He–CO<sub>2</sub> van der Waals Complex.*

A formal citation and DOI will be added after publication.
