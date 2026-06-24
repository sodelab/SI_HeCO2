# Supporting Information for He--CO2 Vibrational Structure Study

This repository contains supporting data for the manuscript:

**Anharmonic Vibrational Structure and Bound States of the He–CO$_2$ van der Waals Complex**

Larry Rodriguez, Michael Natalizio, John Oliver Lobo, Angel Vasquez, and Olaseni Sode

Department of Chemistry and Biochemistry, California State University, Los Angeles

## Overview

This repository provides the numerical data and supporting files associated with the development and analysis of the He--CO2 potential energy surface and the vibrational structure calculations discussed in the manuscript. The data include optimized geometries, potential-energy-surface information, harmonic/VSCF/VCI vibrational results, VCI wavefunction coefficients, effective intermolecular potentials, and bound-state calculation files.

The accompanying Supporting Information document summarizes the dominant VCI wavefunction contributions and the J = 0 bound-state energies most directly relevant to the assignments discussed in the manuscript. This repository contains the more complete numerical data used to generate those summaries.

## Directory Structure

```text
SI_HeCO2/
│
├── README.md
├── LICENSE
├── CITATION.cff
│
├── supporting_information/
│   ├── SI_HeCO2.pdf
│   └── SI_HeCO2.tex
│
├── geometries/
│   ├── cartesian/
│   └── internal_coordinates/
│
├── PES/
│   ├── coefficients/
│   ├── training_set/
│   ├── test_set/
│   └── fitting_statistics/
│
├── vibrational_calculations/
│   ├── harmonic/
│   ├── VSCF/
│   ├── VCI/
│   └── VCI_coefficients/
│
├── bound_state_calculations/
│   ├── effective_potentials/
│   ├── BOUND_inputs/
│   ├── BOUND_outputs/
│   └── bound_state_energies/
│
├── isotopologues/
│   ├── effective_potentials/
│   ├── BOUND_outputs/
│   └── band_origin_shifts/
│
└── scripts/
    ├── extract_vci_coefficients.py
    └── README.md
```

## Repository Contents

### `supporting_information/`

This directory contains the Supporting Information document associated with the manuscript. The PDF version is intended for direct reading, while the LaTeX source is provided for transparency and reproducibility.

### `geometries/`

This directory contains structural data for the He--CO2 complex. Cartesian coordinates are provided in XYZ format, and internal-coordinate files are included where relevant for reproducing the vibrational calculations.

The equilibrium structure corresponds to the T-shaped He--CO2 minimum discussed in the manuscript.

### `PES/`

This directory contains data associated with the flexible-monomer two-body He--CO2 potential energy surface.

The `coefficients/` directory contains the fitted potential-energy-surface coefficients. The `training_set/` and `test_set/` directories contain the electronic-structure data used to fit and validate the potential. The `fitting_statistics/` directory contains summary statistics such as RMSE and maximum errors for the full and low-energy regions of the dataset.

Large data files may be provided in compressed format.

### `vibrational_calculations/`

This directory contains files associated with the harmonic, VSCF, and VCI vibrational calculations.

The `harmonic/` directory contains harmonic frequencies and normal-mode information. The `VSCF/` directory contains VSCF frequency results and relevant output files. The `VCI/` directory contains VCI frequency results, including the states discussed in the manuscript. The `VCI_coefficients/` directory contains the VCI wavefunction expansion coefficients used to assign the dominant vibrational configurations reported in the Supporting Information.

VCI configurations are reported in the mode order:

```text
(nu_b, nu_s, nu_2^o, nu_2^i, nu_1, nu_3)
```

where `nu_b` and `nu_s` are the intermolecular bend and stretch, `nu_2^o` and `nu_2^i` are the out-of-plane and in-plane CO2 bends, `nu_1` is the symmetric stretch, and `nu_3` is the asymmetric stretch.

### `bound_state_calculations/`

This directory contains files associated with the intermolecular bound-state calculations.

The `effective_potentials/` directory contains the vibrational-state-specific effective intermolecular potentials, including the ground-state and nu3-excited surfaces. The `BOUND_inputs/` directory contains representative input files for the BOUND calculations. The `BOUND_outputs/` directory contains raw BOUND output files. The `bound_state_energies/` directory contains processed J = 0 bound-state energy tables.

Energies are reported relative to the appropriate dissociation limit. Excitation energies are reported relative to the lowest bound state on each effective potential.

### `isotopologues/`

This directory contains bound-state and band-origin-shift data for the isotopologues discussed in the manuscript. These files support the isotope-dependent nu3 band-origin shifts reported in the main text.

### `scripts/`

This directory contains analysis scripts used to process output files and generate tables for the Supporting Information. In particular, `extract_vci_coefficients.py` extracts the leading VCI wavefunction contributions from VCI output files and formats them for inspection or direct inclusion in LaTeX tables.

## Notes on Reproducibility

The files in this repository are intended to allow readers to inspect the numerical data supporting the vibrational assignments, potential-energy-surface fitting, and bound-state calculations reported in the manuscript. Some calculations require external software packages, including NITROGEN for the VSCF/VCI calculations and BOUND for the intermolecular bound-state calculations.

Where possible, processed summary files are provided alongside raw output files to make the data easier to inspect without rerunning the full calculations.

## Citation

If you use this data, please cite the associated manuscript:

Rodriguez, L.; Natalizio, M.; Lobo, J. O.; Vasquez, A.; Sode, O. *Anharmonic Vibrational Structure and Bound States of the He–CO$_2$ van der Waals Complex.*

A formal citation and DOI will be added after publication.
