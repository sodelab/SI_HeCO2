"""
units.py - Physical constants and unit conversions for molecular spectroscopy

Uses scipy.constants for CODATA-recommended values.
All conversions assume atomic units (hartree, bohr, m_e) as internal representation.
"""
from __future__ import annotations

import numpy as np
import scipy.constants as sc

# ============================================================================
# Fundamental Constants from scipy.constants
# ============================================================================

# Planck constant (J·s)
HBAR = sc.hbar  # 1.054571817e-34 J·s

# Speed of light (m/s)
SPEED_OF_LIGHT = sc.c  # 299792458 m/s

# Electron mass (kg)
ELECTRON_MASS = sc.m_e  # 9.1093837015e-31 kg

# Atomic mass unit (kg)
AMU = sc.physical_constants['atomic mass constant'][0]  # 1.66053906660e-27 kg

# Avogadro constant (mol^-1)
AVOGADRO = sc.Avogadro  # 6.02214076e23 mol^-1

# Boltzmann constant (J/K)
BOLTZMANN = sc.k  # 1.380649e-23 J/K

# ============================================================================
# Derived Conversion Factors (from scipy.constants)
# ============================================================================

# Hartree energy to various units
EH_TO_J = sc.physical_constants['Hartree energy'][0]  # 4.3597447222071e-18 J
EH_TO_EV = sc.physical_constants['Hartree energy in eV'][0]  # 27.211386245988 eV

# Hartree to wavenumbers (cm^-1)
# E[cm^-1] = E[J] / (h·c) where c is in cm/s
EH_TO_CM1 = EH_TO_J / (sc.h * sc.c * 100)  # ≈ 219474.63 cm^-1
CM1_TO_EH = 1.0 / EH_TO_CM1

# Hartree to kcal/mol
# 1 hartree = EH_TO_J * N_A / 4184 kcal/mol
EH_TO_KCALMOL = EH_TO_J * AVOGADRO / 4184.0  # ≈ 627.509 kcal/mol
KCALMOL_TO_EH = 1.0 / EH_TO_KCALMOL

# kcal/mol to wavenumbers (cm^-1)
# Via hartree: kcal/mol → hartree → cm^-1
KCALMOL_TO_CM1 = KCALMOL_TO_EH * EH_TO_CM1  # ≈ 349.755 cm^-1
CM1_TO_KCALMOL = 1.0 / KCALMOL_TO_CM1

# Atomic mass unit to electron mass
AMU_TO_ME = AMU / ELECTRON_MASS  # ≈ 1822.888 m_e

# Bohr radius (m)
BOHR = sc.physical_constants['Bohr radius'][0]  # 5.29177210903e-11 m

# Angstrom to bohr
ANGSTROM_TO_BOHR = 1e-10 / BOHR  # ≈ 1.88973
BOHR_TO_ANGSTROM = 1.0 / ANGSTROM_TO_BOHR

# ============================================================================
# Conversion Functions
# ============================================================================

def hartree_to_cm1(x_h: np.ndarray | float) -> np.ndarray | float:
    """Convert hartree to wavenumbers (cm^-1)"""
    return x_h * EH_TO_CM1

def cm1_to_hartree(x_cm1: np.ndarray | float) -> np.ndarray | float:
    """Convert wavenumbers (cm^-1) to hartree"""
    return x_cm1 * CM1_TO_EH

def hartree_to_kcalmol(x_h: np.ndarray | float) -> np.ndarray | float:
    """Convert hartree to kcal/mol"""
    return x_h * EH_TO_KCALMOL

def kcalmol_to_hartree(x_kcal: np.ndarray | float) -> np.ndarray | float:
    """Convert kcal/mol to hartree"""
    return x_kcal * KCALMOL_TO_EH

def kcalmol_to_cm1(x_kcal: np.ndarray | float) -> np.ndarray | float:
    """Convert kcal/mol to wavenumbers (cm^-1)"""
    return x_kcal * KCALMOL_TO_CM1

def cm1_to_kcalmol(x_cm1: np.ndarray | float) -> np.ndarray | float:
    """Convert wavenumbers (cm^-1) to kcal/mol"""
    return x_cm1 * CM1_TO_KCALMOL

def hartree_to_ev(x_h: np.ndarray | float) -> np.ndarray | float:
    """Convert hartree to electron volts"""
    return x_h * EH_TO_EV

def ev_to_hartree(x_ev: np.ndarray | float) -> np.ndarray | float:
    """Convert electron volts to hartree"""
    return x_ev / EH_TO_EV

def amu_to_me(m_amu: float) -> float:
    """Convert atomic mass units to electron masses"""
    return m_amu * AMU_TO_ME

def me_to_amu(m_me: float) -> float:
    """Convert electron masses to atomic mass units"""
    return m_me / AMU_TO_ME

def angstrom_to_bohr(r_ang: np.ndarray | float) -> np.ndarray | float:
    """Convert Angstroms to bohr"""
    return r_ang * ANGSTROM_TO_BOHR

def bohr_to_angstrom(r_bohr: np.ndarray | float) -> np.ndarray | float:
    """Convert bohr to Angstroms"""
    return r_bohr * BOHR_TO_ANGSTROM

# ============================================================================
# Quick Reference
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Physical Constants and Unit Conversions")
    print("=" * 70)
    
    print("\nEnergy Conversions (1 hartree =):")
    print(f"  {EH_TO_CM1:.6f} cm⁻¹")
    print(f"  {EH_TO_KCALMOL:.6f} kcal/mol")
    print(f"  {EH_TO_EV:.6f} eV")
    print(f"  {EH_TO_J:.6e} J")
    
    print("\nEnergy Conversions (1 kcal/mol =):")
    print(f"  {KCALMOL_TO_CM1:.6f} cm⁻¹")
    print(f"  {KCALMOL_TO_EH:.6e} hartree")
    
    print("\nEnergy Conversions (1 cm⁻¹ =):")
    print(f"  {CM1_TO_KCALMOL:.6e} kcal/mol")
    print(f"  {CM1_TO_EH:.6e} hartree")
    
    print("\nLength Conversions:")
    print(f"  1 Å = {ANGSTROM_TO_BOHR:.6f} bohr")
    print(f"  1 bohr = {BOHR_TO_ANGSTROM:.6f} Å")
    
    print("\nMass Conversions:")
    print(f"  1 amu = {AMU_TO_ME:.6f} mₑ")
    
    print("\nExample: 10 kcal/mol barrier")
    barrier_kcal = 10.0
    print(f"  {barrier_kcal} kcal/mol = {kcalmol_to_cm1(barrier_kcal):.2f} cm⁻¹")
    print(f"  {barrier_kcal} kcal/mol = {kcalmol_to_hartree(barrier_kcal):.6f} hartree")
    
    print("\nExample: 100 cm⁻¹ vibrational frequency")
    freq_cm1 = 100.0
    print(f"  {freq_cm1} cm⁻¹ = {cm1_to_kcalmol(freq_cm1):.4f} kcal/mol")
    print(f"  {freq_cm1} cm⁻¹ = {cm1_to_hartree(freq_cm1):.6e} hartree")
    
    print("=" * 70)
