#!/usr/bin/env python3
"""
co2_iso_constants.py

Compute:
  1) m(CO2) and reduced mass mu for He–CO2 (amu)
  2) CO2 rotational constant B (cm^-1) at r_eq (Å)

Usage:
  python co2_iso_constants.py 626
  python co2_iso_constants.py 638 --r_eq 1.16221264055
  python co2_iso_constants.py 828 --he 4.002603
"""

from __future__ import annotations
import argparse
import numpy as np

# ---- Physical constants (exact where defined) ----
AMU_KG = 1.66053906660e-27
H = 6.62607015e-34
C = 299792458.0  # m/s
ANGSTROM_M = 1e-10

# ---- Minimal isotopologue mass database (amu) ----
# HITRAN-style 3-digit code: O–C–O for CO2
# 6 = 16O, 7 = 17O, 8 = 18O ; 2 = 12C, 3 = 13C
MASS = {
    "12C": 12.0,
    "13C": 13.003355,
    "16O": 15.994915,
    "18O": 17.999160,
}

ISOTOPOLOGUES = {
    "626": dict(name="12C16O2", mass_C=MASS["12C"], mass_O1=MASS["16O"], mass_O2=MASS["16O"]),
    "636": dict(name="13C16O2", mass_C=MASS["13C"], mass_O1=MASS["16O"], mass_O2=MASS["16O"]),
    "628": dict(name="12C18O2", mass_C=MASS["12C"], mass_O1=MASS["18O"], mass_O2=MASS["18O"]),
    "828": dict(name="12C18O2 (828)", mass_C=MASS["12C"], mass_O1=MASS["18O"], mass_O2=MASS["18O"]),
    "838": dict(name="13C18O2 (838)", mass_C=MASS["13C"], mass_O1=MASS["18O"], mass_O2=MASS["18O"]),
    "638": dict(name="16O13C18O (638)", mass_C=MASS["13C"], mass_O1=MASS["16O"], mass_O2=MASS["18O"]),
}

def he_co2_reduced_mass_amu(m_he: float, m_c: float, m_o1: float, m_o2: float):
    m_co2 = m_c + m_o1 + m_o2
    mu = (m_he * m_co2) / (m_he + m_co2)
    return mu, m_co2

def co2_rotational_constant_cm1(m_c: float, m_o1: float, m_o2: float, r_eq_ang: float) -> float:
    """
    Linear CO2 along z with carbon at 0 and oxygens at ±r_eq.
    Compute moment of inertia about x (or y) through COM, then B = h/(8π^2 I c).
    Returns B in cm^-1.
    """
    # masses in kg
    mC = m_c * AMU_KG
    mO1 = m_o1 * AMU_KG
    mO2 = m_o2 * AMU_KG
    r = r_eq_ang * ANGSTROM_M

    M = mC + mO1 + mO2

    # COM shift for asymmetric O masses: z_com = r*(mO1 - mO2)/M
    z_com = r * (mO1 - mO2) / M

    # positions relative to COM
    zC  = -z_com
    zO1 =  r - z_com
    zO2 = -r - z_com

    # moment of inertia about x (or y)
    I = mC * zC**2 + mO1 * zO1**2 + mO2 * zO2**2  # kg m^2

    # B in 1/m
    B_m_inv = H / (8.0 * np.pi**2 * I * C)
    # convert 1/m -> 1/cm
    return B_m_inv / 100.0

def main():
    ap = argparse.ArgumentParser(description="Compute He–CO2 reduced mass and CO2 rotational constant for a CO2 isotopologue code (e.g. 626).")
    ap.add_argument("iso", help="Isotopologue code (CO2 HITRAN-style O–C–O), e.g. 626, 636, 828, 638")
    ap.add_argument("--r_eq", type=float, default=1.16221264055, help="CO bond length r_eq in Å (default from your PES)")
    ap.add_argument("--he", type=float, default=4.002603, help="He mass in amu (default 4He)")
    args = ap.parse_args()

    if args.iso not in ISOTOPOLOGUES:
        known = ", ".join(sorted(ISOTOPOLOGUES.keys()))
        raise SystemExit(f"Unknown isotopologue '{args.iso}'. Known: {known}")

    iso = ISOTOPOLOGUES[args.iso]
    mC = iso["mass_C"]
    mO1 = iso["mass_O1"]
    mO2 = iso["mass_O2"]
    mHe = args.he

    mu, mco2 = he_co2_reduced_mass_amu(mHe, mC, mO1, mO2)
    B = co2_rotational_constant_cm1(mC, mO1, mO2, args.r_eq)

    print(f"Isotopologue: {args.iso}  ({iso['name']})")
    print(f"r_eq (Å): {args.r_eq:.12f}")
    print("")
    print("Masses (amu):")
    print(f"  m(He)  = {mHe:.9f}")
    print(f"  m(C)   = {mC:.9f}")
    print(f"  m(O1)  = {mO1:.9f}")
    print(f"  m(O2)  = {mO2:.9f}")
    print(f"  m(CO2) = {mco2:.9f}")
    print("")
    print("Derived:")
    print(f"  mu(He–CO2) (amu) = {mu:.12f}")
    print(f"  B(CO2) (cm^-1)   = {B:.12f}")

if __name__ == "__main__":
    main()
