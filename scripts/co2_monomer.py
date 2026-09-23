"""
co2_monomer.py - CO₂ internal coordinates and vibrational wavefunctions

Handles the ν₃ asymmetric stretch normal mode for vibrational averaging.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Tuple

from units import amu_to_me, cm1_to_hartree, angstrom_to_bohr

# ============================================================================
# CO₂ Isotopologue Data
# ============================================================================

@dataclass
class CO2Isotopologue:
    """Physical parameters for CO₂ isotopologues"""
    name: str
    mass_C: float  # amu
    mass_O1: float  # amu
    mass_O2: float  # amu
    r_eq: float    # Equilibrium C-O distance (Å)
    omega_3: float  # ν₃ harmonic frequency (cm⁻¹)
    x_33: float    # Anharmonicity constant (cm⁻¹)
    
    @property
    def mu_3_analytical(self) -> float:
        """
        Effective mass (m_e) for the q3 coordinate you actually use:
            q3 = (Δr1 - Δr2)/sqrt(2)
        with COM-fixed motion (C moves to avoid translation).
        """
        mC = amu_to_me(self.mass_C)
        m1 = amu_to_me(self.mass_O1)
        m2 = amu_to_me(self.mass_O2)

        s = (m1 + m2) / mC
        da = 1.0 / (np.sqrt(2.0) * (1.0 + s))
        db = -s * da

        mu_eff = (m1 + m2) * (da**2) + mC * (db**2)
        return mu_eff


    @property
    def mu_delta(self) -> float:
        """
        Effective mass for the antisymmetric stretch coordinate δ (atomic units),
        where r₁ = r_eq + δ and r₂ = r_eq - δ, with CO₂ center-of-mass fixed.
        
        For COM-fixed linear CO₂:
            μ_δ = m_C × (m_O₁ + m_O₂) / (m_C + m_O₁ + m_O₂)
        
        Returns
        -------
        float
            Reduced mass in atomic units (m_e)
        """
        mC  = amu_to_me(self.mass_C)
        mO1 = amu_to_me(self.mass_O1)
        mO2 = amu_to_me(self.mass_O2)
        return (mC * (mO1 + mO2)) / (mC + mO1 + mO2)

    @property
    def mu_3(self) -> float:
        """
        Effective mass for q₃ = (Δr₁ - Δr₂)/√2 (atomic units).
        
        The asymmetric stretch normal coordinate q₃ is related to δ by:
            q₃ = √2 × δ
        
        Since kinetic energy T = (1/2) μ (dq/dt)², the effective mass is:
            μ_q₃ = μ_δ
        
        (The √2 factor cancels between coordinate and kinetic energy definitions)
        
        Returns
        -------
        float
            Reduced mass in atomic units (m_e)
        """
        #return self.mu_delta
        return self.mu_3_analytical

    def mu_eff_from_cartesian_jacobian(self, dq: float = 1e-4, R_ang: float = 999.0, theta_rad: float = 0.0) -> float:
        """
        Compute effective mass from Cartesian metric tensor g = Σ_i m_i (∂r_i/∂q)².
        
        This numerically verifies the analytical mu_3 formula by computing the
        Jacobian of the q₃ → Cartesian transformation.
        
        Parameters
        ----------
        dq : float
            Finite difference step size (Å)
        R_ang : float
            He-CO₂ distance (Å), use large value for isolated CO₂
        theta_rad : float
            He-CO₂ angle (radians)
        
        Returns
        -------
        float
            Effective mass μ_eff = 1/g (atomic units, m_e)
        """
        import numpy as np
        from units import amu_to_me, angstrom_to_bohr
        
        # Masses in atomic units
        m = np.array([
            amu_to_me(4.002602),      # He (not used for CO₂-only, but kept for generality)
            amu_to_me(self.mass_C),
            amu_to_me(self.mass_O1),
            amu_to_me(self.mass_O2),
        ], dtype=float)

        # Compute Jacobian: ∂r/∂q₃ by finite difference
        c_plus  = q3_to_cartesian(+dq, R_ang, theta_rad, self)
        c_minus = q3_to_cartesian(-dq, R_ang, theta_rad, self)
        
        # Convert to bohr for proper metric calculation
        c_plus_bohr  = angstrom_to_bohr(c_plus)
        c_minus_bohr = angstrom_to_bohr(c_minus)
        dq_bohr = angstrom_to_bohr(dq)
        
        dr_dq = (c_plus_bohr - c_minus_bohr) / (2.0 * dq_bohr)  # (4,3) in bohr
        
        # Metric tensor: g = Σ_i m_i |∂r_i/∂q|²
        metric = np.sum(m * np.sum(dr_dq**2, axis=1))  # scalar
        
        # Effective mass: μ = 1/g
        #return 1.0 / metric
        return metric

    @property
    def energy_levels_cm1(self) -> np.ndarray:
        """
        Vibrational energy levels for ν₃ mode.
        
        G(v) = ω_e(v + 1/2) - x_e(v + 1/2)² (cm⁻¹)
        
        Returns
        -------
        np.ndarray
            Energy levels for v = 0 to 9 (cm⁻¹)
        """
        v = np.arange(10)
        return self.omega_3 * (v + 0.5) - self.x_33 * (v + 0.5)**2


# Standard isotopologues (HITRAN-style codes)
ISOTOPOLOGUES = {
    '626': CO2Isotopologue(
        name='¹²C¹⁶O₂',
        mass_C=12.0,
        mass_O1=15.994915,
        mass_O2=15.994915,
        r_eq=1.16221264055,  # Å (from your PES)
        omega_3=2349.143,     # cm⁻¹ (ν₃ fundamental)  (your value)
        x_33=12.494           # cm⁻¹ (anharmonicity)   (your value)
    ),
    '636': CO2Isotopologue(
        name='¹³C¹⁶O₂',
        mass_C=13.003355,
        mass_O1=15.994915,
        mass_O2=15.994915,
        r_eq=1.16221264055,
        omega_3=2283.489,
        x_33=12.184
    ),
    '628': CO2Isotopologue(
        name='¹²C¹⁸O₂',
        mass_C=12.0,
        mass_O1=15.994915,
        mass_O2=17.999160,
        r_eq=1.16221264055,
        omega_3=2272.355,
        x_33=11.845
    ),

    # ---- Additions ----
    '828': CO2Isotopologue(
        name='¹²C¹⁸O₂ (828)',
        mass_C=12.0,
        mass_O1=17.999160,
        mass_O2=17.999160,
        r_eq=1.16221264055,
        omega_3=2300,   # fill with your chosen reference
        x_33=11.0       # fill with your chosen reference
    ),
    '838': CO2Isotopologue(
        name='¹³C¹⁸O₂ (838)',
        mass_C=13.003355,
        mass_O1=17.999160,
        mass_O2=17.999160,
        r_eq=1.16221264055,
        omega_3=2300,   # fill
        x_33=11.0       # fill
    ),
    '638': CO2Isotopologue(
        name='¹⁶O¹³C¹⁸O (638)',
        mass_C=13.003355,
        mass_O1=15.994915,   # 16O
        mass_O2=17.999160,   # 18O
        r_eq=1.16221264055,
        omega_3=2300,   # fill
        x_33=11.0       # fill
    ),
}


# ============================================================================
# Normal Mode Coordinate Definition
# ============================================================================

def cartesian_to_q3(
    coords_ang: np.ndarray,
    iso: CO2Isotopologue
) -> float:
    """
    Convert Cartesian coordinates to ν₃ asymmetric stretch coordinate.
    
    Asymmetric stretch: q₃ = (Δr₁ - Δr₂) / √2
    where Δr₁ = r(C-O₁) - r_eq, Δr₂ = r(C-O₂) - r_eq
    
    Parameters
    ----------
    coords_ang : np.ndarray, shape (4, 3)
        [He, C, O₁, O₂] in Angstroms
    iso : CO2Isotopologue
        Isotopologue parameters
    
    Returns
    -------
    q3 : float
        Asymmetric stretch coordinate (Å·amu^½), mass-weighted
    """
    C = coords_ang[1]
    O1 = coords_ang[2]
    O2 = coords_ang[3]
    
    # Bond lengths
    r1 = np.linalg.norm(O1 - C)
    r2 = np.linalg.norm(O2 - C)
    
    # Displacements from equilibrium
    Delta_r1 = r1 - iso.r_eq
    Delta_r2 = r2 - iso.r_eq
    
    # Mass-weighted normal coordinate (asymmetric stretch)
    q3 = (Delta_r1 - Delta_r2) / np.sqrt(2)
    
    return q3

def q3_to_cartesian_co2_only(q3_ang: float, iso: CO2Isotopologue) -> np.ndarray:
    # same internal mapping as q3_to_cartesian, but no He row
    Delta_r1 = +q3_ang / np.sqrt(2)
    Delta_r2 = -q3_ang / np.sqrt(2)

    C  = np.array([0.0, 0.0, 0.0])
    O1 = np.array([0.0, 0.0, iso.r_eq + Delta_r1])
    O2 = np.array([0.0, 0.0, -(iso.r_eq + Delta_r2)])

    return np.vstack([C, O1, O2])   # shape (3,3)

# def q3_to_cartesian(
#     q3_ang: float,
#     R_ang: float,
#     theta_rad: float,
#     iso: CO2Isotopologue
# ) -> np.ndarray:
#     """
#     Generate Cartesian coordinates for given q₃ displacement.
    
#     CO₂ geometry:
#         - C at origin
#         - O atoms at ±(r_eq + Δr) along z-axis
#         - Asymmetric stretch: Δr₁ = -Δr₂ = q₃/√2
    
#     Parameters
#     ----------
#     q3_ang : float
#         Asymmetric stretch coordinate (Å)
#     R_ang : float
#         Jacobi R distance (Å)
#     theta_rad : float
#         Jacobi angle (radians)
#     iso : CO2Isotopologue
#         Isotopologue parameters
    
#     Returns
#     -------
#     coords : np.ndarray, shape (4, 3)
#         [He, C, O₁, O₂] in Angstroms
#     """
#     # Asymmetric stretch displacements
#     Delta_r1 = +q3_ang / np.sqrt(2)
#     Delta_r2 = -q3_ang / np.sqrt(2)
    
#     # CO₂ geometry (C at origin, axis along z)
#     C = np.array([0.0, 0.0, 0.0])
#     O1 = np.array([0.0, 0.0, iso.r_eq + Delta_r1])
#     O2 = np.array([0.0, 0.0, -(iso.r_eq + Delta_r2)])
    
#     # He position
#     He = np.array([
#         R_ang * np.sin(theta_rad),
#         0.0,
#         R_ang * np.cos(theta_rad)
#     ])
    
#     return np.vstack([He, C, O1, O2])

def q3_to_cartesian(q3_ang: float, R_ang: float, theta_rad: float, iso) -> np.ndarray:
    import numpy as np

    # bond lengths for asymmetric stretch coordinate q3 = (Δr1 - Δr2)/sqrt(2)
    r1 = iso.r_eq + q3_ang / np.sqrt(2.0)
    r2 = iso.r_eq - q3_ang / np.sqrt(2.0)

    mC  = iso.mass_C
    mO1 = iso.mass_O1
    mO2 = iso.mass_O2
    M = mC + mO1 + mO2

    # Place atoms on z-axis with correct bond lengths, then shift so COM is at 0
    # Coordinates before COM shift:
    #   C  at zC
    #   O1 at zC + r1
    #   O2 at zC - r2
    # Enforce COM=0:
    zC = (mO2 * r2 - mO1 * r1) / M

    C  = np.array([0.0, 0.0, zC])
    O1 = np.array([0.0, 0.0, zC + r1])
    O2 = np.array([0.0, 0.0, zC - r2])

    He = np.array([R_ang * np.sin(theta_rad), 0.0, R_ang * np.cos(theta_rad)])
    return np.vstack([He, C, O1, O2])

# ============================================================================
# 1D Harmonic Oscillator Basis (for vibrational averaging)
# ============================================================================

def harmonic_oscillator_wavefunction(
    q: np.ndarray,
    v: int,
    omega_cm1: float,
    mu_au: float
) -> np.ndarray:
    """
    1D harmonic oscillator wavefunction in mass-weighted coordinates.
    
    ψ_v(q) = N_v × H_v(α·q) × exp(-β·q²/2)
    
    where β = μω in atomic units (ℏ=1), α = √β
    
    Parameters
    ----------
    q : np.ndarray
        Coordinate grid (Å)
    v : int
        Vibrational quantum number
    omega_cm1 : float
        Harmonic frequency (cm⁻¹)
    mu_au : float
        Reduced mass (atomic units, m_e)
    
    Returns
    -------
    psi : np.ndarray
        Wavefunction values, normalized such that ∫|ψ|² dq = 1 in Å
    """
    from scipy.special import hermite, factorial
    from units import ANGSTROM_TO_BOHR
    
    # Convert to atomic units
    omega_au = cm1_to_hartree(omega_cm1)
    q_bohr = angstrom_to_bohr(q)
    
    # Harmonic oscillator parameters
    beta = mu_au * omega_au  # β = μω
    alpha = np.sqrt(beta)     # α = √β
    
    # Hermite polynomial H_v(α·q)
    H_v = hermite(v)
    
    # Normalization constant: N_v = (β/π)^(1/4) / √(2^v × v!)
    N_v = (beta / np.pi)**0.25 / np.sqrt(2**v * factorial(v))
    
    # Wavefunction in atomic units (normalized in bohr^(-1/2))
    psi_bohr = N_v * H_v(alpha * q_bohr) * np.exp(-0.5 * beta * q_bohr**2)
    
    # Convert to Å normalization: ψ_Å = ψ_bohr × √(ANGSTROM_TO_BOHR)
    psi_ang = psi_bohr * np.sqrt(ANGSTROM_TO_BOHR)
    
    return psi_ang


# ============================================================================
# Validation
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("CO₂ Asymmetric Stretch Module Test")
    print("=" * 70)
    
    iso = ISOTOPOLOGUES['626']
    print(f"\nIsotopologue: {iso.name}")
    print(f"  r_eq = {iso.r_eq:.8f} Å")
    print(f"  ω₃ = {iso.omega_3:.3f} cm⁻¹")
    print(f"  x₃₃ = {iso.x_33:.3f} cm⁻¹")
    print(f"  μ₃ = {iso.mu_3:.2f} m_e")
    
    print(f"\nVibrational levels (cm⁻¹):")
    for v in range(5):
        G_v = iso.energy_levels_cm1[v]
        print(f"  v={v}: G(v) = {G_v:.2f}")
    
    # Test coordinate transformation
    print(f"\nCoordinate transformation test:")
    q3_test = 0.05  # Å displacement
    R_test = 7.0    # Å
    theta_test = np.pi / 4
    
    coords = q3_to_cartesian(q3_test, R_test, theta_test, iso)
    q3_back = cartesian_to_q3(coords, iso)
    
    print(f"  q₃ input  = {q3_test:.6f} Å")
    print(f"  q₃ output = {q3_back:.6f} Å")
    print(f"  Error = {abs(q3_test - q3_back):.2e} Å")
    
    # Test wavefunction
    print(f"\nHarmonic oscillator wavefunctions:")
    q_grid = np.linspace(-0.3, 0.3, 200)
    
    for v in [0, 1]:
        psi = harmonic_oscillator_wavefunction(q_grid, v, iso.omega_3, iso.mu_3)
        
        # Check normalization
        dq = q_grid[1] - q_grid[0]
        norm = np.sum(psi**2) * dq
        
        print(f"  v={v}: ∫|ψ|² dq = {norm:.6f} (should be 1.0)")
        print(f"        max|ψ| = {np.max(np.abs(psi)):.4f} Å^-½")
    
    print("\n" + "=" * 70)

    from co2_monomer import ISOTOPOLOGUES

    iso = ISOTOPOLOGUES['626']
    mu_numerical = iso.mu_eff_from_cartesian_jacobian(dq=1e-5)
    print(f"Numerical μ_eff from Jacobian: {mu_numerical:.4f} m_e")
    mu_3 = iso.mu_3_analytical
    print(f"Analytical μ_3: {mu_3:.4f} m_e")
