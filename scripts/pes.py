"""
pes.py - He-CO₂ Potential Energy Surface Interface

Provides V(R, θ) in Jacobi coordinates by wrapping the ab initio PES library.
Handles unit conversions and coordinate transformations internally.
"""
from __future__ import annotations

import numpy as np
import ctypes
from pathlib import Path
from typing import Union

from units import bohr_to_angstrom, cm1_to_hartree, hartree_to_cm1

# ============================================================================
# Constants
# ============================================================================
CM1_PER_KCALMOL = 349.755  # Exact conversion factor
DEFAULT_R_CO = 1.16221264055      # Equilibrium CO₂ C-O distance (Å)

# ============================================================================
# C++ Library Interface (Singleton Pattern)
# ============================================================================
_LIB = None

def _get_lib():
    """Load the C++ PES library (cached after first call)"""
    global _LIB
    if _LIB is None:
        libpath = Path(__file__).resolve().parent / "libCO2He.dylib"
        if not libpath.exists():
            raise FileNotFoundError(f"PES library not found: {libpath}")
        
        _LIB = ctypes.CDLL(str(libpath))
        
        # Define C++ function signatures
        _LIB.getN.argtypes = ()
        _LIB.getN.restype = ctypes.c_int
        
        # He-CO₂ surface (4 atoms, 12 coordinates)
        _LIB.calcSurface_HeCO2.argtypes = (
            ctypes.POINTER(ctypes.c_double),  # coords[12]
            ctypes.POINTER(ctypes.c_int),     # perm (unused)
        )
        _LIB.calcSurface_HeCO2.restype = ctypes.c_double
        
        # CO₂ monomer surface (3 atoms, 9 coordinates)
        _LIB.calcSurface_CO2.argtypes = (
            ctypes.POINTER(ctypes.c_double),  # coords[9]
            ctypes.POINTER(ctypes.c_int),     # perm (unused)
        )
        _LIB.calcSurface_CO2.restype = ctypes.c_double
        
        # Validate library
        n_atoms = _LIB.getN()
        if n_atoms != 4:
            raise RuntimeError(
                f"PES library reports {n_atoms} atoms; expected 4 for He-CO₂"
            )
    
    return _LIB


# ============================================================================
# Coordinate Transformations
# ============================================================================

def jacobi_to_cartesian(
    R_ang: float,
    theta_rad: float,
    rCO_ang: float = DEFAULT_R_CO
) -> np.ndarray:
    """
    Convert Jacobi (R, θ) to Cartesian coordinates.
    
    Geometry:
        - CO₂ axis along z, C at origin
        - He position: (R sin θ, 0, R cos θ)
        - θ = 0: He along +z (linear He-C-O-O)
        - θ = π: He along -z (linear O-O-C-He)
    
    Parameters
    ----------
    R_ang : float
        He-CO₂ center-of-mass distance (Å)
    theta_rad : float
        Jacobi angle (radians, 0 ≤ θ ≤ π)
    rCO_ang : float, optional
        C-O bond length (Å), default = 1.1632 Å
    
    Returns
    -------
    coords : np.ndarray, shape (4, 3)
        Cartesian coordinates [He, C, O₁, O₂] in Å
    """
    # CO₂ geometry (linear, symmetric)
    C  = np.array([0.0, 0.0,  0.0])
    O1 = np.array([0.0, 0.0, +rCO_ang])
    O2 = np.array([0.0, 0.0, -rCO_ang])
    
    # He position in spherical coordinates
    He = np.array([
        R_ang * np.sin(theta_rad),
        0.0,
        R_ang * np.cos(theta_rad)
    ])
    
    return np.vstack([He, C, O1, O2])


# ============================================================================
# PES Evaluation
# ============================================================================

def V_cartesian_HeCO2(coords_angstrom: np.ndarray) -> float:
    """
    Evaluate He-CO₂ PES in Cartesian coordinates.
    
    Parameters
    ----------
    coords_angstrom : np.ndarray, shape (4, 3)
        Cartesian positions [He, C, O₁, O₂] in Å
    
    Returns
    -------
    energy : float
        Potential energy in hartree
    """
    lib = _get_lib()
    
    coords = np.asarray(coords_angstrom, dtype=np.float64)
    if coords.shape != (4, 3):
        raise ValueError(
            f"coords_angstrom must have shape (4, 3); got {coords.shape}"
        )
    
    # Flatten to C++ expected format: double[12]
    x_flat = np.ascontiguousarray(coords.ravel(), dtype=np.float64)
    X = x_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    # PERM argument (required by interface but unused)
    perm = ctypes.c_int(0)
    
    # Call C++ PES (returns cm⁻¹)
    E_cm1 = float(lib.calcSurface_HeCO2(X, ctypes.byref(perm)))
    
    # Convert to hartree
    return cm1_to_hartree(E_cm1)


def V_cartesian_CO2(coords_angstrom: np.ndarray) -> float:
    """
    Evaluate isolated CO₂ PES in Cartesian coordinates.
    
    Parameters
    ----------
    coords_angstrom : np.ndarray, shape (3, 3)
        Cartesian positions [C, O₁, O₂] in Å
    
    Returns
    -------
    energy : float
        Potential energy in hartree
    """
    lib = _get_lib()
    
    coords = np.asarray(coords_angstrom, dtype=np.float64)
    if coords.shape != (3, 3):
        raise ValueError(
            f"coords_angstrom must have shape (3, 3); got {coords.shape}"
        )
    
    # Flatten to C++ expected format: double[9]
    x_flat = np.ascontiguousarray(coords.ravel(), dtype=np.float64)
    X = x_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    # PERM argument (required by interface but unused)
    perm = ctypes.c_int(0)
    
    # Call C++ PES (returns cm⁻¹)
    E_cm1 = float(lib.calcSurface_CO2(X, ctypes.byref(perm)))
    
    # Convert to hartree
    return cm1_to_hartree(E_cm1)


# Keep old name for backward compatibility
def V_cartesian(coords_angstrom: np.ndarray) -> float:
    """
    Evaluate PES in Cartesian coordinates (auto-detects He-CO₂ vs CO₂).
    
    Parameters
    ----------
    coords_angstrom : np.ndarray, shape (4, 3) or (3, 3)
        Cartesian positions in Å
        - (4, 3): [He, C, O₁, O₂] for He-CO₂
        - (3, 3): [C, O₁, O₂] for isolated CO₂
    
    Returns
    -------
    energy : float
        Potential energy in hartree
    """
    coords = np.asarray(coords_angstrom, dtype=np.float64)
    
    if coords.shape == (4, 3):
        return V_cartesian_HeCO2(coords)
    elif coords.shape == (3, 3):
        return V_cartesian_CO2(coords)
    else:
        raise ValueError(
            f"coords must have shape (4,3) for He-CO₂ or (3,3) for CO₂; got {coords.shape}"
        )


def V(
    R_bohr: Union[float, np.ndarray],
    theta_rad: Union[float, np.ndarray],
    rCO_ang: float = DEFAULT_R_CO
) -> Union[float, np.ndarray]:
    """
    Evaluate He-CO₂ PES in Jacobi coordinates.
    
    Parameters
    ----------
    R_bohr : float or array_like
        Jacobi distance in bohr
    theta_rad : float or array_like
        Jacobi angle in radians (0 ≤ θ ≤ π)
    rCO_ang : float, optional
        CO₂ C-O bond length (Å), default = 1.1612 Å
    
    Returns
    -------
    energy : float or np.ndarray
        Potential energy in hartree, same shape as input
    
    Notes
    -----
    - For 2D DVR: pass R as scalar, theta_rad as 1D array
    - For grid evaluation: use np.meshgrid(R_vals, theta_vals)
    """
    # Convert R to Angstrom
    R_ang = bohr_to_angstrom(R_bohr)
    
    # Handle scalar/array inputs uniformly
    theta_rad = np.atleast_1d(theta_rad)
    scalar_input = theta_rad.ndim == 0 or theta_rad.size == 1
    
    # Evaluate PES at each angle
    energies = np.empty_like(theta_rad, dtype=float)
    for i, theta in enumerate(theta_rad.flat):
        coords = jacobi_to_cartesian(float(R_ang), float(theta), rCO_ang)
        energies.flat[i] = V_cartesian_HeCO2(coords)
    
    return float(energies[0]) if scalar_input else energies.reshape(theta_rad.shape)

import numpy as np
from co2_monomer import q3_to_cartesian, q3_to_cartesian_co2_only, CO2Isotopologue, ISOTOPOLOGUES

def V_int(q3_ang: float, R_ang: float, theta_rad: float, iso: CO2Isotopologue) -> float:
    """
    Interaction energy in hartree:
      V_int = V(He+CO2) - V(CO2)
    """
    coords4 = q3_to_cartesian(q3_ang, R_ang, theta_rad, iso)     # (4,3)
    E_tot = V_cartesian_HeCO2(coords4)

    coords3 = q3_to_cartesian_co2_only(q3_ang, iso)              # (3,3)
    E_co2 = V_cartesian_CO2(coords3)

    return E_tot - E_co2

def V_int_cartesian(coords4_ang: np.ndarray) -> float:
    """Interaction only: V(HeCO2) - V(CO2) in hartree."""
    if coords4_ang.shape != (4,3):
        raise ValueError("coords must be (4,3) [He,C,O1,O2] in Å")
    HeCO2 = V_cartesian_HeCO2(coords4_ang)          # hartree
    CO2   = V_cartesian_CO2(coords4_ang[1:4, :])    # hartree
    return HeCO2 - CO2


def sanity_check(iso):
    q3 = 0.04
    theta = np.pi/2  # T-shaped
    for R in [3.0, 6.0, 10.0, 20.0, 50.0]:  # Å
        print(R, V_int(q3, R, theta, iso))

    R=7.0; theta=np.pi/2
    for q3 in [0.02, 0.05, 0.08]:
        vp = V_int(+q3, R, theta, iso)
        vm = V_int(-q3, R, theta, iso)
        print(q3, vp, vm, vp-vm)



# ============================================================================
# Test/Model Potential
# ============================================================================

def V_model(
    R_bohr: Union[float, np.ndarray],
    theta_rad: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Simple model potential for testing (Morse + P₂ anisotropy).
    
    V(R, θ) = D_e[1 - exp(-a(R - R_e))]² + V_anis P₂(cos θ)
    
    Parameters
    ----------
    R_bohr : float or array_like
        Jacobi distance in bohr
    theta_rad : float or array_like
        Jacobi angle in radians
    
    Returns
    -------
    energy : float or np.ndarray
        Model potential in hartree
    """
    # Isotropic part: Morse-like attractive well
    De = 10.0 / 219474.63  # ~10 cm⁻¹ depth in hartree
    Re = 6.5               # Equilibrium distance (bohr)
    a = 0.8                # Width parameter
    
    V_iso = De * (1.0 - np.exp(-a * (R_bohr - Re)))**2 - De
    
    # Anisotropic part: P₂(cos θ) modulation
    P2 = 0.5 * (3 * np.cos(theta_rad)**2 - 1)
    V_anis = 0.2 * De * P2  # 20% anisotropy
    
    return V_iso + V_anis

# -------------------------------------------------------------------------
# Effective surfaces from vibrational averaging (grid-based, NPZ)
# -------------------------------------------------------------------------


import numpy as np
from dataclasses import dataclass
from typing import Optional, Union

from units import angstrom_to_bohr, bohr_to_angstrom


ArrayLike = Union[float, np.ndarray]



def veff_domain():
    """
    Returns (Rmin_bohr, Rmax_bohr, thmin_rad, thmax_rad) for loaded Veff grids.
    """
    if _VEFF0 is None:
        raise RuntimeError("Veff not loaded. Call load_veff_npz(...) first.")
    return (_VEFF0.R_bohr[0], _VEFF0.R_bohr[-1], _VEFF0.theta_rad[0], _VEFF0.theta_rad[-1])



@dataclass
class VeffSurface:
    """
    Holds a grid V(R,theta) and provides bilinear interpolation.
    R is stored in bohr, theta in radians, V in hartree.
    """
    R_bohr: np.ndarray              # shape (nR,)
    theta_rad: np.ndarray           # shape (nT,)
    V: np.ndarray                   # shape (nR,nT)

    def __post_init__(self):
        # Ensure increasing grids
        if not (np.all(np.diff(self.R_bohr) > 0) and np.all(np.diff(self.theta_rad) > 0)):
            raise ValueError("R_bohr and theta_rad grids must be strictly increasing.")
        if self.V.shape != (self.R_bohr.size, self.theta_rad.size):
            raise ValueError("V has wrong shape compared to grids.")

    def interp(self, R_bohr_in: ArrayLike, theta_in: ArrayLike) -> np.ndarray:
        """
        Bilinear interpolation on (R,theta). Supports scalar or array input.
        theta is assumed in [0,pi]. Values outside are clipped.
        R outside grid is clipped to nearest endpoint (safe for long-range tail).
        """
        R = np.asarray(R_bohr_in, dtype=float)
        th = np.asarray(theta_in, dtype=float)

        # Broadcast to common shape
        Rb, thb = np.broadcast_arrays(R, th)

        # Clip into grid domain
        Rb = np.clip(Rb, self.R_bohr[0], self.R_bohr[-1])
        thb = np.clip(thb, self.theta_rad[0], self.theta_rad[-1])

        # Find indices i such that R[i] <= x < R[i+1]
        iR = np.searchsorted(self.R_bohr, Rb, side="right") - 1
        iT = np.searchsorted(self.theta_rad, thb, side="right") - 1

        # Clamp so i+1 is valid
        iR = np.clip(iR, 0, self.R_bohr.size - 2)
        iT = np.clip(iT, 0, self.theta_rad.size - 2)

        R0 = self.R_bohr[iR]
        R1 = self.R_bohr[iR + 1]
        T0 = self.theta_rad[iT]
        T1 = self.theta_rad[iT + 1]

        # Fractional coordinates
        tR = (Rb - R0) / (R1 - R0)
        tT = (thb - T0) / (T1 - T0)

        # Gather the four corners
        V00 = self.V[iR,     iT]
        V10 = self.V[iR + 1, iT]
        V01 = self.V[iR,     iT + 1]
        V11 = self.V[iR + 1, iT + 1]

        # Bilinear interpolation
        V0 = V00 * (1 - tR) + V10 * tR
        V1 = V01 * (1 - tR) + V11 * tR
        Vout = V0 * (1 - tT) + V1 * tT

        return Vout


# Module-level cache
_VEFF0: Optional[VeffSurface] = None
_VEFF1: Optional[VeffSurface] = None



def load_veff_npz(npz_path: str) -> None:
    """
    Load vibrationally averaged effective surfaces from an NPZ created by averaging.py.

    Expected keys:
      R_grid_ang (Å), theta_grid_rad, Veff0 (hartree), Veff1 (hartree)

    After calling this, pes.Veff0 and pes.Veff1 become usable.
    """
    global _VEFF0, _VEFF1

    data = np.load(npz_path)

    # Required arrays
    R_ang = data["R_grid_ang"].astype(float)
    theta = data["theta_grid_rad"].astype(float)
    V0 = data["Veff0"].astype(float)
    V1 = data["Veff1"].astype(float)

    # Ensure theta is increasing (0..pi). If not, sort and reorder columns.
    if np.any(np.diff(theta) <= 0):
        idxT = np.argsort(theta)
        theta = theta[idxT]
        V0 = V0[:, idxT]
        V1 = V1[:, idxT]

    if np.any(np.diff(R_ang) <= 0):
        idxR = np.argsort(R_ang)
        R_ang = R_ang[idxR]
        V0 = V0[idxR, :]
        V1 = V1[idxR, :]

    # Convert R to bohr (solver uses bohr)
    R_bohr = angstrom_to_bohr(R_ang)

    _VEFF0 = VeffSurface(R_bohr=R_bohr, theta_rad=theta, V=V0)
    _VEFF1 = VeffSurface(R_bohr=R_bohr, theta_rad=theta, V=V1)

    print(f"[pes] Loaded Veff surfaces from: {npz_path}")
    print(f"[pes] R range: {R_bohr[0]:.3f}–{R_bohr[-1]:.3f} bohr  (n={R_bohr.size})")
    print(f"[pes] theta range: {theta[0]:.3f}–{theta[-1]:.3f} rad (n={theta.size})")


def Veff0(R_bohr: ArrayLike, theta_rad: ArrayLike) -> np.ndarray:
    """
    Effective surface for v3=0 (hartree). Call load_veff_npz() first.
    """
    if _VEFF0 is None:
        raise RuntimeError("Veff0 not loaded. Call pes.load_veff_npz('file.npz') first.")
    return _VEFF0.interp(R_bohr, theta_rad)


def Veff1(R_bohr: ArrayLike, theta_rad: ArrayLike) -> np.ndarray:
    """
    Effective surface for v3=1 (hartree). Call load_veff_npz() first.
    """
    if _VEFF1 is None:
        raise RuntimeError("Veff1 not loaded. Call pes.load_veff_npz('file.npz') first.")
    return _VEFF1.interp(R_bohr, theta_rad)

# ============================================================================
# Validation
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("He-CO₂ PES Interface Test")
    print("=" * 70)
    
    # Test coordinate transformation
    R_test = 9.0  # bohr
    theta_test = np.pi / 4  # 45 degrees
    
    coords = jacobi_to_cartesian(
        bohr_to_angstrom(R_test),
        theta_test
    )
    print(f"\nJacobi → Cartesian:")
    print(f"  R = {R_test:.3f} bohr, θ = {np.degrees(theta_test):.1f}°")
    print(f"  Coordinates (Å):")
    for i, atom in enumerate(['He', 'C', 'O₁', 'O₂']):
        print(f"    {atom}: [{coords[i, 0]:7.4f}, {coords[i, 1]:7.4f}, {coords[i, 2]:7.4f}]")
    
    # Test He-CO₂ PES evaluation
    try:
        E_abs = V(R_test, theta_test)
        print(f"\n  V_HeCO2(R={R_test:.1f} bohr, θ={np.degrees(theta_test):.1f}°) = {E_abs:.6e} hartree")
        print(f"                                             = {E_abs * 219474.63:.2f} cm⁻¹")
    except FileNotFoundError:
        print("\n  ⚠️  libCO2He.dylib not found; using model potential")
        E_model = V_model(R_test, theta_test)
        print(f"  V_model = {E_model:.6e} hartree = {E_model * 219474.63:.2f} cm⁻¹")
    
    # Test CO₂ monomer PES
    print(f"\nTesting CO₂ monomer surface:")
    try:
        coords_co2 = np.array([
            [0.0, 0.0, 0.0],           # C
            [0.0, 0.0, DEFAULT_R_CO],  # O1
            [0.0, 0.0, -DEFAULT_R_CO]  # O2
        ])
        E_co2 = V_cartesian_CO2(coords_co2)
        print(f"  V_CO2(r_eq) = {E_co2:.6e} hartree = {E_co2 * 219474.63:.2f} cm⁻¹")
    except Exception as e:
        print(f"  ⚠️  CO₂ surface not available: {e}")
    
    # Test array evaluation
    theta_array = np.linspace(0, np.pi, 5)
    E_array = V(R_test, theta_array)
    
    print(f"\nArray evaluation:")
    print(f"  {'θ (deg)':>10} {'V (cm⁻¹)':>12}")
    print("  " + "-" * 23)
    for th, E in zip(theta_array, E_array):
        print(f"  {np.degrees(th):10.1f} {E * 219474.63:12.2f}")
    
    # Test symmetry
    print(f"\nSymmetry test:")
    R = 5.824
    th = np.linspace(0, np.pi, 401)
    a = V(R, th)
    b = V(R, np.pi - th)
    print(f"  sym err (cm⁻¹): {hartree_to_cm1(np.max(np.abs(a - b))):.2e}")
    
    print("\n" + "=" * 70)

    coords4 = q3_to_cartesian(0.0, R_ang=999.0, theta_rad=0.7, iso=ISOTOPOLOGUES['626'])   # [He,C,O,O]
    E4 = V_cartesian_HeCO2(coords4)     # hartree
    E3 = V_cartesian_CO2(coords4[1:])   # hartree
    print((E4 - E3) * 219474.6313705)   # should be ~0 cm^-1 for all q3
    print(E4,E3)
    print("hello")


    sanity_check(ISOTOPOLOGUES['626'])




