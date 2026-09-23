"""
averaging_nitrogen.py - Vibrational averaging of He–CO2 interaction potential
over NITROGEN VSCF coordinate q'_6 defined by COORD_TRANS=NC_TRANS.

Builds:
  Veff0(R,theta) = <psi_v0 | V_int(R,theta,q'_6) | psi_v0>
  Veff1(R,theta) = <psi_v1 | V_int(R,theta,q'_6) | psi_v1>

Assumes:
- You created q3_v0.npz and q3_v1.npz using extract_vscf_target.py
  (keys: grid, psi).
- The VSCF job used COORD=ZMAT and COORD_TRANS=NC_TRANS.
- V_int_cartesian expects a 4-atom Cartesian geometry from q3_to_cartesian-like
  routine (we provide one here), and returns Hartree (your old behavior).

Key changes vs your old averaging.py:
- weights are |psi|^2 (already normalized in DVR basis)
- geometries come from inverting NC_TRANS to get ZMAT bonds R13,R14 at each q'_6
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from typing import Tuple

from units import angstrom_to_bohr
from pes import V_int_cartesian  # keep using your existing PES evaluator


# -----------------------------------------------------------------------------
# Grids (unchanged)
# -----------------------------------------------------------------------------

def make_R_grid_ang_eff(R_min=2.0, R_max=25.0) -> np.ndarray:
    R1 = np.arange(R_min, 7.0 + 1e-12, 0.05)
    R2 = np.arange(7.0 + 0.10, 11.0 + 1e-12, 0.10)
    R3 = np.arange(11.0 + 0.25, 20.0 + 1e-12, 0.25)
    R4 = np.arange(20.0 + 0.50, R_max + 1e-12, 0.50)
    return np.concatenate([R1, R2, R3, R4])

def make_theta_grid_rad_gl(n_theta=48) -> np.ndarray:
    x, _w = leggauss(n_theta)
    theta = np.arccos(x)
    theta = np.sort(theta)
    return theta


# -----------------------------------------------------------------------------
# NITROGEN VSCF helpers
# -----------------------------------------------------------------------------

def load_vscf_npz(npz_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load (grid, psi) from q3_v0.npz / q3_v1.npz made by extract_vscf_target.py
    """
    d = np.load(npz_path)
    q = np.array(d["grid"], dtype=float)
    psi = np.array(d["psi"], dtype=float)
    return q, psi

def vscf_weights(psi: np.ndarray) -> np.ndarray:
    """
    NITROGEN SBWF columns are normalized in DVR basis (you saw sum|psi|^2 = 1).
    So weights for diagonal operator expectation are simply |psi|^2.
    """
    w = np.abs(psi) ** 2
    w /= w.sum()
    return w

def read_nc_trans(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read NC_TRANS as T-type transform: q' = T q + t.
    Each line has [T_row ... t_i].
    """
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            rows.append([float(x) for x in s.split()])
    A = np.array(rows, dtype=float)
    T = A[:, :-1]
    t = A[:, -1]
    if T.shape != (6, 6) or t.shape != (6,):
        raise ValueError(f"Expected 6x6 T and length-6 t from NC_TRANS; got {T.shape}, {t.shape}")
    return T, t

def zmat_from_q6(q6: float, Tinv: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    For your VSCF job: REF1 = 0 0 0 0 0 0, and you sample coordinate 6.
    So q' = (0,0,0,0,0,q6). Then q = T^{-1}(q' - t).

    q order matches ZMAT variables:
      (R12, R13, A213, R14, A214, D3214)

    We only use R13 and R14 for CO2 internal distortion.
    """
    qprime = np.zeros(6, dtype=float)
    qprime[5] = float(q6)
    q = Tinv @ (qprime - t)
    return q


# -----------------------------------------------------------------------------
# Geometry builder (replacement for old q3_to_cartesian)
# -----------------------------------------------------------------------------

def _get_masses_co2(iso):
    """
    Best-effort: try to pull masses from your iso object.
    Fallback to 12C16O16O if not present.
    """
    # common patterns people use
    for attr in ("masses", "mass", "atom_masses"):
        if hasattr(iso, attr):
            m = getattr(iso, attr)
            if isinstance(m, (list, tuple, np.ndarray)) and len(m) >= 3:
                # Expect [C, O, O] or similar
                return float(m[0]), float(m[1]), float(m[2])
    # fallback: 12C, 16O
    return 12.0, 15.999, 15.999

def co2_cartesian_from_bonds(r13_ang: float, r14_ang: float, iso) -> np.ndarray:
    """
    Build linear CO2 along z with COM at origin.
    C at zC, O's at zC +/- bond lengths.

    Returns (3,3) in Angstrom in order [C, O2, O3].
    """
    mC, mO2, mO3 = _get_masses_co2(iso)

    # Put CO2 on z-axis; choose zC to make COM=0
    # COM z = (mC*zC + mO2*(zC+r13) + mO3*(zC-r14)) / (mC+mO2+mO3) = 0
    zC = -(mO2 * r13_ang - mO3 * r14_ang) / (mC + mO2 + mO3)

    C  = np.array([0.0, 0.0, zC])
    O2 = np.array([0.0, 0.0, zC + r13_ang])
    O3 = np.array([0.0, 0.0, zC - r14_ang])
    return np.vstack([C, O2, O3])

def he_position_from_Rtheta(R_ang: float, theta_rad: float) -> np.ndarray:
    """
    Place He at Jacobi (R,theta) relative to CO2 COM at origin,
    with CO2 axis along +z and phi=0 (in xz-plane).
    """
    return np.array([R_ang * np.sin(theta_rad), 0.0, R_ang * np.cos(theta_rad)])

def q6_to_cartesian_bohr(q6: float, R_ang: float, theta_rad: float, iso,
                         Tinv: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Build 4-atom geometry [He, C, O2, O3] in BOHR suitable for V_int_cartesian.

    Steps:
      - invert NC_TRANS at this q6 -> primitive ZMAT coords q
      - extract R13,R14 bond lengths (Angstrom)
      - build linear CO2 geometry with COM at origin
      - place He at Jacobi (R,theta)
      - convert Angstrom -> Bohr
    """
    q = zmat_from_q6(q6, Tinv, t)
    r13 = float(q[1])  # R13
    r14 = float(q[3])  # R14

    He = he_position_from_Rtheta(R_ang, theta_rad)
    C_O2_O3 = co2_cartesian_from_bonds(r13, r14, iso)  # (3,3) Angstrom
    C, O2, O3 = C_O2_O3

    X_ang = np.vstack([He, C, O2, O3])
    return angstrom_to_bohr(X_ang)  # shape (4,3) in Bohr

def q6_to_cartesian_angstrom(q6: float, R_ang: float, theta_rad: float, iso,
                             Tinv: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Build 4-atom geometry [He, C, O2, O3] in BOHR suitable for V_int_cartesian.

    Steps:
      - invert NC_TRANS at this q6 -> primitive ZMAT coords q
      - extract R13,R14 bond lengths (Angstrom)
      - build linear CO2 geometry with COM at origin
      - place He at Jacobi (R,theta)
      - convert Angstrom -> Bohr
    """
    q = zmat_from_q6(q6, Tinv, t)
    r13 = float(q[1])  # R13
    r14 = float(q[3])  # R14

    He = he_position_from_Rtheta(R_ang, theta_rad)
    C_O2_O3 = co2_cartesian_from_bonds(r13, r14, iso)  # (3,3) Angstrom
    C, O2, O3 = C_O2_O3

    X_ang = np.vstack([He, C, O2, O3])
    return X_ang  # shape (4,3) in Bohr


# -----------------------------------------------------------------------------
# Main averaging routine (same structure, new weights + geometry mapping)
# -----------------------------------------------------------------------------

def build_veff_grid_from_vscf(
    q_grid: np.ndarray,
    psi0: np.ndarray,
    psi1: np.ndarray,
    iso,
    R_grid_ang: np.ndarray,
    theta_grid_rad: np.ndarray,
    nctrans_path: str = "NC_TRANS",
    progress_every: int = 10
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Veff0 and Veff1 on (R,theta) grid using NITROGEN VSCF wavefunctions.
    Energies are whatever V_int_cartesian returns (your old code: Hartree).
    """
    if len(q_grid) != len(psi0) or len(q_grid) != len(psi1):
        raise ValueError("q_grid, psi0, psi1 must have the same length.")

    w0 = vscf_weights(psi0)
    w1 = vscf_weights(psi1)

    # Read/invert NC_TRANS once
    T, t = read_nc_trans(nctrans_path)
    Tinv = np.linalg.inv(T)

    nR = len(R_grid_ang)
    nT = len(theta_grid_rad)
    Veff0 = np.zeros((nR, nT), dtype=float)
    Veff1 = np.zeros((nR, nT), dtype=float)

    # Loop
    for iR, R in enumerate(R_grid_ang):
        for it, th in enumerate(theta_grid_rad):
            # Evaluate V_int on the q-grid (q is q'_6 values from NPZ)
            Vq = np.array(
                [V_int_cartesian(q6_to_cartesian_angstrom(q6, R, th, iso, Tinv, t))
                 for q6 in q_grid],
                dtype=float
            )

            #X_bohr = q6_to_cartesian_bohr(0, R, th, iso, Tinv, t)   # what you're feeding now
            #E_bohr = V_int_cartesian(X_bohr)

            #X_ang = X_bohr / 1.8897259886
            #E_ang = V_int_cartesian(X_ang)

            #print("E(bohr coords passed) =", E_bohr)
            #print("E(angstrom coords passed) =", E_ang)

            #print(q6_to_cartesian_bohr(-4, R, th, iso, Tinv,t))
            #print(V_int_cartesian(q6_to_cartesian_bohr(0, R, th, iso, Tinv, t)))

            #print(Vq)
            #exit()
            Veff0[iR, it] = np.sum(w0 * Vq)
            Veff1[iR, it] = np.sum(w1 * Vq)

        if progress_every and (iR % progress_every == 0):
            print(f"Completed R[{iR}/{nR-1}] = {R:.3f} Å")

    return Veff0, Veff1


# -----------------------------------------------------------------------------
# Export (unchanged)
# -----------------------------------------------------------------------------

def write_grid_dat(filename: str,
                   R_grid_ang: np.ndarray,
                   theta_grid_rad: np.ndarray,
                   Veff: np.ndarray) -> None:
    with open(filename, "w") as f:
        f.write("# R_ang  theta_rad  Veff_hartree\n")
        for iR, R in enumerate(R_grid_ang):
            for it, th in enumerate(theta_grid_rad):
                f.write(f"{R:12.6f} {th:12.8f} {Veff[iR, it]: .16e}\n")


def save_npz(filename: str,
             q_grid: np.ndarray,
             psi0: np.ndarray,
             psi1: np.ndarray,
             R_grid_ang: np.ndarray,
             theta_grid_rad: np.ndarray,
             Veff0: np.ndarray,
             Veff1: np.ndarray,
             iso_name: str = "") -> None:
    np.savez(
        filename,
        q_grid=q_grid,
        psi0=psi0,
        psi1=psi1,
        R_grid_ang=R_grid_ang,
        theta_grid_rad=theta_grid_rad,
        Veff0=Veff0,
        Veff1=Veff1,
        iso_name=iso_name
    )

if __name__ == "__main__":
    import argparse
    from co2_monomer import ISOTOPOLOGUES

    ap = argparse.ArgumentParser(
        description="Build Veff(R,theta) grids from NITROGEN VSCF modals + NC_TRANS"
    )
    ap.add_argument(
        "--iso",
        required=True,
        help="CO2 isotopologue key, e.g. 626, 636, 828",
    )
    ap.add_argument(
        "--q0",
        default=None,
        help="NPZ for v3=0 modal (defaults to q3_<iso>_v0.npz if present, else q3_v0.npz)",
    )
    ap.add_argument(
        "--q1",
        default=None,
        help="NPZ for v3=1 modal (defaults to q3_<iso>_v1.npz if present, else q3_v1.npz)",
    )
    ap.add_argument("--nctrans", default="NC_TRANS", help="NC_TRANS file path")
    ap.add_argument("--Rmin", type=float, default=1.50)
    ap.add_argument("--Rmax", type=float, default=25.0)
    ap.add_argument("--ntheta", type=int, default=48)
    ap.add_argument("--progress", type=int, default=10)
    args = ap.parse_args()

    if args.iso not in ISOTOPOLOGUES:
        raise KeyError(f"--iso {args.iso} not in ISOTOPOLOGUES keys: {list(ISOTOPOLOGUES.keys())}")

    iso = ISOTOPOLOGUES[args.iso]

    # Default filenames that include isotopologue label
    default_q0 = f"q3_{args.iso}_v0.npz"
    default_q1 = f"q3_{args.iso}_v1.npz"

    # Backward-compatible fallback names
    q0_path = args.q0 or (default_q0 if __import__("os").path.exists(default_q0) else "q3_v0.npz")
    q1_path = args.q1 or (default_q1 if __import__("os").path.exists(default_q1) else "q3_v1.npz")

    # Load NITROGEN VSCF wavefunctions (from your extracted NPZ files)
    q0, psi0 = load_vscf_npz(q0_path)
    q1, psi1 = load_vscf_npz(q1_path)

    # sanity: grids should match
    if not np.allclose(q0, q1):
        raise ValueError(f"{q0_path} and {q1_path} have different grids; fix before averaging.")
    q_grid = q0

    R_grid_ang = make_R_grid_ang_eff(R_min=args.Rmin, R_max=args.Rmax)
    theta_grid_rad = make_theta_grid_rad_gl(n_theta=args.ntheta)

    print(f"iso={args.iso}  nq={len(q_grid)}  nR={len(R_grid_ang)}  ntheta={len(theta_grid_rad)}")

    Veff0, Veff1 = build_veff_grid_from_vscf(
        q_grid, psi0, psi1, iso, R_grid_ang, theta_grid_rad,
        nctrans_path=args.nctrans,
        progress_every=args.progress
    )

    # Output names tagged by isotopologue
    out_npz = f"Veff_grids_from_NITROGEN_{args.iso}.npz"
    out_v0 = f"Veff0_grid_{args.iso}.dat"
    out_v1 = f"Veff1_grid_{args.iso}.dat"

    save_npz(out_npz, q_grid, psi0, psi1,
             R_grid_ang, theta_grid_rad, Veff0, Veff1,
             iso_name=getattr(iso, "name", args.iso))
    write_grid_dat(out_v0, R_grid_ang, theta_grid_rad, Veff0)
    write_grid_dat(out_v1, R_grid_ang, theta_grid_rad, Veff1)

    print("Wrote:")
    print(" ", out_npz)
    print(" ", out_v0)
    print(" ", out_v1)

