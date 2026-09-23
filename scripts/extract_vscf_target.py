#!/usr/bin/env python3
"""
Read NITROGEN VSCF_MODES_TARGETxxx.dat written by VCI_saveVSCFTargettoDisk()
(as in the C code you provided).

Usage:
  python extract_vscf_target.py VSCF_MODES_TARGET000.dat --mode 2 --col 0
  python extract_vscf_target.py VSCF_MODES_TARGET000.dat --summary
"""

from __future__ import annotations
import argparse
import os
import struct
from dataclasses import dataclass
import numpy as np


@dataclass
class ModeBlock:
    n: int
    grid: np.ndarray     # (n,)
    VI: np.ndarray       # (n,)
    GI: np.ndarray       # (n,)
    UI: np.ndarray       # (n,)
    SCFE: np.ndarray     # (n,)
    SBWF: np.ndarray     # (n,n) column-order in file
    DSBWF: np.ndarray    # (n,n)
    D1OP: np.ndarray     # (n,n)
    D2OP: np.ndarray     # (n,n)


@dataclass
class VSCFTargetFile:
    Nm: int
    target: list[int]
    modes: list[ModeBlock]
    RVSCF: int
    J: int | None = None
    RE: np.ndarray | None = None     # (NJ,)
    RWF: np.ndarray | None = None    # (NJ,NJ)


def _read_exact(f, nbytes: int) -> bytes:
    b = f.read(nbytes)
    if len(b) != nbytes:
        raise EOFError(f"Unexpected EOF at pos {f.tell()} (wanted {nbytes}, got {len(b)})")
    return b


def _read_int32(f, endian: str) -> int:
    return struct.unpack(endian + "i", _read_exact(f, 4))[0]


def _read_uint16(f, endian: str) -> int:
    return struct.unpack(endian + "H", _read_exact(f, 2))[0]


def _read_f64_array(f, n: int, endian: str) -> np.ndarray:
    # np dtype uses '<f8' or '>f8'
    dt = (endian + "f8")
    return np.frombuffer(_read_exact(f, 8 * n), dtype=dt, count=n).copy()


def _read_f64_mat(f, n: int, endian: str) -> np.ndarray:
    # File is "column order" => Fortran-order reshape to get columns as eigenvectors
    arr = _read_f64_array(f, n * n, endian)
    return arr.reshape((n, n), order="F")


def read_vscf_modes_target(path: str, endian: str = "<") -> VSCFTargetFile:
    """
    Parse VSCF_MODES_TARGETxxx.dat per the provided C I/O functions.
    endian: '<' little-endian (expected), '>' big-endian
    """
    fsize = os.path.getsize(path)
    with open(path, "rb") as f:
        Nm = _read_int32(f, endian)
        if not (1 <= Nm <= 500):
            raise ValueError(f"Unreasonable Nm={Nm}. Try endian='>' or check file.")

        target = list(struct.unpack(endian + f"{Nm}i", _read_exact(f, 4 * Nm)))

        modes: list[ModeBlock] = []
        for i in range(Nm):
            n = _read_uint16(f, endian)
            if not (2 <= n <= 20000):
                raise ValueError(f"Unreasonable n={n} for mode {i+1} at pos {f.tell()} (file size {fsize})")

            grid = _read_f64_array(f, n, endian)
            VI   = _read_f64_array(f, n, endian)
            GI   = _read_f64_array(f, n, endian)
            UI   = _read_f64_array(f, n, endian)
            SCFE = _read_f64_array(f, n, endian)

            SBWF  = _read_f64_mat(f, n, endian)
            DSBWF = _read_f64_mat(f, n, endian)
            D1OP  = _read_f64_mat(f, n, endian)
            D2OP  = _read_f64_mat(f, n, endian)

            modes.append(ModeBlock(n, grid, VI, GI, UI, SCFE, SBWF, DSBWF, D1OP, D2OP))

        RVSCF = _read_int32(f, endian)

        out = VSCFTargetFile(Nm=Nm, target=target, modes=modes, RVSCF=RVSCF)

        if RVSCF:
            J = _read_int32(f, endian)
            NJ = 2 * J + 1
            RE  = _read_f64_array(f, NJ, endian)
            RWF = _read_f64_array(f, NJ * NJ, endian).reshape((NJ, NJ), order="F")
            out.J = J
            out.RE = RE
            out.RWF = RWF

        # Optional sanity: ensure we're at EOF (or only a few trailing bytes)
        leftover = f.read()
        if len(leftover) not in (0,):
            # Not fatal, but useful warning
            pass

        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", help="VSCF_MODES_TARGETxxx.dat")
    ap.add_argument("--endian", choices=["<", ">"], default="<", help="Byte order (default: little-endian)")
    ap.add_argument("--summary", action="store_true", help="Print summary of modes and sizes")
    ap.add_argument("--mode", type=int, default=None, help="1-based mode index to extract")
    ap.add_argument("--col", type=int, default=0, help="0-based column of SBWF (modal quantum number)")
    ap.add_argument("--out", default=None, help="Output .npz filename (optional)")
    args = ap.parse_args()

    data = read_vscf_modes_target(args.file, endian=args.endian)

    if args.summary or args.mode is None:
        print(f"Nm = {data.Nm}")
        print(f"target indices = {data.target}")
        for i, m in enumerate(data.modes, start=1):
            print(f"mode {i}: n={m.n}, grid[min,max]=({m.grid.min():.6g},{m.grid.max():.6g}), "
                  f"SCFE[0:3]={m.SCFE[:3]}")
        print(f"RVSCF = {data.RVSCF}")
        if data.RVSCF:
            print(f"J = {data.J}, NJ = {2*data.J+1}")

    if args.mode is not None:
        i = args.mode - 1
        if i < 0 or i >= data.Nm:
            raise SystemExit(f"--mode must be 1..{data.Nm}")
        m = data.modes[i]
        col = args.col
        if col < 0 or col >= m.n:
            raise SystemExit(f"--col must be 0..{m.n-1} for this mode (n={m.n})")

        psi = m.SBWF[:, col]  # modal eigenfunction on DVR grid
        # Not all DVRs have trivial weights; but this is still a useful sanity check:
        norm2 = float(np.sum(np.abs(psi) ** 2))
        print(f"Extracted mode {args.mode} col {col}: psi shape={psi.shape}, sum|psi|^2={norm2:.6g}")

        if args.out:
            np.savez(
                args.out,
                Nm=data.Nm,
                target=np.array(data.target, dtype=np.int32),
                mode=args.mode,
                col=col,
                grid=m.grid,
                psi=psi,
                VI=m.VI, GI=m.GI, UI=m.UI,
                SCFE=m.SCFE,
            )
            print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()

