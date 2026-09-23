import numpy as np
import struct

"""
def read_vscf_modes(fname):
    with open(fname, "rb") as f:
        def read_int():
            return struct.unpack("<i", f.read(4))[0]
        def read_uint16():
            return struct.unpack("<H", f.read(2))[0]
        def read_doubles(n):
            return np.frombuffer(f.read(8*n), dtype="<f8", count=n)

        Nm = read_int()
        mode_indices = [read_int() for _ in range(Nm)]  # "Target mode index list"
        modes = []

        for _ in range(Nm):
            ni = read_uint16()
            q  = read_doubles(ni)
            Vi = read_doubles(ni)
            Gi = read_doubles(ni)
            Ui = read_doubles(ni)
            Ei = read_doubles(ni)         # "VSCF energies for mode i"
            vec = read_doubles(ni*ni)     # eigenvector array, column order

            # Column order => reshape with Fortran order so columns are eigenvectors
            evecs = vec.reshape((ni, ni), order="F")
            modes.append(dict(ni=ni, q=q, Vi=Vi, Gi=Gi, Ui=Ui, Ei=Ei, evecs=evecs))

        RVSCF_flag = read_int()
        return dict(Nm=Nm, mode_indices=mode_indices, modes=modes, RVSCF_flag=RVSCF_flag)
"""
import numpy as np
import os, struct

def read_vscf_modes(fname, endian="<"):
    fsize = os.path.getsize(fname)

    def u(fmt, f):
        n = struct.calcsize(fmt)
        b = f.read(n)
        if len(b) != n:
            raise EOFError(f"EOF at pos {f.tell()} reading {fmt} (wanted {n}, got {len(b)})")
        return struct.unpack(fmt, b)

    with open(fname, "rb") as f:
        (Nm,) = u(endian+"i", f)
        mode_indices = list(u(endian+f"{Nm}i", f))

        modes = []
        for im in range(Nm):
            (ni,) = u(endian+"H", f)

            # *** critical: align to 8-byte boundary before reading doubles ***
            pad = (-f.tell()) % 8
            if pad:
                f.read(pad)

            need = 8*(5*ni + ni*ni)  # Table 2.7 logical content :contentReference[oaicite:1]{index=1}
            if f.tell() + need > fsize:
                raise ValueError(f"Not enough bytes for mode {im+1}: ni={ni}, pos={f.tell()}, need={need}, size={fsize}")

            q  = np.frombuffer(f.read(8*ni), dtype=endian+"f8", count=ni)
            Vi = np.frombuffer(f.read(8*ni), dtype=endian+"f8", count=ni)
            Gi = np.frombuffer(f.read(8*ni), dtype=endian+"f8", count=ni)
            Ui = np.frombuffer(f.read(8*ni), dtype=endian+"f8", count=ni)
            Ei = np.frombuffer(f.read(8*ni), dtype=endian+"f8", count=ni)

            vec = np.frombuffer(f.read(8*ni*ni), dtype=endian+"f8", count=ni*ni)
            evecs = vec.reshape((ni, ni), order="F")  # column order :contentReference[oaicite:2]{index=2}

            modes.append(dict(ni=ni, q=q, Ei=Ei, evecs=evecs))

        (rvscf_flag,) = u(endian+"i", f)
        return dict(Nm=Nm, mode_indices=mode_indices, modes=modes, rvscf_flag=rvscf_flag)


if __name__=="__main__":

    data = read_vscf_modes("VSCF_MODES_TARGET000.dat")
    q3_grid = data["modes"][idx]["q"]
    psi_v0 = data["modes"][idx]["evecs"][:, 0]   # ground modal
#    psi_v1 = data["modes"][idx]["evecs"][:, 1]   # 1-quantum excited modal

