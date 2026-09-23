#!/usr/bin/env python3
import argparse
import numpy as np

def read_nc_trans(path: str):
    """
    Read NC_TRANS as T-type transformation: q' = T q + t.
    File format: each line is [T_row ... t_i]. :contentReference[oaicite:6]{index=6}
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
    return T, t

def zmat_heco2_xyz(R12, R13, A213, R14, A214, D3214):
    """
    ZMAT:
      1 He
      2 C1  1 R12
      3 O2  2 R13  1 A213
      4 O3  2 R14  1 A214  3 D3214
    Angles assumed radians; distances Angstrom.
    Body-fixed convention matches standard ZMAT: 1 at origin, 2 on +z, 3 in xz-plane with x>0. :contentReference[oaicite:7]{index=7}
    """
    # Atom 1: He at origin
    He = np.array([0.0, 0.0, 0.0])

    # Atom 2: C on +z
    C = np.array([0.0, 0.0, R12])

    # Atom 3: O2, bonded to C, angle O2-C-He = A213, in xz-plane (y=0)
    # Place relative to C with polar angle measured from vector C->He (which points -z from C).
    # We'll construct using a generic "place atom with distance+angle+dihedral" routine by faking a dihedral.
    def place_atom(p_b, p_a, p_d, r, theta, phi):
        """
        Place new atom i with:
          |i - b| = r
          angle(i-b-a) = theta
          dihedral(i-b-a-d) = phi
        """
        b = p_b; a = p_a; d = p_d
        e1 = (b - a)
        e1 /= np.linalg.norm(e1)
        n = np.cross(d - a, e1)
        n /= np.linalg.norm(n)
        e2 = n
        e3 = np.cross(e1, e2)
        return b + r * (-np.cos(theta)*e1 + np.sin(theta)*(np.cos(phi)*e3 + np.sin(phi)*e2))

    # For atom 3, choose a dummy dihedral reference so it lies in xz-plane with +x.
    # Let d_ref be +x direction from He.
    d_ref = He + np.array([1.0, 0.0, 0.0])
    O2 = place_atom(C, He, d_ref, R13, A213, 0.0)

    # Atom 4: O3, bonded to C, angle O3-C-He = A214, dihedral O3-C-He-O2 = D3214
    O3 = place_atom(C, He, O2, R14, A214, D3214)

    return np.vstack([He, C, O2, O3])

def write_xyz(frames, symbols, outpath):
    with open(outpath, "w") as f:
        for k, X in enumerate(frames):
            f.write(f"{len(symbols)}\n")
            f.write(f"frame {k}\n")
            for s, (x,y,z) in zip(symbols, X):
                f.write(f"{s:2s} {x: .10f} {y: .10f} {z: .10f}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True, help="NPZ containing DVR grid (key: grid)")
    ap.add_argument("--nctrans", default="NC_TRANS", help="NC_TRANS file")
    ap.add_argument("--coord", type=int, default=6, help="Which transformed coordinate q'_i to scan (1..6)")
    ap.add_argument("--out", default="scan.xyz")
    ap.add_argument("--symbols", default="He C O O", help="Element symbols in ZMAT order")
    ap.add_argument("--print_q0", action="store_true", help="Print primitive ZMAT coords at q'=0")
    args = ap.parse_args()

    qgrid = np.load(args.npz)["grid"]
    T, t = read_nc_trans(args.nctrans)

    if T.shape[0] != 6 or T.shape[1] != 6:
        raise SystemExit(f"Expected 6x6 T in NC_TRANS, got {T.shape}")

    # q = T^{-1}(q' - t)
    Tinv = np.linalg.inv(T)

    # optional: primitive coordinates at q'=0
    if args.print_q0:
        qprime0 = np.zeros(6)
        q0 = Tinv @ (qprime0 - t)
        print("Primitive ZMAT coords at q'=0 (R12,R13,A213,R14,A214,D3214):")
        print(q0)

    # Build XYZ frames by scanning q'_coord over the DVR grid
    frames = []
    coord_idx = args.coord - 1
    for val in qgrid:
        qprime = np.zeros(6)
        qprime[coord_idx] = float(val)
        q = Tinv @ (qprime - t)
        X = zmat_heco2_xyz(*q)
        frames.append(X)

    symbols = args.symbols.split()
    write_xyz(frames, symbols, args.out)
    print(f"Wrote {len(frames)} frames to {args.out}")

if __name__ == "__main__":
    main()
