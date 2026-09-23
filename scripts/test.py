import numpy as np

def read_xyz_frames(path):
    frames = []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            n = int(line.strip())
            f.readline()  # comment
            xyz = []
            for _ in range(n):
                parts = f.readline().split()
                xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
            frames.append(np.array(xyz))
    return frames

frames = read_xyz_frames("mode6_scan.xyz")

R12 = []
A213 = []
A214 = []
CO1 = []
CO2 = []

for X in frames:
    He, C, O2, O3 = X
    R12.append(np.linalg.norm(C - He))
    CO1.append(np.linalg.norm(O2 - C))
    CO2.append(np.linalg.norm(O3 - C))

print("He–C R12 range:", min(R12), max(R12))
print("C–O ranges:", (min(CO1), max(CO1)), (min(CO2), max(CO2)))
print("delta(CO): range", min(np.array(CO1)-np.array(CO2)), max(np.array(CO1)-np.array(CO2)))
