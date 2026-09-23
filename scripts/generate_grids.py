import numpy as np
from extract_vscf_target import read_vscf_modes_target  # from the script you used

def get_q3_modal_probs(vscf_file="VSCF_MODES_TARGET000.dat", q3_mode=6, v=0):
    """
    q3_mode is 1-based; v is 0-based modal quantum number (0=ground, 1=1-quantum, ...).
    Returns: q_grid, prob (unnormalized if DVR weights != 1; normalized by simple sum).
    """
    data = read_vscf_modes_target(vscf_file, endian="<")
    m = data.modes[q3_mode - 1]
    psi = m.SBWF[:, v]
    prob = np.abs(psi)**2
    prob = prob / prob.sum()  # simple normalization; replace with weighted norm if you have DVR weights
    return m.grid.copy(), prob

# Example:
q3, p0 = get_q3_modal_probs("VSCF_MODES_TARGET000.dat", q3_mode=6, v=0)
q3, p1 = get_q3_modal_probs("VSCF_MODES_TARGET000.dat", q3_mode=6, v=1)

print(q3,p0)
