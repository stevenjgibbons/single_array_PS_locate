#!/usr/bin/env python3

import sys
import numpy as np

# ----------------------------------------------------------------------
# Your existing interpolation functions
# ----------------------------------------------------------------------

def gfdcfd(X, XARR):
    XARR = np.array(XARR, dtype=float)
    NNDS = len(XARR)

    # Check uniqueness
    if len(np.unique(XARR)) != NNDS:
        raise ValueError("Grid points in XARR are not distinct.")

    # Compute h_i = x_i - X
    h = XARR - X

    # Build matrix
    COEFM = np.zeros((NNDS, NNDS), dtype=float)
    for nder in range(NNDS):
        for inode in range(NNDS):
            if nder == 0:
                COEFM[inode, nder] = 1.0
            else:
                COEFM[inode, nder] = COEFM[inode, nder - 1] * (h[inode] / nder)

    # Invert
    return np.linalg.inv(COEFM)

def interp_1d(x, xarr, yarr, n=6):
    # exact grid hit → return exact value
    if x in xarr:
        return yarr[np.where(xarr == x)[0][0]]

    if x < xarr[0] or x > xarr[-1]:
        return None

    idx = np.searchsorted(xarr, x)
    half = n // 2
    i0 = max(0, idx - half)
    i1 = i0 + n
    if i1 > len(xarr):
        i1 = len(xarr)
        i0 = i1 - n

    xa = xarr[i0:i1]
    ya = yarr[i0:i1]

    M = gfdcfd(x, xa)
    w = M[0, :]
    return float(np.dot(w, ya))

# ----------------------------------------------------------------------
# Read traveltime table
# ----------------------------------------------------------------------

def read_diff_tt(filename):
    """
    Reads a traveltime table of format:
    P1  S1  (S1-P1)  distance_deg  depth_km

    Returns arrays: ttp, tts, tts_minus_ttp, distdeg
    """
    data = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            data.append([float(p) for p in parts[:4]])

    arr = np.array(data)
    ttp = arr[:, 0]
    tts = arr[:, 1]
    tts_minus_ttp = arr[:, 2]
    distdeg = arr[:, 3]

    return ttp, tts, tts_minus_ttp, distdeg

# ----------------------------------------------------------------------
# Compute distance and travel times from S-P value
# ----------------------------------------------------------------------

def dist_from_table(S_minus_P, ttp, tts, tts_minus_ttp, distdeg):
    """
    1. Interpolate Distance from S_minus_P
    2. Interpolate P and S times from Distance
    Returns: (distance, P_time, S_time)
    """

    # Step 1: Get distance from S-P
    dist = interp_1d(S_minus_P, tts_minus_ttp, distdeg)
    if dist is None:
        raise ValueError("S_minus_P is out of table range.")

    # Step 2: Interpolate times from distance
    P_time = interp_1d(dist, distdeg, ttp)
    S_time = interp_1d(dist, distdeg, tts)

    return dist, P_time, S_time

# ----------------------------------------------------------------------
# Standalone program
# ----------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print("Usage:  python prog.py  tablefile.txt  S_minus_P")
        sys.exit(1)

    filename = sys.argv[1]
    S_minus_P = float(sys.argv[2])

    ttp, tts, tts_minus_ttp, distdeg = read_diff_tt(filename)
    dist, P_time, S_time = dist_from_table(
        S_minus_P, ttp, tts, tts_minus_ttp, distdeg
    )

    print(f"P_travel_time = {P_time:.6f}")
    print(f"S_travel_time = {S_time:.6f}")
    print(f"S_minus_P     = {S_minus_P:.6f}")
    print(f"Distance_deg  = {dist:.6f}")

if __name__ == "__main__":
    main()
