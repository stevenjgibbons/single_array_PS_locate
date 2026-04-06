import numpy as np
from datetime import datetime, timedelta

# --------------------------------------------------------------
# Parse ISO timestamps
# --------------------------------------------------------------

def _parse_iso8601(s):
    return datetime.fromisoformat(s)

# --------------------------------------------------------------
# Read arrival-time distribution
# --------------------------------------------------------------

def arrtime_distr_file_read(filename):
    """
    Reads a polygonal probability distribution for arrival times.

    First line:
        t_lo  t_hi      (ISO-8601 timestamps)

    Following lines:
        coefficients defining polygon heights
    """
    with open(filename, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    # --- Time limits ---
    t_lo_str, t_hi_str = lines[0].split()
    t_lo = _parse_iso8601(t_lo_str)
    t_hi = _parse_iso8601(t_hi_str)

    dt = (t_hi - t_lo).total_seconds()
    if dt <= 0:
        raise ValueError("Upper time must be later than lower time.")

    # --- Coefficients ---
    coeffs = np.array([float(x) for x in lines[1:]], dtype=float)
    N = len(coeffs)

    if N == 0:
        raise ValueError("No coefficients found.")

    if N == 1:
        coeffs = np.ones(2)
        N = 2

    # --- Edges in seconds from t_lo ---
    edges_sec = np.linspace(0.0, dt, N)

    # --- Normalize polygonal PDF ---
    bin_w = dt / (N - 1)
    area = 0.5 * bin_w * np.sum(coeffs[:-1] + coeffs[1:])
    if area <= 0:
        raise ValueError("Distribution area must be positive.")

    pdf = coeffs / area

    # --- Build CDF ---
    cdf = np.zeros(N)
    for i in range(1, N):
        cdf[i] = cdf[i-1] + 0.5 * (pdf[i] + pdf[i-1]) * bin_w

    cdf /= cdf[-1]

    return {
        "t_lo": t_lo,
        "t_hi": t_hi,
        "edges_sec": edges_sec,
        "pdf": pdf,
        "cdf": cdf,
        "N": N,
        "dt": dt
    }

# --------------------------------------------------------------
# Sample from polygon distribution
# --------------------------------------------------------------

def arrtime_distr_pick(u, dist):
    """
    u ∈ [0,1] returns a sampled datetime
    """
    cdf = dist["cdf"]
    edges = dist["edges_sec"]
    pdf = dist["pdf"]
    N = dist["N"]
    t_lo = dist["t_lo"]

    if u <= 0:
        return t_lo
    if u >= 1:
        return dist["t_hi"]

    # --- Locate segment ---
    idx = np.searchsorted(cdf, u) - 1
    if idx < 0: idx = 0
    if idx >= N - 1: idx = N - 2

    x0 = edges[idx]
    x1 = edges[idx + 1]
    y0 = pdf[idx]
    y1 = pdf[idx + 1]
    w = x1 - x0

    C0 = cdf[idx]
    Ulocal = u - C0

    # --- Flat segment ---
    if abs(y1 - y0) < 1e-12:
        frac = Ulocal / (y0 * w)
        return t_lo + timedelta(seconds=x0 + frac*w)

    # --- Solve quadratic ---
    a = 0.5 * (y1 - y0)
    b = y0
    A = a
    B = b
    C = -Ulocal / w

    disc = max(0.0, B*B - 4*A*C)
    t = (-B + np.sqrt(disc)) / (2*A)
    t = max(0.0, min(1.0, t))

    x = x0 + t * w
    return t_lo + timedelta(seconds=x)
