import numpy as np

# --------------------------------------------------------------
# Read and normalize backazimuth distribution file
# --------------------------------------------------------------

def azi_distr_file_read(filename):
    with open(filename, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    # First line: azimuth limits
    lo, hi = map(float, lines[0].split())

    # Handle wrap-around
    span = (hi - lo) % 360.0
    if span == 0:
        raise ValueError("Upper and lower azimuth limits cannot be equal.")

    # Remaining lines = coefficients
    coeffs = np.array([float(x) for x in lines[1:]], dtype=float)
    N = len(coeffs)

    if N == 0:
        raise ValueError("File must contain coefficients after the first line.")

    # If only 1 coefficient -> uniform
    if N == 1:
        coeffs = np.ones(2)
        N = 2

    # Create azimuth bin edges
    edges = lo + np.linspace(0, span, N)

    # Compute piecewise-linear PDF area for normalization
    bin_width = span / (N - 1)
    area = 0.5 * bin_width * np.sum(coeffs[:-1] + coeffs[1:])
    if area <= 0:
        raise ValueError("Invalid distribution: total area must be positive.")

    pdf = coeffs / area  # normalized heights

    # Build CDF from area under piecewise-linear segments
    cdf = np.zeros(N)
    for i in range(1, N):
        cdf[i] = cdf[i-1] + 0.5 * (pdf[i] + pdf[i-1]) * bin_width

    # Normalize (may be tiny floating drift)
    cdf /= cdf[-1]

    return {
        "az_lo": lo,
        "az_hi": hi,
        "span": span,
        "edges_deg": edges,
        "pdf": pdf,
        "cdf": cdf,
        "N": N
    }

# --------------------------------------------------------------
# Inverse sampling from polygonal CDF
# --------------------------------------------------------------

def azi_distr_pick(u, dist):
    """
    u = uniform random number in [0,1]
    dist = dictionary returned by azi_distr_file_read
    """
    cdf = dist["cdf"]
    edges = dist["edges_deg"]
    pdf = dist["pdf"]
    N = dist["N"]
    span = dist["span"]

    # Exact 0 → lower limit
    if u <= 0:
        return edges[0]
    if u >= 1:
        return edges[-1]

    # Locate which segment the value falls in
    idx = np.searchsorted(cdf, u) - 1
    if idx < 0:
        idx = 0
    if idx >= N - 1:
        idx = N - 2

    # Segment info
    x0 = edges[idx]
    x1 = edges[idx + 1]
    y0 = pdf[idx]
    y1 = pdf[idx + 1]
    w = x1 - x0

    # Local CDF in this segment:
    # CDF(x) = CDF0 + ∫y dx = CDF0 + (y0 + y)/2 * dx
    CDF0 = cdf[idx]
    CDF1 = cdf[idx + 1]
    U = (u - CDF0)

    # Solve quadratic for linear PDF segment
    # PDF(x) = y0 + (y1 - y0)*(t),   t = (x - x0)/w
    # dCDF/dt = 0.5*w*(2*y0 + (y1 - y0)t)
    # but easier: use numerical local inversion

    # If y0 == y1 simple linear
    if abs(y1 - y0) < 1e-12:
        frac = U / ((y0) * w)
        return x0 + frac * w

    # Else solve quadratic:
    a = 0.5 * (y1 - y0)
    b = y0
    # Integrated: CDF increment = w*(a*t^2 + b*t)
    # Solve a t^2 + b t - U/w = 0
    A = a
    B = b
    C = -U / w

    disc = B*B - 4*A*C
    if disc < 0:
        disc = 0
    t = (-B + np.sqrt(disc)) / (2*A)

    if t < 0:
        t = 0
    if t > 1:
        t = 1

    return x0 + t * w
