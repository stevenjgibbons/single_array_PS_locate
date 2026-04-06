#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# single_array_PS_locate.py
#
# Fully standalone Monte‑Carlo single‑array P–S locator
#
# Author: (generate name)
# Date: 2026-04-06
#
# --------------------------------------------------------------
# PART 1 — Imports, Time Utilities, Travel‑Time Table Reader
# --------------------------------------------------------------

import sys
import math
import numpy as np
from datetime import datetime, timedelta

# ==============================================================
# TIME UTILITIES
# ==============================================================

def parse_iso8601(s):
    """Parse ISO timestamp like 2020-06-22T09:20:01.500"""
    try:
        return datetime.fromisoformat(s)
    except Exception:
        raise ValueError(f"Could not parse timestamp: {s}")

def datetime_to_epoch_seconds(dt):
    """Convert datetime to seconds since Unix epoch."""
    epoch = datetime(1970,1,1)
    return (dt - epoch).total_seconds()

# ==============================================================
# TRAVEL‑TIME TABLE READER
# ==============================================================

def read_diff_tt(filename):
    """
    Reads table with columns:
    P1_travel_time S1_travel_time (S1-P1) distance_deg depth_km
    Returns four numpy arrays: ttp, tts, tts_minus_ttp, distdeg
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
            vals = [float(p) for p in parts[:4]]
            data.append(vals)

    arr = np.array(data, dtype=float)
    ttp = arr[:,0]
    tts = arr[:,1]
    tts_minus_ttp = arr[:,2]
    distdeg = arr[:,3]

    return ttp, tts, tts_minus_ttp, distdeg

# ==============================================================
# GENERAL INTERPOLATION (your gfdcfd + interp_1d)
# ==============================================================

def gfdcfd(X, XARR):
    """Generate finite difference interpolation weights."""
    XARR = np.array(XARR, dtype=float)
    NNDS = len(XARR)

    if len(np.unique(XARR)) != NNDS:
        raise ValueError("Grid points in XARR are not distinct.")

    h = XARR - X
    COEFM = np.zeros((NNDS, NNDS), dtype=float)

    for nder in range(NNDS):
        for inode in range(NNDS):
            if nder == 0:
                COEFM[inode, nder] = 1.0
            else:
                COEFM[inode, nder] = COEFM[inode, nder-1] * (h[inode] / nder)

    return np.linalg.inv(COEFM)

def interp_1d(x, xarr, yarr, n=6):
    """Finite‑difference interpolation."""
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
    w = M[0,:]
    return float(np.dot(w, ya))

# ==============================================================
# DISTANCE & TRAVEL‑TIME FROM S–P
# ==============================================================

def dist_from_table(S_minus_P, ttp, tts, tts_minus_ttp, distdeg):
    """Given S–P, interpolate distance, P_time, S_time."""
    dist = interp_1d(S_minus_P, tts_minus_ttp, distdeg)
    if dist is None:
        raise ValueError("S_minus_P out of travel‑time table range.")

    P_time = interp_1d(dist, distdeg, ttp)
    S_time = interp_1d(dist, distdeg, tts)

    return dist, P_time, S_time

# --------------------------------------------------------------
# PART 2 — Backazimuth & Arrival-Time Distributions
# --------------------------------------------------------------

# ==============================================================
# BACKAZIMUTH DISTRIBUTION
# ==============================================================

def azi_distr_file_read(filename):
    """
    Reads backazimuth distribution file of form:

        az_lo  az_hi        (two floats)
        c1
        c2
        ...
        cN

    If only 1 coefficient -> uniform.
    If N coefficients -> piecewise-linear polygon with N-1 segments.
    """

    with open(filename, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    # first line: azimuth range
    lo_str, hi_str = lines[0].split()
    az_lo = float(lo_str)
    az_hi = float(hi_str)

    # clockwise wrap-around: convert to span in degrees
    span = (az_hi - az_lo) % 360.0
    if span <= 0:
        raise ValueError("Backazimuth span must be > 0.")

    # coefficients
    coeffs = np.array([float(x) for x in lines[1:]], dtype=float)
    N = len(coeffs)
    if N == 0:
        raise ValueError("Backazimuth distribution requires coefficients.")

    if N == 1:
        coeffs = np.ones(2)
        N = 2

    # edges equally spaced in degrees
    edges_deg = az_lo + np.linspace(0, span, N)

    # piecewise-linear normalization
    bin_w = span / (N - 1)
    area = 0.5 * bin_w * np.sum(coeffs[:-1] + coeffs[1:])
    if area <= 0:
        raise ValueError("Backazimuth polygon has zero or negative area.")

    pdf = coeffs / area

    # build CDF
    cdf = np.zeros(N)
    for i in range(1, N):
        cdf[i] = cdf[i-1] + 0.5 * (pdf[i] + pdf[i-1]) * bin_w

    cdf /= cdf[-1]

    return {
        "az_lo": az_lo,
        "az_hi": az_hi,
        "span": span,
        "edges_deg": edges_deg,
        "pdf": pdf,
        "cdf": cdf,
        "N": N
    }


def azi_distr_pick(u, dist):
    """Inverse-CDF sampler for backazimuth distribution."""

    cdf = dist["cdf"]
    edges = dist["edges_deg"]
    pdf = dist["pdf"]
    N = dist["N"]
    span = dist["span"]

    if u <= 0:
        return edges[0]
    if u >= 1:
        return edges[-1]

    idx = np.searchsorted(cdf, u) - 1
    if idx < 0:
        idx = 0
    if idx >= N - 1:
        idx = N - 2

    x0 = edges[idx]
    x1 = edges[idx+1]
    y0 = pdf[idx]
    y1 = pdf[idx+1]
    w = x1 - x0

    C0 = cdf[idx]
    Ulocal = u - C0

    # flat
    if abs(y1 - y0) < 1e-12:
        frac = Ulocal / (y0 * w)
        return x0 + frac * w

    # quadratic solve
    a = 0.5 * (y1 - y0)
    b = y0
    A = a
    B = b
    C = -Ulocal / w

    disc = max(0.0, B*B - 4*A*C)
    t = (-B + math.sqrt(disc)) / (2*A)
    t = max(0.0, min(1.0, t))

    return x0 + t * w


# ==============================================================
# ARRIVAL-TIME DISTRIBUTION
# ==============================================================

def arrtime_distr_file_read(filename):
    """
    File format:

      t_lo   t_hi        (ISO-8601 timestamps)
      c1
      c2
      ...
      cN

    Returns dictionary containing:
      t_lo, t_hi (datetime)
      edges_sec (float seconds from t_lo)
      pdf, cdf
    """

    with open(filename, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    # first line: time range
    t_lo_str, t_hi_str = lines[0].split()
    t_lo = parse_iso8601(t_lo_str)
    t_hi = parse_iso8601(t_hi_str)
    dt = (t_hi - t_lo).total_seconds()
    if dt <= 0:
        raise ValueError("Arrival-time upper bound must be after lower bound.")

    coeffs = np.array([float(x) for x in lines[1:]], dtype=float)
    N = len(coeffs)
    if N == 0:
        raise ValueError("Arrival-time distribution requires coefficients.")
    if N == 1:
        coeffs = np.ones(2)
        N = 2

    edges_sec = np.linspace(0.0, dt, N)
    bin_w = dt / (N - 1)

    area = 0.5 * bin_w * np.sum(coeffs[:-1] + coeffs[1:])
    if area <= 0:
        raise ValueError("Arrival-time polygon has zero or negative area.")

    pdf = coeffs / area

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


def arrtime_distr_pick(u, dist):
    """Sample a datetime from a polygonal arrival-time distribution."""

    cdf = dist["cdf"]
    edges = dist["edges_sec"]
    pdf = dist["pdf"]
    N = dist["N"]
    t_lo = dist["t_lo"]

    if u <= 0:
        return t_lo
    if u >= 1:
        return dist["t_hi"]

    idx = np.searchsorted(cdf, u) - 1
    if idx < 0:
        idx = 0
    if idx >= N - 1:
        idx = N - 2

    x0 = edges[idx]
    x1 = edges[idx+1]
    y0 = pdf[idx]
    y1 = pdf[idx+1]
    w = x1 - x0

    C0 = cdf[idx]
    Ulocal = u - C0

    if abs(y1 - y0) < 1e-12:
        frac = Ulocal / (y0 * w)
        return t_lo + timedelta(seconds=(x0 + frac * w))

    a = 0.5 * (y1 - y0)
    b = y0
    A = a
    B = b
    Cc = -Ulocal / w

    disc = max(0.0, B*B - 4*A*Cc)
    t = (-B + math.sqrt(disc)) / (2*A)
    t = max(0.0, min(1.0, t))

    x = x0 + t * w
    return t_lo + timedelta(seconds=x)


# --------------------------------------------------------------
# PART 3 — Geodesy + Core Monte‑Carlo Event Locator
# --------------------------------------------------------------

# ==============================================================
# GEODESY — GREAT-CIRCLE FORWARD PROJECTION
# ==============================================================

def project_from_station(statlat, statlon, dist_deg, backazimuth_deg):
    """
    Forward great-circle projection.

    Inputs:
      statlat, statlon       degrees
      dist_deg               great-circle distance (degrees)
      backazimuth_deg        Convention A:
                              0 = North, 90 = East, 180 = South, 270 = West
                              Event lies outward along this direction.

    Returns:
      evlat, evlon           degrees
    """

    # Convert degrees → radians
    lat1 = math.radians(statlat)
    lon1 = math.radians(statlon)
    baz  = math.radians(backazimuth_deg)

    # Convert angular distance to radians
    ang_dist = math.radians(dist_deg)

    # Spherical forward projection
    sin_lat2 = math.sin(lat1)*math.cos(ang_dist) + \
               math.cos(lat1)*math.sin(ang_dist)*math.cos(baz)
    lat2 = math.asin(sin_lat2)

    y = math.sin(baz)*math.sin(ang_dist)*math.cos(lat1)
    x = math.cos(ang_dist) - math.sin(lat1)*sin_lat2
    lon2 = lon1 + math.atan2(y, x)

    # Normalize longitude
    lon2 = (lon2 + math.pi) % (2*math.pi) - math.pi

    return math.degrees(lat2), math.degrees(lon2)


# ==============================================================
# GEOCENTRIC DISTANCE → KM
# ==============================================================

def dist_km_from_deg(dist_deg):
    """Convert geodesic angular distance (deg) to km.
       Uses Earth's mean radius 6371 km."""
    return dist_deg * (math.pi/180.0) * 6371.0


# ==============================================================
# CORE MONTE-CARLO EVENT LOCATION LOOP
# ==============================================================

def run_monte_carlo(
        N,
        statlat, statlon,
        OMPcorrectP1, OMPcorrectS1,
        ttp, tts, tts_minus_ttp, distdeg,
        bazidist,
        P1_arrtime_dist,
        S1_arrtime_dist):

    """
    Perform N independent location trials.

    Returns dictionaries containing arrays:
       evlon[], evlat[], evOrigTime_dt[], evOrig_secs[],
       distdeg_arr[], distkm_arr[],
       backaz[], P1TimeMeasured[], S1TimeMeasured[]
    """

    # Outputs
    evlons = np.zeros(N)
    evlats = np.zeros(N)
    evOrigTimes_dt = [None]*N
    evOrigTimes_sec = np.zeros(N)
    distdeg_arr = np.zeros(N)
    distkm_arr = np.zeros(N)
    backaz_arr = np.zeros(N)
    P1_meas_list = [None]*N
    S1_meas_list = [None]*N

    for i in range(N):

        # ----------------------------------------------------------
        # (1) Sample backazimuth
        # ----------------------------------------------------------
        u_bazi = np.random.rand()
        backaz = azi_distr_pick(u_bazi, bazidist)
        backaz_arr[i] = backaz

        # ----------------------------------------------------------
        # (2) Sample measured P1 arrival-time
        # ----------------------------------------------------------
        u_p = np.random.rand()
        P1_meas = arrtime_distr_pick(u_p, P1_arrtime_dist)
        P1_meas_list[i] = P1_meas

        # ----------------------------------------------------------
        # (3) Sample measured S1 arrival-time
        # ----------------------------------------------------------
        u_s = np.random.rand()
        S1_meas = arrtime_distr_pick(u_s, S1_arrtime_dist)
        S1_meas_list[i] = S1_meas

        # ----------------------------------------------------------
        # (4) Correct P1
        # ----------------------------------------------------------
        P1_corr = P1_meas - timedelta(seconds=OMPcorrectP1)

        # ----------------------------------------------------------
        # (5) Correct S1
        # ----------------------------------------------------------
        S1_corr = S1_meas - timedelta(seconds=OMPcorrectS1)

        # ----------------------------------------------------------
        # (6) S_minus_P
        # ----------------------------------------------------------
        S_minus_P = (S1_corr - P1_corr).total_seconds()

        # ----------------------------------------------------------
        # (7) Sanity check
        # ----------------------------------------------------------
        if S_minus_P <= 0 or S_minus_P > 300:
            print("WARNING: S_minus_P out of expected range.")
            print("P1_corr =", P1_corr.isoformat())
            print("S1_corr =", S1_corr.isoformat())
            # Continue, but skip this iteration or clamp?
            # Here we choose to skip:
            continue

        # ----------------------------------------------------------
        # (8) Distance and theoretical P,S travel times
        # ----------------------------------------------------------
        try:
            dist_in_deg, P_th, S_th = dist_from_table(
                S_minus_P, ttp, tts, tts_minus_ttp, distdeg)
        except Exception as e:
            print(f"Interpolation failure at trial {i}: {e}")
            continue

        distdeg_arr[i] = dist_in_deg
        distkm_arr[i] = dist_km_from_deg(dist_in_deg)

        # ----------------------------------------------------------
        # (9) Project event location from station
        # ----------------------------------------------------------
        evlat, evlon = project_from_station(
            statlat, statlon, dist_in_deg, backaz)

        evlats[i] = evlat
        evlons[i] = evlon

        # ----------------------------------------------------------
        # (10) Origin time = corrected P1 time - theoretical P travel time
        # ----------------------------------------------------------
        evOrig_dt = P1_corr - timedelta(seconds=P_th)
        evOrigTimes_dt[i] = evOrig_dt
        evOrigTimes_sec[i] = datetime_to_epoch_seconds(evOrig_dt)

    # end Monte‑Carlo loop

    return {
        "evlon": evlons,
        "evlat": evlats,
        "evOrig_dt": evOrigTimes_dt,
        "evOrig_sec": evOrigTimes_sec,
        "distdeg": distdeg_arr,
        "distkm": distkm_arr,
        "backaz": backaz_arr,
        "P1_meas": P1_meas_list,
        "S1_meas": S1_meas_list
    }

# --------------------------------------------------------------
# PART 4 — Output File Writing + Statistics + 95% Uncertainty Ellipse
# --------------------------------------------------------------

# ==============================================================
# WRITE OUTPUT FILE (ONE LINE PER SAMPLE)
# ==============================================================

def format_datetime_for_output(dt):
    """
    Return datetime in format:
    yyyy-mm-ddThh.mm.ss.ssss
    """
    return dt.strftime("%Y-%m-%dT%H.%M.%S.%f")[:-2]  # keep 4 decimals


def write_output_file(outfile, results, statlat, statlon):
    """
    Writes one line per Monte-Carlo sample with formatted fields:

      i, backazimuth, P1TimeMeasured, S1TimeMeasured,
      distdegrees, distkm, evlon, evlat,
      evOrigTime (formatted), evOrigTime (epoch seconds)
    """

    with open(outfile, "w") as f:

        for i in range(len(results["evlon"])):

            evlon = results["evlon"][i]
            evlat = results["evlat"][i]
            backaz = results["backaz"][i]
            distdeg = results["distdeg"][i]
            distkm = results["distkm"][i]
            P1_meas = results["P1_meas"][i]
            S1_meas = results["S1_meas"][i]
            evOrig_dt = results["evOrig_dt"][i]
            evOrig_sec = results["evOrig_sec"][i]

            # skip blanks (e.g., interpolation failures)
            if P1_meas is None or S1_meas is None or evOrig_dt is None:
                continue

            # format time strings
            P1_str = format_datetime_for_output(P1_meas)
            S1_str = format_datetime_for_output(S1_meas)
            evOrig_str = format_datetime_for_output(evOrig_dt)

            f.write(
                f"{i:7d} "
                f"{backaz:7.3f} "
                f"{P1_str} "
                f"{S1_str} "
                f"{distdeg:7.4f} "
                f"{distkm:8.3f} "
                f"{evlon:9.4f} "
                f"{evlat:8.4f} "
                f"{evOrig_str} "
                f"{evOrig_sec:12.4f} "
                f"{distdeg:7.4f} "
                f"{distkm:8.3f}\n"
            )


# ==============================================================
# STATISTICS
# ==============================================================

def compute_basic_statistics(arr):
    """Return (median, mean, stddev) for a 1D array."""
    arr = np.array(arr, dtype=float)
    med = float(np.median(arr))
    mean = float(np.mean(arr))
    std = float(np.std(arr))
    return med, mean, std


# ==============================================================
# 95% CONFIDENCE ELLIPSE FROM COVARIANCE
# ==============================================================

def compute_uncertainty_ellipse(lons, lats):
    """
    Given arrays of longitudes & latitudes, compute:

      - Covariance matrix
      - Eigenvalues/eigenvectors
      - 95% confidence ellipse:
           semi-major axis
           semi-minor axis
           strike/azimuth of major axis

    Uses 2D Gaussian approximation with chi-square scaling:
      χ²(2 df, 95%) = 5.991
    """

    X = np.column_stack((lons, lats))

    # Covariance matrix
    C = np.cov(X, rowvar=False)

    # Eigen decomposition (sorted largest → smallest)
    eigvals, eigvecs = np.linalg.eigh(C)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Major/minor axis lengths (95% conf)
    chi2_95 = 5.991
    major_axis = math.sqrt(eigvals[0] * chi2_95)
    minor_axis = math.sqrt(eigvals[1] * chi2_95)

    # Orientation of major axis
    vx, vy = eigvecs[:,0]
    angle = math.degrees(math.atan2(vy, vx))  # degrees from East CCW

    return major_axis, minor_axis, angle


# ==============================================================
# WRITE STATISTICS OUTPUT
# ==============================================================

def write_statistics(stats_outfile, results):
    """
    Write median, mean, stddev for:
      evlon, evlat, evOrigTime
    plus 95% confidence ellipse for (lon, lat).
    """

    evlon = results["evlon"]
    evlat = results["evlat"]
    evOrig = results["evOrig_sec"]

    # Basic statistics
    med_lon, mean_lon, std_lon = compute_basic_statistics(evlon)
    med_lat, mean_lat, std_lat = compute_basic_statistics(evlat)
    med_ot,  mean_ot,  std_ot  = compute_basic_statistics(evOrig)

    # Uncertainty ellipse
    major, minor, angle = compute_uncertainty_ellipse(evlon, evlat)

    with open(stats_outfile, "w") as f:
        f.write("STATISTICS OF EVENT LOCATIONS\n")
        f.write("---------------------------------------\n")
        f.write(f"median(evlon)      = {med_lon:.6f}\n")
        f.write(f"median(evlat)      = {med_lat:.6f}\n")
        f.write(f"median(evOrigTime) = {med_ot:.6f}  (epoch seconds)\n\n")

        f.write(f"mean(evlon)        = {mean_lon:.6f}\n")
        f.write(f"mean(evlat)        = {mean_lat:.6f}\n")
        f.write(f"mean(evOrigTime)   = {mean_ot:.6f}\n\n")

        f.write(f"std(evlon)         = {std_lon:.6f}\n")
        f.write(f"std(evlat)         = {std_lat:.6f}\n")
        f.write(f"std(evOrigTime)    = {std_ot:.6f}\n\n")

        f.write("95% CONFIDENCE ELLIPSE (Gaussian)\n")
        f.write("---------------------------------------\n")
        f.write(f"Semi-major axis (deg): {major:.6f}\n")
        f.write(f"Semi-minor axis (deg): {minor:.6f}\n")
        f.write(f"Ellipse azimuth (deg): {angle:.3f} (CCW from East)\n")

# --------------------------------------------------------------
# PART 5 — Command Line Interface + Main Program
# --------------------------------------------------------------

def parse_arguments(argv):
    """
    Expected arguments:

    statlat             (required)
    statlon             (required)
    OMPcorrectP1        (default = 0.0)
    OMPcorrectS1        (default = 0.0)
    TraveltimeTableFile
    BaziDistFile
    P1ArrtimeDistFile
    S1ArrtimeDistFile
    N                   (default = 10000)
    outfile             (default = single_array_PS_locate_output.txt)
    stats_outfile       (default = single_array_PS_locate_stats_output.txt)
    """

    if len(argv) < 9:
        print("Usage:")
        print("  python single_array_PS_locate.py "
              "statlat statlon OMPcorrectP1 OMPcorrectS1 "
              "TraveltimeTableFile BaziDistFile P1ArrtimeDistFile "
              "S1ArrtimeDistFile [N=10000] "
              "[outfile=single_array_PS_locate_output.txt] "
              "[stats_outfile=single_array_PS_locate_stats_output.txt]")
        sys.exit(1)

    statlat = float(argv[1])
    statlon = float(argv[2])
    OMPcorrectP1 = float(argv[3])
    OMPcorrectS1 = float(argv[4])

    TraveltimeTableFile = argv[5]
    BaziDistFile        = argv[6]
    P1ArrFile           = argv[7]
    S1ArrFile           = argv[8]

    # defaults
    if len(argv) > 9:
        N = int(argv[9])
    else:
        N = 10000

    if len(argv) > 10:
        outfile = argv[10]
    else:
        outfile = "single_array_PS_locate_output.txt"

    if len(argv) > 11:
        stats_outfile = argv[11]
    else:
        stats_outfile = "single_array_PS_locate_stats_output.txt"

    return (statlat, statlon,
            OMPcorrectP1, OMPcorrectS1,
            TraveltimeTableFile, BaziDistFile,
            P1ArrFile, S1ArrFile,
            N, outfile, stats_outfile)


def main(argv):

    # --------------------------------------------------------------
    # Parse arguments
    # --------------------------------------------------------------
    (statlat, statlon,
     OMPcorrectP1, OMPcorrectS1,
     TravTable, BaziFile,
     P1ArrFile, S1ArrFile,
     N, outfile, stats_outfile) = parse_arguments(argv)

    print("\n--- Single Array P–S Monte Carlo Locator ---")
    print(f"Station: lat={statlat}, lon={statlon}")
    print(f"N = {N}")
    print()

    # --------------------------------------------------------------
    # Read travel‑time table
    # --------------------------------------------------------------
    print("Reading travel‑time table...")
    ttp, tts, tts_minus_ttp, distdeg = read_diff_tt(TravTable)

    # --------------------------------------------------------------
    # Read backazimuth distribution
    # --------------------------------------------------------------
    print("Reading backazimuth distribution...")
    bazidist = azi_distr_file_read(BaziFile)

    # --------------------------------------------------------------
    # Read arrival-time distributions
    # --------------------------------------------------------------
    print("Reading P1 arrival-time distribution...")
    P1_arr = arrtime_distr_file_read(P1ArrFile)

    print("Reading S1 arrival-time distribution...")
    S1_arr = arrtime_distr_file_read(S1ArrFile)

    # --------------------------------------------------------------
    # Run Monte‑Carlo location
    # --------------------------------------------------------------
    print("Running Monte‑Carlo event location...")
    results = run_monte_carlo(
        N,
        statlat, statlon,
        OMPcorrectP1, OMPcorrectS1,
        ttp, tts, tts_minus_ttp, distdeg,
        bazidist,
        P1_arr,
        S1_arr
    )

    # --------------------------------------------------------------
    # Write full output file
    # --------------------------------------------------------------
    print(f"Writing detailed output: {outfile}")
    write_output_file(outfile, results, statlat, statlon)

    # --------------------------------------------------------------
    # Write summary statistics
    # --------------------------------------------------------------
    print(f"Writing statistics file: {stats_outfile}")
    write_statistics(stats_outfile, results)

    print("\nDone.\n")


if __name__ == "__main__":
    main(sys.argv)
