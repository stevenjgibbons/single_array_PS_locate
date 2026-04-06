#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from datetime import timedelta
from arrtime_distr import arrtime_distr_file_read, arrtime_distr_pick


def main():
    if len(sys.argv) != 3:
        print("Usage: python test_arrtime_distr.py file.txt N")
        sys.exit(1)

    filename = sys.argv[1]
    N = int(sys.argv[2])

    # Read distribution
    dist = arrtime_distr_file_read(filename)

    # --- Generate samples ---
    us = np.random.rand(N)
    samples_dt = [arrtime_distr_pick(u, dist) for u in us]

    # Convert to seconds for histogram
    secs = np.array([(s - dist["t_lo"]).total_seconds() for s in samples_dt])

    # --- Histogram normalized by total counts (sum = 1) ---
    counts, bins = np.histogram(secs, bins=40)
    total = np.sum(counts)
    counts_norm = counts / total
    centers_sec = 0.5 * (bins[1:] + bins[:-1])
    bin_width = bins[1] - bins[0]

    # Convert bin centers back to datetime for plotting
    centers_dt = [dist["t_lo"] + timedelta(seconds=float(x)) for x in centers_sec]
    centers_dt_num = mdates.date2num(centers_dt)

    # --- Polygon PDF scaled to probability per bin ---
    edges_sec = dist["edges_sec"]
    pdf = dist["pdf"]
    pdf_scaled = pdf * bin_width

    edges_dt = [dist["t_lo"] + timedelta(seconds=float(e)) for e in edges_sec]
    edges_dt_num = mdates.date2num(edges_dt)

    # --- Plot ---
    plt.figure(figsize=(12,5))

    # Normalized histogram (counts / N)
    plt.bar(centers_dt_num, counts_norm, 
            width=(mdates.date2num(dist["t_lo"] + timedelta(seconds=bin_width)) -
                   mdates.date2num(dist["t_lo"])),
            alpha=0.6, label="Normalized histogram")

    # Scaled polygon PDF
    plt.plot(edges_dt_num, pdf_scaled, 'r-o', linewidth=2,
             label="Scaled polygon PDF")

    # Beautiful time labels
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S.%f"))
    plt.gcf().autofmt_xdate()

    plt.xlabel("Arrival Time")
    plt.ylabel("Probability per bin (counts / N)")
    plt.title(f"Arrival Time Distribution — {N} samples")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
