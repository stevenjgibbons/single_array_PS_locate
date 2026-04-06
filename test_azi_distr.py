import sys
import numpy as np
import matplotlib.pyplot as plt

from azi_distr import azi_distr_file_read, azi_distr_pick

def main():
    if len(sys.argv) != 3:
        print("Usage: python test_azi_distr.py file.txt  N")
        sys.exit(1)

    filename = sys.argv[1]
    N = int(sys.argv[2])

    dist = azi_distr_file_read(filename)

    # Generate N random samples
    us = np.random.rand(N)
    samples = np.array([azi_distr_pick(u, dist) for u in us])

    # Plot histogram
    plt.figure(figsize=(8,5))
    plt.hist(samples, bins=40, density=True, alpha=0.6, label="Sample histogram")

    # Plot the input polygon for comparison
    centers = dist["edges_deg"]
    plt.plot(centers, dist["pdf"], 'r-o', label="Input polygon PDF")

    plt.xlabel("Backazimuth (deg)")
    plt.ylabel("Probability density")
    plt.title(f"Histogram of {N} sampled backazimuths")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
