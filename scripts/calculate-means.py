#!/usr/bin/env python3

import numpy as np

# Read in files.
for infile in [
    # "output/times-baseline.txt",
    # "output/times-baseline-numba.txt",
    "output/times-squash.txt",
    "output/times-squash-numba.txt",
    "output/times-linkancestors.txt"
]:
    times = np.loadtxt(infile)
    mn = np.mean(times)
    sd = np.std(times)
    
    with open("output/mean-times.txt", "a") as f:
        print(f"{infile}: mean {mn}, sd {sd}", file=f)