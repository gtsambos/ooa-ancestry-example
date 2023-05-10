#!/usr/bin/env python3

import tskit
import argparse
import time
import os

# Constants
census_time = 100 
outfile = "output/times-baseline.txt"

# Import the correct input
parser = argparse.ArgumentParser()
parser.add_argument(
    '--rep', type=int, dest='rep', nargs=1,
    help="The replicate number.")
args = parser.parse_args()
rep = args.rep[0]

ts = tskit.load(f"output/ts-{rep}.trees")

# Run tree-by-tree method 
l = []
r = []
a = []
c = []
time_start = time.time()
for t in ts.trees():
    # Get the set of census nodes in this tree.
    for n in t.nodes():
        if t.time(n) == census_time:
            samps = t.samples(n)
            for s in samps:
                c.append(s)
                a.append(n)
                l.append(t.interval.left)
                r.append(t.interval.right)
time_end = time.time()

with open(outfile, "a") as f:
    print(time_end - time_start, file=f)

# make the tmp directory if it doesn't exist
if not os.path.exists("tmp"):
    os.makedirs("tmp")

with open(f"tmp/base-{rep}.txt", "w") as f:
    print(time_end - time_start, file=f)