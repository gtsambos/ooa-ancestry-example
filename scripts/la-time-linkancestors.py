#!/usr/bin/env python3

import tskit
import argparse
import time


# Constants
census_time = 2093.5
outfile = "output/times-linkancestors.txt"

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
ts.tables.link_ancestors(
    samples=ts.samples(),
    ancestors=[u.id for u in ts.nodes() if u.time == census_time])
time_end = time.time()

with open(outfile, "a") as f:
    print(time_end - time_start, file=f)

with open(f"tmp/la-{rep}.txt", "w") as f:
    print(time_end - time_start, file=f)