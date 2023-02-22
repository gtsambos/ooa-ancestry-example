#!/usr/bin/env python3

import tskit
import argparse
import time

# Constants
census_time = 2093.5
outfile = "output/times-squash.txt"

# Import the correct input
parser = argparse.ArgumentParser()
parser.add_argument(
    '--rep', type=int, dest='rep', nargs=1,
    help="The replicate number.")
args = parser.parse_args()
rep = args.rep[0]

ts = tskit.load(f"output/ts-{rep}.trees")

# Run tree-by-tree method 
la_output = {}
for s in ts.samples():
    la_output[s] = {
        'right_bp' : [],
        'ancestor' : []
    }

time_start = time.time()
for t in ts.trees():
    # Get the set of census nodes in this tree.
    for n in t.nodes():
        if t.time(n) == census_time:
            samps = t.samples(n)
            for s in samps:
                prev_a = la_output[s]['ancestor']
                if len(prev_a) > 0 and prev_a[-1] == n:
                    la_output[s]['right_bp'][-1] = t.interval.right
                else:
                    la_output[s]['ancestor'].append(n)
                    la_output[s]['right_bp'].append(t.interval.right)
time_end = time.time()

with open(outfile, "a") as f:
    print(time_end - time_start, file=f)

with open(f"tmp/sq-{rep}", "w") as f:
    print(time_end - time_start, file=f)