#!/usr/bin/env python3

import tskit
import argparse
import time
import numba

# Constants
census_time = 100 
outfile = "output/times-squash-numba.txt"

# Import the correct input
parser = argparse.ArgumentParser()
parser.add_argument(
    '--rep', type=int, dest='rep', nargs=1,
    help="The replicate number.")
args = parser.parse_args()
rep = args.rep[0]

ts = tskit.load(f"output/ts-{rep}.trees")

# Define numba recursive function with parent_array.
@numba.jit
def find_census_ancestor(sample, parent_arr, times, censustime):
    v = sample
    while times[v] != censustime:
        v = parent_arr[v]
        if v == -1:
            break
    return v

# Get the set of all node times.
node_times = ts.tables.nodes.time

# Run tree-by-tree method 
la_output = {}
sa = ts.samples()
for s in sa:
    la_output[s] = {
        'right_bp' : [],
        'ancestor' : []
    }

time_start = time.time()
for t in ts.trees():
    pa = t.parent_array
    for s in sa:
        ca = find_census_ancestor(s, pa, node_times, census_time)

        # Squash the output.
        prev_a = la_output[s]['ancestor']
        if len(prev_a) > 0 and prev_a[-1] == ca:
            la_output[s]['right_bp'][-1] = t.interval.right
        else:
            la_output[s]['ancestor'].append(ca)
            la_output[s]['right_bp'].append(t.interval.right)
time_end = time.time()

with open(outfile, "a") as f:
    print(time_end - time_start, file=f)

with open(f"tmp/sq_nm-{rep}.txt", "w") as f:
    print(time_end - time_start, file=f)