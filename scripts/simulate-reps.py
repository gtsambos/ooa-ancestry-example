#!/usr/bin/env python3

import demes
import stdpopsim, msprime
import tskit
import time
import argparse
import os

# Import the replicate number
parser = argparse.ArgumentParser()
parser.add_argument(
    '--rep', type=int, dest='rep', nargs=1,
    help="The replicate number.")
args = parser.parse_args()
rep = args.rep[0]

# Make the output directory if it doesn't exist.
if not os.path.exists("output"):
    os.makedirs("output")
outfile = "output/times-sim.txt"

# Import the American admixture model from stdpopsim and convert to demes format.
species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model('AmericanAdmixture_4B11')
graph = model.model.to_demes()
census_times = [100]
num_inds = 1000

# Convert to msprime format and add a census time just before OOA
msp_dem = msprime.Demography.from_demes(graph)
msp_dem.add_census(time=census_times[0])
msp_dem.sort_events()

# Simulate with msprime.
print("Simulating with msprime...")

time_start = time.time()
ts = msprime.sim_ancestry(
    samples={"ADMIX": num_inds},
    demography=msp_dem,
    random_seed=1016,
    sequence_length=248956422, # stdpopsim default for chr1
    recombination_rate=1.15235e-08 # stdpopsim default for chr1
)
time_end = time.time()

# Print runtime to file.
with open(outfile, "a") as f:
    print(time_end - time_start, file=f)

# Output tree file.
ts.dump(f"output/ts-{rep}.trees")

# print("Total simulation time: ", time_end - time_start)
