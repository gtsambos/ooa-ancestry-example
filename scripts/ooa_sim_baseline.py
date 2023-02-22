#!/usr/bin/env python3

import demes
# import demesdraw
import tspop
import stdpopsim
import msprime, tskit
import matplotlib.pyplot as plt
import time

# Import the OOA_archaic model from stdpopsim and convert to demes format.
species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model('OutOfAfricaArchaicAdmixture_5R19')
graph = model.model.to_demes()
census_time = 2093.5
num_inds = 100

# Convert to msprime format and add a census time just before OOA
msp_dem = msprime.Demography.from_demes(graph)
msp_dem.add_census(time=census_time)
msp_dem.sort_events()

# Simulate with msprime.
print("Simulating with msprime...")
time_start = time.time()
ts = msprime.sim_ancestry(
    samples={"YRI": num_inds, "CEU": num_inds, "CHB" : num_inds},
    demography=msp_dem,
    random_seed=1016,
    sequence_length=248956422, # stdpopsim default for chr22
    recombination_rate=1.15235e-08 # stdpopsim default for chr1
)
time_end = time.time()
print("Total simulation time: ", time_end - time_start)

# # Save the ts in case you want to look at later.
# ts.dump("ooa.trees")

# Check there are no coalescences more recent than OOA.
for t in ts.trees():
    if t.time(t.root) < census_time:
        print("Warning: some MRCAs post-date the census time!")
        break

# Run tree-by-tree method 
print("Timing naive tree-by-tree method")
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
print("Time running naive-tree-by-tree method", time_end - time_start)
