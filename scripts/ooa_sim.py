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
census_time = 500
# census_time = 50
num_inds = 50

# # Uncomment to plot
# _, ax = plt.subplots(1, 1, figsize=(4.5, 4), tight_layout=True)
# ax = demesdraw.tubes(graph, ax=ax, log_time=True)

# Convert to msprime format and add a census time just before OOA
msp_dem = msprime.Demography.from_demes(graph)
msp_dem.add_census(time=census_time)
# msp_dem.add_census(time=census_time2)
msp_dem.sort_events()

# Simulate with msprime.
print("Simulating with msprime...")
time_start = time.time()
ts = msprime.sim_ancestry(
    samples={"YRI": num_inds, "CEU": num_inds, "CHB" : num_inds},
    demography=msp_dem,
    random_seed=1016,
    sequence_length=1e6,
    # sequence_length=51304566, # stdpopsim default for chr22
    recombination_rate=1.4445e-08, # stdpopsim default for chr1
    end_time=501
)
time_end = time.time()
print("Total simulation time: ", time_end - time_start)

# Save the ts in case you want to look at later.
ts.dump("ooa_50inds_500c-stop501.trees")

# # Check there are no coalescences more recent than OOA.
# for t in ts.trees():
#     if t.time(t.root) < census_time:
#         print("Warning: some MRCAs post-date the census time!")
#         break

# # Run link-ancestors
# print("Timing link_ancestors (ancient census)...")
# time_start = time.time()
# ts.tables.link_ancestors(
#     samples=ts.samples(),
#     ancestors=[u.id for u in ts.nodes() if u.time == census_time])
# time_end = time.time()
# print("Time running link_ancestors (ancient census): ", time_end - time_start)

# # Run link-ancestors
# print("Timing link_ancestors (recent census)...")
# time_start = time.time()
# ts.tables.link_ancestors(
#     samples=ts.samples(),
#     ancestors=[u.id for u in ts.nodes() if u.time == census_time2])
# time_end = time.time()
# print("Time running link_ancestors (recent census): ", time_end - time_start)


# # make a PopAncestry object
# print("Timing get_pop_ancestry...")
# time_start = time.time()
# pa = tspop.get_pop_ancestry(ts, census_time=census_time)
# time_end = time.time()
# print("Time running get_pop_ancestry: ", time_end - time_start)

# # Plot the first CEU individual.
# pair = [ts.samples(3)[0], ts.samples(3)[1]]
# pa.plot_karyotypes(
#     sample_pair= pair,
#     title="Hominin ancestry in an Homo Sapiens individual",
#     length_in_Mb=True,
#     pop_labels=['HomSap', 'HomNea', 'ArcAfr'],
#     height=4,
#     width=13)

