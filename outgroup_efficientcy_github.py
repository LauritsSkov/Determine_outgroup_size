import msprime
import demes
import demesdraw
import sys
import matplotlib.pyplot as plt
import numpy as np
import argparse

# -----------------------------------------------------------------------------------------------------
# Parameters for demography (plot with demes)
# -----------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-demography", metavar='',help="File with demography", type=str, default = 'Demography_a.yaml')
parser.add_argument("-iterations", metavar='',help="Number of iterations", type=int, default = 100)
parser.add_argument("-chrom_size", metavar='',help="chromosome size", type=int, default = 50_000)
parser.add_argument("-outgroup_size", metavar='',help="n individual per outgroup.", type=int, default = 100)
parser.add_argument("-outfile", metavar='',help="outplot name", type=str, default = 'test.pdf')

args = parser.parse_args()



# Parameters
CHROM_SIZE = args.chrom_size
n_windows = int(CHROM_SIZE/1000)
mutation_rate = 1.45e-8 
recombination_rate = 1.45e-8 

# Plot demography
MAX_Y = 60_000 
MODEL_NAME = args.demography
graph = demes.load(f"{MODEL_NAME}")
demography = msprime.Demography.from_demes(graph)

INGROUP_SPLIT = None
for event in demography.events:
    if type(event) == msprime.demography.PopulationSplit:
        if event.derived[0] == 'ingroup':
            INGROUP_SPLIT = event.time


OUTGROUPS = []
COLORS = {}

for pop in demography.populations:
    if pop.name.startswith('out'):
        OUTGROUPS.append(pop.name)  

    if not (pop.name.startswith('out') or pop.name.startswith('ingroup')):
        sys.exit('Population names must start be EITHER out1, out2... or ingroup')

    COLORS[pop.name] = "#999999"

w = demesdraw.utils.separation_heuristic(graph)
x_positions  = {}
for index, outgroup in enumerate(OUTGROUPS):
    x_positions[outgroup] = index * w
x_positions['ingroup'] = len(OUTGROUPS) * w

plt.figure(figsize=(5,15))
fig, ax = plt.subplots(1, 2)
demesdraw.tubes(graph, ax=ax[0], seed=1, colours = COLORS, max_time=MAX_Y, positions=x_positions)
ax[0].set_ylim(0, MAX_Y)
ax[0].set_title(f'Demography')


Outgroup_sizes = []
percent_removed = []
emissions = []


print( 'n_outgroup', 'removed', 'ingroup', 'total snps', 'removed_snps', 'genome size', 'obs_mean', 'obs_variance', sep = '\t')

for outgroup_size in range(1, args.outgroup_size):

    SAMPLES = []
    SAMPLES.append(msprime.SampleSet(1, population='ingroup')) 
    for outgroup in OUTGROUPS:
        SAMPLES.append(msprime.SampleSet(outgroup_size, population=outgroup)) 
   
    removed = 0
    ingroup = 0

    obs = np.zeros(args.iterations * n_windows )

    for iteration in range(args.iterations):

        offset = iteration * n_windows

        ts = msprime.sim_ancestry(
            samples= SAMPLES, 
            demography=demography,
            sequence_length=CHROM_SIZE,
            recombination_rate=recombination_rate,
            random_seed = outgroup_size + iteration + 1)

        # mutations
        mts = msprime.sim_mutations(ts, rate=mutation_rate)
        for var in mts.variants():

            if np.sum(var.genotypes[0]) > 0:
                if np.sum(var.genotypes[2:]) < 1:
                    if var.site.mutations[0].time > INGROUP_SPLIT:
                        removed += 1
                    else:
                        ingroup += 1

                    position = int(var.site.position)
                    rounded_pos = int((position - position % 1000)/1000)
                    obs[offset + rounded_pos] += 1
                

    Outgroup_sizes.append(outgroup_size * len(OUTGROUPS))

    fraction = ingroup/(removed + ingroup) if removed + ingroup > 0 else 0
    percent_removed.append(round(fraction * 100,2))

    obs_mean = np.mean(obs)
    obs_var = np.var(obs)
    obs_ratio = obs_var / obs_mean

    print(outgroup_size * len(OUTGROUPS), removed, ingroup, (removed + ingroup), round(fraction * 100,2), CHROM_SIZE * args.iterations, obs_mean, obs_var, sep = '\t')
    single_emission = round(obs_mean, 4)


# set grid
ax[1].set_yticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax[1].grid(True, axis='y', linestyle='--', alpha=0.7)
ax[1].set_title(f'% removed')

# plot amount of missing SNPs
ax[1].plot(Outgroup_sizes, percent_removed,  color='black', linestyle='-', label=f"emission={single_emission}")

# legend
ax[1].legend(loc='upper left', bbox_to_anchor=(0, 1), frameon=False)
ax[1].set_ylim(0, 100)

plt.savefig(args.outfile) 


