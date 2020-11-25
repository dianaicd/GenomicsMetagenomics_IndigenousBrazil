
# %%

import msprime
import tskit

# %%
# Simulating trees
tree_sequence = msprime.simulate(sample_size=6, Ne=1000)
tree = tree_sequence.first()
print(tree.draw(format="unicode"))


u = 2
while u != tskit.NULL:
    print("node {}: time = {}".format(u, tree.time(u)))
    u = tree.parent(u)

print(tree.branch_length(6))    
print(tree.total_branch_length)

# %%
# Recombination





# %%
# Multiple chromosomes
# recombination rate. Sikora et al. set it at 1.25e-8

rho = 1e-8
# positions in the recombination map
positions = [0, 1e8-1, 1e8, 2e8-1]
# recombination rate specified at a given position
rates = [rho, 0.5, rho, 0]

# The maximum number of non-recombining loci in the underlying simulation.
num_loci = int(positions[-1])

recombination_map = msprime.RecombinationMap(
    positions=positions, rates=rates, num_loci=num_loci)

tree_sequence = msprime.simulate(
    sample_size = 20, Ne = 100, recombination_map=recombination_map,
    model="dtwf",
    mutation_rate=2e-8)

# %%
# Variants
for variant in tree_sequence.variants():
     print(
        variant.site.id, variant.site.position,
        variant.alleles, variant.genotypes, sep="\t")
# %%
