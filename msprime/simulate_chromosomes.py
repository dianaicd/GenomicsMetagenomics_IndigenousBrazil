# %%
import msprime, tskit, random, getopt, sys, numpy as np

# # %%
# chr = 19
# test_path = '/Users/dcruz/Projects/Botocudos/Files/test/msprime/'
# positions_path = test_path + 'positions_19.bim'
# recomb_map_path = '/Users/dcruz/Projects/Botocudos/Files/HapMap_genetic_maps/HapmapII_GRCh37_RecombinationHotspots/genetic_map_GRCh37_chr19.txt'
# vcf_path_out = test_path + 'chr19.vcf'

# %%
def load_recomb_map_chr(recomb_map_path):
    recomb_map = msprime.RecombinationMap.read_hapmap(recomb_map_path)
    return(recomb_map)
# %%

# %%
def load_positions(positions_path):
    with open(positions_path, 'r') as positions_fd:
        positions = {int(line.split()[1].replace('\n', '')):True for line in positions_fd.readlines()}
    return positions

# %%
def simulate_data_chr(Ne, sample_size, recombination_map, 
    mutation_rate=1.25e-8, model='dtwf', num_replicates=1):
    tree_sequence = msprime.simulate(
        sample_size=sample_size, 
        Ne=Ne, 
        recombination_map=recombination_map,
        #model=model,
        mutation_rate=mutation_rate,
        num_replicates=num_replicates)

    return tree_sequence

# %%
def adjust_number_sites(tree_sequence, positions):
    total_num_positions = len(positions.keys())
    num_old_positions = tree_sequence.get_num_sites()
    num_missing_positions = int(total_num_positions - num_old_positions)
    if num_missing_positions > 0:
        old_positions = [ np.round(site.position) for site in tree_sequence.sites() ]
        new_tree_sequence = add_variants(tree_sequence = tree_sequence,
                                        num_missing_positions = num_missing_positions,
                                        old_positions = old_positions)
    elif num_missing_positions < 0:
        old_positions = [ site.position for site in tree_sequence.sites() ]
        new_tree_sequence = remove_variants(tree_sequence = tree_sequence, 
                                            old_positions = old_positions, 
                                            total_num_positions = total_num_positions)
    else:
        new_tree_sequence = tree_sequence
    return new_tree_sequence

# %%
def add_variants(tree_sequence, num_missing_positions, old_positions):
    tables = tree_sequence.dump_tables()

    possible_positions = set(range(1, int(tree_sequence.get_sequence_length()))).difference(old_positions)
    new_positions = random.sample(possible_positions, num_missing_positions)

    for position in new_positions:
        site_id = tables.sites.add_row(position = position, ancestral_state = '0')
        internal_nodes = list(set(tables.edges.child))
        # add mutations on the leaves
        for node_id in range(internal_nodes[-1] + 1, tables.nodes.num_rows):
            tables.mutations.add_row(site = site_id, node = node_id, derived_state = '1')

    tables.sort()
    new_tree_sequence = tables.tree_sequence()
    
    return new_tree_sequence

# %%
def remove_variants(tree_sequence, old_positions, total_num_positions):
    new_positions = random.sample(old_positions, total_num_positions)

    tables = tree_sequence.dump_tables()
    sites = tables.sites.copy()
    mutations = tables.mutations.copy()
    tables.sites.clear()
    tables.mutations.clear()
    
    old_node = {mut.site:mut.node for mut in mutations}

    for old_site_id,site in enumerate(sites, 0):
        if site.position in new_positions:
            site_id = tables.sites.add_row(position = site.position, ancestral_state = '1')
            tables.mutations.add_row(site = site_id, node = old_node[old_site_id], derived_state = '1')
    tables.sort()
    new_tree_sequence = tables.tree_sequence()

    return new_tree_sequence

# %%
def write_chromosomes(tree_sequence, contig_id, vcf_path_out, ploidy = 2):
    n_dip_indv = int(tree_sequence.num_samples / 2)
    indv_names = ["tsk_{i}indv".format(str(i)) for i in range(n_dip_indv)]

    with open(vcf_path_out, "w") as vcf_file: 
        tskit.ALLELES_01 = ('1', '2')
        tree_sequence.write_vcf(
            output = vcf_file, 
            contig_id = contig_id, 
            ploidy = ploidy,
            individual_names=indv_names)

# %%
def change_ref_alt(tree_sequence):
    tables = tree_sequence.dump_tables()

    sites = tables.sites.copy()
    mutations = tables.mutations.copy()

    tables.sites.clear()
    tables.mutations.clear()

    for site in sites:
        site_id = tables.sites.add_row(position = site.position, ancestral_state = '1')

    for mut in mutations:
        tables.mutations.add_row(site = mut.site, node = mut.node, derived_state = '2')

    new_tree_sequence = tables.tree_sequence()

    return new_tree_sequence

# %%


def main():

    print('ARGV      :', sys.argv[1:])

    options, remainder = getopt.getopt(sys.argv[1:], 
                                    'r:p:N:s:i:n:o:', ['recombination-map=',
                                                    'positions=',
                                                    'Ne=',
                                                    'sample-size=',
                                                    'chromosome=',
                                                    'num-replicates=',
                                                    'vcf-out='
                                                    ])
    print('OPTIONS   :', options)

    # %%
    num_replicates = 1
    for opt, arg in options:
        if opt in ('-r', '--recombination-map'):
            recomb_map_path = arg
        elif opt in ('-p', '--positions'):
            positions_path = arg
        elif opt in ('-N', '--Ne'):
            Ne = int(arg)
        elif opt in ('-s', '--sample-size'):
            sample_size = np.int(arg)
        elif opt in ('-i', '--chromosome'):
            chr = int(arg)
        elif opt in ('-n', '--num-replicates'):
            num_replicates = int(arg)
        elif opt in ('-o', '--vcf-out'):
            vcf_path_out = arg


    my_recomb_map = load_recomb_map_chr(recomb_map_path )
    simulated_raw = simulate_data_chr(Ne = Ne, 
                                        sample_size = sample_size, 
                                        recombination_map = my_recomb_map,
                                        num_replicates=num_replicates)
    # print(simulated_raw.num_sites)
    S = np.zeros(num_replicates)
    for j, tree_sequence in enumerate(simulated_raw):
        S[j] = tree_sequence.num_sites
    print(np.mean(S))
    # print("              mean              variance")
    # print("Observed      {}\t\t{}".format(np.mean(S), np.var(S)))

    # positions = load_positions(positions_path)
    # simulated_full = adjust_number_sites(simulated_raw, positions)
    # simulated_full_recoded = change_ref_alt(tree_sequence = simulated_full)

    # write_chromosomes(tree_sequence = simulated_full_recoded,
    #                     contig_id = chr, 
                        # vcf_path_out = vcf_path_out)


# %%
if __name__ == "__main__":
    main()

