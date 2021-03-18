#!/usr/bin/python

# Jan 29 2020.
# Script to select reads that mapped against a human virus in the DIAMOND BLAST like output.
# The selected reads are written to another file in a DIAMOND tax format way.
# The script receives a file with the human virus taxids and a file with DIAMOND output

# Example run:
# python count_reads.py HSapiensVirus_taxids MN00010_alg.txt MN00010_Human_Virus_hits.csv

# %%
import sys
from itertools import groupby
import csv

# %%
# sys.argv = [
#     __file__,
#     '/Users/yami_ommar/axiom_virome/Boto_virome/true_complete/DIAMOND_output/HSapiensVirus_taxids',
#     '/Users/yami_ommar/axiom_virome/Boto_virome/MN0003_OneHit_BLAST6.txt'
# ]

# %%
# Read human virus tax ids, save them in a set.
human_viruses = open(sys.argv[1], 'r')
lines = human_viruses.readlines()
lines = [lines[i].rstrip('\n') for i in range(len(lines))]
taxid_hs_virus = set(lines)
human_viruses.close()

# %%
# Read DIAMOND's output.
diamond_output = open(sys.argv[2], 'r')

# Save first line in a data frame.
line = diamond_output.readline()
read_hits = [line.rsplit(sep="\t")]
# The first read, in the first hit.
current_read = read_hits[0]

# To save best human virus hits.
read_best_hit = []


# %%
# Function to favour human virus hits in case of ties.
# It uses groupby from itertools.

def human_virus_tie(hits, viruses):
    # Select the e-values
    e_values = [x[13] for x in hits]
    # Count repetitions of values, store them.
    counts = [len(list(group)) for key, group in groupby(e_values)]
    # Select taxids
    tax_ids = [x[2] for x in hits]

    # Check if there's a tie.
    if counts[0] > 1:
        print("We have a tie")
        # Subset the tied hits.
        tied_subset = tax_ids[0:counts[0]]
        # Ask if they are human viruses.
        for i in range(1, counts[0]):
            if tied_subset[i] in viruses:
                print("There is a human virus with the same score, returning it")
                return i
                # [i for i in range(1, counts[0]) if tied_subset[i] in human_viruses]


# %%
while True:
    line = diamond_output.readline()
    if line:
        # Create a list.
        current_line = line.rsplit(sep="\t")
        # Float to compare e-values.
        current_line[13] = float(current_line[13])

        # Group hits involving the same read.
        if current_read == current_line[0]:
            read_hits.append(current_line)
        # No more hits with the same read? Then, select the best one of the group.
        else:
            # Sort the list of lists by e-value.
            read_hits.sort(key=lambda x: x[13])
            # Check if the best hit is against a human virus, save it if that's the case
            if read_hits[0][2] in taxid_hs_virus:
                read_hits[0][14] = read_hits[0][14].rstrip('\n')
                read_best_hit.append(read_hits[0])
            # Check for ties with human viruses
            else:
                index = human_virus_tie(read_hits, taxid_hs_virus)
                if index:
                    read_hits[index][14] = read_hits[0][14].rstrip('\n')
                    read_best_hit.append(read_hits[index])
            # Update current read, and group of hits
            current_read = current_line[0]
            read_hits = [current_line]
    # For last line read
    else:
        read_hits.sort(key=lambda x: x[13])
        if read_hits[0][2] in taxid_hs_virus:
            read_hits[0][14] = read_hits[0][14].rstrip('\n')
            read_best_hit.append(read_hits[0])
        # Check for ties with human viruses
        else:
            index = human_virus_tie(read_hits, taxid_hs_virus)
            if index:
                read_hits[0][14] = read_hits[0][14].rstrip('\n')
                read_best_hit.append(read_hits[index])
        break

diamond_output.close()

print(read_best_hit)

with open(sys.argv[3], "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(read_best_hit)


