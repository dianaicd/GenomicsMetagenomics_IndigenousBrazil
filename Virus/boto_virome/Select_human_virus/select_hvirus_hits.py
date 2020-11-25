#!/usr/bin/python

# February 11 2019
# Script to select the hits against human viruses
# This come from the count_reads.py script

# Example run:
# python count_reads.py HSapiensVirus_taxids MN00010_alg.txt MN00010_Human_Virus_hits.csv
# %%
import sys
from itertools import groupby
import csv


# %%
# Function to select hits against human viruses. It receives the human virus taxids (in a set) and the alignment file
# (file name).
def select_human_virus(viruses, alg_file):
    # Open DIAMOND's output alignment file.
    diamond_file = open(alg_file, 'r')
    # List to save the human hits.
    human_virus_hits = []
    # Read alignment file
    while True:
        line = diamond_file.readline()
        if line:
            line = line.rsplit(sep="\t")
            # Check if the virus is a human virus
            if line[2] in viruses:
                line[14] = line[14].rstrip('\n')
                human_virus_hits.append([line])
        else:
            break
    diamond_file.close()
    return human_virus_hits


# %%
# Read human virus tax ids, save them in a set.
human_viruses = open(sys.argv[1], 'r')
lines = human_viruses.readlines()
lines = [lines[i].rstrip('\n') for i in range(len(lines))]
taxid_hs_virus = set(lines)
human_viruses.close()

# %%

human_v_hits_test = select_human_virus(taxid_hs_virus, sys.argv[2])

# %%
with open(sys.argv[3], "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(human_v_hits_test)
