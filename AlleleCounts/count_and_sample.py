#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 23:36:21 2019

@author: dcruz
"""
# %%
import re
import numpy as np
import numpy.ma as ma
import time
import sys
# %%
path_mpileup = sys.argv[1]
path_out_counts = sys.argv[2]
path_out_sampled = sys.argv[3]
path_sites = sys.argv[4]

#%%
#path_mpileup = "/Users/dcruz/Projects/Botocudos/Files/test/test_mpileup.txt"
#path_out_counts = "/Users/dcruz/Projects/Botocudos/Files/test/test_counts.gz"
#path_out_sampled = "/Users/dcruz/Projects/Botocudos/Files/test/test_sampled.gz"
#path_sites = "/Users/dcruz/Projects/Botocudos/Files/test/mpileup.sites.txt"

# %%
base_column = {0:"A", 1:"C", 2:"G", 3:"T"}
refalt = {}
j = 0
# %%
# Parse the amazing mpileup format
def parse_line(pileup, nInd):
    # remove missing, beginning, end
    #pattern = "\*|\^.|\$"
    j += 1
    print(j)
    pattern = "\^."
    parsed = re.sub(pattern, "", pileup)
    pattern = re.compile(r"\+\d+|\-\d+")
    matches = [int(s) for s in re.findall(pattern, parsed)]
    # remove indels
    for p in matches:
        pattern = re.compile("[+|-]" + str(abs(p)) + ".{" + str(abs(p)) +"}")
        parsed = re.sub(pattern, "", parsed)

    # Count matches
    allele_counts = []
    for ind in range(0, nInd):
        i = 3*ind + 4
        #print("Value of i:")

        # for base in ref, alt
        pos = parsed.split("\t")[0] + "_" + parsed.split("\t")[1]
        string = parsed.split("\t")[i]
        for base in refalt[pos]:
            pattern = base_column[base]
            count = len(tuple(re.finditer(pattern, string, flags = re.I)))
            allele_counts.append(count)

    return(allele_counts)


def add_key(line):
    refalt[line.split()[0]] = [int(line.split()[1]), int(line.split()[2])]
# %%
# Get reference and alternative alleles
with open(path_sites, 'r') as sites:
    [add_key(line) for line in sites.readlines()]

#print(refalt.keys())
#print("we have "+str(len(refalt.keys()))+" keys. For example " + refalt.keys()[0])

# %%

nInd = int((len(open(path_mpileup, 'r').readline().split("\t")) -3 ) / 3)
print(nInd)
#%%

print("Parsing and counting bases.")
start = time.time()
with open(path_mpileup, "r") as file:
    counts = np.array([parse_line(line, nInd) for line in file.readlines()])
end = time.time()
print(end - start)
np.savetxt(fname = path_out_counts, X = counts, fmt = "%1.f")

# %%
# Dimensions
print("Doing other things")
start = time.time()
nSites, nInd = counts.shape

index_ref = range(0, nInd, 2)
index_alt = range(1, nInd, 2)
nInd = int(nInd/2)

counts_ref = ma.array(counts[:, index_ref])
counts_alt = ma.array(counts[:, index_alt])
del counts

# %%
# Find missing data
missing_data = np.logical_and(counts_ref==0, counts_alt==0)

# %%
# Find sites where there is only one allele
only_one = np.logical_xor(counts_ref, counts_alt)
lonely_alleles = counts_ref[np.where(only_one)] > 0

# %%
# Mask array where data are missing or there is only one allele
to_mask = np.logical_or(only_one, missing_data)
counts_ref[to_mask] = ma.masked
counts_alt[to_mask] = ma.masked

# %%
# Calculate base frequencies

freqs = counts_ref/(counts_ref + counts_alt)

alt_is_major_allele = np.where(freqs < 0.5)
freqs[alt_is_major_allele] = 1 - freqs[alt_is_major_allele]

# %%
# Sample
probs = np.random.sample(size = missing_data.shape)

sampled = probs < freqs
sampled[alt_is_major_allele] = np.logical_not(sampled[alt_is_major_allele])
# %%
del probs
del freqs
del alt_is_major_allele
# %%
alleles = np.empty(sampled.shape)
alleles[:] = np.nan
alleles[only_one] = lonely_alleles
alleles[np.where(sampled == True)] = 0
alleles[np.where(sampled == False)] = 1
end = time.time()
print(end - start)
np.savetxt(fname = path_out_sampled, X = alleles, fmt = "%1.f")