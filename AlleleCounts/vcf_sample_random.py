#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 19:39:24 2019

@author: dcruz
"""

# %%
import numpy as np
import re
import time
import numpy.ma as ma
import sys
# %%
vcf_path = sys.argv[1]
#vcf_path = "/Users/dcruz/Projects/Botocudos/Files/test/88ind_head.vcf"
counts_path = sys.argv[2]
#counts_path = "/Users/dcruz/Projects/Botocudos/Files/test/88ind_counts.txt"
path_out_sampled = sys.argv[3]
# path_out_sampled = "/Users/dcruz/Projects/Botocudos/Files/test/88ind_sampled.txt"
# %%
# Parse VCF format to allele counts

def parse_genos(line):

    pattern = "/"
    line = "\t".join(line.split()[9:])
    parsed_line = re.sub(pattern, "\t", line)

    pattern = "\."
    parsed_line = re.sub(pattern, "nan", parsed_line)
    pattern = "$"
    parsed_line = re.sub(pattern, "\n", parsed_line)
    parsed_line.split("\t")

    return(parsed_line)

# %%
# Skip header and get nInd
counts = open(counts_path, "w")
with open(vcf_path, 'r') as vcf:
    while re.match("#", vcf.readline()):
                   continue
    [counts.write(parse_genos(line)) for line in vcf.readlines()]

counts.close()
# %%

counts = np.loadtxt(counts_path)

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
heterozygote = np.zeros(freqs.shape)+0.5

alt_is_major_allele = np.where(ma.less(freqs,heterozygote ))
freqs[alt_is_major_allele] = 1 - freqs[alt_is_major_allele]

# %%
# Sample
probs = np.random.sample(size = missing_data.shape)

sampled = ma.less(probs, freqs)
sampled[alt_is_major_allele] = np.logical_not(sampled[alt_is_major_allele])
# %%
del probs
del freqs
del alt_is_major_allele
# %%
alleles = np.empty(sampled.shape)
alleles[:] = np.nan
alleles[only_one] = lonely_alleles
alleles[np.where(ma.equal(sampled, True))] = 0
alleles[np.where(ma.equal(sampled, False))] = 1
end = time.time()
print(end - start)
np.savetxt(fname = path_out_sampled, X = alleles, fmt = "%1.f")
