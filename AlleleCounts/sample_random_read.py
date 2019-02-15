#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 14:27:38 2019

@author: Diana I. Cruz Davalos
"""

# %%
# Import modules

import numpy as np

# %%
# Load counts of reference and alternative alleles

counts_path = "/Users/dcruz/Projects/Botocudos/Files/test/counts.gz"
output_name = "/Users/dcruz/Projects/Botocudos/Files/test/sampled.txt"
counts = np.loadtxt(counts_path)

# %%
# Dimensions
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

    # Remember to switch
#alt_is_major_allele = ma.array(np.ones(counts_ref.shape))
#alt_is_major_allele[np.where(freqs > 0.5)] = 0
#alt_is_major_allele[missing_data] = ma.masked

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
alleles[np.where(sampled == True)] = 1
alleles[np.where(sampled == False)] = 0

np.savetxt(fname = output_name, X = alleles, fmt = "%1.f")