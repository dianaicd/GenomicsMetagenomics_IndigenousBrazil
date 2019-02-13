#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 18:15:09 2018

@author: Diana I. Cruz Davalos
"""

import numpy as np
import scipy.spatial as sc
import pandas as pd
import sys
# Use genotype likelihoods to compute the distance matrix


# %%<
beagle_name = sys.argv[1]
expected_name = sys.argv[2]
dist_name = sys.argv[3]
# %%
print("Input file: " + beagle_name)
# "/Users/dcruz/Projects/Botocudos/Files/test/head.beagle"

beagle = np.loadtxt(beagle_name, dtype="str")
# Forget about the 3 first columns (position, ref, alt)
beagle = beagle[1:, 3:].astype(np.float16)


# %%
def calculate_expected(index, beagle):
    # Calculate the expected number of alleles
    expected = beagle[:, index]*1 + \
    beagle[:, index + 1]*2 + beagle[:, index + 2]*3

    return expected


# %%
def missing_to_na(snp, index, beagle):
    # values of 0.33333,0.333333, 0.3333333 are considered as
    # missing data
    # missing data are encoded as 0
    if(beagle[snp, index] == beagle[snp, index + 1] and
       beagle[snp, index] == beagle[snp, index + 2]):
        return np.array([0,0,0])
    else:
        return beagle[snp, index:index + 3]

# %%
def distance_ind(ind1, ind2):
    # no missing value
    is_not_na = np.logical_and(ind1,ind2)
    # calculate distance ignoring missing values
    #dist = np.sum(is_not_na * np.power((ind1 - ind2),2))/np.sum(is_not_na)
    dist = np.sum(is_not_na * np.abs((ind1 - ind2)))/np.sum(is_not_na)
    return(dist)

# %%
(nSNP, nInd) = beagle.shape
print("Number of individuals: " + str(nInd/3))
print("Number of markers: " + str(nSNP))
beagle = np.vstack([np.hstack([missing_to_na(snp, index, beagle)
                    for index in range(0, nInd, 3)])
                   for snp in range(0, nSNP)])

expected_alleles = np.array([calculate_expected(i, beagle)
                            for i in range(0, nInd, 3)]).T
# %%
nInd = int(nInd/3)
final_dist = np.zeros(shape = (nInd, nInd))
for i in range(0, nInd):
    for j in range(0, nInd):
        if i < j:
            final_dist[i][j] = distance_ind(expected_alleles[:,i], expected_alleles[:, j])

#%%
# %%

print("Expected alleles file: " + expected_name)
np.savetxt(fname=expected_name, fmt='%.6f',
           X=expected_alleles, delimiter='\t')

print("Distance file: " + dist_name)
np.savetxt(fname=dist_name, fmt='%.6f',
           X=final_dist, delimiter='\t')