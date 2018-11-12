#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 18:15:09 2018

@author: Diana I. Cruz Davalos
"""
import numpy as np
import sys
# Use genotype likelihoods to compute the distance matrix


# %%
beagle_name = sys.argv[1]
expected_name = sys.argv[2]

# %%
print("Input file: " + beagle_name)
# "/Users/dcruz/Projects/Botocudos/Files/test/head.beagle"

beagle = np.loadtxt(beagle_name, dtype="str")
beagle = beagle[1:, 3:].astype(np.float16)


# %%
def calculate_expected(index, beagle):
    expected = beagle[:, index]*0
    + beagle[:, index + 1]*1
    + beagle[:, index + 2]*2

    return expected


# %%
def missing_to_na(snp, index, beagle):
    if(beagle[snp, index] == beagle[snp, index + 1] and
       beagle[snp, index] == beagle[snp, index + 2]):
        return np.array([np.nan, np.nan, np.nan])
    else:
        return beagle[snp, index:index + 3]


# %%
nInd = beagle.shape[1]
nSNP = beagle.shape[0]
print("Number of individuals: " + str(nInd/3))
print("Number of markers: " + str(nSNP))
beagle = np.vstack([np.hstack([missing_to_na(snp, index, beagle)
                    for index in range(0, nInd - 1, 3)])
                   for snp in range(0, nSNP)])

expected_alleles = np.array([calculate_expected(i, beagle)
                            for i in range(0, nInd-4, 3)]).T
# %%

print("Output file: " + expected_name)
np.savetxt(fname=expected_name, fmt='%.6f',
           X=expected_alleles, delimiter='\t')
