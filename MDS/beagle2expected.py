#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 18:15:09 2018

@author: Diana I. Cruz Davalos
"""
import numpy as np
import sys
# Use genotype likelihoods to compute the distance matrix


#%%

beagle_name = sys.argv[1]
print("Input file: " + beagle_name)
#"/Users/dcruz/Projects/Botocudos/Files/test/head.beagle"
expected_name = sys.argv[2]
beagle = np.loadtxt(beagle_name, dtype = "str")
beagle = beagle[1:, 3:].astype(np.float16)
#%% 
def calculate_expected(index, beagle):
    expected = beagle[:,index]*np.int(0) + beagle[:,index + 1]*1 + beagle[:, index  + 2]*2
    return expected
#%%
nInd = beagle.shape[1] - 1
print("Number of individuals: " + str(nInd + 1))
expected_alleles = np.array([calculate_expected(i, beagle) for i in range(0,nInd, 3)])
#%%
print("Output file: " + expected_name)
np.savetxt(fname = expected_name, fmt = '%.6f',
           X = expected_alleles, delimiter='\t')  