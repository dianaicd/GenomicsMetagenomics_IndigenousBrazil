#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 13:58:59 2019

@author: Diana I. Cruz Davalos
"""

import numpy as np
import numpy.ma as ma
import sys
# Use genotype likelihoods to compute the distance matrix


# %%
counts_name = sys.argv[1]
dist_name = sys.argv[2]
# %%
print("Input file: " + counts_name)
# "/Users/dcruz/Projects/Botocudos/Files/test/head.beagle"

counts = np.loadtxt(counts_name)

# %%
def distance_ind(ind1, ind2):
    # mask missing values
    ind1 = ma.masked_invalid(ind1)
    ind2 = ma.masked_invalid(ind2)
    # calculate dstance ignoring missing values
    #dist = np.sum(is_not_na * np.power((ind1 - ind2),2))/np.sum(is_not_na)
    diff = ma.abs(ind1 - ind2)
    dist = ma.sum(diff)/diff.count()
    return(dist)

# %%
(nSNP, nInd) = counts.shape

print("Number of individuals: " + str(nInd))
print("Number of markers: " + str(nSNP))

# %%

final_dist = np.zeros(shape = (nInd, nInd))
for i in range(0, nInd):
    for j in range(0, nInd):
        if i < j:
            final_dist[i][j] = distance_ind(counts[:,i], counts[:, j])

# %%
print("Distance file: " + dist_name)
np.savetxt(fname=dist_name, fmt='%.6f',
           X=final_dist, delimiter='\t')