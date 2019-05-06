#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 00:21:42 2019

@author: dcruz

Script to calculate pi as in Skoglund et al. 2014 (Science)
"""
# %%
import numpy as np
import numpy.ma as ma
import sys, getopt
from astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats
# %%

path_sampled = "/Users/dcruz/Projects/Botocudos/Files/test/Botocudos.counts.sampled.txt"
path_sites = "/Users/dcruz/Projects/Botocudos/Files/test/Nigeria_B_Yoruba-3.refalt"
basename = "/Users/dcruz/Projects/Botocudos/Files/test/out"
blockSize = 5e6
autosomes = False

print('ARGV      :', sys.argv[1:])

options, remainder = getopt.getopt(sys.argv[1:], 'c:s:o:ab:', ['counts=',
                                                         'sites=',
                                                         'output=',
                                                         'autosomes_only',
                                                         'block_size='])
print('OPTIONS   :', options)

for opt, arg in options:
    if opt in ('-c', '--counts'):
        path_sampled = arg
    elif opt in ('-s', '--sites'):
        path_sites = arg
    elif opt in ('-o', '--output'):
        basename = arg
    elif opt in ('-a', '--autosomes_only'):
        autosomes = True
    elif opt in ('-b', '--block_size'):
        blockSize = arg


#%%
sampled = ma.array(np.loadtxt(path_sampled))
sampled[np.isnan(sampled)] = ma.masked

# %%
#

def add_sites(line):
    line = line.replace("_", "\t")

    if line.split()[0] == "X":
        chr=23
    elif line.split()[0] == "Y":
        chr=24
    else:
        chr = int(line.split()[0])

    pos = int(line.split()[1])

    return((chr,pos))

with open(path_sites, 'r') as file:
    sites = np.array(
            [add_sites(line) for line in file.readlines()]
            )


# %%
def pontus_pi(mismatches):
    PI = ma.sum(mismatches)/ma.count(mismatches)

    return(PI)

# %%
# Define blocks that will be jackknifed (usually 5 Mb)

def define_blocks(bllockSize, chr, sites):

    positions = sites[np.where(sites[:,0] == chr),1]
    blockLimits = []
    i = 0

    while positions[0][i] < positions[0][-1]:
        blockLimits.append(i)#positions[0][i])

        i = np.where(positions[0] - positions[0][blockLimits[-1]] > blockSize)[0][0]

        if positions[0][-1] - positions[0][i] < blockSize:
            blockLimits.append(i)#positions[0][i])
            i = -1
    return(np.array(blockLimits))

# %%
# Get sites with data for at least 2 individuals
# and sample two alleles
def sample_2ind(Block, sites):
    counts = np.array([ma.count(Block[x,]) for x in range(0,Block.shape[0])])
    index2alleles = ma.where(np.greater_equal(counts, 2))[0]

    counts = counts[index2alleles]
    Block = Block[index2alleles, :].reshape(counts.shape[0], Block.shape[1])

    probs = np.random.sample(size = counts.shape)
    index2sampled = np.array(probs*counts, dtype = np.int)
    firstChr = []
    for i in range(0, Block.shape[0]):
        firstChr.append(Block[i,~Block[i,:].mask][index2sampled[i]])
        Block[i,~Block[i,:].mask].mask[index2sampled[i]] = True

    counts = counts - 1
    probs = np.random.sample(size = counts.shape)
    index2sampled = np.array(probs*counts, dtype = np.int)
    secondChr = [Block[i,~Block[i,:].mask][index2sampled[i]]
                    for i in range(0, Block.shape[0])]

    twoChrs = np.array((firstChr, secondChr), dtype=np.int)

    return(twoChrs, sites[index2alleles])
#%%

twoChrs, newSites = sample_2ind(sampled, sites)

# %%
blockGenome = np.hstack([define_blocks(blockSize, chr, newSites)
                        for chr in np.unique(newSites[:, 0])])
nChr = len(np.unique(newSites[:, 0]))
pi_genome = np.zeros(shape=(blockGenome.shape[0] - nChr,4))

j = 0

for chr in np.unique(newSites[:, 0]):
    blockLimits = define_blocks(blockSize, chr, newSites)
    print(chr)
    currentChr = twoChrs[:, np.where(newSites[:,0] == chr)[0]]
    #print(pi_genome)
    for i in range(0, blockLimits.shape[0]-1):
        pi_genome[j,0] = chr
        pi_genome[j,1] = blockLimits[i]
        pi_genome[j,2] = blockLimits[i+1]
        pi_genome[j,3] = pontus_pi(np.abs(currentChr[0,blockLimits[i]:blockLimits[i+1]] -
                            currentChr[1,blockLimits[i]:blockLimits[i+1]])
            )
        j += 1
# %%
np.savetxt(fname = basename+".pi.twochrs.txt", X=twoChrs.T, delimiter="\t",
    fmt = ['%i', '%i'])

np.savetxt(fname = basename+".pi.txt", X=pi_genome, delimiter="\t",
           fmt = ['%i','%i','%i','%f'], header = 'chr\tstart\tend\tpi')
# %%

resamples = jackknife_resampling(pi_genome[:,3])
np.savetxt(fname = basename+".pi.resampled.txt", X=resamples, delimiter="\t")
# %%

test_statistic = np.mean
jackKnifed = np.hstack((jackknife_stats(pi_genome[:,3], test_statistic, 0.95))).reshape((1,5))

np.savetxt(fname = basename+'.pi.stats.txt', X = jackKnifed,
           delimiter = "\t", fmt = '%f',
           header = 'estimate\tbias\tstderror\tconfint1\tconfint2')
# %%
#print("Estimate:", estimate)
#
#print("BIAS:", bias)
#
#print("STDERROR:",stderr)
#
#print("CONFINTERVAL:", conf_interval)












