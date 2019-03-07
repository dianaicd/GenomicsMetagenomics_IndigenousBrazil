#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 23:55:13 2018

@author: dcruz
Missing data are coded as 9
"""

import pysam
import numpy as np
import random
import re
import sys

# %%
print "This is the name of the script: ", sys.argv[0]
print "Number of arguments: ", len(sys.argv)
print "The arguments are: " , str(sys.argv)

# %%
# Function to verify that a base:
# has a base qualiity >= baseQ_thres
# has the same state as the reference or alternative allele
# am I missing other filters?
# In the future, we might want to account for damage

def meets_req(base, baseQ, ref, alt, baseQ_thres):
    if baseQ >= baseQ_thres and (base == ref or base == alt):
        return True
    else:
        return False


# %%
# Function to sample a random base aligned to a position in the genome
# From BAM to a single base
# position in format "chromosome_position"
def sample_allele_bam(position, samfile, ref, alt, baseQ_thres=30):
    [Chr, pos] = re.split("_", position)
    pos = int(pos)
    bases = []
    n_reads = samfile.count(contig=Chr, start=pos, end=pos + 1)

    # Iterate over the position
    for pileupcolumn in samfile.pileup(Chr, pos, pos + 1):
        if pileupcolumn.pos == pos:
            # Iterate over the reads
            for pileupread in pileupcolumn.pileups:
                # Control for base quality
                Q = pileupread.alignment.query_qualities[pileupread.query_position]
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if  meets_req(base = base, baseQ = Q, ref = ref,
                              alt = alt, baseQ_thres = baseQ_thres):
                    # Debug: Are we looking at bases at the right position?
                    #print( pileupread.alignment.positions[pileupread.query_position])
                    if base == ref:
                        bases.append(0)
                    else:
                        bases.append(1)


    if n_reads == 0 or len(bases) == 0:
        return 9
    else:
        selected = random.randint(0, len(bases)-1)
        return bases[selected]


#%%
# Function to sample a random allele from a VCF file (one position in the genome)
# Returns one allele (0, 1 or None) indicating (reference, alternative, missing data)

def sample_allele_vcf(Chr,
                      vcf, anc_bases=[], recode = False, bases_ind = [],
                      positions = {}
                      ):

    # This goes over the rows with data of a vcf (skipping header)
    index_pos = -1
    for record in vcf.fetch(contig = Chr):
        index_pos += 1
        index_ind = 0
        pos = str(Chr) + "_" + str(record.pos)
        # Skip sites where there is not info in the ancestral sample
        if pos not in positions:
            index_pos -= 1
            continue

        for ind in record.samples.itervalues():
            #print(ind.name)
            # 0/0 means it has two reference alleles
            # 1/1 means it has two alternative alleles

            selected = random.randint(0,1)
            base = ind.values()[0][selected]

            # If by chance we have 1/. or 0/., try to sample the non-missing allele
            if base == None:
                selected = abs(selected - 1)
                base = ind.values()[0][selected]
                if base == None:
                    base = 9

            # Recode according to ancestral state
            if recode:
                #print(base)
                #print(positions[pos])
                if base == positions[pos]:
                    base = 0
                else:
                    base = 1

            bases_ind[index_pos, index_ind] = base
            index_ind += 1

    return(bases_ind)

# %%
# Define blocks of block_size
# Breakpoints indicate start and end of the block
# Default is 5 megabases (5e6)

def find_breakpoints(block_size, contig, positions):

    pattern = str(contig) + "_"
    sites = np.array([int(re.split("_", pos)[1])
    for pos in positions  if re.match(pattern, pos) ])

    i = 0
    breakpoints = []

    while i < len(positions) - 1:
        index = np.where(sites - sites[i] >= block_size)

        if(np.shape(index)[1] > 1):
            breakpoints.append((i, index[0][0]))
            i = index[0][1]
        else:
            breakpoints.append((i, len(positions)))
            i = len(positions)

    return breakpoints

# %%

def abba_baba(samples, index_anc, breakpoints, full_data, positions, contig):
    nInd = len(samples)
    nBlocks = len(breakpoints)
    abbababa = np.zeros(shape = (nBlocks, 2*(nInd-1)*(nInd-2)*(nInd-3) + 3), dtype = int)
    l = 3
    for h3 in range(0, nInd):

        if h3 == index_anc:
            continue
        for h2 in range(0, nInd):
            if h2 == index_anc or h2 == h3:
                continue
            #l+=2           #print(h2)
            for h1 in range(0, nInd) :
                if h1 == h3 or h1 == h2 or h1 == index_anc:
                    continue
                #print(h1)

                for k in range(0, nBlocks):
                    i = breakpoints[k][0]
                    j = breakpoints[k][1]

                    abbababa[k][0] = contig
                    abbababa[k][1] = re.split("_", positions[breakpoints[k][0]])[1]

                    abbababa[k][2] = re.split("_", positions[breakpoints[k][1] - 1])[1]
                    counts = full_data[i:j, h1] + 2*full_data[i:j, h2] + 3*full_data[i:j, h3]

                    # Number of ABBA
                    abbababa[k][l] = sum(counts == 5)
                    # Number of BABA
                    abbababa[k][l+1] = sum(counts == 4)
                #print(l)
                l += 2

    return abbababa


# %%
# path to vcf.gz file
vcf_path = sys.argv[1]
# file with paths to bam files
bam_list = open(sys.argv[2], "r")
# basename
out = sys.argv[3] + ".abbababa"
logfile = sys.argv[3] + ".log"
logfile = open(logfile, "w")
# Ancestral
anc = sys.argv[4]
# Contig
Chr = sys.argv[5]

#%%
# VCF file
# do not forget to have the corresponding index:
# bgzip test.vcf
# bcftools index test.vcf.gz
Chr = "1"
vcf_anc = pysam.VariantFile(vcf_path)
vcf_anc.subset_samples([anc])

vcf = pysam.VariantFile(vcf_path)

bam_paths = []

for line in bam_list:
    bam_paths.append(line.split())

logfile.write("starting")
#%%
# Obtain positions and ancestral allele
# Save samples' names

chroms = []
positions = {}
ref_alt = {}
samples_names = []

i = 0
for name in vcf.header.samples:
    samples_names.append(name)
    if(name == anc):
        index_anc = i
    i += 1

for contig in vcf_anc.header.contigs:
    chroms.append(contig)
    for rec in vcf_anc.fetch(contig = contig):
        positions[str(contig) + "_" + str(rec.pos)] = 1
        ref_alt[str(contig) + "_" + str(rec.pos)] = rec.alleles

anc_bases = np.empty(shape = (len(positions), 1), dtype = int)
anc_bases = sample_allele_vcf(Chr = chroms[0], vcf = vcf_anc,
                              bases_ind = anc_bases,
                              positions = positions)

logfile.write("ancestral bases obtained")

#%%
# Build dictionary of positions and bases in ancestral state
# remove sites with missing data
i = 0
to_remove = []
for key in positions:
    if(anc_bases[i][0] != 9):
        positions[key] = anc_bases[i][0]
    else:
        to_remove.append(key)
    i += 1

for key in to_remove:
    positions.pop(key)

logfile.write("positions dictionary built")
#%%
# Recode alleles as 0 for ancestral, 1 for derived and 3 for missing
nInd = len(vcf.header.samples)
bases_ind = np.empty(shape = (len(positions), nInd), dtype = int)
bases_ind = sample_allele_vcf(Chr = chroms[0], vcf = vcf,
                            positions = positions, anc_bases = anc_bases,
                             recode = True, bases_ind = bases_ind)
vcf.close()

logfile.write("Alleles recoded")
#%%
# Sample an allele at random for each of the positions in the panel
# for the bam file
# Merge sampled alleles to panel
# samfile = target BAM

baseQ = 30
#bam_paths = ["/Users/dcruz/Projects/Botocudos/Files/test/test_1.bam"]

for path in bam_paths:
    samfile = pysam.AlignmentFile(path[0], "rb")

    sampled = np.array([sample_allele_bam(position=pos, samfile=samfile,
                             ref=ref_alt[pos][0],
                             alt=ref_alt[pos][1])
        for pos in positions],
            ndmin = 2, dtype = int).T

    samfile.close()
    logfile.write("Alleles sampled from bam file")

    # Merge sampled alleles to panel

    full_data = np.hstack([bases_ind, sampled])
    samples_names.append(re.split("\.", re.split("/", path[0])[-1])[0])

    #del sampled

# %%
# Define blocks of block_size
# Breakpoints indicate start and end of the block
# Default is 5 megabases (5e6)

breakpoints = find_breakpoints(block_size = 5e6, contig = Chr,
                               positions = positions)

logfile.write("breakpoints obtained")

# %%
#tmp = full_data
#tmp_pos = positions
#tmp_samples = samples_names
#samples = list(range(0, 5))
#positions = list(positions.keys())
#index_anc = 0
#full_data = full_data[:, 0:5]

#%%
result = abba_baba(samples = samples_names, index_anc = index_anc,
                   breakpoints = breakpoints, full_data = full_data,
                   positions = list(positions.keys()), contig = Chr)

print("result's dimensions")
print(result.shape)
np.savetxt(fname=out, 
           X=result, delimiter='\t', fmt = "%.0f")
logfile.close()

