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
import sys, getopt
import gzip
#%%
path_mpileup = False
path_out_counts = "counts.txt"
path_out_sampled = "sampled.txt"
path_sites = "sites.refalt"
ped = False
all_mutations = False

path_mpileup = "/Users/dcruz/Projects/Botocudos/Files/test/bamfile.list_19.mpileup"
path_out_counts = "/Users/dcruz/Projects/Botocudos/Files/test/test_counts.gz"
path_out_sampled = "/Users/dcruz/Projects/Botocudos/Files/test/test_sampled.gz"
path_sites = "/Users/dcruz/Projects/Botocudos/Files/test/MAF_NoLaz_HOP_MinMyc_19.refalt"

print('ARGV      :', sys.argv[1:])

options, remainder = getopt.getopt(sys.argv[1:], 'm:c:s:rpa', ['mpileup=',
                                                         'counts=',
                                                         'sampled=',
                                                         'refalt=',
                                                         'ped',
                                                         'allmutations'])
print('OPTIONS   :', options)

for opt, arg in options:
    if opt in ('-m', '--mpileup'):
        path_mpileup = arg
    elif opt in ('-c', '--counts'):
        path_out_counts = arg
    elif opt in ('-s', '--sampled'):
        path_out_sampled = arg
    elif opt in ('-r', '--refalt'):
        path_sites = arg
    elif opt in ('-p', '--ped'):
        ped = True
    elif opt in ('-a', '--allmutations'):
        all_mutations = True

# %%
base_column = {0:"A", 1:"C", 2:"G", 3:"T"}
refalt = {}
# %%
# Parse the amazing mpileup format
def parse_line(pileup, nInd):
#%%
    # remove missing, beginning, end
    #pattern = "\*|\^.|\$"
    # Get olnly the columns with the bases
    pos = pileup.split("\t")[0] + "_" + pileup.split("\t")[1]
    pileup = "\t".join(pileup.split("\t")[4:(nInd*3+8):3])
    #%%
    pattern = "\^."
    parsed = re.sub(pattern, "", pileup)
    pattern = re.compile(r"\+\d+|\-\d+")
    matches = [int(s) for s in re.findall(pattern, parsed)]
#%%
    # remove indels
    for p in matches:
        pattern = re.compile("[+|-]" + str(abs(p)) + ".{" + str(abs(p)) +"}")
        parsed = re.sub(pattern, "", parsed)
#%%
    # Count matches
    allele_counts = []
        
    for ind in range(0, nInd):
        # for base in ref, alt
        string = parsed.split("\t")[ind]
        for base in refalt[pos]:
            pattern = base_column[base]
            count = len(tuple(re.finditer(pattern, string, flags = re.I)))
            allele_counts.append(count)
        allele_counts.append(len(string))
#%%
    return(allele_counts)

# %%
# Build dictionary with position as key
def add_key(line):
    refalt[line.split()[0]] = [int(line.split()[1]), int(line.split()[2])]
    return(line.split()[0])

# %%
# Parse 0 and 1 to nucleotides
def int2nucleotide(line, nucleotides):
    # in ped, do it diploid
    
    ref = base_column[nucleotides[0]] + " " + base_column[nucleotides[0]]
    alt = base_column[nucleotides[1]] + " " + base_column[nucleotides[1]]
    lineParsed = re.sub("0.", ref, line)
    lineParsed = re.sub("1.", alt, lineParsed)
    lineParsed = re.sub("nan", "0 0", lineParsed)
    lineParsed = re.sub("\[", "", lineParsed)
    lineParsed = re.sub("\]", "", lineParsed)
    lineParsed = re.sub("\n", "", lineParsed)
    # Remove extra spaces
    lineParsed = re.sub("\s{2,}", " ", lineParsed)
    lineParsed = re.sub("^\s", "", lineParsed)
    lineParsed = re.sub("\s$", "", lineParsed)
    
    return(lineParsed+"\n")

# %%
#print("Parsing and counting bases.")
if path_mpileup:
# Get reference and alternative alleles
    with open(path_sites, 'r') as sites:
        positions = [add_key(line) for line in sites.readlines()]
    nInd = int((len(open(path_mpileup, 'r').readline().split("\t")) -3 ) / 3)

    #counts_all_mut = np.zeros((len(positions), nInd))
    
    start = time.time()
    with open(path_mpileup, "r") as file:
        counts = np.array([parse_line(line, nInd) 
                                        for line in file.readlines()])
    end = time.time()
    #print(end - start)
    np.savetxt(fname = path_out_counts, X = counts, fmt = "%1.f")

# %%
# Dimensions
#print("Subsetting reference and alternative alleles.")
counts = np.loadtxt(path_out_counts)
start = time.time()
nSites, nInd = counts.shape

# First column is reference, second is alternative
# if asking to draw an allele using the frequencies of all observed states,
# then there is a third column per individual with total counts
# Regardless of what was asked, index every 3 columns
#if all_mutations:
#    nColumns = 3
#else:
#    nColumns = 2
nColumns = 3
    
index_ref = range(0, nInd, nColumns)
index_alt = range(1, nInd, nColumns)

if all_mutations:
    index_total_counts = range(2, nInd, nColumns)
    all_counts = counts[:, index_total_counts]
    print(all_counts.shape)

nInd = int(nInd/nColumns)

counts_ref = ma.array(counts[:, index_ref])
counts_alt = ma.array(counts[:, index_alt])
del counts

# %%
# Find missing data
#print("Finding missing data")
missing_data = np.logical_and(np.equal(counts_ref, 0),
                              np.equal(counts_alt, 0))

# Mask array where data are missing 
#print("Mask sites where data are missing or have only one allele.")
#to_mask = np.logical_or(only_one, missing_data)
counts_ref[missing_data] = ma.masked
counts_alt[missing_data] = ma.masked

# %%
# Calculate base frequencies
#print("Calculate base frequencies")
if not all_mutations:
    all_counts = counts_ref + counts_alt
# Frequencies of the reference allele
freq = counts_ref/all_counts#(counts_ref + counts_alt)
alt_is_major_allele = ma.where(freq < 0.5)
freq[alt_is_major_allele] = 1 - freq[alt_is_major_allele]

#print("frequencies for ind1:")
x = freq[freq[:,0].mask == False, 0]

#freq_ref[alt_is_major_allele] = 1 - freq_ref[alt_is_major_allele]

# %%
# Sample
#print("Sampling random alleles")
probs = np.random.sample(size = missing_data.shape)

sampled = ma.less_equal(probs, freq)
#sampled[alt_is_major_allele] = ma.greater(probs[alt_is_major_allele], 
#       freq_ref[alt_is_major_allele])
# Where alternative allele is at higher frequency,
# switch the selected allele
sampled[alt_is_major_allele] = ma.logical_not(sampled[alt_is_major_allele])
# %%
del probs
del freq
del alt_is_major_allele
# %%
# True if observing only alternative alleles


alleles = np.empty(sampled.shape)
alleles[:] = np.nan
# Must assign True to cells with reference alleles
#alleles[only_one] = np.greater_equal(counts_ref[only_one], 0)

alleles[np.where(sampled == True)] = 0
alleles[np.where(sampled == False)] = 1

if ped:
    bases = [refalt[pos] for pos in positions]
    #alleles = 
    
    with gzip.open(path_out_sampled, "wt") as file:
        file.writelines([int2nucleotide(np.array2string(alleles[i,:]), bases[i]) 
                for i in range(0, alleles.shape[0])])
    file.close()
else:
    np.savetxt(fname = path_out_sampled, X = alleles, fmt = "%1.f")
    
end = time.time()
print(end - start)
