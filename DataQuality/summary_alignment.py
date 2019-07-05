#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 11:10:26 2019

@author: dcruz
"""

import pysam,sys,getopt,math

# %%
mito = 'MT'
Xchr = 'X'
Ychr = 'Y'
lib_type = "Single-end"
options, remainder = getopt.getopt(sys.argv[1:],
                                   's:o:m:t:l:x:y:',
                                   ['sample=', 'output=', 'mitochondrial=',
                                    'target=', 'lib_type=', 'x_chrom=',
                                    'y_chrom='])
for opt, arg in options:
    if opt in ('-o', '--output'):
        output = arg
    elif opt in ('-s', '--sample'):
        sample = arg
        target = sample
    elif opt in ('-m', '--mitochondrial'):
        mito = arg
    elif opt in ('-t', '--target'):
        target = arg
    elif opt in ('-l', '--lib_type'):
        library_type = arg
    elif opt in ('-x', '--x_chrom'):
        Xchr = arg
    elif opt in ('y', '--y_chrom'):
        Ychr = arg

# %%
def print_line(file):
    myLine = "\t".join([sample, target, library,
                        sex, str(Ry), str(Ry_confint),
                        lib_type, str(seq_reads_se), str(seq_trash_se),
                        str(seq_trash_se_frac), str(seq_retained_reads),
                        str(seq_retained_nts), str(seq_retained_length),
                        str(hits_raw_endogenous),
                        str(hits_raw_frac_endogenous), str(hits_clonality_endogenous),
                        str(hits_unique_endogenous), str(hits_unique_frac_endogenous),
                        str(hits_coverage_endogenous), str(hits_length_endogenous),
                        #hits_raw_mitochondrial, hits_raw_frac_mitochondrial,
                        #hits_clonality_mitochondrial,
                        str(hits_unique_mitochondrial), str(hits_unique_frac_mitochondrial),
                        str(hits_coverage_mitochondrial), str(hits_length_mitochondrial)#,
                        #hits_raw_nuclear,
                        #hits_raw_frac_nuclear, hits_clonality_nuclear,
                        #hits_unique_nuclear, hits_unique_frac_nuclear,
                        #hits_coverage_nuclear, hits_length_nuclear
                        ]) + "\n"
    file.writelines(myLine)
    #print(myLine)
# %%

path = sample + "/" + sample + ".bam"

myBam = pysam.AlignmentFile(path, "rb")
# %%
# Dict with libraries and their IDs
readGroup = myBam.header.as_dict()['RG']
myLibs = {}

for rg in readGroup:
    if rg['LB'] in myLibs:

        myLibs[rg['LB']].append(rg['ID'])
    else:
        myLibs[rg['LB']] = [rg['ID']]

# %% Read in settings from AdapterRemoval

reads = {}
trashed = {}
retained = {}
nucleotides = {}

for lib in myLibs.keys():
    path = sample + "/" + lib + "/" + lib + ".settings"
    with open(path, 'r') as file:
        file.readline()
        for line in file.readlines():
            if lib not in reads:
                reads[lib] = int(line.split()[-1])
            elif lib not in trashed:
                trashed[lib] = int(line.split()[-1])
            elif lib not in retained:
                retained[lib] = int(line.split()[-1])
            elif lib not in nucleotides:
                nucleotides[lib] = int(line.split()[-1])

# %%
# Read in rmdup info
path = sample + "/" + sample + "_rmdup.stats"
rmdup = {}
with open(path) as file:
    for line in file.readlines():
        rmdup[line.split()[0].split("/")[-1]] = int(line.split()[1])

# %%
# Read in genomecov info

mitoCov = {}
coverage = {}
basesMito = {}
bases = {}
for lib in myLibs.keys():
    path = sample + "/" + lib + "/" + lib + ".genomecov"
    with open(path, 'r') as file:
        for line in file.readlines():
            sequence, depth, times, length = line.split()[0:4]
            depth = int(depth)
            times = int(times)
            length = int(length)
            if sequence == mito:
                if lib not in basesMito:
                    basesMito[lib] = depth*times
                    mitoLength = length
                else:
                    basesMito[lib] += depth*times
            elif sequence == "genome":
                if lib not in bases:
                    totalLength = length
                    bases[lib] = depth*times
                else:
                    bases[lib] += depth*times
    coverage[lib] = bases[lib]/totalLength
    mitoCov[lib] = basesMito[lib]/mitoLength

# %%
# Get number of hits and read length

hits = {}
lengths = {}

mitoHits = {}
mitoLengths = {}

for lib in myLibs.keys():
    path = sample + "/" + lib + "/" + lib + ".length"
    with open(path, 'r') as file:
        sumLen = 0
        for line in file.readlines():
            l,times = line.split()[0:2]
            l = int(l)
            times = int(times)
            if lib not in hits:
                hits[lib] = times
            else:
                hits[lib] += times
            sumLen += times*l
    path = sample + "/" + lib + "/" + lib + ".mito.length"
    with open(path, 'r') as file:
        sumMitoLen = 0
        for line in file.readlines():
            l,times = line.split()[0:2]
            l = int(l)
            times = int(times)
            if lib not in mitoHits:
                mitoHits[lib] = times
            else:
                mitoHits[lib] += times
            sumMitoLen += times*l
    lengths[lib] = sumLen/hits[lib]
    mitoLengths[lib] = sumMitoLen/mitoHits[lib]

# %%
# Determine sex

x_counts = {}
y_counts = {}
ry = {}
ry_confint = {}
mySex = {}
maleLim = 0.075
femLim = 0.016

def determine_sex(ry, confint):
    sex='NA'
    if (ry < femLim) and (ry > maleLim):
        sex = 'Not Assigned'
    elif ry == 0:
        sex = 'consistent with XX'
    elif ry + confint < femLim:
        sex = 'XX'
    elif ry - confint > maleLim:
        sex = 'XY'
    elif ry - confint > femLim and ry + confint > maleLim:
        sex = 'consistent with XY but not XX'
    elif ry - confint < femLim and ry + confint < maleLim:
        sex = 'consistent with XX but not XY'
    else:
        sex = 'Not Assigned'
    return(sex)


for lib in myLibs.keys():
    path = sample + "/" + lib + "/" + lib + "_idxstats.txt"
    with open(path, "r") as file:
        for line in file.readlines():
            chr = line.split()[0]
            if chr == Xchr:
                x_counts[lib] = int(line.split()[2])
            elif chr == Ychr:
                y_counts[lib] = int(line.split()[2])
    ry[lib] = y_counts[lib]/(x_counts[lib] + y_counts[lib])
    ry_confint[lib] = 1.96*math.sqrt((ry[lib]*(1-ry[lib]))/(x_counts[lib] + y_counts[lib]))
    mySex[lib] = determine_sex(ry[lib], ry_confint[lib])

x_counts[sample] = sum(x_counts.values())
y_counts[sample] = sum(y_counts.values())
ry[sample] = y_counts[sample]/(x_counts[sample] + y_counts[sample])
ry_confint[sample] = 1.96*math.sqrt((ry[sample]*(1-ry[sample]))/(x_counts[sample] + y_counts[sample]))
mySex[sample] = determine_sex(ry[sample], ry_confint[sample])

# %%

#sample, target, library, lib_type
library = "*"

# sex, Ry, Ry_confint
sex = mySex[sample]
Ry = ry[sample]
Ry_confint = ry_confint[sample]

# seq_reads_se, seq_trash_se, seq_trash_se_frac,
seq_reads_se = sum(reads.values())
seq_trash_se = sum(trashed.values())
seq_trash_se_frac = seq_trash_se / seq_reads_se

#seq_retained_reads, seq_retained_nts, seq_retained_length,
seq_retained_reads = sum(retained.values())
seq_retained_nts = sum(nucleotides.values())
seq_retained_length = seq_retained_nts / seq_retained_reads

# hits_raw_endogenous, hits_raw_frac_endogenous, hits_clonality_endogenous,
hits_raw_endogenous = sum(rmdup.values()) + sum(hits.values())
hits_raw_frac_endogenous = hits_raw_endogenous / seq_retained_reads
hits_clonality_endogenous = sum(rmdup.values())

# hits_unique_endogenous, hits_unique_frac_endogenous,
# hits_coverage_endogenous, hits_length_endogenous,
hits_unique_endogenous = sum(hits.values())
hits_unique_frac_endogenous = hits_unique_endogenous/seq_retained_reads
hits_coverage_endogenous = sum(bases.values())/totalLength
hits_length_endogenous = sum(bases.values())/sum(hits.values())

# hits_unique_mitochondrial, hits_unique_frac_mitochondrial,
# hits_coverage_mitochondrial, hits_length_mitochondrial
hits_unique_mitochondrial = sum(mitoHits.values())
hits_unique_frac_mitochondrial = hits_unique_mitochondrial/seq_retained_reads
hits_coverage_mitochondrial = sum(basesMito.values())/mitoLength
hits_length_mitochondrial = sum(basesMito.values())/sum(mitoHits.values())

# %%
header = "\t".join(["sample", "target", "library",
                    "sex", "Ry", "Ry_confint","lib_type",
                    "seq_reads_se", "seq_trash_se", "seq_trash_se_frac",
                    "seq_retained_reads", "seq_retained_nts",
                    "seq_retained_length", "hits_raw_endogenous",
                    "hits_raw_frac_endogenous", "hits_clonality_endogenous",
                    "hits_unique_endogenous", "hits_unique_frac_endogenous",
                    "hits_coverage_endogenous", "hits_length_endogenous",
                    "hits_unique_mitochondrial", "hits_unique_frac_mitochondrial",
                    "hits_coverage_mitochondrial", "hits_length_mitochondrial"]) + "\n"

with open(output, 'w') as file:
    file.writelines(header)
    print_line(file)

    for lib in myLibs.keys():
        library = lib
        # sex, Ry, Ry_confint
        sex = mySex[lib]
        Ry = ry[lib]
        Ry_confint = ry_confint[lib]

        # seq_reads_se, seq_trash_se, seq_trash_se_frac,
        seq_reads_se = reads[lib]
        seq_trash_se = trashed[lib]
        seq_trash_se_frac = seq_trash_se / seq_reads_se

        #seq_retained_reads, seq_retained_nts, seq_retained_length,
        seq_retained_reads = retained[lib]
        seq_retained_nts = nucleotides[lib]
        seq_retained_length = seq_retained_nts / seq_retained_reads

        # hits_raw_endogenous, hits_raw_frac_endogenous, hits_clonality_endogenous,
        hits_raw_endogenous = rmdup[lib] + hits[lib]
        hits_raw_frac_endogenous = hits_raw_endogenous / seq_retained_reads
        hits_clonality_endogenous = rmdup[lib]

        # hits_unique_endogenous, hits_unique_frac_endogenous,
        # hits_coverage_endogenous, hits_length_endogenous,
        hits_unique_endogenous = hits[lib]
        hits_unique_frac_endogenous = hits_unique_endogenous/seq_retained_reads
        hits_coverage_endogenous = coverage[lib]
        hits_length_endogenous = lengths[lib]

        # hits_unique_mitochondrial, hits_unique_frac_mitochondrial,
        # hits_coverage_mitochondrial, hits_length_mitochondrial
        hits_unique_mitochondrial = mitoHits[lib]
        hits_unique_frac_mitochondrial = hits_unique_mitochondrial/seq_retained_reads
        hits_coverage_mitochondrial = mitoCov[lib]
        hits_length_mitochondrial = mitoLengths[lib]

        print_line(file)

file.close()