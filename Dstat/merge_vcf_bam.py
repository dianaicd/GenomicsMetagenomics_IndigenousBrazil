#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 16:34:31 2018

@author: Diana I. Cruz Davalos

This script will merge bam samples to a VCF panel
"""

import sys
import getopt, pysam, re, random

# %%
# Parse arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],"htv:b:",["vcf=","bam="])
except getopt.GetoptError:
    print('sys.argv[0] -i <inputfile> -o <outputfile>')
    sys.exit(2)
print(opts)
print(args)
for opt, arg in opts:
    # path to vcf.gz file
    if opt in ("--vcf", "-v"):
        vcf_path = arg
    # file with paths to bam files
    elif opt in ("--bam", "-b"):
        bam_list = open(arg)
    # basename
    elif opt in ("--output", "-o"):
        out = sys.argv[3] + ".vcf"
        logfile = sys.argv[3] + ".log"
        logfile = open(logfile, "w")
        new_vcf = open(out, "w")
    elif opt == '-t':
        vcf_path = "/Users/dcruz/Projects/Botocudos/Files/test/Maanasa_mask1_flip.vcf.gz"
        bam_list = open("/Users/dcruz/Projects/Botocudos/Files/test/bamlist")
        logfile = open("test.log", "w")
        new_vcf = open("out.vcf", "w")
        

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
                Q = pileupread.alignment.query_qualities[
                        pileupread.query_position]
                base = pileupread.alignment.query_sequence[
                        pileupread.query_position]
                if  meets_req(base=base, baseQ=Q, ref=ref,
                              alt=alt, baseQ_thres=baseQ_thres):
                    # Debug: Are we looking at bases at the right position?
                    if base == ref:
                        bases.append("0/.")
                    else:
                        bases.append("1/.")

    if n_reads == 0 or len(bases) == 0:
        return "./."
    else:
        selected = random.randint(0, len(bases)-1)

    answer = re.sub("1", "1/.", str(bases[selected]))
    answer = re.sub("0", "0/.", str(bases[selected]))

    return answer


# %%
# Function to sample a random allele from a VCF file
# (one position in the genome)
# Returns one allele (0, 1 or None) indicating
# (reference, alternative, missing data)

def sample_allele_vcf(Chr,
                      vcf, anc_bases=[], recode=False, bases_ind=[],
                      positions={}
                      ):

    # This goes over the rows with data of a vcf (skipping header)
    index_pos = -1
    for record in vcf.fetch(contig=Chr):
        index_pos += 1
        index_ind = 0
        pos = str(Chr) + "_" + str(record.pos)
        # Skip sites where there is not info in the ancestral sample
        if pos not in positions:
            index_pos -= 1
            continue

        for ind in record.samples.itervalues():
            # print(ind.name)
            # 0/0 means it has two reference alleles
            # 1/1 means it has two alternative alleles

            selected = random.randint(0, 1)
            base = ind.values()[0][selected]

            # If by chance we have 1/. or 0/., try to
            # sample the non-missing allele
            if base is None:
                selected = abs(selected - 1)
                base = ind.values()[0][selected]
                if base is None:
                    base = 9

            # Recode according to ancestral state
            if recode:
                if base == positions[pos]:
                    base = 0
                else:
                    base = 1

            bases_ind[index_pos, index_ind] = base
            index_ind += 1

    return(bases_ind)


# %%
# VCF file
# do not forget to have the corresponding index:
# bgzip test.vcf
# bcftools index test.vcf.gz

vcf = pysam.VariantFile(vcf_path)

# %%
# Add samples to header

baseQ = 30
bam_paths = []

for line in bam_list:
    bam_paths.append(line)
    print(line)
header_copy = vcf.header.copy()

# Add samples' names to VCF header
for path in bam_paths:
    new = re.split("\.", re.split("/", path)[-1])[0]
    print(new)
    header_copy.add_sample(new)

new_vcf.write(str(header_copy))

# %%
# Print each vcf record with the new samples' genotypes
for contig in vcf.header.contigs:
    logfile.write(str(contig))
    i = 0
    for rec in vcf.fetch(contig=contig):
        i += 1
        if i == 0:
            logfile.write("10,000 lines written")
            i = 0
        pos = str(contig) + "_" + str(rec.pos)
        # VCF record to print
        txt = str(rec.copy())
        txt = re.sub("\n", "", txt)

        for path in bam_paths:
            samfile = pysam.AlignmentFile(path.strip(), "rb")
            sampled = sample_allele_bam(position=pos, samfile=samfile,
                                        ref=rec.alleles[0],
                                        alt=rec.alleles[1])
            txt = txt + "\t" + sampled
            samfile.close()

        txt = txt + "\n"
        new_vcf.write(txt)

vcf.close()
new_vcf.close()
