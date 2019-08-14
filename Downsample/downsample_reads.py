#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 13:58:43 2019

@author: Diana I. Cruz Davalos
"""
#%%
# Trying to downsample reads from a bam file
# It is compulsory to have the idxstats by samtools

import pysam
import numpy as np
import datetime
import sys, subprocess, getopt, os.path

# %%
options, remainder = getopt.getopt(sys.argv[1:],
                                   'b:o:n:f:i:',
                                   ['bam=', 'output=', 
                                   'num_reads=', 'frac_reads=','idxstats='])
                                
bam_in = False
bam_out = False 
f_reads = False

for opt, arg in options:
    if opt in ('-b', '--bam'):
        bam_in = arg
        if not os.path.isfile(bam_in):
            print("Input bam file does not exist.")
            sys.exist()
        idxstats = bam_in.replace(".bam", "_idxstats.txt")
    elif opt in ('-o', '--output'):
        bam_out = arg
    elif opt in ('-n', '--num_reads'):
        n_reads = int(arg)
    elif opt in ('-f', '--frac_reads'):
        f_reads = arg
    elif opt in ('-i', '--idxstats'):
        idxstats = arg

# %%
# Check if idxstats files exist
if not os.path.isfile(idxstats):
    bai = bam_in + ".bai"
    if not os.path.isfile(bai):
        print("You do not have the index file. I will samtools index for you, poor student.")
        command = "samtools index " + bam_in
        subprocess.Popen(command, shell=True, 
                                  stderr = open("err_index.txt", 'w'))
                                  
    command = "samtools idxstats " + bam_in
    f = open(idxstats, 'w')
    subprocess.Popen(command, shell=True, 
                                  stderr = open("err_idxstats.txt", 'w'),
                                  stdout = f)
    f.close()
    

# %%
total_reads = 0
with open(idxstats, 'r') as file:
    for line in file.readlines():
        total_reads += int(line.split('\t')[2])

file.close()

#%%
my_msg = bam_in + " has " + str(total_reads) + " reads."
print(my_msg)
my_msg = "I will downsample to " + str(n_reads) + " reads."
print(my_msg)
selected_lines = np.zeros(total_reads,dtype=int)
chosen = np.random.choice(a=total_reads, size = n_reads, replace = False)
selected_lines[chosen] = 1

# %%
i = 0
a = datetime.datetime.now()

if not bam_out:
    bam_out = "dummy_downsampled.bam"

samfile = pysam.AlignmentFile(bam_in, "rb")
downsampled = pysam.AlignmentFile(bam_out, "wb", template=samfile)
samfile.close()

samfile = pysam.AlignmentFile("-", "rb")

for read in samfile:
    if(selected_lines[i]):
        downsampled.write(read)
    i+=1
samfile.close()
