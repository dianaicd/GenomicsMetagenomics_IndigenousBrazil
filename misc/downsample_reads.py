#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 13:58:43 2019

@author: Diana I. Cruz Davalos
"""
#%%

import pysam
import numpy as np
import datetime
import sys, os, subprocess 
# %%

# file with 1529173 reads
path = sys.argv[1] 
#path = "/Users/dcruz/Projects/Botocudos/Files/test/test_1.bam"

down_path = sys.argv[2] 
#down_path = "/Users/dcruz/Projects/Botocudos/Files/test/test_downsampled.bam"

total_selected = int(sys.argv[3])
#total_selected = 100000

command = "samtools idxstats " + path + " | awk -F '\t' '{s+=$3+$4}END{print s}'"
total_size = int(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().rstrip())

#total_size = 1529173
#print(total_selected) 
#%%
selected_sites = np.zeros(total_size,dtype=int)
position1 = np.random.choice(a=total_size,size = total_selected, replace = False)
selected_sites[position1] = 1

#%%

samfile = pysam.AlignmentFile(path, "rb")
total_downsampled = sum(selected_sites)

print(total_downsampled)
# %%
i = 0
a = datetime.datetime.now()

downsampled = pysam.AlignmentFile(down_path, "wb", template=samfile)
for i,read in enumerate(samfile):
    if(selected_sites[i]):
        downsampled.write(read)
    i = i +1 
    
print(datetime.datetime.now()-a)

downsampled.close()
samfile.close()
#%%

