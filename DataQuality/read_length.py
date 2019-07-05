#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:35:40 2019

@author: dcruz
"""

import pysam
import sys, getopt
# %%
#path = "/Users/dcruz/Projects/Botocudos/Files/test/MN1943/MN1943.bam"

# File from stream
path = "-"
options, remainder = getopt.getopt(sys.argv[1:],
                                   'o:b:', ['output=', 'bam'])
for opt, arg in options:
    if opt in ('-o', '--output'):
        output = arg
    elif opt in ('-b', 'bam'):
        path = arg

# %%
infile = pysam.AlignmentFile(path, "rb")

#%%
myLengths = {}

for s in infile:
    l = s.infer_read_length()
    if l in myLengths:
        myLengths[l] += 1
    else:
        myLengths[l] = 1

# %%
with open(output, "w") as file:
    for key,value in myLengths.items():
        file.write("{0}\t{1}\n".format(key,value))
    file.close()