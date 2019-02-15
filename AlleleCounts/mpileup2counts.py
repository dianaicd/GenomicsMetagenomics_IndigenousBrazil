#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 23:36:21 2019

@author: dcruz
"""
# %%
import re

# %%
path = "/Users/dcruz/Projects/Botocudos/Files/test/test_mpileup.txt"
pileup = "*ggG$tta^cC***t*+3gTAtta-5TaCCTAG"
file = open(path, 'r')

# %%

def parse_line(pileup, nInd):
    # remove missing, beginning, end
    #pattern = "\*|\^.|\$"
    pattern = "\^."
    parsed = re.sub(pattern, "", pileup)
    pattern = re.compile(r"\+\d+|\-\d+")
    matches = [int(s) for s in re.findall(pattern, parsed)]
    # remove indels
    for p in matches:
        pattern = re.compile("[+|-]" + str(abs(p)) + ".{" + str(abs(p)) +"}")
        paarsed = re.sub(pattern, "", parsed)
    
    # Count matches
    allele_counts = []
    for ind in range(0, nInd):
        i = 3*ind + 4
        for base in ["A", "C", "G", "T"]:
            count = len(tuple(re.finditer(base, parsed.split("\t")[i], flags = re.I)))
            allele_counts.append(count)

    return(allele_counts)

# %%
pileup = file.readline()
pileup = parse_line(pileup, 2)
#nInd = (len(pileup.split("\t")) - 3 ) / 3
#print(nInd)

parsed_file = [parse_line(line, 2) for line in file.readlines()]

