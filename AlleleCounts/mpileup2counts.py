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

def parse_line(pileup):
    # remove missing, beginning, end
    #pattern = "\*|\^.|\$"
    pattern = "\^."
    result = re.sub(pattern, "", pileup)
    pattern = re.compile(r"\+\d+|\-\d+")
    matches = [int(s) for s in re.findall(pattern, result)]
    # remove indels
    for p in matches:
        pattern = re.compile("[+|-]" + str(abs(p)) + ".{" + str(abs(p)) +"}")
        result = re.sub(pattern, "", result)

    return(result)

# %%
pileup = file.readline()
pileup = parse_line(pileup)
nInd = (len(pileup.split("\t")) - 3 ) / 3
print(nInd)

parsed_file = [parse_line(line) for line in file.readlines()]

