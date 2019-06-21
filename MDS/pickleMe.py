#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 11:35:00 2019

@author: Diana I. Cruz Davalos
"""

import pickle

#%%
delta_path = ""

print('ARGV      :', sys.argv[1:])

options, remainder = getopt.getopt(sys.argv[1:], 'd:', ['distance='])
print('OPTIONS   :', options)

for opt, arg in options:
    if opt in ('-d', '--distance'):
        delta_path = arg
        delta = delta_path + ".pckl"

# %%
with open(delta_path, 'wb') as f:
	pickle.dump(delta, f)