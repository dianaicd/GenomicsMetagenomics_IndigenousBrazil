import sys,getopt
import numpy as np 
from io import StringIO

options, remainder = getopt.getopt(sys.argv[1:], 'i:l:d:', ['idxstats=',
                                                         'length=',
                                                         'depth='])

for opt, arg in options:
    if opt in ('-i', '--idxstats'):
        idxstats = arg
    elif opt in ('-l', '--length'):
        length = arg
    elif opt in ('-d', '--depth'):
        depth = float(arg)

# %%
with open(idxstats, 'r') as file:
    genome_len = np.nansum([np.genfromtxt(StringIO(line))[1] for line in file.readlines()])
# Number of mapped bases
mapped_len = np.genfromtxt(length, delimiter = '\t')
total_reads = np.nansum(mapped_len[:,1])
total_bases = np.nansum(mapped_len[:,0]*total_reads)
coverage = total_bases/genome_len

reads_to_sample = int(depth*total_reads/coverage)
        #np.savetxt(fname = output, X = reads_to_sample, fmt = "%.5f %d")
print(reads_to_sample)