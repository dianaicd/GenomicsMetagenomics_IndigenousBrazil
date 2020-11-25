# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import multiprocessing as mp


print(f"CPU Count: {mp.cpu_count()}")

# %%
# Define some paths and prefixes

# The path before the .ind, .snp, .geno
target_ind_prefix               = './v37.2.1240K_HumanOrigins/Jan2020'
ind_id                          = 'MN00056'
chromosomes                     = range(1, 23)
haplotypes_prefix               = './all1240/chr'
metadata_path                   = './all1240/meta_df_all.csv'
folder_out                      = './test_roh/'
prefix_out                      = ''
emission_model                  = 'haploid'
preprocessing_model             = 'EigenstratPacked'
postprocessing_model            = 'Standard'
num_processes                   = 6
delete_raw_posterior            = False
print_output                    = True 
save_ROHs                       = True
save_full_posterior             = False
n_diploid_ref_to_use            = 2504
pops_to_exclude                 = []
load_readcounts                 = True
pick_random_allele              = True
trans_rate_into_ROH             = 1
trans_rate_out_ROH              = 20
trans_rate_within_ROH           = 300
error_rate                      = 0.01
error_rate_in_reference         = 0.0
cutoff_posterior                = 0.999
max_gap_to_merge                = 0
roh_min_length                  = 0.01
use_logfile                     = False
combine_output_all_chrs         = False
appendix_results_per_ind        = '_roh_full.csv'

# %%
### If wanting to use local version and not  pip installed version
#sys.path.append("./package/") # Append local Hapsburg Folder

# %%
# Set path

path = "/Users/dcruz/Projects/Botocudos/Files/Heterozygosity/2020_07_08/hapROH/"  
os.chdir(path)  
print(f"Set path to: {os.getcwd()}") 

# %%
from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind  # Need this import

# %%
# Test calling ROH on single chromosome from Individual with full output printed

hapsb_ind(iid           = ind_id, 
          chs           = chromosomes, 
          path_targets  = target_ind_prefix, 
          h5_path1000g  = haplotypes_prefix, 
          meta_path_ref = metadata_path, 
          folder_out    = folder_out, 
          prefix_out    = prefix_out, 
          e_model       = emission_model, 
          p_model       = preprocessing_model, 
          post_model    = postprocessing_model, 
          processes     = num_processes,
          delete        = delete_raw_posterior,
          output        = print_output, 
          save          = save_ROHs, 
          save_fp       = save_full_posterior, 
          n_ref         = n_diploid_ref_to_use, 
          exclude_pops  = pops_to_exclude,
          readcounts    = load_readcounts,
          random_allele = pick_random_allele, 
          roh_in        = trans_rate_into_ROH, 
          roh_out       = trans_rate_out_ROH, 
          roh_jump      = trans_rate_within_ROH, 
          e_rate        = error_rate, 
          e_rate_ref    = error_rate_in_reference, 
          cutoff_post   = cutoff_posterior, 
          max_gap       = max_gap_to_merge, 
          roh_min_l     = roh_min_length, 
          logfile       = use_logfile, 
          combine       = combine_output_all_chrs, 
          file_result   = appendix_results_per_ind)

# %%
# Example run of whole individual (all chromosomes) with output to logfile
"""

This example runs a whole individual in parallel, 
with the output send to a logfile

Attention: hapROH uses quite a bit of RAM, 
for an Eigenstrat it can be spiking up to 6gb for a long chromosome - 
allocate memory accoringly
when running multiple processes, or set that number lower!

"""

chromosomes             = range(1, 23)
load_readcounts         = False
use_logfile             = True
combine_output_all_chrs = True

hapsb_ind(iid               = ind_id, 
          chs               = chromosomes, 
          processes         = num_processes, 
          path_targets      = target_ind_prefix, 
          h5_path1000g      = haplotypes_prefix, 
          meta_path_ref     = metadata_path, 
          folder_out        = folder_out, 
          prefix_out        = prefix_out, 
          e_model           = emission_model, 
          p_model           = preprocessing_model, 
          n_ref             = n_diploid_ref_to_use,
          random_allele     = pick_random_allele,
          readcounts        = load_readcounts,
          delete            = delete_raw_posterior,
          logfile           = use_logfile,
          combine           = combine_output_all_chrs)

# %%
# Run multiple individuals, all chromosomes

### This are all individuals with 400k SNPs covered
iids = ['I1178', 'I0644', 'I1160', 'I1152', 'I1168', 'I1166', 'I1170', 'I1165', 'I1182']
load_readcounts = False

for iid in iids:
    print(f"Doing Individual: {iid}")
    hapsb_ind(iid           = iid, 
              chs           = chromosomes, 
              processes     = num_processes, 
              path_targets  = target_ind_prefix, 
              h5_path1000g  = haplotypes_prefix, 
              meta_path_ref = metadata_path, 
              folder_out    = folder_out,
              prefix_out    = prefix_out, 
              e_model       = emission_model,
              p_model       = preprocessing_model,
              n_ref         = n_diploid_ref_to_use,
              random_allele = pick_random_allele, 
              readcounts    = load_readcounts,
              delete        = delete_raw_posterior, 
              logfile       = use_logfile, 
              combine       = combine_output_all_chrs)

# %%
#  Postprocess Results into one results.csv (copying in Meta Data)
"""
Take indivdiual output .csvs and combine into one big results .csv
Merging of output gaps, and ROH>x cM happens here
"""
from hapsburg.PackagesSupport.pp_individual_roh_csvs import pp_individual_roh



# %%
# Create example metafile
iids = ['I1178', 'I0644', 'I1160', 'I1152', 'I1168', 'I1166', 'I1170', 'I1165', 'I1182']
df = pd.DataFrame({"iid":iids})

df["age"] = 5950
df["clst"] = "Israel_C"
df["lat"] = 32.974167
df["lon"] = 35.331389

df.to_csv("./Data/ExampleData/meta_blank.csv", 
          sep=",", index=False)

# %%
# Combine into output file
%%time
### Postprocess the two Individuals from above and combine into one results .csv
iids = ['I1178', 'I0644', 'I1160', 'I1152', 
        'I1168', 'I1166', 'I1170', 
        'I1165', 'I1182']

df1 = pp_individual_roh(iids, 
                        meta_path="./Data/ExampleData/meta_blank.csv", 
                        base_folder="./Empirical/Eigenstrat/Example/",
                        save_path="./Empirical/Eigenstrat/Example/combined_roh05.csv", 
                        output=False, 
                        min_cm=[4, 8, 12, 20], 
                        snp_cm=50, 
                        gap=0.5, 
                        min_len1=2.0, 
                        min_len2=4.0)

