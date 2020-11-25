#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:58:42 2020

@author: dcruz
"""


import os as os
import pandas as pd

path = "/Users/dcruz/Projects/Botocudos/Files/Heterozygosity/2020_07_08/hapROH/"  

os.chdir(path)  
print(f"Set path to: {os.getcwd()}") 

# %%

folder_results = './test_roh/MN00056/chr22/'
savepath = ''
plot_latent_states = True
num_reads = 1
plot_groundtruth = False
cm_lim = []
readcount = False
figsize = (10, 4)
min_cm = 1

# %%
# Plot the Posterior along one Chromosome. 
"""
Needs the posterior to be saved,
(so in hapsb_ind use delete=False to keep it)
Set folder to the chromosome output you want to plot
"""

from hapsburg.figures.plot_posterior import plot_posterior_cm

title = "Chromosome "

plot_posterior_cm(folder        = folder_results, 
                  savepath      = savepath, 
                  empirical     = plot_latent_states,
                  m             = num_reads, 
                  cm_lim        = cm_lim,
                  groundtruth   = plot_groundtruth,
                  min_cm        = min_cm,
                  readcount     = readcount, 
                  figsize       = figsize,
                  title         = title
                  )

# %%
# Zoom in
cm_lim = [ 0, 50 ]

plot_posterior_cm(folder        = folder_results, 
                  savepath      = savepath, 
                  empirical     = plot_latent_states, 
                  m             = num_reads, 
                  cm_lim        = cm_lim,
                  groundtruth   = plot_groundtruth,
                  min_cm        = min_cm,
                  readcount     = readcount, 
                  figsize       = figsize, 
                  title         = title
                  )
# %%
# Plot karyotipe with ROHs 
from hapsburg.figures.plot_individual_roh import plot_roh_individual

prefix_output = ""
min_cm = 4
folder_results = './test_roh/'


plot_roh_individual(iid         ="MN00056", 
                    folder      = folder_results, 
                    prefix_out  = prefix_output, 
                    min_cm      = min_cm, 
                    plot_bad    = False)  # MA89

# %%
# Do histogram of the length distribution and theoretical expectations
from hapsburg.figures.plot_individual_roh import plot_pde_indivdiual

min_cm = 3
bin_width_cm = 4
plotlim = [4, 60]
plot_pde_indivdiual(iid         ="MN00056", 
                    min_cm      = min_cm, 
                    bw_cm       = bin_width_cm, 
                    kde_plot    = False, 
                    plotlim     = plotlim, 
                    prefix_out  = prefix_output,
                    savepath    = savepath, 
                    folder      = folder_results)

# %%
# Plot summary over multiple individuals
from hapsburg.figures.plot_bars import plot_legend_only, plot_panel_row, prepare_dfs_plot

figsize = (4,3.5)

plot_legend_only(savepath   = savepath, 
                 figsize    = figsize)

#%% 
# plot an empirical dataset

df1 = pd.read_csv("./Empirical/Eigenstrat/Example/combined_roh05.csv", sep='\t')
plot_dfs, cols = prepare_dfs_plot(df1, cms=[4, 8, 12, 20])

# ./figures_test/freilich20_bars.pdf
plot_panel_row(plot_dfs, 
               savepath ="", 
               wspace   =0.1, 
               r_title  =20, 
               leg_pos  =-2, 
               ylim     =[0,750], 
               figsize  =(8,3.5))

# %%
# plot timelines

from hapsburg.figures.plot_timelines import plot_map_time, extract_pop, prep_label

# one example, using a metafile that has the information
pop = "Rome"
df1 = pd.read_csv("./Empirical/roh_all_inds_final_v42.csv", sep="\t")
df_t = extract_pop(df1, age_range=[0,12200], pop=pop) # Cut out the right samples
label = prep_label(df_t, pop)

plot_map_time(df_t, 
              figsize=(3.6 , 1.1), 
              crs_m=[28, 63, -11, 38], 
              width_ratios=(8, 20),
              height_ratios=[15, 1],
              hspace=0.06, 
              wspace=0.015,
              s_tl=12, 
              ec="k", 
              lw=0.09, 
              x_lim_tl=[-500, 12200], 
              vrange_m=[0,12200], 
              y_lim_tl=[0,220], 
              fsl=5,
              fs=5, 
              fs_leg=5,
              leg_loc_tl="", 
              title_tl="",
              map_title=label, 
              title_loc=(0.2,0.01), 
              cm=4,
              cm1=8, 
              frac=0, 
              lgth_s=[1500,1500], 
              bottomrow=False,
              rightcol=True, 
              lw_fit=1.5, 
              ticks=[158.74, 83.82, 21.84], 
              tick_l=[f"$2N_e=250$", f"$2N_e=500$", f"$2N_e=2000$"], 
              width_t=0.6,
              length_t=2, 
              lbl_pad_time=5,
              lbl_pad_age=0, 
              xl_pad=1.5, 
              yl_pad=1, 
              widths=800, 
              alpha_vio=0.4, 
              savepath="")