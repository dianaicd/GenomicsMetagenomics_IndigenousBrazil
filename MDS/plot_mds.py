#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 11:04:03 2019

@author: Diana I. Cruz Davalos
"""
#%%
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
# %%
dist_path = "/Users/dcruz/Projects/Botocudos/Files/MDS/2019_02_05/Maanasa_mask1_flip_24ind.dist"

# %%
#botocudos = "/Users/dcruz/Projects/Botocudos/Files/MDS/2019_02_05/bam.filelist"
botocudos = "/Users/dcruz/Projects/Botocudos/Files/Panels/Maanasa_pop.txt"
f = open(botocudos)
label_given = []

for x in f:
    label_given.append(x.strip())
# %%
#label_given = panel
labelUnique,index = np.unique(label_given,return_inverse=True)
colors = cm.rainbow(np.linspace(0, 1, len(labelUnique)))
color_value = colors[index]

# %%
dist = np.loadtxt(dist_path)
#dist = dist[0:23, 0:23]
dist = dist[0:2598, 0:2598]
dist = dist + dist.T

#%%


plt.imshow(dist, vmin = 0, vmax = 1)
plt.colorbar()

#%%

empty = np.where(np.isnan(dist))
full = np.zeros(dist.shape)
full[empty] = 1

plt.imshow(full)

# %%
# Remove nan
individualToRemove = (5, 7, 9, 15, 21, 22)
sub_dist = np.delete(dist, individualToRemove, axis = 0)
sub_dist = np.delete(sub_dist, individualToRemove, axis = 1)

plt.imshow(sub_dist)

# %%
empty = np.where(np.isnan(sub_dist))
full = np.zeros(sub_dist.shape)
full[empty] = 1

plt.imshow(full)

#%%

lambdas,vect = calc_mds(dist)

# %%
plt.figure(figsize = (10.,10.))
plt.scatter(vect[:,1], vect[:,2], c = color_value)
plt.legend(labelUnique,scatterpoints=1)

#%%
plt.figure(figsize = (10.,10.))
for index_label, label in enumerate(labelUnique):
	position = np.where(np.array(label_given) == label)
	plt.scatter(vect[position,0], vect[position,1], c = colors[index_label],
             label = label )
#plt.legend(label_given)

# %%
# plot eigenvalues
 plot(-sort(-lambdas), "o")

 # %%
 # (should be) Marchenko Pastur

  hist(-sort(-lambdas))