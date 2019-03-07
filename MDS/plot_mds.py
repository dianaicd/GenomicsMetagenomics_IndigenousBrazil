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
import operator
import re
# %%
# Read in distance file and labels
dist_path = "/Users/dcruz/Projects/Botocudos/Files/MDS/2019_02_05/Botocudos_fromVictor.Maanasa_americas_reheaded_filtered.dist"

#botocudos = "/Users/dcruz/Projects/Botocudos/Files/MDS/2019_02_05/bam.filelist"
#botocudos = "/Users/dcruz/Projects/Botocudos/Files/Panels/Maanasa_pop.txt"
#botocudos = "/Users/dcruz/Projects/Botocudos/Files/MDS/2019_02_05/Maanasa_americas.inds"
botocudos = "/Users/dcruz/Projects/Botocudos/Files/MDS/2019_02_05/Botocudos_fromVictor_Maanasa.ind"

f = open(botocudos)
label_given = []

for x in f:
    label_given.append(x.strip())

#%%
# Colors for the plot
label_given = [re.sub("MN.*", "Botocudos", s) for s in label_given]
label_given = [re.sub("LS.*", "Lagoa Santa", s) for s in label_given]

labelUnique,index = np.unique(label_given,return_inverse=True)
colors = cm.rainbow(np.linspace(0, 1, len(labelUnique)))
color_value = colors[index]

# %%
dist = np.loadtxt(dist_path)

# Bug in previous version
if not dist.shape[0] == dist.shape[1]: 
    dist = dist[0:dist.shape[0], 0:dist.shape[0]]


#%%
plt.imshow(dist, vmin = 0, vmax = 1)
plt.colorbar()

# %%
# Remove nan

empty = np.where(np.isnan(dist))

unique,counts = np.unique(empty, return_counts = True)
occur = dict(zip(unique, counts))
individualToRemove = []

toremove = sorted(occur.items(), key=operator.itemgetter(1), reverse=True)
for key,value in toremove:

    if(np.sum(np.isnan(dist[:, key]))>0):
        dist[key, :] = 0
        individualToRemove.append(key)

for ind in sorted(individualToRemove, reverse = True):
    del label_given[ind]

#%%
individualToRemove = sort(individualToRemove)

sub_dist = np.delete(dist, individualToRemove, axis = 1)
sub_dist = np.delete(sub_dist, individualToRemove, axis = 0)


#%%
full = np.zeros(dist.shape)
full[empty] = 1

plt.imshow(full)
plt.colorbar()
#%%
plt.imshow(dist)
plt.colorbar()
#%%
plt.imshow(sub_dist)
plt.colorbar()
# %%
sub_dist = sub_dist + sub_dist.T
#%%

lambdas,vect = calc_mds(sub_dist)

# %%
plt.figure(figsize = (10.,10.))
plt.scatter(vect[:,0], vect[:,1], c = color_value)
plt.legend(labels=labelUnique,scatterpoints=1)

#%%
plt.figure(figsize = (10.,10.))
for index_label, label in enumerate(labelUnique):
	position = np.where(np.array(label_given) == label)
	plt.scatter(vect[position,0], vect[position,1], c = colors[index_label],
             label = label)
position = np.where(np.array(label_given) == "Karitiana")
plt.scatter(vect[position,0], vect[position,1], c = "green",
             label = label )
position = np.where(np.array(label_given) == "Surui")
plt.scatter(vect[position,0], vect[position,1], c = "pink",
             label = label )
position = np.where(np.array(label_given) == "Bajo")
plt.scatter(vect[position,0], vect[position,1], c = "red",
             label = label )

plt.legend(labelUnique)
# %%
# plot eigenvalues
plt.plot((lambdas), "o")

 # %%
 # (should be) Marchenko Pastur
plt.hist(-sort(-lambdas))