#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 11:04:03 2019

@author: Diana I. Cruz Davalos
"""
#%%
import matplotlib.cm as cm

# %%
label_given = panel
labelUnique,index = np.unique(label_given,return_inverse=True)
colors = cm.rainbow(np.linspace(0, 1, len(labelUnique)))
color_value = colors[index]
# %%
plt.figure(figsize = (10.,10.))
plt.scatter(vect[:,1], vect[:,2], c = color_value)
plt.legend(labelUnique,scatterpoints=1)

#%%
plt.figure(figsize = (10.,10.))
for index_label, label in enumerate(labelUnique):
	position = np.where(np.array(label_given) == label)
	plt.scatter(vect[position,0], vect[position,1], c = colors[index_label],label = label )
plt.legend()