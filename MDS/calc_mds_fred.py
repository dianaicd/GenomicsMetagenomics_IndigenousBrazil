#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 10:26:35 2019

@author: Diana I. Cruz Davalos
"""
# %%
def calc_mds(delta) -> None:
    """Take the distance matrix delta and return the coordinate
	and the eigenvalues from the mds
    """

    N = len(delta)

    a = -delta**2/2

    at0 = np.nansum(a, 0)/N
    att = np.nansum(a)/N**2
    
    one = np.ones((N,))
    b = a - np.outer(at0, one) - np.outer(one, at0) + att
    print(b)
    plt.imshow(b)
    lambdas, vecs = np.linalg.eigh(b)
    print(lambdas)
    print(vecs)
    print(b.shape)
    N = len(lambdas)

    coordinates = np.hstack((lambdas.reshape((N, 1)), vecs.T)).copy()
    coordinates = sorted(coordinates, key=lambda x: x[0], reverse=True)

    for i, v in enumerate(coordinates.copy()):
        coordinates[i] = np.sqrt(v[0])*v[1:]

    coordinates = np.array(coordinates)
    coordinates = coordinates.T
    return(lambdas, coordinates)