#!/usr/bin/env python

import numpy as np
import MDAnalysis as md

## Loads in the x,y coordinates for each lipid and the corresponding order parameter
## Then divide it unto a grid, where you specify the number of bins



X = np.loadtxt('X_data_script2.txt') -10
Y = np.loadtxt('Y_data_script2.txt') -10
A = np.load('Order_all_ava_script2.npy')


u = md.Universe('gromacs_files/6-prod-nvt.gro', 'gromacs_files/concat_nanodisc_1195frames_fit_2.xtc')
x_dim, y_dim, z_dim = u.dimensions[:3]

nbins = 10

x_steps = np.linspace(0, x_dim, num=nbins+1)
y_steps = np.linspace(0, y_dim, num=nbins+1)



nframes = 1195
grid = np.zeros([nbins,nbins, nframes])

for f in range(nframes):
    for i in range(nbins):
        for j in range(nbins):
            x = X[f]
            y = Y[f]
            bool_idx = [(x>=x_steps[i]) & (x<x_steps[i+1]) & (y>=y_steps[j]) & (y<y_steps[j+1])]
            idx = np.where(bool_idx[0] ==True)
            ava_bin = np.average(A[idx,f])
            if np.isnan(ava_bin):
                grid[i,j,f] = False
            else:
                grid[i,j,f] = ava_bin

np.save('Order_Grid_{0:d}nbins_test.npy'.format(nbins), grid)
