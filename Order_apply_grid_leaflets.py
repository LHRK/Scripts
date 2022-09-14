#!/usr/bin/env python

import numpy as np
import MDAnalysis as md

## Same as Order_apply_grid.py, excepts does it for each leaflet
## Loads in the x,y coordinates for each lipid and the corresponding order parameter
## Then divide it unto a grid, where you specify the number of bins

X_top = np.loadtxt('X_data_script_top.txt') 
Y_top = np.loadtxt('Y_data_script_top.txt') 
A_top = np.load('Order_all_ava_script_top.npy')


X_bot = np.loadtxt('X_data_script_bot.txt') 
Y_bot = np.loadtxt('Y_data_script_bot.txt') 
A_bot = np.load('Order_all_ava_script_bot.npy')

u = md.Universe('gromacs_files/6-prod-nvt.gro', 'gromacs_files/concat_nanodisc_1195frames_fit_2.xtc')
x_dim, y_dim, z_dim = u.dimensions[:3]

nbins = 10

x_steps = np.linspace(0, x_dim, num=nbins+1)
y_steps = np.linspace(0, y_dim, num=nbins+1)


def grid(nbins, X, Y, A):
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
    return grid

grid_top = grid(nbins, X_top, Y_top, A_top)
grid_bot = grid(nbins, X_bot, Y_bot, A_bot)

grid_ava = grid_top + grid_bot / 2

np.save('Order_Grid_{0:d}nbins_top_test.npy'.format(nbins), grid_top)
np.save('Order_Grid_{0:d}nbins_bot_test.npy'.format(nbins), grid_bot)
np.save('Order_Grid_{0:d}nbins_ava_test.npy'.format(nbins), grid_ava)
