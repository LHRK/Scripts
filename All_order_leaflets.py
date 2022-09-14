#!/usr/bin/env python 
from MDAnalysis.analysis.leaflet import LeafletFinder
import numpy as np
import MDAnalysis as md
import matplotlib.pyplot as plt 
#from progressbar import *

## Same script as All_order.py, excepts it does it separtly for each leaflet


def Sn2_AA_per_lipid(lipid, resid, u):
    """Works for getting the weighted Scd for the Sn1 tail out for a specific lipid and timeframe"""
    carbons = u.select_atoms("resname {0:s} and name C3* and not name C3 and not name C31".format(lipid))
    lenght = len(np.unique(carbons.names))

    C3_order = []
    normal = [0,0,1]
    for i in range(2, lenght+2): #Loop over the carbon tail (C2* tail)
        #print (i)
        cp = u.select_atoms("resname {1:s} and name C3{0:d} and resid {2:d}".format(i, lipid, resid), updating=True)
        hx = u.select_atoms("resname {1:s} and name H{0:d}X and resid {2:d}".format(i, lipid, resid), updating=True)
        hy = u.select_atoms("resname {1:s} and name H{0:d}Y and resid {2:d}".format(i, lipid, resid), updating=True)
        hz = u.select_atoms("resname {1:s} and name H{0:d}Z and resid {2:d}".format(i, lipid, resid), updating=True)
        summ = 0 #Order sum
        nh = 0 #normalizing factor
        nres = cp.names.shape[0] #number of carbons
        cpx = cp.positions[:,0] #Positions for carbon
        cpy = cp.positions[:,1]
        cpz = cp.positions[:,2]

        if hx.names.shape[0] != 0:
            hxx = cpx - hx.positions[:,0] # the vector describing C-H
            hxy = cpy - hx.positions[:,1]
            hxz = cpz - hx.positions[:,2]
            for idx, x in enumerate(hxx): #Calculate order for HX atoms
                norm2 = hxx[idx]**2 + hxy[idx]**2 + hxz[idx]**2
                summ += hxz[idx]**2 / norm2
            nh += nres
        if hy.names.shape[0] != 0: #Calculate order for HY atoms
            hyx = cpx - hy.positions[:,0]
            hyy = cpy - hy.positions[:,1]
            hyz = cpz - hy.positions[:,2]
            for idx, x in enumerate(hyx):
                norm2 = hyx[idx]**2 + hyy[idx]**2 + hyz[idx]**2
                summ += hyz[idx]**2 / norm2
            nh += nres
        if hz.names.shape[0] != 0: #Calculate order for HZ atoms
            hzx = cpx - hz.positions[:,0]
            hzy = cpy - hz.positions[:,1]
            hzz = cpz - hz.positions[:,2]
            for idx, x in enumerate(hzx):
                norm2 = hzx[idx]**2 + hzy[idx]**2 + hzz[idx]**2
                summ += hzz[idx]**2 / norm2
            nh += nres
        order = -1.5 * (summ/(nh)) + 0.5
        C3_order.append(order)
    return C3_order


def Sn1_AA_per_lipid(lipid, resid, u):
    """Works for getting the weighted Scd for the Sn1 tail out for a specific lipid and timeframe"""
    carbons = u.select_atoms("resname {0:s} and name C2* and not name C2 and not name C21".format(lipid))
    lenght = len(np.unique(carbons.names))

    C2_order = []
    normal = [0,0,1]
    for i in range(2, lenght+2): #Loop over the carbon tail (C2* tail)
        #print (i)
        cp = u.select_atoms("resname {1:s} and name C2{0:d} and resid {2:d}".format(i, lipid, resid), updating=True)
        hx = u.select_atoms("resname {1:s} and name H{0:d}R and resid {2:d}".format(i, lipid, resid), updating=True)
        hy = u.select_atoms("resname {1:s} and name H{0:d}S and resid {2:d}".format(i, lipid, resid), updating=True)
        hz = u.select_atoms("resname {1:s} and name H{0:d}T and resid {2:d}".format(i, lipid, resid), updating=True)
        h9 = u.select_atoms("resname {1:s} and name H{0:d}1 and resid {2:d}".format(i, lipid, resid), updating=True)
        summ = 0 #Order sum
        nh = 0 #normalizing factor
        nres = cp.names.shape[0] #number of carbons
        cpx = cp.positions[:,0] #Positions for carbon
        cpy = cp.positions[:,1]
        cpz = cp.positions[:,2]

        if hx.names.shape[0] != 0:
            hxx = cpx - hx.positions[:,0] # the vector describing C-H
            hxy = cpy - hx.positions[:,1]
            hxz = cpz - hx.positions[:,2]
            for idx, x in enumerate(hxx): #Calculate order for HX atoms
                norm2 = hxx[idx]**2 + hxy[idx]**2 + hxz[idx]**2
                summ += hxz[idx]**2 / norm2
            nh += nres
        if hy.names.shape[0] != 0: #Calculate order for HY atoms
            hyx = cpx - hy.positions[:,0]
            hyy = cpy - hy.positions[:,1]
            hyz = cpz - hy.positions[:,2]
            for idx, x in enumerate(hyx):
                norm2 = hyx[idx]**2 + hyy[idx]**2 + hyz[idx]**2
                summ += hyz[idx]**2 / norm2
            nh += nres
        if h9.names.shape[0] != 0: #Calculate order for HY atoms
            h9x = cpx - h9.positions[:,0]
            h9y = cpy - h9.positions[:,1]
            h9z = cpz - h9.positions[:,2]
            for idx, x in enumerate(h9x):
                norm2 = h9x[idx]**2 + h9y[idx]**2 + h9z[idx]**2
                summ += h9z[idx]**2 / norm2
            nh += nres
        if hz.names.shape[0] != 0: #Calculate order for HZ atoms
            hzx = cpx - hz.positions[:,0]
            hzy = cpy - hz.positions[:,1]
            hzz = cpz - hz.positions[:,2]
            for idx, x in enumerate(hzx):
                norm2 = hzx[idx]**2 + hzy[idx]**2 + hzz[idx]**2
                summ += hzz[idx]**2 / norm2
            nh += nres
        order = -1.5 * (summ/(nh)) + 0.5
        C2_order.append(order)
    return C2_order


def order (u, nframes, sel):
    nresid  = u.select_atoms('resname DMPC and name P and resid {0:s}'.format(sel)).resids.shape[0]
    resids  = u.select_atoms('resname DMPC and name P and resid {0:s}'.format(sel)).resids

    nresid = len(resids)
    Output            = np.zeros([13, nresid, nframes])
    Output_ava        = np.zeros([nresid, nframes])
    X = []
    Y = []
    #Looping over each frame, then each grid bin, selecting the lipids in each bin. 
    #For each bin the order parameter is calculated for each tail of the DMPC lipid
    #An averaged over the two tails, the number of lipids in each bin is made for each frame and saved to the output array 
    for idx, ts in enumerate(u.trajectory[1:]):
    	X_loc = []
    	Y_loc = []
    	frame = u.trajectory.frame
    	for ndx, r in enumerate(resids):
    		sel = u.select_atoms('resname DMPC and resid {0:d} and name P'.format(r))
    		x,y,z = sel.positions[0]
    		X_loc.append(x)
    		Y_loc.append(y)
    		 
    		C2_order  = np.array(Sn1_AA_per_lipid('DMPC', r, u)) 
    		C3_order  = np.array(Sn2_AA_per_lipid('DMPC', r, u))
    		C_profile = np.average(np.vstack((C2_order, C3_order)), axis=0)
    		C_order   = np.average(C_profile)
    		#print ('Order parameters calculated and assigned to array')
    		Output[:,ndx,idx]= C_profile #Take the weighted averaged afterwards
    		Output_ava[ndx,idx] = C_order
    	X.append(X_loc)
    	Y.append(Y_loc)
    return X, Y, Output , Output_ava

weight = np.loadtxt('clustering/QT_weights.txt', comments='#')[:,1]
u = md.Universe('gromacs_files/6-prod-nvt.gro', 'gromacs_files/concat_nanodisc_1195frames_fit_2.xtc')
print ('Trajectory: gromacs_files/concat_nanodisc_1195frames_fit_cluster2.xtc')
print ('Topology: gromacs_files/6-prod-nvt.gro')
print ('Universe loaded')
nframes = len(u.trajectory[1:]) #Skipping the gro file


L = LeafletFinder(u, 'resname DMPC and name P')
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)
top_leaf = " ".join(['{0:d}'.format(i) for i in leaflet0.resids ])
bot_leaf = " ".join(['{0:d}'.format(i) for i in leaflet1.resids ])


X_top, Y_top, Output_top, Output_ava_top = order (u, nframes, top_leaf)
X_bot, Y_bot, Output_bot, Output_ava_bot = order (u, nframes, bot_leaf)

np.save('Order_all_profile_script_top.npy', Output_top)
np.save('Order_all_ava_script_top.npy', Output_ava_top)
np.savetxt('X_data_script_top.txt', np.array(X_top))
np.savetxt('Y_data_script_top.txt', np.array(Y_top))


np.save('Order_all_profile_script_bot.npy', Output_bot)
np.save('Order_all_ava_script_bot.npy', Output_ava_bot)
np.savetxt('X_data_script_bot.txt', np.array(X_bot))
np.savetxt('Y_data_script_bot.txt', np.array(Y_bot))

print ('Calculation done! Now saving to files')
