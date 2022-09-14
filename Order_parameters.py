#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as md


def Order_sn2_AA (u, lipid):
	"""Calculate the order parameters for the sn2 carbon tail, for all-atom.
	Takes the universe (MDAnalysis) and the lipid resname as input.
	Returns a list of the averaged order parameters for each carbon,
	averaged over lipids and frames.
	WARNING: It does not reproduce the correct order parameter for double bonds (unsaturation)"""
	nframes = u.trajectory.n_frames
	carbons = u.select_atoms("resname {0:s} and name C2* and not name C2 and not name C21".format(lipid))
	lenght = len(np.unique(carbons.names))
	
	C2_order = []
	normal = [0,0,1]
	for i in range(2, lenght+2): #Loop over the carbon tail (C2* tail)
	    cp = u.select_atoms("resname {1:s} and name C2{0:d}".format(i, lipid), updating=True)
	    hx = u.select_atoms("resname {1:s} and name H{0:d}R".format(i, lipid), updating=True)
	    hy = u.select_atoms("resname {1:s} and name H{0:d}S".format(i, lipid), updating=True)
	    hz = u.select_atoms("resname {1:s} and name H{0:d}T".format(i, lipid), updating=True)
	    h9 = u.select_atoms("resname {1:s} and name H{0:d}1".format(i, lipid), updating=True)
	    summ = 0 #Order sum
	    nh = 0 #normalizing factor
	    nres = cp.names.shape[0] #number of carbons 
	    frames = len(u.trajectory[::1])
	    for ts in u.trajectory[::1]: #Loop over trajectory
	        cpx = cp.positions[:,0] #Positions for carbon
	        cpy = cp.positions[:,1]
	        cpz = cp.positions[:,2]
	
	        if hx.names.shape[0] != 0:
	            hxx = cpx - hx.positions[:,0] # the vector describing C-H
	            hxy = cpy - hx.positions[:,1]
	            hxz = cpz - hx.positions[:,2]
	            for idx, x in enumerate(hxx): #Calculate order for HX atoms
	                # P2 = 0.5∙(3∙(𝑐𝑜𝑠)^2 〈𝜃〉−1)
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


def Order_sn1_AA (u, lipid):
	"""Calculate the order parameters for the sn1 carbon tail, for all-atom.
	Takes the universe (MDAnalysis) and the lipid resname as input.
	Returns a list of the averaged order parameters for each carbon,
	averaged over lipids and frames.
	WARNING: It does not reproduce the correct order parameter for double bonds (unsaturation)"""
	nframes = u.trajectory.n_frames
	carbons = u.select_atoms("resname {0:s} and name C3* and not name C3 and not name C31".format(lipid))
	lenght = len(np.unique(carbons.names))
	
	C3_order = []
	normal = [0,0,1]
	for i in range(2, lenght+2): #Loop over the carbon tail (C3* tail)
	    #print (i)
	    cp = u.select_atoms("resname {1:s} and name C3{0:d}".format(i, lipid), updating=True)
	    hx = u.select_atoms("resname {1:s} and name H{0:d}X".format(i, lipid), updating=True)
	    hy = u.select_atoms("resname {1:s} and name H{0:d}Y".format(i, lipid), updating=True)
	    hz = u.select_atoms("resname {1:s} and name H{0:d}Z".format(i, lipid), updating=True)
	    summ = 0 #Order sum
	    nh = 0 #normalizing factor
	    nres = cp.names.shape[0] #number of carbons 
	    frames = len(u.trajectory[::1])
	    for ts in u.trajectory[::1]: #Loop over trajectory
	        cpx = cp.positions[:,0] #Positions for carbon
	        cpy = cp.positions[:,1]
	        cpz = cp.positions[:,2]
	        if hx.names.shape[0] != 0:
	            hxx = cpx - hx.positions[:,0] # the vector describing C-H
	            hxy = cpy - hx.positions[:,1]
	            hxz = cpz - hx.positions[:,2]
	            for idx, x in enumerate(hxx): #Calculate order for HX atoms
	                # P2 = 0.5∙(3∙(𝑐𝑜𝑠)^2 〈𝜃〉−1)
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



def P2_order_per_bond (u, resids, b=""):
    """Function for calculating the averaged order per bond for each resid.
    The function only works for Coarse Grained, since it is the second rank P2 order 
    parameter, which is calculated. 
    The function takes a list of resids to loop over.
    The list can contain either one or more resids. 
    you can spcify the bonds, otherwise they are by default
    specific for POPC tails only in Martini.
    The normal is assumed to be [0,0,1]"""
    #bonds = [ ['GL1', 'C1A'], ['C1A', 'D2A'], ['D2A', 'C3A'], ['C3A', 'C4A'], ['GL2','C1B'], ['C1B','C2B'], ['C2B','C3B'], ['C3B','C4B'] ]
    if len(b) == 0:
        bonds = [ ['C1A', 'D2A'], ['D2A', 'C3A'], ['C3A', 'C4A'], ['C1B','C2B'], ['C2B','C3B'], ['C3B','C4B'] ]
    else:
        bonds = b
    normal = [0,0,1]
    current_order_parameters= []
    
    if len(resids) == 1:
        for bond in bonds:
            order_parameter = 0
            a = u.select_atoms('resid {0:d} and name {1:s}'.format(resids[0],bond[0])).positions
            b = u.select_atoms('resid {0:d} and name {1:s}'.format(resids[0],bond[1])).positions
            v = a-b
            vector = v[0]
            norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
            projection = vector[0]*normal[0] + vector[1]*normal[1] + vector[2]*normal[2]
            order_parameter += projection**2/norm2
            current_order_parameters.append(0.5*(3.0*(order_parameter/len(resids)) - 1.0))
            
    elif len(resids) > 1:
        for bond in bonds:
            order_parameter = 0
            for r in resids:
                a = u.select_atoms('resid {0:d} and name {1:s}'.format(r,bond[0])).positions
                b = u.select_atoms('resid {0:d} and name {1:s}'.format(r,bond[1])).positions
                v = a-b
                vector = v[0]
                norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
                projection = vector[0]*normal[0] + vector[1]*normal[1] + vector[2]*normal[2]
                order_parameter += projection**2/norm2
            current_order_parameters.append(0.5*(3.0*(order_parameter/len(resids)) - 1.0))
    return current_order_parameters



def select_all_lipids (u):
    """Selects all the lipids in the system, 
    assuming there is only protein and lipids in the universe.
    The universe needs to be loaded with the tpr file.
    The function returns both the membrane selection
    and the unique resids and resnames of the lipids."""
    memb = u.select_atoms('all and not protein')
    memb_resids = np.unique(memb.resids)
    lipid_resnames = np.unique(memb.resnames)
    return memb, memb_resids, lipid_resnames


def extract_bonds (lipid_selection):
    """The function extracts and append
    the bonds for the selection. 
    The output is a nested list with the bond_names"""
    bond_names =  []
    for b in lipid_selection.bonds:
        bond_names.append(b.atoms.names)
    B = [a for i, a in enumerate(bond_names) if not any(all(c in h for c in a) for h in bond_names[:i])]
    return B


def get_normal (selection):
    """Calculate the normal of the selection,
    by taking the smallest principal axis
    from the moment of inertia."""
    p1, p2, p3 = selection.principal_axes()
    return p3

def divide_bonds_into_tails (bonds):
    """The function divides the bond names into tail A and tail B"""
    Tail_A = []
    Tail_B = []
    for b in bonds:
        if b[0].find('A')>0 or b[1].find('A')>0:
            Tail_A.append(b)
        elif b[0].find('B')>0 or b[1].find('B')>0:
            Tail_B.append(b)
    return Tail_A, Tail_B
