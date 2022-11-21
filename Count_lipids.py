#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import string
import MDAnalysis as md
import numpy as np
import matplotlib.pyplot as plt
from numba import jit


# In[2]:


systems=["5XEZ_LB1", "5XEZ_LB2", "5XEZ_LB3", "5XEZ_LB4", "5XEZ_LB_semi-complex", "5XEZ_LB_semi-complex2", "5YQZ_LB_semi-complex", "5YQZ_LB_semi-complex2"]


# In[3]:


lipid_dir = {"5XEZ_LB1":"POPC, POPS",
            "5XEZ_LB2":"POPC, POPE",
            "5XEZ_LB3":"POPC, POPA",
            "5XEZ_LB4":"POPC, POPG",
            "5XEZ_LB_semi-complex":"POPC, POPS, POPE, POPA, POPG",
            "5XEZ_LB_semi-complex2":"POPC, POPS, POPE, POPA, POPG, POP2, CHOL",
            "5YQZ_LB_semi-complex":"POPC, POPS, POPE, POPA, POPG",
            "5YQZ_LB_semi-complex2":"POPC, POPS, POPE, POPA, POPG, POP2, CHOL"}


# # Frist shell

# In[4]:


@jit
def count_lipids (u, lipid, stride=10):
    '''Counts the number of lipids within a shell of 0.7 nm of the TMD'''
    tot_frames = len(u.trajectory[::stride])
    O = np.zeros([tot_frames])
    
    for idx, ts in enumerate(u.trajectory[::stride]):
        sel = u.select_atoms('resname {0:s} and global around 7 resid 111-394'.format(lipid))
        num_contacts = np.unique(sel.resids).shape[0]
        O[idx] = num_contacts
    return O


# In[7]:


for idx, s in enumerate(systems):
    lipids = lipid_dir[s]
    for r in range(3):
        gro = 'GRO/{0:s}_{1:d}_nowat.gro'.format(s,r)
        xtc = 'XTC/{0:s}_{1:d}_nowat_mol_center_fit.xtc'.format(s,r)
        u = md.Universe(gro,xtc)
        for lipid in lipids.split(','):
            print ('lipid:', lipid)
            llipid = lipid.translate({ord(c): None for c in string.whitespace})
            O = count_lipids(u, lipid)
            np.savetxt('COUNTS/{0:s}_{1:d}_number_of_{2:s}_within_0.7-nm.txt'.format(s,r,llipid), O)


# In[8]:


@jit
def count_lipids_1 (u, lipid, stride=10):
    '''Counts the number of lipids within a shell of 1 nm of the TMD'''
    tot_frames = len(u.trajectory[::stride])
    O = np.zeros([tot_frames])
    
    for idx, ts in enumerate(u.trajectory[::stride]):
        sel = u.select_atoms('resname {0:s} and global around 10 resid 111-394'.format(lipid))
        num_contacts = np.unique(sel.resids).shape[0]
        O[idx] = num_contacts
    return O


# In[10]:


for idx, s in enumerate(systems):
    lipids = lipid_dir[s]
    for r in range(3):
        gro = 'GRO/{0:s}_{1:d}_nowat.gro'.format(s,r)
        xtc = 'XTC/{0:s}_{1:d}_nowat_mol_center_fit.xtc'.format(s,r)
        u = md.Universe(gro,xtc)
        for lipid in lipids.split(','):
            print ('lipid:', lipid)
            llipid = lipid.translate({ord(c): None for c in string.whitespace})
            O = count_lipids_1(u, lipid)
            np.savetxt('COUNTS/{0:s}_{1:d}_number_of_{2:s}_within_1-nm.txt'.format(s,r,llipid), O)


# # Second shell

# In[11]:


@jit
def count_lipids_1 (u, lipid, stride=10):
    '''Counts the number of lipids within a shell of 1.2 nm of the TMD'''
    tot_frames = len(u.trajectory[::stride])
    O = np.zeros([tot_frames])
    
    for idx, ts in enumerate(u.trajectory[::stride]):
        sel = u.select_atoms('resname {0:s} and global around 12 resid 111-394'.format(lipid))
        num_contacts = np.unique(sel.resids).shape[0]
        O[idx] = num_contacts
    return O


# In[12]:


for idx, s in enumerate(systems):
    lipids = lipid_dir[s]
    for r in range(3):
        gro = 'GRO/{0:s}_{1:d}_nowat.gro'.format(s,r)
        xtc = 'XTC/{0:s}_{1:d}_nowat_mol_center_fit.xtc'.format(s,r)
        u = md.Universe(gro,xtc)
        for lipid in lipids.split(','):
            print ('lipid:', lipid)
            llipid = lipid.translate({ord(c): None for c in string.whitespace})
            O = count_lipids_1(u, lipid)
            np.savetxt('COUNTS/{0:s}_{1:d}_number_of_{2:s}_within_1.2-nm.txt'.format(s,r,llipid), O)


# In[ ]:




