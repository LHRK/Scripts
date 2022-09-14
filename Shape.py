#!/usr/bin/env python
# coding: utf-8



import MDAnalysis as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from matplotlib.ticker import FormatStrFormatter


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)



plt.rcParams["figure.figsize"] = cm2inch(8.5,4)
plt.rcParams.update({'font.size':8})




def get_shape_quantity(u, stride=10):
    #prot = u.select_atoms('name CA', updating=True)
    prot = u.select_atoms('name BB', updating=True)
    aspher = []
    for ts in u.trajectory[::stride]:
        aspher.append(prot.asphericity())
    return aspher


# In[6]:


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


# In[7]:


def area_leaflet(u, stride=10):
    A = []
    for ts in u.trajectory[::stride]:
        points = u.select_atoms('resname DMPC and name P').positions
        hull = ConvexHull(points)
        A.append(hull.area)
    return A


# In[8]:


def plot_convexhull (hull, points):
    plt.plot(points[hull.vertices,0], sel[hull.vertices,1], 'r--', lw=2)
    plt.plot(points[hull.vertices[0],0], sel[hull.vertices[0],1], 'ro')
    return


# In[9]:


systems_M2 = ['1D1_69', '1E3D1_134', '2N2_302', 'NW9_58', 'NW11_73', 'NW13_151']
systems_M3 = ['1D1_69', '1E3D1_134', '2N2_302', 'NW9_57', 'NW11_72', 'NW13_150']


# In[10]:


def load_replicas (sys, M=2):
    u1 = md.Universe('ANALYSIS_Build_Martini{0:d}/GRO/{1:s}_v1_nowat.gro'.format(M,sys), 
                     'ANALYSIS_Build_Martini{0:d}/XTC/{1:s}_v1_nowat_fit.xtc'.format(M, sys))
    u2 = md.Universe('ANALYSIS_Build_Martini{0:d}/GRO/{1:s}_v2_nowat.gro'.format(M,sys), 
                     'ANALYSIS_Build_Martini{0:d}/XTC/{1:s}_v2_nowat_fit.xtc'.format(M, sys))
    u3 = md.Universe('ANALYSIS_Build_Martini{0:d}/GRO/{1:s}_v3_nowat.gro'.format(M,sys), 
                     'ANALYSIS_Build_Martini{0:d}/XTC/{1:s}_v3_nowat_fit.xtc'.format(M, sys))
    return u1, u2, u3


# In[19]:


def Shape_quantity_main (u1,u2,u3):
    a1 = get_shape_quantity(u1)
    a2 = get_shape_quantity(u2)
    a3 = get_shape_quantity(u3)
    return a1, a2, a3





