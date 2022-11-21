#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import string
import matplotlib.pyplot as plt


# In[2]:


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
plt.rcParams["figure.figsize"] = cm2inch(16,8)
plt.rcParams.update({'font.size':10})
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


# In[3]:


systems=["FA_no_G-prot","FA_G-prot","IA","PA"]
begin_time_ps=["15000000","15000000","15000000","15000000"]
end_time_ps=["-1","-1","-1","-1"]


# In[4]:


lipid_dir = {"FA_G-prot":"POPC, POPS, POPA, POPG, POPE, CHOL, DPGM3, DPSM, PAP6",
		"FA_no_G-prot":"POPC, POPS, POPA, POPG, POPE, CHOL, DPGM3, DPSM, PAP6",
		"IA":"POPC, POPS, POPA, POPG, POPE, CHOL, DPGM3, DPSM, PAP6",
		"PA":"POPC, POPS, POPA, POPG, POPE, CHOL, DPGM3, DPSM, PAP6"}


# In[ ]:


#ndx_cmd= "for i in `ls GRO/*[0-3].gro`; do base=`basename $i .gro`; ./Make_ndx_rdf.py -n ${base}; done"
#os.system(ndx_cmd)


# ### Calculating RDF for all the lipids in the different systems

# #### only using the entire trajectory

# In[ ]:


for idx, s in enumerate(systems):
    for r in range(3): 
        print ("Replica", r)
        #tpr = "GRO/{0:s}_{1:d}_nowat.gro".format(s,r)
        tpr = "ORIGO_TPR/{0:s}_{1:d}.tpr".format(s,r)
        ndx = "NDX/{0:s}_{1:d}_nowat_rdf.ndx".format(s,r)
        xtc = "XTC/{0:s}_{1:d}_nowat_res_center_fit.xtc".format(s,r)
        b = begin_time_ps[idx]
        e = end_time_ps[idx]
        lipids = lipid_dir[s]
        sel = [ "{0:s}".format(l) for l in lipids.split(',') ]
        output_base = "RDF/{0:s}_{1:d}_nowat".format(s,r)
        
        for sele in sel:
            #print ("lipid {0:s}".format(sele))
            output_xvg = output_base+"_{0:s}_rdf.xvg".format(sele.split('&')[0]).translate({ord(c): None for c in string.whitespace})
            if os.path.exists(output_xvg):
                print ("{0:s} exists".format(output_xvg))
                continue
            else:
                lipid_sel = sele.translate({ord(c): None for c in string.whitespace})
                print ('lipid selection', lipid_sel)
                rdf_cmd = "gmx_s rdf -f {0:s} -s {1:s} -n {2:s} -o {3:s} -dt 1000 -ref Protein -sel {4:s} -b {5:s}".format(xtc, tpr, ndx, output_xvg, lipid_sel,b)
                #rdf_cmd = "gmx rdf -f {0:s} -s {1:s} -n {2:s} -o {3:s} -dt 10 -ref Protein -sel {4:s} -b {5:s} -e {6:s}".format(xtc, tpr, ndx, output_xvg, lipid_sel, b, e)
                print (rdf_cmd)
                os.system(rdf_cmd)





