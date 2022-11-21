#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import os
import string


# # Notebook for running the densmap calculations of the different systems
## The density is calculated using only the headgroups of the different lipids

# In[4]:


systems = ['IA','PA', 'FA_G-prot','FA_no_G-prot']

unique_lipids = ['POPC','POPS','POPE','POPA','POPG','CHOL','DPGM3', 'DPSM', 'PAP6']

all_lipids_sele = " ".join(['{0:s}'.format(i) for i in unique_lipids ])

lipid_dir = {"FA_G-prot":"POPC,POPS,POPA,POPG,POPE,CHOL,DPGM3,DPSM,PAP6",
            "FA_no_G-prot":"POPC,POPS,POPA,POPG,POPE,CHOL,DPGM3,DPSM,PAP6",
            "IA":"POPC,POPS,POPA,POPG,POPE,CHOL,DPGM3,DPSM,PAP6",
            "PA":"POPC,POPS,POPA,POPG,POPE,CHOL,DPGM3,DPSM,PAP6"}

begin_time_ps=["19000000"]




## FITTET XTC

r = 0
for idx, s in enumerate(systems):
    tpr = "ORIGO_TPR/{0:s}_{1:d}.tpr".format(s,r)
    ndx = "NDX/{0:s}_{1:d}_nowat_rdf_chains.ndx".format(s,r)
    xtc = "XTC/{0:s}_all_nowat_res_center_fit_last5_DENSMAP.xtc".format(s)
    lipids = lipid_dir[s]
    sel = [ "{0:s}".format(l) for l in lipids.split(',') ]
    output_base = "DENSMAP/FIT/COUNT/{0:s}_all_nowat".format(s)
    for sele in sel:
        print ("lipid {0:s}".format(sele))
        sele=sele.strip(' ')
        #output_xpm = output_base+"_{0:s}_densmap.xpm".format(sele).translate({ord(c): None for c in string.whitespace})
        #output_dat = output_base+"_{0:s}_densmap.dat".format(sele).translate({ord(c): None for c in string.whitespace})
        #output_eps = output_base+"_{0:s}_densmap.eps".format(sele).translate({ord(c): None for c in string.whitespace})
        output_xpm = output_base+"_{0:s}_ONE_REF_densmap.xpm".format(sele)
        output_dat = output_base+"_{0:s}_ONE_REF_densmap.dat".format(sele)
        output_eps = output_base+"_{0:s}_ONE_REF_densmap.eps".format(sele)
        if os.path.exists(output_dat):
            print ("{0:s} exists".format(output_dat))
            continue
        else:
            grep_cmd = "grep '\[' "+ndx+"| gawk '{print NR-1, $2}'"+"| grep '{0:s}' |".format(sele)+"gawk '{print $1}' > sel"
            print (grep_cmd)
            os.system(grep_cmd)
            #densmap_cmd1 = "cat sel | gmx_s densmap -f {0:s} -s {1:s} -n {2:s} -o {3:s} -dt 10 -b 19000000 -unit count".format(xtc, tpr, ndx, output_xpm)
            #print (densmap_cmd1)
            #os.system(densmap_cmd1)
            densmap_cmd2 = "cat sel | gmx_s densmap -f {0:s} -s {1:s} -n {2:s} -od {3:s} -dt 10 -b 19000000 -unit count".format(xtc, tpr, ndx, output_dat)
            print (densmap_cmd2)
            os.system(densmap_cmd2)


## # Protein Density
#
#
### FITTET XTC
#
#for idx, s in enumerate(systems):
#    tpr = "ORIGO_TPR/{0:s}_{1:d}.tpr".format(s,r)
#    ndx = "NDX/{0:s}_{1:d}_nowat_rdf_chains.ndx".format(s,r)
#    xtc = "XTC/{0:s}_{1:d}_all_nowat_res_center_fit_last5.xtc".format(s,r)
#    output_base = "DENSMAP/FIT/Protein/{0:s}_all_nowat".format(s)
#    output_xpm = output_base+"_{0:s}_densmap.xpm".format(sele)
#    output_dat = output_base+"_{0:s}_densmap.dat".format(sele)
#    output_eps = output_base+"_{0:s}_densmap.eps".format(sele)
#    if os.path.exists(output_dat):
#        print ("{0:s} exists".format(output_dat))
#        continue
#    else:
#        grep_cmd = "echo 4 > sel"
#        print (grep_cmd)
#        os.system(grep_cmd)
#        densmap_cmd1 = "cat sel | gmx densmap -f {0:s} -s {1:s} -n {2:s} -o {3:s} -dt 10 -b 19000000 -unit count".format(xtc, tpr, ndx, output_xpm)
#        print (densmap_cmd1)
#        os.system(densmap_cmd1)
#        densmap_cmd2 = "cat sel | gmx densmap -f {0:s} -s {1:s} -n {2:s} -od {3:s} -dt 10 -b 19000000 -unit count".format(xtc, tpr, ndx, output_dat)
#        print (densmap_cmd2)
#        os.system(densmap_cmd2)






