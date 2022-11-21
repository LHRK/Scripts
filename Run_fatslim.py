#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os



systems=["FA_G-prot"]
begin_time_ps=["15000000"]
end_time_ps=["18000000"]


# THICKNESS

for idx, s in enumerate(systems):
    for r in range(3):
        b = begin_time_ps[idx]
        e = end_time_ps[idx]
        tpr = "TPR/{0:s}_{1:d}_nowat.tpr".format(s,r)
        gro = "GRO/{0:s}_{1:d}_nowat.gro".format(s,r)
        ndx = "NDX/{0:s}_{1:d}_nowat_selection.ndx".format(s,r)
        xtc = "XTC/{0:s}_{1:d}_nowat_res_center_fit.xtc".format(s,r)
        output = "THICKNESS/{0:s}_{1:d}_nowat_thickness.xvg".format(s,r)
        output_raw = "THICKNESS/RAW/{0:s}_{1:d}_nowat_thickness.csv".format(s,r)
        if os.path.exists(output):
                print ("{0:s} exists".format(output))
                continue
        else:
                cmd = "fatslim thickness -c {0:s} -t {1:s} -n {2:s} --plot-thickness {3:s} --export-thickness-raw {4:s} --hg-group Headgroup_POPC --interacting-group GCGR --thickness-cutoff 10.0 -b 1000000".format(gro, xtc, ndx, output, output_raw)
                os.system(cmd)



# APL

for idx, s in enumerate(systems):
    for r in range(3):
        b = begin_time_ps[idx]
        e = end_time_ps[idx]
        tpr = "TPR/{0:s}_{1:d}_nowat.tpr".format(s,r)
        gro = "GRO/{0:s}_{1:d}_nowat.gro".format(s,r)
        ndx = "NDX/{0:s}_{1:d}_nowat_selection.ndx".format(s,r)
        xtc = "XTC/{0:s}_{1:d}_nowat_res_center_fit.xtc".format(s,r)
        output = "APL/{0:s}_{1:d}_nowat_apl.xvg".format(s,r)
        output_raw = "APL/RAW/{0:s}_{1:d}_nowat_apl.csv".format(s,r)
        if os.path.exists(output):
                print ("{0:s} exists".format(output))
                continue
        else:
                cmd = "fatslim apl -c {0:s} -t {1:s} -n {2:s} --plot-apl {3:s} --export-apl-raw {4:s} --hg-group Headgroups --interacting-group GCGR -b 1000000 --apl-by-type".format(gro, xtc, ndx, output, output_raw)
                os.system(cmd)

