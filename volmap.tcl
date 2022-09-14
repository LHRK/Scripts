!/usr/bin/tclsh


set SYS [lindex $argv 0]
set lipid [lindex $argv 1]

mol new GRO/${SYS}.gro
mol addfile XTC/${SYS}_res_center_densmap_fit.xtc first 0 last -1 step 1 waitfor -1 0 

volmap occupancy [atomselect 0 "resname ${lipid}"] -res 1.0 -allframes -combine avg -o VOLMAP/${SYS}_${lipid}_last_0ns_volmap.dx

exit

