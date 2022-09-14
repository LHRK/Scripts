#!/usr/local/bin/vmd
# For visualizing the different TMs as colored in the RMSF plots in the notebook.
# run it by: vmd -e this_script.tcl gro_file
# The systems without glycans, the starting residue is numbered 26, so you just need to renumber thoses gro files beforehand.
# There is a script in the GRO directory to do so


set OUTDIR "/home/au447022/Documents/GU/PIP_AA/GRENDELS/IA_POPC_PIP2_16_18/ANALYSIS"


color Display Background white
axes location off
display rendermode GLSL
display ambientocclusion on
display shadows on
display cuedensity 0.09
display aodirect 0.5
display projection Orthographic
display depthcue off

#color change rgb 27 0.520000 0.270000 0.500000
#color change rgb 7  0.000000 0.320000 0.160000
material change opacity AOEdgy 0.500000

#TM0
# Have changed color from 21 to 11
mol addrep 0
mol modselect 1 0 "protein and resid 26 to 132"
mol modcolor 1 0 ColorID 10  
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 1 0 AOEdgy

#TM1
mol addrep 0
mol modselect 2 0 "protein and resid 133 to 166"
mol modcolor 2 0 ColorID 25
mol modstyle 2 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 2 0 AOChalky

#TM2
mol addrep 0
mol modselect 3 0 "protein and resid 171 to 199"
mol modcolor 3 0 ColorID 0
mol modstyle 3 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 3 0 AOChalky

#TM3
mol addrep 0
mol modselect 4 0 "protein and resid 219 to 255"
mol modcolor 4 0 ColorID 12
mol modstyle 4 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 4 0 AOChalky

#TM4
mol addrep 0
mol modselect 5 0 "protein and resid 261 to 290"
mol modcolor 5 0 ColorID 7
mol modstyle 5 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 5 0 AOChalky


#TM5
mol addrep 0
mol modselect 6 0 "protein and resid 301 to 334"
mol modcolor 6 0 ColorID 9
mol modstyle 6 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 6 0 AOChalky

#TM6
mol addrep 0
mol modselect 7 0 "protein and resid 341 to 371"
mol modcolor 7 0 ColorID 1
mol modstyle 7 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 7 0 AOChalky

#TM7
mol addrep 0
mol modselect 8 0 "protein and resid 376 to 404"
mol modcolor 8 0 ColorID 4
mol modstyle 8 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 8 0 AOChalky

#TM8
mol addrep 0
mol modselect 9 0 "protein and resid 405 to 419"
mol modcolor 9 0 ColorID 3
mol modstyle 9 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 9 0 AOChalky


#Protein
mol addrep 0
mol modselect 10 0 protein
mol modstyle 10 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor 10 0 ColorID 2
mol modstyle 10 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 10 0 AOChalky


#TM0
# Have changed color from 21 to 11
mol addrep 0
mol modselect 11 0 "protein and resid 26 to 132"
mol modcolor 11 0 ColorID 10  
mol modstyle 11 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 11 0 AOEdgy

#TM1
mol addrep 0
mol modselect 12 0 "protein and resid 133 to 166"
mol modcolor 12 0 ColorID 25
mol modstyle 12 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 12 0 AOChalky

#TM2
mol addrep 0
mol modselect 13 0 "protein and resid 171 to 199"
mol modcolor 13 0 ColorID 0
mol modstyle 13 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 13 0 AOChalky

#TM3
mol addrep 0
mol modselect 14 0 "protein and resid 219 to 255"
mol modcolor 14 0 ColorID 12
mol modstyle 14 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 14 0 AOChalky

#TM4
mol addrep 0
mol modselect 15 0 "protein and resid 261 to 290"
mol modcolor 15 0 ColorID 7
mol modstyle 15 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 15 0 AOChalky


#TM5
mol addrep 0
mol modselect 16 0 "protein and resid 301 to 334"
mol modcolor 16 0 ColorID 9
mol modstyle 16 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 16 0 AOChalky

#TM6
mol addrep 0
mol modselect 17 0 "protein and resid 341 to 371"
mol modcolor 17 0 ColorID 1
mol modstyle 17 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 17 0 AOChalky

#TM7
mol addrep 0
mol modselect 18 0 "protein and resid 376 to 404"
mol modcolor 18 0 ColorID 4
mol modstyle 18 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 18 0 AOChalky

#TM8
mol addrep 0
mol modselect 19 0 "protein and resid 405 to 419"
mol modcolor 19 0 ColorID 3
mol modstyle 19 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modmaterial 10 0 AOChalky

mol delrep 0 0

set cent "{{1 0 0 -57.3869} {0 1 0 -59.2546} {0 0 1 -70.759} {0 0 0 1}} {{-0.703063 0.70111 0.118862 0} {-0.0377366 -0.203703 0.978308 0} {0.710111 0.683323 0.169668 0} {0 0 0 1}} {{0.0205545 0 0 0} {0 0.0205545 0 0} {0 0 0.0205545 0} {0 0 0 1}} {{1 0 0 0.01} {0 1 0 -0.26} {0 0 1 0} {0 0 0 1}}"

molinfo top set {center_matrix rotate_matrix scale_matrix global_matrix} $cent

material change opacity AOEdgy 0.500000
material change diffuse AOEdgy 0.890000
material change specular AOEdgy 0.200000
material change mirror AOEdgy 0.000000
material change outline AOEdgy 0.820000
material change outlinewidth AOEdgy 0.930000


#TM1
color change rgb 25 [expr 202/256.] [expr 178/265.] [expr 214/256.]  
#TM2
color change rgb 0 [expr 31/256.] [expr 120/256.] [expr 180/256.]
#TM3
color change rgb 12 [expr 178/256.] [expr 223/256.] [expr 138/256.]
#TM4
color change rgb 7 [expr 51/256.] [expr 160/256.] [expr 44/256.]
#TM5
color change rgb 9 [expr 251/256.] [expr 154/256.] [expr 163/256.]
#TM6
color change rgb 1 [expr 227/256.] [expr 26/256.] [expr 28/256.]
#TM7
color change rgb 4 [expr 253/256.] [expr 191/256.] [expr 111/256.]
#H8
color change rgb 3 [expr 255/256.] [expr 127/256.] [expr 0/256.]


#render Tachyon Sideview_1_$system "/mnt/software/vmd/1.9.3/tachyon_LINUXAMD64" -aasamples 12 -rescale_lights 1 -add_skylight 1 %s -format TARGA -res 2000 2000 -o $OUTDIR/%s.tga
#
#rotate y by 180
#
#render Tachyon Sideview_2_$system "/mnt/software/vmd/1.9.3/tachyon_LINUXAMD64" -aasamples 12 -rescale_lights 1 -add_skylight 1 %s -format TARGA -res 2000 2000 -o $OUTDIR/%s.tga
#

### VOLMAP

#mol mol new $VOLMAP type {dx}
#mol modstyle 0 1 Isosurface 0.250000 0 0 0 1 1
#mol modmaterial 0 1 AOChalky

molinfo top set {center_matrix rotate_matrix scale_matrix global_matrix} $cent
molinfo 0 set {center_matrix rotate_matrix scale_matrix global_matrix} $cent

#render Tachyon Sideview_1_$system "/mnt/software/vmd/1.9.3/tachyon_LINUXAMD64" -aasamples 12 -rescale_lights 1 -add_skylight 1 %s -format TARGA -res 2000 2000 -o $OUTDIR/%s.tga
#
#rotate y by 180
#
#render Tachyon Sideview_2_$system "/mnt/software/vmd/1.9.3/tachyon_LINUXAMD64" -aasamples 12 -rescale_lights 1 -add_skylight 1 %s -format TARGA -res 2000 2000 -o $OUTDIR/%s.tga

#exit
