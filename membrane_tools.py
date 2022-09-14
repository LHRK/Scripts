#!/bin/python
import numpy as np

#Get list of lipids found within your system
def make_ndx(arr, title="selection", filename="dummy"):
    """The following function takes a array or list and writes it to an ndx file"""
    with open("{}.ndx".format(filename), "w+") as f:
        f.write("[ {} ]\n".format(title))
        lst=[]
        print("WORKING")
        for num, i in enumerate(arr):
            print(i)
            f.write("{0:8d}".format(np.array(i).astype(int)))
            lst.append(i)
            if len(lst)%15==0:
                f.write("\n")
                lst=[]
        f.write("\n")
        
def read_ndx(filename, DIR):
    """reads an index file and returns a numpy array"""
    data_all=[]
    v="%s/%s.ndx" %(filename,DIR)
    with open(v) as f:
        for num, line in enumerate(f.readlines()[1:]):
            data=np.array(line.strip().split(' '))
            data_all=np.append(data_all,data)
    return data_all.astype(int)
        
def lipid_list_unique(u):
    """Reads from universe and detects all unique lipid types. A list with the different lipid names present is returned"""
    resnames = [a.resname for a in u.select_atoms("not (protein or resname ION W WF CL- NA+)")]
    u_res = np.unique(resnames)
    return u_res

#Makes index file when MDAnalysis cannot be used
def append_make_ndx(arr, title="selection", filename="dummy", DIR="DIR"):
    """A input array or list is appended to an excisting index file"""
    with open("{0}/{1}.ndx".format(DIR,filename), "a") as f:
        f.write("[ {} ]\n".format(title))
        lst=[]
        print("WORKING")
        for num, i in enumerate(arr):
            f.write("{0:6d}".format(i))
            lst.append(i)
            if len(lst)%15==0:
                f.write("\n")
                lst=[]
        f.write("\n")

#MAKES UPPER AND LOWER LEAFLETS .NDX FOR LIPID TYPES FOUND IN THE SELECTION
def ndx_upper_lower_leaflet(select, u, select_leaflet="name PO4 GM1", name="index", OUT_DIR="/home"):
        """Index files are made for all lipid types including only a single bead. When possible to destinguish
        both upper and lower leaflet lipids are selected"""
        filename="{0}/{1}.ndx".format(OUT_DIR, name)
        print(filename)
        up_low={0:"upper", 1:"lower"}
        try:
            both_leaflets=u.select_atoms(select).write(filename="{0}/{1}.ndx".format(OUT_DIR, name), file_format="ndx")
        except:
            print("WARNING: NO ATOMS WERE FOUND IN SELECTION")
            return
        #try:
        #    cutoff, N = optimize_cutoff(u,select_leaflet, pbc=True)
        #    print(cutoff)
        #except: 
        #    print("WARNING: NO CUTOFF WAS FOUND")
        #    return
        leaflets=LeafletFinder(u, select_leaflet, cutoff=14.5, pbc=True, sparse=None)
        if len(leaflets.components)>2:
            print("WARNING: THE LEAFLET FINDER FINDS MORE THAN 2 LEAFLETS")
        for groupindex in range(len(leaflets.components)):
            resnames = [a.resname for a in leaflets.groups(groupindex).select_atoms(select)]
            u_res = np.unique(resnames)
            selection=leaflets.groups(groupindex).select_atoms(select)
            if not selection:
                continue
            indices=[a.index+1 for a in selection]
            filename="{0}/{1}_{2}".format(OUT_DIR,name, up_low[groupindex])
            print(filename)
            make_ndx(indices, title="%s_%i" %(name,groupindex), filename="{0}/{1}_{2}".format(OUT_DIR,name, up_low[groupindex]))
            print(u_res, groupindex)

#Given a priority atom name list and lipid resname a selection is returned of the specific lipid and beads of highest
#priority

def select_prio(resname,u,priority=["PO4", "GL1", "AM1"]):
    """Givwn a resname a priority selection is made only containing a single bead for the same lipid type"""
    selection=u.select_atoms("resname %s and name PO4 GM1 AM1 ROH GL1" %resname)
    if len(selection)!=selection.n_residues:
        beads=np.unique(selection.names)
        print("%s contains multiple of the listed bead types:" %resname)
        print("The bead with the highest priority is selected")
        bead_type=priority[np.amin([priority.index(v) for v in beads])]
        print("Which is %s" %bead_type)
    else:
        beads=np.unique(selection.names)
        bead_type=beads[0]
    new_selection=u.select_atoms("resname %s and name %s" %(resname, bead_type))
    return new_selection, bead_type
              
def calculate_density(universe, selection,delta=1.0):
    """The density is calculated for a selection. The MDAnalysis density function returns a 3D tensor,
    which is subsequently flatterned. Only the 2D grid is returned"""
    minv,maxv=min_max_coord(universe)
    gridcenter=(maxv-minv)/2
    site_density = density_from_Universe(universe, delta=delta,
                                     atomselection=selection, update_selection=True, gridcenter=gridcenter,xdim=maxv[0]-minv[0], ydim=maxv[1]-minv[1], zdim=maxv[2]-minv[2])
    arr=site_density.grid
    flatten=np.zeros([np.shape(arr)[0],np.shape(arr)[1]])
    #flatten_all=np.zeros([int(np.shape(arr)[0]/2),int(np.shape(arr)[1]/2)])
    for k, i in enumerate(range(np.shape(arr)[0])):
        for v,j in enumerate(range(np.shape(arr)[1])):
                   flatten[k,v]=np.sum(arr[i,j,:])
    return flatten
    #for i in range(np.shape(flatten)[0]):
    #    for j in range(np.shape(flatten)[1]):
    #               flatten_all[i%np.shape(flatten_all)[0],j%np.shape(flatten_all)[1]]+=flatten[i,j]
    #return flatten, flatten_all
def average_cent_geo(universe, selection):
    """The aberage center og geometry is calculated from a given selection"""
    v=[]
    for ts in universe.trajectory:
        prot=universe.select_atoms(selection)
        v.append(prot.center_of_geometry())
    center=np.average(np.array(v), axis=0)
    return center

def min_max_coord(universe):
    """The minimum an maximum coordinates are returned from the first frame or single structure present in
    the universe"""
    #maybe calculate avereage if single frame does not work
    sel1=universe.select_atoms("all")
    min=np.amin(sel1.atoms.positions, axis=0)
    max=np.amax(sel1.atoms.positions, axis=0)
    return min, max
        
def sub_matrix(universe, density, atoms_in_prot):
    """This function is only relavant if one has 4 proteins fixed in the same membrane"""
    mini, maxi = min_max_coord(universe)
    cent_ndx_x=[]
    cent_ndx_y=[]
    min_x=[]
    min_y=[]
    start_ndx=0
    end_ndx=atoms_in_prot
    cent_geo=[]
    num_prot=4
    #COG for all proteins are collected
    for i in range(num_prot):
        val=average_cent_geo(universe,"bynum %i:%i" %(start_ndx+1, end_ndx))
        cent_geo.append(val)
        start_ndx+=atoms_in_prot
        end_ndx=atoms_in_prot+start_ndx
    #The COG coordinates are converted to indices in the density array
    for i in range(num_prot):
        cog_x=int((cent_geo[i][0]-mini[0])/(maxi[0]-mini[0])*np.shape(density)[0])
        cog_y=int((cent_geo[i][1]-mini[1])/(maxi[1]-mini[1])*np.shape(density)[1])
        cent_ndx_x.append(cog_x)
        cent_ndx_y.append(cog_y)
        min_x.append(cog_x-0)
        min_x.append(np.shape(density)[0]-cog_x)
        min_y.append(cog_y-0)
        min_y.append(np.shape(density)[0]-cog_y)
    min_xslice=np.amin(np.array(min_x))
    min_yslice=np.amin(np.array(min_y))
    #The density array is sliced relative to the COG indices
    dens=np.zeros([min_xslice*2,min_yslice*2])
    for i in range(num_prot):
        add_dens=density[cent_ndx_x[i]-min_xslice:cent_ndx_x[i]+min_xslice,cent_ndx_y[i]-min_yslice:cent_ndx_y[i]+min_yslice]
        shape_x=np.shape(add_dens)[0]
        shape_y=np.shape(add_dens)[1]
        dens[:shape_x,:shape_y]+=add_dens
    return dens

def superimpose(name1, outname):
    """A .tga picture of ones protein is superimposed onto the density map"""
    import svgwrite
    O = svgwrite.drawing.Drawing(filename="%s" %outname,size=("200mm","200mm"))
    bottom=svgwrite.image.Image(href="%s" %name1, insert=("0mm","0mm"), size=("200mm","200mm"))
    top=svgwrite.image.Image(href="/home/au341226/LLNL/4_proteins_complex_brain/hDAT/MD/hDAT.pov.tga",insert=("48mm","48mm"), size=("120mm","120mm"))
    O.add(bottom)
    O.add(top)
    O.save()
    
