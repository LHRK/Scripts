#!/usr/bin/env python

'''
import functions you need
''' 
import numpy
import scipy
import sys
import string
import MDAnalysis
import os,shutil,math,timeit,multiprocessing
import cPickle as pickle
from multiprocessing import Value,Lock
import MDAnalysis.analysis.contacts as contacts
import MDAnalysis.analysis.distances as distances
import scipy.spatial.distance as scipydistance

'''
Option class
 
borrowed from martinize.py - this gives each Option object a type that the function is, e.g. 
boolean or string), a number (the number of arguments expected for the option), the value 
of the option (can be a default, or specified by the user), and the description of the 
option (what the option is for, how to use it, when is it best used, etc.)
'''
class Option:
	def __init__(self,func=str,num=1,default=None,description=""):
		self.func        = func
		self.num         = num
		self.value       = default
		self.description = description
	def __nonzero__(self):
		if self.func == bool:
		    return self.value != False
		return bool(self.value)
	def __str__(self):
		return self.value and str(self.value) or ""
	def setvalue(self,v):
		if len(v) == 1:
		    self.value = self.func(v[0])
		else:
		    self.value = [ self.func(i) for i in v ]

'''
options

Structure borrowed from martinize.py - list of all the options names, as well as their respective 
option objects (with type, number of arguments, default values, and explanations).
'''
options=[
	("-g",          Option(str,  1,None,                        "Input file (gro)")),
	("-x",          Option(str,  1,None,                        "Input file (xtc)")),
	("-s",          Option(int,  1,1,                           "Number of frames to skip (default=1)")),
	("-b",          Option(int,  1,0,                           "Starting frame (default=0")),
	("-e",          Option(int,  1,-1,                          "Last frame for calculations (default=-1)")),
	("-z",          Option(float,1,140,                         "Z cutoff for determining which lipids are on each leaflet (default=140; in A)")),
	("-numprot",    Option(int,  1,1,                           "Number of proteins in your system (default=1)")),
	("-protsel",    Option(str,  1,"name@BB",                   "Residue selection of protein in your system, separated by @'s, an \n\t\t\t\t\t example is resnum@72:643 for selecting residues 72-643 (default=name@BB)")),
	("-cutoff",     Option(float,1,20,                          "Cutoff for determining the depletion/enrichment index (default=20; in A)")),
	("-bilCentInd", Option(int,  1,52637,                       "Atom number of the bilayer middle, i.e. all lipids < this atom are \n\t\t\t\t\t on the top leaflet (default=52637)")),
	("-lastLipInd", Option(int,  1,99093,                       "Atom number of the last lipid particle in the bilayer, i.e. all lipids < this atom\n\t\t\t\t\t and > -bilCentInd are on the bottom leaflet (default=99093)")),
	("-outfile",    Option(str,  1,"DepletementEnrichment.xvg", "Output file (default: DepletionEnrichment.xvg)")),
	("-h",          Option(bool, 0,False,                       "Display this help screen and exit"))
]

'''
help

prints the help screen for this script and exits.  Will list the option tag, the type it is/takes,
the number of arguments it takes, and the description of what it is/what it does.  
'''
def help():
	print
	print "\t\t\t      ==> lipid-PM-counter-fix-v3.0.py <=="
	print
	print "\t\t\t    Written by Helgi Ingolfsson and Ruo-Xu Gu"
	print
	print "\t\t\t           Modified by Amanda Buyan"
	print
	print "\tDESCRIPTION"
	print "\t--------------------------------------------------------------------"
	print
	print "\t lipid-PM-counter-fix-v3.0.py calculates the depletion-enrichment "
	print "\t index for each lipid in your system, for as many proteins that you "
	print "\t have in your system.  Currently supports the following lipids:"
	print
	print "\t 	    1.  Phosphatidylcholines      (POPC,DOPC,PEPC,PAPC,DAPC,PUPC) "
	print "\t 	    2.  Phosphatidylethanolamines (POPE,DOPE,PIPE,PQPE,PAPE,DAPE,PUPE,DUPE) "
	print "\t 	    3.  Sphingomyelins            (DPSM,DBSM,DXSM,POSM,PGSM,PNSM,BNSM,XNSM) "
	print "\t 	    4.  Ceramides                 (DPCE,DXCE,PNCE,XNCE) "
	print "\t 	    5.  Glycolipids               (DPG1,DXG1,PNG1,XNG1,DPG3,DXG3,PNG3,XNG3) "
	print "\t 	    6.  Phosphatidic Acids        (POPA,PIPA,PAPA,PUPA) "
	print "\t 	    7.  Lysophosphatidylcholines  (PPC,OPC,IPC,APC,UPC) "
	print "\t 	    8.  Diacylglycerols           (PODG,PIDG,PADG,PUDG) "
	print "\t 	    9.  Phosphatidylserines       (POPS,PIPS,PQPS,PAPS,DAPS,PUPS,DUPS) "
	print "\t 	    10. Phosphatidylinositols     (POPI,PIPI,PAPI,PUPI) "
	print "\t 	    11. PIPs (1-3 Phosphates)     (POP1,POP2,POP3) "
	print "\t 	    12. Cholesterol "
	print
	print "\t To calculate the lipid depletion/enrichment index, the following "
	print "\t formula is used: "
	print
	print "\t            DEI(L)x = Ratio(L)x/Ratio(L)bulk "
	print
	print "\t where "
	print
	print "\t            Ratio(L)x =     (number of lipids)x "
	print "\t                        -------------------------- "
	print "\t                        (total number of lipids)x "
	print
	print "\t         Ratio(L)bulk =       total number L "
	print "\t                        ------------------------ "
	print "\t                         total number of lipids "
	print
	print "\t L represents one type of lipid, and x represents the protein in "
	print "\t your system.  To calculate this, the following variables are  "
	print "\t available:"
	print
	print "\t     Option   Type  Num Args   Description"
	print "\t--------------------------------------------------------------------"
	print
	for item in options:
		if type(item) == str:
			print item
	for item in options:
		if type(item) != str:
			print "\t%11s  %5s     %s       %s"%(item[0],item[1].func.__name__,item[1].num,item[1].description)
	print
	sys.exit()

'''
option_parser

creates a library of options, which the function then returns for use by the program.
args: arguments from the options
options: options library that is returned
'''
def option_parser(args,options):

	# Check whether there is a request for help
	if '-h' in args or '--help' in args:
		help()

	# Convert the option list to a dictionary, discarding all comments
	options = dict([i for i in options if not type(i) == str])

	# get the arguments in the args
	while args:
		ar = args.pop(0)
		if ar=="-eigen" or ar=="-pick" or ar=="-mindist":
			options[ar].setvalue([True])
		else:
			options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])

	#return options
	return options

'''
update_progress

updates progress bar for calculations

progress: decimal between 0 and 1 describing what percentages the calculations are at
'''
def update_progress(progress):
        barLength=50
        status=""
        if isinstance(progress,int):
                progress=float(progress)
        if not isinstance(progress,float):
                progress=0
                status="error: progress var must be float\r\n"
        if progress<0:
                progress=0
                status="Halt...\r\n"
        if progress >=1:
                progress=1
                status="Done...\r\n"
        block=int(round(barLength*progress))
        text="\r\t\tPercent: [{0}] {1:.1f}% {2}".format("#"*(block)+" "*((barLength-block)),progress*100,status) #{1}%
        sys.stdout.write(text)
        sys.stdout.flush()

'''
getlindex

Get lipid list index
'''
def getlindex(atomgroup,lipidDictionary):

	atomindices=atomgroup.atoms.indices

	cIndex = numpy.full((len(atomindices),1),-1).astype(int)
	clayer = numpy.full((len(atomindices),1),-1).astype(int)

	for i,entry in enumerate(atomindices):
		if entry < options["-bilCentInd"].value:
			clayer[i] = 0
		elif entry < options["-lastLipInd"].value:
			clayer[i] = 1
		else:
			clayer[i] = -1
		lipid=str(atomgroup[i].resname) + "_" + str(clayer[i][0])
		try:
			cIndex[i] = lipidDictionary[lipid]
		except KeyError:
			cIndex[i] -1
	return cIndex

'''
getSelection

Return lipid selection string
'''
def getSelection(lipid):
	cbead = ""
	if lipid[3] == 0:
		cbead = "GL1"
	elif lipid[3] == 1:
		cbead = "AM1"
	elif lipid[3] == 2:
		cbead = "ROH"
	else:
		cbead = ""
	cIndex = ""
	if lipid[1] == 0:
		cIndex = "1:"+str(options["-bilCentInd"].value)
	elif lipid[1] == 1:
		cIndex = str(options["-bilCentInd"].value)+":"+str(options["-lastLipInd"].value)
	elif lipid[1] == 2:
		cIndex = "1:"+str(options["-lastLipInd"].value)
	else:
		cIndex = "" 
	return "resname " + lipid[0] + " and name " + cbead + " and bynum " + cIndex

'''
get_counts

determines the number and identity of lipids within the pretermined cutoff of the protein.  Function
returns the process number, the total number of lipids around the protein, and the total number of each
different type of lipid around the protein.
'''
def get_counts(indexlist,lipidDictionary,indexlistcount,protein_selection_str,lipidList):

	# create a universe object for this function, and selections for protein and relevant headgroups 
	syst=MDAnalysis.Universe(options["-g"].value,options["-x"].value)
	headgroups = syst.select_atoms("name ROH GL1 AM1")
	protein = syst.select_atoms(protein_selection_str)

        # figure out number of frames
        num_frames_local=int(indexlist[1])-int(indexlist[0])+1

	# keep track of frames processed locally 
        num_frames_processed_locally=0

	# allocate these arrays
	plTotal =  numpy.zeros((num_frames_local,options["-numprot"].value), dtype=int)
	plCounts = numpy.zeros((num_frames_local,len(lipidList),options["-numprot"].value), dtype=int)

	# loop over the number of frames for this part of the trajectory 
	for i,ts in enumerate(syst.trajectory[int(indexlist[0]):int(indexlist[1])+1]):

		# loop over all frames and print things (change this)
		for j in range(options["-numprot"].value):
			coords_protein=protein.atoms.positions
			coords_headgroups=headgroups.atoms.positions
			pairs=MDAnalysis.lib.distances.capped_distance(coords_protein,coords_headgroups,options["-cutoff"].value,return_distances=False)
			indices=list((set(pairs[:,1])))
			tempSel=headgroups[indices]
			plTotal[i,j]+=len(tempSel.atoms)
			curlIndices = getlindex(tempSel,lipidDictionary)
			if curlIndices.any < 0:
				for k,index in enumerate(curlIndices):
					if index < 0:
						print("WARNING: found %s that is not in lipid list, go add it" % tempSel[k].resname)
						sys.exit()
			else:
				for curlIndex in curlIndices:
					plCounts[i,curlIndex,j] += 1

		# update the number of frames processed
		num_frames_processed_locally+=1

		# update progress       
		if num_frames_processed.value==num_frames.value-1:
			update_progress(1)
		else:
			if num_frames_processed_locally==10:
				update_num_frames_processed(num_frames_processed_locally)
				update_progress(float(num_frames_processed.value)/num_frames.value)
				num_frames_processed_locally=0
			elif ts.frame==int(indexlist[1]):
				update_num_frames_processed(num_frames_processed_locally)
				update_progress(float(num_frames_processed.value)/num_frames.value)
				num_frames_processed_locally=0
			else:
				pass

	# return desired values
	return indexlistcount,plTotal,plCounts	

'''
update_num_frames

update the global total number of frames for the workers

x: total number of frames (from main function)
'''
def update_num_frames(x):
        global num_frames
        num_frames.value=x

'''
update_num_frames_processed

update the global number of frames that have been calculated through by the workers.
'''
def update_num_frames_processed(num):
        global num_frames_processed
        global lock
        lock.acquire()
        num_frames_processed.value+=num
        lock.release()

'''
declare global variables here
'''

# declare counting variables and lock for keeping track of frames
num_frames_processed=Value('i',1)
num_frames=Value('i',1)
lock=Lock()

# get arguments
args=sys.argv[1:]
options=option_parser(args,options)

'''
main function

Calculates the depletement/enrichment index of all lipids in your system (this is also parallelised)
'''
def main():

	# variable lipidList
	lipidList = []

	# initialise indices for grouping
	oldlipidindex=0
	lipidindex=0

	# saturation index groups
	# fully saturated: PPC, DPSM, DBSM, DXSM, DPCE, DXCE, DPG1, DXG1, DPG3, DXG3
	# poly-unsaturated: DAPC, DAPE, DUPE, DAPS, DUPS, APC, UPC
	fully_saturated=[]
	poly_unsaturated=[]
	neither=[]

	### make the list of all available lipids
	# 0 ==> Name of lipid
	# 1 ==> Leaflet (Upper=0,Lower=1,Both=2)
	# 2 ==> Lipid Count 
	# 3 ==> Bead to use as representative (0=GL1, 1=AM1, 2=ROH)

	# PCs (Phosphatidylcholines)
	upperPCindices=[]
	lowerPCindices=[]
	lipidList.append(["POPX", 0, 0, 0])
	lipidList.append(["POPC", 0, 0, 0])
	lipidList.append(["POPC", 1, 0, 0])
	lipidList.append(["DOPC", 0, 0, 0])
	lipidList.append(["DOPC", 1, 0, 0])
	lipidList.append(["PIPX", 0, 0, 0])
	lipidList.append(["PIPC", 0, 0, 0])
	lipidList.append(["PIPC", 1, 0, 0])
	lipidList.append(["PEPC", 0, 0, 0])
	lipidList.append(["PEPC", 1, 0, 0])
	lipidList.append(["PAPC", 0, 0, 0])
	lipidList.append(["PAPC", 1, 0, 0])
	lipidList.append(["DAPC", 0, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	lipidList.append(["DAPC", 1, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	lipidList.append(["PUPC", 0, 0, 0])
	lipidList.append(["PUPC", 1, 0, 0])
	lipidindex+=len(lipidList)
	for i in range(oldlipidindex,lipidindex):
		if lipidList[i][1]==0:
			upperPCindices.append(i)
		else:
			lowerPCindices.append(i)

	# PEs (Phosphatidylethanolamines)
	upperPEindices=[]
	lowerPEindices=[]
	lipidList.append(["POPE", 0, 0, 0])
	lipidList.append(["POPE", 1, 0, 0])
	lipidList.append(["DOPE", 0, 0, 0])
	lipidList.append(["DOPE", 1, 0, 0])
	lipidList.append(["PIPE", 0, 0, 0])
	lipidList.append(["PIPE", 1, 0, 0])
	lipidList.append(["PQPE", 0, 0, 0])
	lipidList.append(["PQPE", 1, 0, 0])
	lipidList.append(["PAPE", 0, 0, 0])
	lipidList.append(["PAPE", 1, 0, 0])
	lipidList.append(["DAPE", 0, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	lipidList.append(["DAPE", 1, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	lipidList.append(["PUPE", 0, 0, 0])
	lipidList.append(["PUPE", 1, 0, 0])
	lipidList.append(["DUPE", 0, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	lipidList.append(["DUPE", 1, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		if lipidList[i][1]==0:
			upperPEindices.append(i)
		else:
			lowerPEindices.append(i)

	# PSs (Phosphatidylserines)
	PSindices=[]
	lipidList.append(["POPS", 1, 0, 0]) 
	lipidList.append(["PIPS", 1, 0, 0])
	lipidList.append(["PQPS", 1, 0, 0])
	lipidList.append(["PAPS", 1, 0, 0])
	lipidList.append(["DAPS", 1, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	lipidList.append(["PUPS", 1, 0, 0])
	lipidList.append(["DUPS", 1, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		PSindices.append(i)

	# PAs (Phosphatidic Acids)
	PAindices=[]
	lipidList.append(["POPA", 1, 0, 0]) 
	lipidList.append(["PIPA", 1, 0, 0])
	lipidList.append(["PAPA", 1, 0, 0])
	lipidList.append(["PUPA", 1, 0, 0])
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		PAindices.append(i)

	# DAGs (Diacylglycerols)
	# Calculate inner/outer together, as these flip flop
	DAGindices=[]
	lipidList.append(["PODG", 2, 0, 0])
	lipidList.append(["PIDG", 2, 0, 0])
	lipidList.append(["PADG", 2, 0, 0])
	lipidList.append(["PUDG", 2, 0, 0])
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		DAGindices.append(i)
 
	# LPCs - Lysophosphatidylcholines
	LPCindices=[]
	lipidList.append(["PPC", 0, 0, 0])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["OPC", 0, 0, 0])
	lipidList.append(["IPC", 0, 0, 0])
	lipidList.append(["APC", 0, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	lipidList.append(["UPC", 0, 0, 0])
	poly_unsaturated.append(len(lipidList)-1)
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		LPCindices.append(i)

	# SMs - Sphingomyelins
	upperSMindices=[]
	lowerSMindices=[]
	lipidList.append(["DPSM", 0, 0, 1]) 
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DPSM", 1, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DBSM", 0, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DBSM", 1, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DXSM", 0, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DXSM", 1, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["POSM", 0, 0, 1])
	lipidList.append(["POSM", 1, 0, 1])
	lipidList.append(["PGSM", 0, 0, 1])
	lipidList.append(["PGSM", 1, 0, 1])
	lipidList.append(["PNSM", 0, 0, 1])
	lipidList.append(["PNSM", 1, 0, 1])
	lipidList.append(["BNSM", 0, 0, 1])
	lipidList.append(["BNSM", 1, 0, 1])
	lipidList.append(["XNSM", 0, 0, 1])
	lipidList.append(["XNSM", 1, 0, 1])
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		if lipidList[i][1]==0:
			upperSMindices.append(i)
		else:
			lowerSMindices.append(i)

	# CEs - Ceramides
	# Calculate inner/outer together, as these flip flop
	CEindices=[]
	lipidList.append(["DPCE", 2, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DXCE", 2, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["PNCE", 2, 0, 1])
	lipidList.append(["XNCE", 2, 0, 1])
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		CEindices.append(i)

	# PIs - Phosphatidylinositols
	PIindices=[]
	lipidList.append(["POPI", 1, 0, 0])
	lipidList.append(["PIPI", 1, 0, 0])
	lipidList.append(["PAPI", 1, 0, 0])
	lipidList.append(["PUPI", 1, 0, 0])
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		PIindices.append(i)

	# PIPs (1-3 Phosphates)
	PIPindices=[]
	lipidList.append(["POP1", 1, 0, 0]) 
	lipidList.append(["POP2", 1, 0, 0])
	lipidList.append(["POP3", 1, 0, 0])
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		PIPindices.append(i)

	# Glycolipids - Skipping large heads
	GMindices=[]
	GM1indices=[]
	GM3indices=[]
	lipidList.append(["DPG1", 0, 0, 1])
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DXG1", 0, 0, 1]) 
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["PNG1", 0, 0, 1]) 
	lipidList.append(["XNG1", 0, 0, 1]) 
	lipidList.append(["DPG3", 0, 0, 1]) 
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["DXG3", 0, 0, 1]) 
	fully_saturated.append(len(lipidList)-1)
	lipidList.append(["PNG3", 0, 0, 1]) 
	lipidList.append(["XNG3", 0, 0, 1]) 
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		if i < oldlipidindex+(lipidindex-oldlipidindex)/2:
			GM1indices.append(i)
		else:
			GM3indices.append(i)
		GMindices.append(i)

	# CHOL - Calculate inner/outer together, as it tends to flip flop
	CHOLindices=[]
	lipidList.append(["CHOL", 2, 0, 2])
	oldlipidindex=lipidindex
	interval=len(lipidList)-lipidindex
	lipidindex+=interval
	for i in range(oldlipidindex,lipidindex):
		CHOLindices.append(i)

	# create a list of all lipids in list (minus cholesterol)
	for i in range(0,len(lipidList)-1):
		neither.append(i)

	# remove all entries that are fully saturated
	for entry in fully_saturated:
		if entry in neither:
			neither.remove(entry)

	# remove all entries that are poly-unsaturated
	for entry in poly_unsaturated:
		if entry in neither:
			neither.remove(entry)

	# Lipids groups
	groupList = []
	groupList.append(["PCupper" , 0, 0, upperPCindices])
	groupList.append(["PClower" , 1, 0, lowerPCindices])
	groupList.append(["PEupper" , 0, 0, upperPEindices])
	groupList.append(["PElower" , 1, 0, lowerPEindices])
	groupList.append(["PS"  , 1, 0, PSindices])
	groupList.append(["PA"  , 1, 0, PAindices])
	groupList.append(["DAG" , 2, 0, DAGindices])
	groupList.append(["LPC" , 0, 0, LPCindices])
	groupList.append(["SMupper" , 0, 0, upperSMindices])
	groupList.append(["SMlower" , 1, 0, lowerSMindices])
	groupList.append(["CER" , 2, 0, CEindices])
	groupList.append(["PI"  , 1, 0, PIindices])
	groupList.append(["PIPs", 1, 0, PIPindices])
	groupList.append(["GM"  , 0, 0, GMindices])
	groupList.append(["GM1" , 0, 0, GM1indices])
	groupList.append(["GM3" , 0, 0, GM3indices])
	groupList.append(["fully_saturated" , 2, 0, fully_saturated])
	groupList.append(["poly_unsaturated" , 2, 0, poly_unsaturated])
	groupList.append(["neither_fully_saturated_nor_polyunsaturated" , 2, 0, neither])

	# create your system - do I need this?
	syst=MDAnalysis.Universe(options["-g"].value,options["-x"].value)

	# check if the value of -e is -1; if it is, change it to length of trajectory
	if options["-e"].value == -1:
		options["-e"].value == len(syst.trajectory)

	# get number of frames you will process
        update_num_frames(options["-e"].value-options["-b"].value)

	# get protein selection
	s = " "
	protein_selection_str=s.join(options["-protsel"].value.split('@'))

	# Create a dictionary with all lipid names/indexes
	# change this so it points to separate leaflets
	lipidDictionary = {}
	for i, lipid in enumerate(lipidList):
	    # comment this out; possibly then get splitting of leaflets?
	    if lipid[1] == 2: # If lipid in both leaflets make inner / outer key point to same index
		lipidDictionary[lipid[0] + "_0"] = i
		lipidDictionary[lipid[0] + "_1"] = i
	    else:
		lipidDictionary[lipid[0] + "_" + str(lipid[1])] = i

	# Get total number of all lipids and how many of each lipid 
	totalLipids = 0
	for i in range(len(lipidList)):
		tempSel = syst.select_atoms(getSelection(lipidList[i]))
		lipidList[i][2] = len(tempSel)
		totalLipids += lipidList[i][2]

        # determine number of available CPUs (can include other options later), leave 2 free so 
        # the memory isn't overloaded
        available_cores = multiprocessing.cpu_count()-2

	# preallocate arrays for multiprocessing
	plTotalTemp = [ [] for x in range(available_cores) ]
	plCountsTemp = [ [] for x in range(available_cores) ]

        # distribute indices of lipids evenly for each core
        frames_list=[]
        for i in range(options["-b"].value,options["-e"].value):
                frames_list.append(i)
        index_array=numpy.asarray(frames_list)
        list_of_index_arrays_preprocessed=numpy.array_split(index_array,available_cores)
        list_of_index_arrays_for_distribution=[]
        for array in list_of_index_arrays_preprocessed:
                temparray=[]
                temparray.append(array[0])
                temparray.append(array[-1])
                list_of_index_arrays_for_distribution.append(temparray)

	# initiate the number of frames processed to 0
	num_frames_processed.value=0

	# activate a pool of workers 
	pool = multiprocessing.Pool()

	# log results from pool so I can concatenate them 
	def log_results_from_pool(listofresults):
		index=int(listofresults[0])
		plTotalTemp[index]=listofresults[1]
		plCountsTemp[index]=listofresults[2]

	# iterate through list of indices and give them contact tasks to do
	indexlistcount=0
	for indexlist in list_of_index_arrays_for_distribution:
		pool.apply_async(get_counts,args=(indexlist,lipidDictionary,indexlistcount,protein_selection_str,lipidList),callback=log_results_from_pool)
		indexlistcount+=1

	# ensure that parent process waits for all of its children to finish
	pool.close()
	pool.join()

	# concatenate arrays from the child processes into one big array
	plTotal = numpy.concatenate(plTotalTemp,axis=0)
	plCounts = numpy.concatenate(plCountsTemp,axis=0)

	# Get stats for proteins - this is if you have multiple proteins, though it will work for one
	pTotalLipids = numpy.zeros(options["-numprot"].value)
	for i in range(len(pTotalLipids)):
		pTotalLipids[i] = numpy.sum(numpy.average(plCounts[:,:,i], 0))
	avgProteinLipids =  numpy.average(numpy.sum(numpy.average(plCounts[:,:,:], 0),0))

	# Output all data for each lipid separately
	output = open(options["-outfile"].value, 'w')
	output.write("# name | leaflet | total lipid count | lipid counts around protein 1..X | average count around protein | SD count around protein | fraction around protein | enrichment around protein\n")
	for i in range(len(lipidList)):
		# print the following in one row:
		# col0 -> name,
		# col1 -> upper=0 / lower=1 / both=2
		# col2 -> total lipid count 
		# col3,..,X -> lipid counts around protein 1,..,X
		# colX+1 -> average count around protein
		# colX+2 -> sd count around protein
		# colX+3 -> fraction around protein
		# colX+4 -> enrichment around protein 
		tempAvgCount = numpy.average(plCounts[:,i,:])
		tempSDCount = numpy.std(plCounts[:,i,:])  # WARNING this is SD of all frames and all proteins together  
		tempPfraction = tempAvgCount / float(avgProteinLipids)
		tempEnrichment =  tempPfraction/(lipidList[i][2]/float(totalLipids))
		outString = "%s %i %i" % (lipidList[i][0], lipidList[i][1], lipidList[i][2])
		for j in range(len(pTotalLipids)):
			outString += " %.3f" % (numpy.average(plCounts[:,i,j]))
		outString += " %.3f %.3f %.3f %.3f" % (tempAvgCount, tempSDCount, tempPfraction, tempEnrichment)
		output.write(outString + "\n")

	# convert lipidList to numpy array
	lipidNumpyList = numpy.array(lipidList)

	# loop over grouped lipids and write an output string
	for i in range(len(groupList)):
		groupLipidCount = 0
		for j in lipidNumpyList[groupList[i][3],2]:
			groupLipidCount += int(j)

		# get the average and standard deviation of each group in the group list 
		tempAvgCount = numpy.average(numpy.sum(plCounts[:,groupList[i][3],:],1))
		tempSDCount = numpy.std(numpy.average(numpy.sum(plCounts[:,groupList[i][3],:],1),0))  # now getting SD between the proteins only

		# get the protein fraction and then the enrichment
		tempPfraction = tempAvgCount/float(avgProteinLipids)
		tempEnrichment =  tempPfraction/(groupLipidCount/float(totalLipids))
		outString = "%s %s %i" % (groupList[i][0], groupList[i][1], groupLipidCount)
		for j in range(len(pTotalLipids)):
			outString += " %.3f" % ( numpy.average(numpy.sum(plCounts[:,groupList[i][3],j],1)) )

		# write the average number of lipids, standard deviation, protein fraction and enrichment to a file
		outString += " %.3f %.3f %.3f %.3f" % (tempAvgCount, tempSDCount, tempPfraction, tempEnrichment)
		output.write(outString + "\n")

	# close output file
	output.close()

	# Print summary
	outString = "Total frames = %i, Total lipids = %i, pdist = %i \n" % ((options["-e"].value-options["-b"].value), totalLipids, options["-cutoff"].value)
	for i in range(len(pTotalLipids)):
		outString += " protein %i total lipids %.3f \n" % (i+1, pTotalLipids[i])
	print(outString)

# main function
if __name__=="__main__":
	main()
