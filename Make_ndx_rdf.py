#!/usr/bin/env python

import numpy as np
import argparse
import MDAnalysis as md
from MDAnalysis.analysis.leaflet import LeafletFinder as LF


parser = argparse.ArgumentParser(description='Script for making index files')
parser.add_argument('-n', dest='name', action='store', type=str, help='-n for name of the system (without nowat, its the raw gro file)')

args = parser.parse_args()
name = args.name

u = md.Universe('GRO/{0:s}.gro'.format(name))


L = LF(u, 'name P*', cutoff=20.0)
leaflet0 = L.groups(0).residues.atoms
leaflet1 = L.groups(1).residues.atoms

memb = u.select_atoms('resname POP*')

lipids = np.unique(u.select_atoms('resname POP*').resnames)

sele = []

for l in lipids:
	sele.append(u.select_atoms('resname {0:s} and name PO4'.format(l)))

prot = u.select_atoms('name BB SC1 SC3 SC2 SC4')
BB   = u.select_atoms('name BB')
All_headgroups = u.select_atoms('name P*')
sel_resnames = " ".join(['{0:s}'.format(i) for i in np.unique(All_headgroups.resnames)])
chol = u.select_atoms("resname CHOL")
pap = u.select_atoms("resname PAP6")
gm3 = u.select_atoms("resname DPG3")
sm  = u.select_atoms("resname DPSM")

with md.selections.gromacs.SelectionWriter('NDX/{0:s}_rdf.ndx'.format(name), mode='w') as ndx:
	ndx.write(chol, name='CHOL')
	ndx.write(pap, name='PAP6')
	ndx.write(gm3, name='DPGM3')
	ndx.write(sm, name='DPSM')
	ndx.write(prot, name='Protein')
	for idx, s in enumerate(sele):
		print (lipids[idx])
		ndx.write(s, name=lipids[idx])


