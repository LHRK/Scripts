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
protein = u.select_atoms('name BB SC1 SC2 SC3 SC4')

All_headgroups = u.select_atoms('name P*')

with md.selections.gromacs.SelectionWriter('NDX/{0:s}_thickness.ndx'.format(name), mode='w') as ndx:
	ndx.write(All_headgroups, name='headgroups')
	ndx.write(protein, name='protein')


