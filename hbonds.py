#!/usr/bin/python

import time
from numpy import *
import MDAnalysis
import MDAnalysis.analysis.hbonds as hb
from base.tools import status

def hbonds(grofile,trajfile,**kwargs):

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']

	#---compute SASA
	uni = MDAnalysis.Universe(grofile,trajfile)
	occupancies,timesteps,counts_per_timestep,startres_idxs = get_hbonds(uni)

	#---pack
	attrs,results = {},{}
	results['hbonds_occupancies']=array(occupancies)
	results['hbonds_timesteps']=array(timesteps)
	results['hbonds_counts_per_timestep']=array(counts_per_timestep)
	results['startres_idxs']=array(startres_idxs)
	return results,attrs

def get_hbonds(universe, distance=3.2, angle=150.0):
	hbond_distance=distance
	hbond_angle=angle
	hbonds=hb.HydrogenBondAnalysis(universe, 'protein', 'protein',
								   update_selection1=True,
								   update_selection2=True,
								   detect_hydrogens='distance',
								   distance = hbond_distance,
								   angle = hbond_angle,
								   distance_type="heavy")
	hbonds.run()
	hbonds=hbonds
	hbonds_occupancy=hbonds.count_by_type()
	hbonds_occupancy_list=hbonds_occupancy.tolist()
	hbonds_by_timestep=hbonds.timesteps_by_type()
	hbonds_by_timestep_list=hbonds_by_timestep.tolist()
	hbonds_per_timestep=hbonds.count_by_time()
	hbonds_per_timestep_list=hbonds_per_timestep.tolist()

	#do a naive splitting of chains based on indexes
	startres_idxs=[res.resindex for ii,res in enumerate(universe.residues) if
				   res.resid-universe.residues[ii-1].resid!=1]
	return hbonds_occupancy_list,hbonds_by_timestep_list,hbonds_per_timestep_list,startres_idxs
