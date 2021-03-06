#!/usr/bin/python

import time,sys
from numpy import *
import MDAnalysis
from base.tools import status
import Bio.PDB as bp
from common_plot import *

master_dssp_location = '/usr/bin/mkdssp'

def sasa(grofile,trajfile,**kwargs):

	"""
	LIPID ABSTRACTOR
	Reduce a bilayer simulation to a set of points.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']

	#---compute SASA
	dssp_in = get_sasa(grofile,trajfile)

	#---pack
	attrs,results = {},{}

	for chain in dssp_in.keys():
		resids = array(sorted(dssp_in[chain].keys()))
		chain_ids = array([chain]*len(dssp_in[chain].keys()))
		resnames = array([dssp_in[chain][k]['resname'] for k in resids])
		dssp = array([dssp_in[chain][k]['dssp'] for k in resids])
		rel_sasa = array([dssp_in[chain][k]['rel_sasa'] for k in resids])
		abs_sasa = array([dssp_in[chain][k]['abs_sasa'] for k in resids])
		phi = array([dssp_in[chain][k]['phi'] for k in resids])
		psi = array([dssp_in[chain][k]['psi'] for k in resids])
		time = array([dssp_in[chain][k]['time'] for k in resids])

		if 'resid' in results.keys() and 'dssp' in results.keys():
			results['dssp'] = append(results['dssp'],dssp,axis=0)
			results['rel_sasa'] = append(results['rel_sasa'],rel_sasa,axis=0)
			results['abs_sasa'] = append(results['abs_sasa'],abs_sasa,axis=0)
			results['phi'] = append(results['phi'],phi,axis=0)
			results['psi'] = append(results['psi'],psi,axis=0)
			results['resid'] = append(results['resid'],resids,axis=0)
			results['chain_id'] = append(results['chain_id'],chain_ids,axis=0)
			results['resname'] = append(results['resname'],resnames,axis=0)
			results['time'] = append(results['time'],resnames,axis=0)
		else:
			results['dssp'] = dssp
			results['rel_sasa'] = rel_sasa
			results['abs_sasa'] = abs_sasa
			results['phi'] = phi
			results['psi'] = psi
			results['resid'] = resids
			results['chain_id'] = chain_ids
			results['resname'] = resnames
			results['time'] = time
	return results,attrs

def get_sasa(topology, trajectory, dssp_loc=master_dssp_location,skip=None):

	"""
	This function currently only works with one or two chains, because I am lazy.
	"""

	dssp_loc = dssp_loc
	DSSP={'A':{}}
	universe = MDAnalysis.Universe(topology, trajectory)

	#set the chain name here. this will only work for MDAnalysis 0.16
	chain_name=universe.add_Segment(segid='A')
	universe.residues[...].segments=chain_name

	protein=universe.select_atoms("protein")
	diff_res=[]
	#this attempt to identify chain breaks will only work if the resids
	#... in the chains are not numbered consecutively
	for i in range(len(protein.resnums)):
		if protein.resnums[i]-protein.resnums[i-1]<0 and i!=0:
			diff_res.append(i)
	if len(diff_res)>=1:
		chain_sep=diff_res.pop(0)
		chain_end=len(protein.resnums)
		bchain=protein[chain_sep:chain_end]
		bchain.set_segids('B')
		DSSP['B']={}

	for ts in universe.trajectory:
		if skip:
			universe.trajectory.skip=skip
		sys.stdout.flush()
		sys.stdout.write('\rsasa [step {0}]  '.format(
			universe.trajectory.frame))
		writer=MDAnalysis.Writer("tmp.pdb")
		writer.write(protein)
		writer.close()
		parser=bp.PDBParser()
		structure=parser.get_structure('tmp','tmp.pdb')
		dssp=bp.DSSP(structure[0],'tmp.pdb',dssp_loc)
		for key in dssp.keys():
			if 0:
				resobj=dssp[key][0]
				resname=dssp[key][0].resname
				residx=resobj.id[1]
				chain=key[0]
				secondary_structure=resobj.xtra['SS_DSSP']
				rel_sasa=resobj.xtra['EXP_DSSP_RASA']
				abs_sasa=resobj.xtra['EXP_DSSP_ASA']
				phi=resobj.xtra['PHI_DSSP']
				psi=resobj.xtra['PSI_DSSP']
			resobj=dssp[key]
			resname=residue_codes_reverse[resobj[1]]
			residx=key[1][1]
			chain=key[0]
			secondary_structure=resobj[2]
			rel_sasa=resobj[3]
			abs_sasa=resobj[3]*dssp.residue_max_acc[resname]
			phi=resobj[4]
			psi=resobj[5]
			if residx in DSSP[chain] and DSSP[chain][residx]['resname']==resname:
				DSSP[chain][residx]['dssp'].append(secondary_structure)
				DSSP[chain][residx]['rel_sasa'].append(rel_sasa)
				DSSP[chain][residx]['abs_sasa'].append(abs_sasa)
				DSSP[chain][residx]['phi'].append(phi)
				DSSP[chain][residx]['psi'].append(psi)
				DSSP[chain][residx]['time'].append(ts.time)
			else:
				DSSP[chain][residx]={'dssp':[secondary_structure],'phi':[phi],'time':[ts.time],
									 'psi':[psi],'rel_sasa':[rel_sasa],'chain':chain,
									 'abs_sasa':[abs_sasa],'resname':resname}
	return DSSP
