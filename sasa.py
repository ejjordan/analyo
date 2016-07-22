#!/usr/bin/python

import time,sys
from numpy import *
import MDAnalysis
from base.tools import status
import Bio.PDB as bp

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
	resids = array(sorted(dssp_in.keys()))
        resnames = array([dssp_in[k]['resname'] for k in resids])
        dssp = array([dssp_in[k]['dssp'] for k in resids])
	rel_sasa = array([dssp_in[k]['rel_sasa'] for k in resids])
	abs_sasa = array([dssp_in[k]['abs_sasa'] for k in resids])
	phi = array([dssp_in[k]['phi'] for k in resids])
	psi = array([dssp_in[k]['psi'] for k in resids])

	#---pack
	attrs,results = {},{}
	results['dssp'] = dssp
	results['rel_sasa'] = rel_sasa
	results['abs_sasa'] = abs_sasa
	results['phi'] = phi
	results['psi'] = psi
	results['resid'] = resids
	results['resname'] = resnames
	return results,attrs	

def get_sasa(topology, trajectory, dssp_loc=master_dssp_location,skip=None):
        dssp_loc = dssp_loc
        DSSP={}
        universe = MDAnalysis.Universe(topology, trajectory)
        protein=universe.select_atoms("protein")
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
                        resobj=dssp[key][0]
                        resname=dssp[key][0].resname
                        residx=resobj.id[1]
                        #secondary_structure=dssp[key][1]
                        #sasa=dssp[key][2]
                        secondary_structure=resobj.xtra['SS_DSSP']
                        rel_sasa=resobj.xtra['EXP_DSSP_RASA']
                        abs_sasa=resobj.xtra['EXP_DSSP_ASA']
                        phi=resobj.xtra['PHI_DSSP']
                        psi=resobj.xtra['PSI_DSSP']
                        if residx in DSSP and DSSP[residx]['resname']==resname:
                                DSSP[residx]['dssp'].append(secondary_structure)
                                DSSP[residx]['rel_sasa'].append(rel_sasa)
                                DSSP[residx]['abs_sasa'].append(abs_sasa)
                                DSSP[residx]['phi'].append(phi)
                                DSSP[residx]['psi'].append(psi)
                        else:
                                DSSP[residx]={'dssp':[secondary_structure],'phi':[phi],
                                              'psi':[psi],'rel_sasa':[rel_sasa],
                                              'abs_sasa':[abs_sasa],'resname':resname}
        return DSSP
