#!/usr/bin/python

import time
from numpy import *
import MDAnalysis
from base.tools import status
from MDAnalysis.analysis.rms import RMSD

def protein_rmsd(grofile,trajfile,**kwargs):

	"""
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
		
	#---prepare universe	
	slice_name = kwargs['calc']['slice_name']
	group = kwargs['calc']['group']
	grofile,trajfile = [work.slice(sn)[slice_name][group][i] for i in ['gro','xtc']]
	uni = MDAnalysis.Universe(work.postdir+grofile,work.postdir+trajfile)
	nframes = len(uni.trajectory)
	protein = uni.select_atoms('protein and name CA')

	#---reference frame
	uni.trajectory[0]
	r0 = protein.coordinates()
	r0 -= mean(r0,axis=0)

	#---collect coordinates
	nframes = len(uni.trajectory)
	coords = []
	for fr in range(0,nframes):
		uni.trajectory[fr]
		r1 = protein.coordinates()
		coords.append(r1)

#def kinase_domains(kinase_domain_file='all_kinase.csv'):
        

def rmsd_multi(grofile, trajfile,
               active_struct='/home/joe/alk/analysis/active_wt_alk.pdb',
               inactive_struct='/home/joe/alk/analysis/inactive_wt_alk.pdb',
               selections='all', **kwargs):
                
	#---unpack
	sn = kwargs['sn']
	slice_name = kwargs['slice_name']
	work = kwargs['workspace']
	result = {}
	#result['rmsds'] = {}

	#---prepare universe	
	slice_name = kwargs['calc']['slice_name']
	group = kwargs['calc']['group']
	grofile,trajfile = [work.slice(sn)[slice_name][group][i] for i in ['gro','xtc']]
	uni = MDAnalysis.Universe(work.postdir+grofile,work.postdir+trajfile)
	nframes = len(uni.trajectory)
	protein = uni.select_atoms('protein and name CA')
        
        #---reference structures
        act_ref=MDAnalysis.Universe(active_struct)
        inact_ref=MDAnalysis.Universe(inactive_struct)

	header = "Protein_name,kd_start,kd_end,ploop_start,ploop_end,alphac_start,alphac_end"+\
		",catloop_start,catloop_end,activation_start,activation_end"
	header = header.split(',')
        alk='ALK,1116,1392,1123,1128,1158,1176,1243,1254,1268,1291'.split(',')

	alphac_start = int(alk[header.index('alphac_start')])
	alphac_end = int(alk[header.index('alphac_end')])
	aloop_start = int(alk[header.index('activation_start')])
	aloop_end = int(alk[header.index('activation_end')])

	if selections=='ach' or selections=='all':
		align_sel='name CA and not (resid %s-%s)'%(alphac_start,alphac_end)
		rmsd_sel='name CA and resid %s-%s'%(alphac_start,alphac_end)
		mod_act='ach_act';mod_inact='ach_inact'
		
		act_rmsd=RMSD(uni,act_ref,select=align_sel,groupselections=[rmsd_sel])
		act_rmsd.run()
		inact_rmsd=RMSD(uni,inact_ref,select=align_sel,groupselections=[rmsd_sel])
		inact_rmsd.run()
                
		result[mod_act] = act_rmsd.rmsd
		result[mod_inact] = inact_rmsd.rmsd

	if selections=='aloop' or selections=='all':
		align_sel='name CA and not resid %s-%s'%(aloop_start,aloop_end)
		rmsd_sel='name CA and resid %s-%s'%(aloop_start,aloop_end)
		mod_act='aloop_act';mod_inact='aloop_inact'
		
		act_rmsd=RMSD(uni,act_ref,select=align_sel,groupselections=[rmsd_sel])
		act_rmsd.run()
		inact_rmsd=RMSD(uni,inact_ref,select=align_sel,groupselections=[rmsd_sel])
		inact_rmsd.run()
		
		result[mod_act] = act_rmsd.rmsd
		result[mod_inact] = inact_rmsd.rmsd
			
	if selections=='ach_aloop' or selections=='all':
		align_sel='name CA and not (resid %s-%s or resid %s-%s)'%(
			alphac_start,alphac_end,aloop_start,aloop_end)
		rmsd_sel='name CA and (resid %s-%s or resid %s-%s)'%(
			alphac_start,alphac_end,aloop_start,aloop_end)
		mod_act='ach_aloop_act';mod_inact='ach_aloop_inact'

		act_rmsd=RMSD(uni,act_ref,select=align_sel,groupselections=[rmsd_sel])
		act_rmsd.run()
		inact_rmsd=RMSD(uni,inact_ref,select=align_sel,groupselections=[rmsd_sel])
		inact_rmsd.run()
		
		result[mod_act] = act_rmsd.rmsd
		result[mod_inact] = inact_rmsd.rmsd
			
	if selections=='CA' or selections=='all':
		align_sel='name CA'
		rmsd_sel='name CA'
		mod_act='CA_act';mod_inact='CA_inact'
		
		act_rmsd=RMSD(uni,act_ref,select=align_sel,groupselections=[rmsd_sel])
		act_rmsd.run()
		inact_rmsd=RMSD(uni,inact_ref,select=align_sel,groupselections=[rmsd_sel])
		inact_rmsd.run()
		
		result[mod_act] = act_rmsd.rmsd
		result[mod_inact] = inact_rmsd.rmsd
			
	#---pack
	attrs = {}
	
	return result,attrs	

#---TRY AGAIN

def compute_rmsd(coords,align_selection=None,rmsd_selection=None):

	"""
	Return RMSDs and rotations for a single trajectory.
	"""

	if type(align_selection)==type(None): alignsel = slice(None,None)
	else: alignsel = array(align_selection)
	if type(rmsd_selection)==type(None): rmsdsel = slice(None,None)
	else: rmsdsel = array(rmsd_selection)
	#---simple RMSD code
	rmsds = []
	r0 = coords[0]
	r0 -= mean(r0,axis=0)
	nframes = len(coords)
	for fr in range(0,nframes):
		status('RMSD',i=fr,looplen=nframes)
		r1 = coords[fr]
		r1 -= mean(r1,axis=0)
		#---computation of RMSD validated against VMD but no reflection
		U,s,Vt = linalg.svd(dot(r0[alignsel].T,r1[alignsel]))
		signer = identity(3)
		signer[2,2] = sign(linalg.det(dot(Vt.T,U)))
		RM = dot(dot(U,signer),Vt)
		rmsds.append(sqrt(mean(sum((r0[rmsdsel].T-dot(RM,r1[rmsdsel].T))**2,axis=0))))
	return rmsds
		
def collect_coordinates(uni,sel):

	"""
	Given a universe and a selection, return the XYZ coordinates for each frame.
	"""

	#---collect coordinates
	nframes = len(uni.trajectory)
	coords = []
	for fr in range(0,nframes):
		uni.trajectory[fr]
		r1 = sel.coordinates()
		coords.append(r1)
	return array(coords)

		
def rmsd_multiX2(grofile,trajfile,*args,**kwargs):

	"""
	???
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	results,attrs = {},{}
	if sn == 'simulation-v00117': return results,attrs

	#---TEMPORARY HARD_CODED DATA
	active_struct='/home/keerthana/erbb2/her2_active.pdb'
	inactive_struct='/home/keerthana/erbb2/winners/her2_inact_renum.pdb'
	header = "Protein_name,kd_start,kd_end,ploop_start,ploop_end,alphac_start,alphac_end"+\
		",catloop_start,catloop_end,activation_start,activation_end"
	header = header.split(',')

	#---prepare universe	
	grofile,trajfile = [work.slices[sn]['current']['protein'][i] for i in ['gro','xtc']]
	uni = MDAnalysis.Universe(work.postdir+grofile,work.postdir+trajfile)

	#---numbering issues
	sel_traj = uni.select_atoms('name CA')
	badnumber = uni.select_atoms('all').resids[0]==1
	
	erbb2 = "ERBB2,720,987,727,732,760,778,839,850,861,884".split(',')
	alphac_start = int(erbb2[header.index('alphac_start')]) - (720 if badnumber else 0)
	alphac_end = int(erbb2[header.index('alphac_end')]) - (720 if badnumber else 0)
	catloop_start = int(erbb2[header.index('catloop_start')]) - (720 if badnumber else 0)
	catloop_end = int(erbb2[header.index('catloop_end')]) - (720 if badnumber else 0)
	comparisons = {
		'ach':[(alphac_start,alphac_end)],
		'catloop':[(catloop_start,catloop_end)],
		'ach-catloop':[(alphac_start,alphac_end),(catloop_start,catloop_end)],
		'ca':None,
		}

	act_ref = MDAnalysis.Universe(active_struct)
	inact_ref = MDAnalysis.Universe(inactive_struct)
	
	for key,startstops in comparisons.items():
		if startstops == None: selection = 'name CA'
		else: selection = '('+' or '.join(['(resid %s-%s)'%s for s in startstops])+')'

		sel_active = act_ref.select_atoms('name CA')
		sel_active_body = act_ref.select_atoms('name CA and not '+selection)
		sel_active_loop = act_ref.select_atoms('name CA and '+selection)
		sel_inactive = inact_ref.select_atoms('name CA')
		sel_inactive_body = inact_ref.select_atoms('name CA and not '+selection)
		sel_inactive_loop = inact_ref.select_atoms('name CA and '+selection)
		
		sel_traj = uni.select_atoms('name CA')
		sel_traj_body = uni.select_atoms('name CA and not %s'%selection)
		sel_traj_loop = uni.select_atoms('name CA and %s'%selection)
	
		active_reference = sel_active.coordinates()
		inactive_reference = sel_inactive.coordinates()
		coords = collect_coordinates(uni,sel_traj)
		#---beware hack below
		coords_active = concatenate(([sel_active.coordinates()],coords[:,:267]))
		coords_inactive = concatenate(([sel_inactive.coordinates()],coords[:,:]))
		#---decrement resids by 1 to match the CA coordinates
		align_selection = sel_traj_body.resids-1-(720 if not badnumber else 0)
		rmsd_selection = sel_traj_loop.resids-1-(720 if not badnumber else 0)
		if selection == 'name CA':
			align_selection = sel_traj.resids[:267]-1-(720 if not badnumber else 0)
			rmsd_selection = sel_traj.resids[:267]-1-(720 if not badnumber else 0)
		rmsd_inactive = compute_rmsd(coords_inactive,
			align_selection=align_selection,rmsd_selection=rmsd_selection)
		#---beware hack continues
		align_selection = array([i for i in align_selection if i<(267 if badnumber else 987)])
		rmsd_active = compute_rmsd(coords_active,
			align_selection=align_selection,rmsd_selection=rmsd_selection)
		results['rmsd_active_%s'%key] = array(rmsd_active)
		results['rmsd_inactive_%s'%key] = array(rmsd_inactive)
	return results,attrs
	
