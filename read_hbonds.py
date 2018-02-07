#!/usr/bin/env python

import numpy as np
import json

def read_hbonds(**kwargs):

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	data=kwargs['upstream']['hbonds']
	divy=kwargs['calc']['specs']['divy']

	#this hack guesses the chain breaks based on residue numbering                                      
	#if residue[i]-residue[i-1]!=1 then assume a chain break                                            
	if 'startres_idxs' in data:
		startres_idxs=data['startres_idxs']
		alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ' #max number of chains=26                              
		if len(startres_idxs)>len(alphabet):
			raise Exception("you can't have more than 26 chains")
		chain_locs=[[alphabet[i],startres_idxs[i],startres_idxs[i+1]-1] for i in range(len(startres_idxs)-1)]
		#for last chain residx just needs to be a number                                            
		chain_locs.append([alphabet[len(startres_idxs)-1],startres_idxs[len(startres_idxs)-1],99999])
	else: chain_locs=None
	if 'kcat' in work.meta[sn]: active=work.meta[sn]['kcat']
	else: active=None
	if 'name' in work.meta[sn]: name=work.meta[sn]['name']
	else: name=None

	#---compute
	bonds=hbonds_info(data,timesteps=True,
					  donor_reslist=None,
					  acceptor_reslist=None,
					  chain_locs=chain_locs,
					  divy=divy)

	#---pack
	attrs,results = {},{}
	attrs={'divy':divy, 'active':active, 'name':name}
	results['bonds']=np.string_(json.dumps(bonds))
	return results,attrs

def hbonds_info(data,timesteps=False,donor_restype=None,acceptor_restype=None,
				donor_reslist=None,acceptor_reslist=None,divy=False,chain_locs=None):
	hbonds={}
	main_donor='N'
	main_acceptor='O'
	if timesteps:
		raw_hbonds=data['hbonds_timesteps']
		num_steps=len(set([i[11] for i in raw_hbonds]))
	else: raw_hbonds=data['hbonds_occupancies']
	if divy: divies={}
	for hbond in raw_hbonds:
		store_hbond=True
		donor_restype=hbond[4];donor_residx=int(hbond[5])
		acceptor_restype=hbond[8];acceptor_residx=int(hbond[9])
		if donor_restype and donor_restype not in donor_restype:
			store_hbond=False
		if acceptor_restype and acceptor_restype not in acceptor_restype:
			store_hbond=False
		if donor_reslist and int(donor_residx) not in donor_reslist:
			store_hbond=False
		if acceptor_reslist and int(acceptor_residx) not in acceptor_reslist:
			store_hbond=False
		if store_hbond:
			ts=hbond[11];donor_ha=hbond[6];acceptor_ha=hbond[10]
			if divy:
				if donor_ha == main_donor:
					donor_loc='main'
				else:
					donor_loc='side'
				if acceptor_ha == main_acceptor:
					acceptor_loc='main'
				else:
					acceptor_loc='side'
				if chain_locs:
					donor_chain,acceptor_chain=get_chains(donor_residx,
														  acceptor_residx,chain_locs)
					label='{0} {1} {2} {3} {4} {5} {6} {7}'.format(donor_restype,
																   donor_residx,
																   donor_loc,
																   donor_chain,
																   acceptor_restype,
																   acceptor_residx,
																   acceptor_loc,
																   acceptor_chain)
				else:
					label='{0} {1} {2} {3} {4} {5}'.format(donor_restype,
														   donor_residx,
														   donor_loc,
														   acceptor_restype,
														   acceptor_residx,
														   acceptor_loc)
				if label in divies.keys():
					if timesteps:
						time=divies[label]
						time.append(ts)
						divies[label]=time
					else: divies[label]+=float(ts)
				else:
					if timesteps:
						divies[label]=[ts]
					else: divies[label]=float(ts)
			if not divy:
				if chain_locs:
					donor_chain,acceptor_chain=get_chains(donor_residx,
														  acceptor_residx,chain_locs)
					label='{0} {1} {2} {3} {4} {5} {6} {7}'.format(donor_restype,
																   donor_residx,
																   donor_ha,
																   donor_chain,
																   acceptor_restype,
																   acceptor_residx,
																   acceptor_ha,
																   acceptor_chain)

				else:
					label='{0} {1} {2} {3} {4} {5}'.format(donor_restype,
														   donor_residx,
														   donor_ha,
														   acceptor_restype,
														   acceptor_residx,
														   acceptor_ha)
				if timesteps:
					if label in hbonds.keys():
						time=hbonds[label]['times']
						time.append(ts)
						hbonds[label]['times']=time
					else:
						if chain_locs:
							donor_chain,acceptor_chain=get_chains(donor_residx,
																  acceptor_residx,chain_locs)
							hbonds[label]={'donor_restype':donor_restype,
										   'donor_residx':donor_residx,
										   'donor_HA':donor_ha,
										   'donor_chain':donor_chain,
										   'acceptor_restype':acceptor_restype,
										   'acceptor_residx':acceptor_residx,
										   'acceptor_HA':acceptor_ha,
										   'acceptor_chain':acceptor_chain,
										   'times':[ts]}
						else:
							hbonds[label]={'donor_restype':donor_restype,
										   'donor_residx':donor_residx,
										   'donor_HA':donor_ha,
										   'acceptor_restype':acceptor_restype,
										   'acceptor_residx':acceptor_residx,
										   'acceptor_HA':acceptor_ha,
										   'times':[ts]}
				if not timesteps:
					if label in hbonds.keys():
						hbonds[label]['occupancy']+=float(ts)
					else:
						if chain_locs:
							donor_chain,acceptor_chain=get_chains(donor_residx,
																  acceptor_residx,chain_locs)
							hbonds[label]={'donor_restype':donor_restype,
										   'donor_residx':donor_residx,
										   'donor_HA':donor_ha,
										   'donor_chain':donor_chain,
										   'acceptor_restype':acceptor_restype,
										   'acceptor_residx':acceptor_residx,
										   'acceptor_HA':acceptor_ha,
										   'acceptor_chain':acceptor_chain,
										   'occupancy':float(ts)}

						else:
							hbonds[label]={'donor_restype':donor_restype,
										   'donor_residx':donor_residx,
										   'donor_HA':donor_ha,
										   'acceptor_restype':acceptor_restype,
										   'acceptor_residx':acceptor_residx,
										   'acceptor_HA':acceptor_ha,
										   'occupancy':float(ts)}

	if divy:
		for key,val in divies.items():
			infos=key.split(' ')
			if timesteps:
				if chain_locs:
					donor_chain,acceptor_chain=get_chains(donor_residx,
														  acceptor_residx,chain_locs)
					hbonds[key]={'donor_restype':infos[0],
								 'donor_residx':infos[1],
								 'donor_loc':infos[2],
								 'donor_chain':infos[3],
								 'acceptor_restype':infos[4],
								 'acceptor_residx':infos[5],
								 'acceptor_loc':infos[6],
								 'acceptor_chain':infos[7],
								 'times':val}
				else:
					hbonds[key]={'donor_restype':infos[0],
								 'donor_residx':infos[1],
								 'donor_loc':infos[2],
								 'acceptor_restype':infos[3],
								 'acceptor_residx':infos[4],
								 'acceptor_loc':infos[5],
								 'times':val}
			if not timesteps:
				if chain_locs:
					donor_chain,acceptor_chain=get_chains(donor_residx,
														  acceptor_residx,chain_locs)
					hbonds[key]={'donor_restype':infos[0],
								 'donor_residx':infos[1],
								 'donor_loc':infos[2],
								 'donor_chain':infos[3],
								 'acceptor_restype':infos[4],
								 'acceptor_residx':infos[5],
								 'acceptor_loc':infos[6],
								 'acceptor_chain':infos[7],
								 'occupancy':float(val)}
				else:
					hbonds[key]={'donor_restype':infos[0],
								 'donor_residx':infos[1],
								 'donor_loc':infos[2],
								 'acceptor_restype':infos[3],
								 'acceptor_residx':infos[4],
								 'acceptor_loc':infos[5],
								 'occupancy':float(val)}
	for hb in hbonds.keys():
		if timesteps:
			hbonds[hb]['occupancy']=float('{:0.5f}'.format(len(
				hbonds[hb]['times'])/float(num_steps)))
		else:
			hbonds[hb]['occupancy']=float('{:0.5f}'.format(hbonds[hb]['occupancy']))
	return hbonds

def get_chains(donor_residx,acceptor_residx,chain_locs):
	donor_chain=[chain[0] for chain in chain_locs if
				 donor_residx>=chain[1] and donor_residx<=chain[2]]
	acceptor_chain=[chain[0] for chain in chain_locs if
					acceptor_residx>=chain[1] and acceptor_residx<=chain[2]]
	return donor_chain[0], acceptor_chain[0]
