#!/usr/bin/env python

import numpy as np
import re
import json

def hbonds_timesteps(work,read_hbonds,hbond_keys,donor_restype=None,acceptor_restype=None,
                     donor_reslist=None,acceptor_reslist=None,divy=False):
    hbonds={}
    main_donor='N'
    main_acceptor='O'
    for sn in hbond_keys:
        raw_hbonds=json.loads(str(read_hbonds[sn]['data']['bonds']))
        hbonds[sn]={}
        if divy: divies={}
        for hbond in raw_hbonds:
            store_hbond=True
            donor_residx=int(raw_hbonds[hbond]['donor_residx'])
            acceptor_residx=int(raw_hbonds[hbond]['acceptor_residx'])
            if donor_restype and raw_hbonds[hbond]['donor_restype'] not in donor_restype:
                store_hbond=False
            if acceptor_restype and raw_hbonds[hbond]['acceptor_restype'] not in acceptor_restype:
                store_hbond=False
            if donor_reslist and donor_residx not in donor_reslist:
                store_hbond=False
            if acceptor_reslist and acceptor_residx not in acceptor_reslist:
                store_hbond=False
            if store_hbond:
                hbonds[sn][hbond]=raw_hbonds[hbond]
    for name in hbonds: hbonds[name]['kcat']=work.meta[name]['kcat']
    return hbonds

def timesteps_discard(hbonds_timesteps,start_time,stop_time,timestep=20):
    total_time=float(stop_time-start_time)/timestep
    for sn in hbonds_timesteps:
        for hbond in hbonds_timesteps[sn]:
            if hbond=='kcat': continue
            timesteps=[]
            for time in hbonds_timesteps[sn][hbond]['times']:
                if float(time)>start_time and float(time)<stop_time:
                    timesteps.append(time)
            occupancy=len(timesteps)/total_time
            hbonds_timesteps[sn][hbond]['times']=timesteps
            hbonds_timesteps[sn][hbond]['occupancy']=occupancy

def intradomain_bonds_discard(hbonds_timesteps,discard_residues,print_bonds=False):
    trash_bonds={}
    for sn in hbonds_timesteps:
        for hbond in hbonds_timesteps[sn]:
            if hbond=='kcat': continue
            donor_residx=int(hbonds_timesteps[sn][hbond]['donor_residx'])
            acceptor_residx=int(hbonds_timesteps[sn][hbond]['acceptor_residx'])
            if donor_residx in discard_residues and acceptor_residx in discard_residues:
                trash_bonds[hbond]=1
    if print_bonds:
        print trash_bonds.keys()
    for bond in trash_bonds.keys():
        for sn in hbonds_timesteps:
            if bond in hbonds_timesteps[sn].keys():
                hbonds_timesteps[sn].pop(bond)

def combine_hbonds(unpacked_hbonds,sorted_keys,num_replicates=2,timesteps=False,divy=False):
    # this will only work with a sorted list of keys that contain num_replicate number of
    # replicates sequentially
    if len(sorted_keys)%num_replicates!=0:
        print "[Error] number of keys is not a multiple of number of replicates; returning."
        return
    combined_hbonds={};count=0
    for i in range(len(sorted_keys)/num_replicates):
        bonds={};key_list=[]
        for j in range(num_replicates):
            name=sorted_keys[count]
            count+=1
            key_list.append(name)
            for label in unpacked_hbonds[name]:
                if label=='kcat': continue
                if label in bonds: bonds[label].append(unpacked_hbonds[name][label]['occupancy'])
                else: bonds[label]=[unpacked_hbonds[name][label]['occupancy']]
        name=min(key_list,key=len)
        combined_hbonds[name]={}
        for label,occupancies in bonds.items():
            if divy:
                infos=label.split(' ')
                combined_hbonds[name][label]={'donor_restype':infos[0], 'donor_residx':infos[1],
                                              'donor_loc':infos[2], 'acceptor_restype':infos[3],
                                              'acceptor_residx':infos[4], 'acceptor_loc':infos[5],
                                              'mean':float('{:0.5f}'.format(np.mean(occupancies))),
                                              'std':float('{:0.5f}'.format(np.std(occupancies))),
                                              'max':max(occupancies), 'min':min(occupancies)}
            else:
                combined_hbonds[name][label]={'donor_restype':infos[0], 'donor_residx':infos[1],
                                              'donor_HA':infos[2], 'acceptor_restype':infos[3],
                                              'acceptor_residx':infos[4], 'acceptor_HA':infos[5],
                                              'std':float('{:0.5f}'.format(np.std(occupancies))),
                                              'max':max(occupancies), 'min':min(occupancies)}
    for name in combined_hbonds:
        combined_hbonds[name]['name']=' '.join(name.split('_'))
        combined_hbonds[name]['kcat']=unpacked_hbonds[name]['kcat']
    return combined_hbonds

def occupancy_diff(hbonds,reference=None,threshold=0.4):
    altered_donors={}
    bond_list=list(set([label for name in hbonds for label in hbonds[name]]))
    for sn in hbonds:
        if sn==reference: continue
        altered_donors[sn]={'delta':0, 'name':hbonds[sn]['name'], 'kcat':hbonds[sn]['kcat'],
                            'bonds':[]}
        for bond in bond_list:
            if bond=='name' or bond=='kcat': continue
            if bond in hbonds[sn].keys() or bond in hbonds[reference].keys():
                if bond not in hbonds[reference].keys(): hbonds[reference][bond]={'mean':0}
                if bond not in hbonds[sn].keys(): hbonds[sn][bond]={'mean':0}
                delta_wt=float(hbonds[sn][bond]['mean'])-float(hbonds[reference][bond]['mean'])
                if abs(delta_wt)>threshold:
                    altered_donors[sn]['delta']+=delta_wt
                    altered_donors[sn]['bonds'].append(bond)
    means=[[k,v['delta']] for k,v in altered_donors.items()]
    sorted_keys=[i[0] for i in sorted(means,key=lambda x: x[1],reverse=True)]
    return altered_donors,sorted_keys

def label_maker(deltas, kcat_cut=3, name_list=None,val_type='deltas'):

    """
    This function takes a 'delta' object and a kcat cut-off and returns activation (or
    non-activation) labels. If supplied a name_list the labels will be returned in the
    specified order.
    """

    if not name_list: name_list=[deltas[sn]['name'] for sn in deltas]
    if val_type=='deltas': kcats=[deltas[name]['kcat'] for name in name_list]
    else: kcats=values
    labels=[]
    for kcat in kcats:
        if kcat=='X': labels.append('maybe')
        elif kcat=='WT': labels.append('wt')
        elif float(kcat)>=kcat_cut:
            labels.append(True)
        else:
            labels.append(False)
    return labels
