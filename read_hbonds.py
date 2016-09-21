#!/usr/bin/python -i
#-*- coding: utf-8 -*-

#---plot prep
execfile('./omni/base/header.py')
import numpy as np
import json

#---settings
plotname = 'hbonds'

#---load the upstream data
#data,calc = plotload(plotname,work)

def read_hbonds(**kwargs):


	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
        data=kwargs['upstream']['hbonds']
        divy=kwargs['calc']['specs']['divy']
        active=work.meta[sn]['active']
        name=work.meta[sn]['name']

        #---compute
        bonds=hbonds_info(data,timesteps=True,donor_reslist=None,
                         acceptor_reslist=None,divy=divy)

        #---pack
	attrs,results = {},{}
        attrs={'divy':divy, 'active':active, 'name':name}
        results['bonds']=np.string_(json.dumps(bonds))
        return results,attrs

def hbonds_info(data,timesteps=False,donor_restype=None,acceptor_restype=None,
                           donor_reslist=None,acceptor_reslist=None,divy=False):
    hbonds={}
    main_donor='N'
    main_acceptor='O'
    if timesteps: 
        raw_hbonds=data['hbonds_timesteps'][0]
        num_steps=len(set([i[9] for i in raw_hbonds]))
    else: raw_hbonds=data['hbonds_occupancies']
    if divy: divies={}
    for hbond in raw_hbonds:
        store_hbond=True
        if donor_restype and hbond[2] not in donor_restype:
            store_hbond=False               
        if acceptor_restype and hbond[6] not in acceptor_restype:
            store_hbond=False
        if donor_reslist and int(hbond[3]) not in donor_reslist:
            store_hbond=False
        if acceptor_reslist and int(hbond[7]) not in acceptor_reslist:
            store_hbond=False
        if store_hbond:
            if divy:
                if hbond[4] == main_donor:
                    donor_loc='main'
                else:
                    donor_loc='side'
                if hbond[8] == main_acceptor:
                    acceptor_loc='main'
                else:
                    acceptor_loc='side'
                label='{0} {1} {2} {3} {4} {5}'.format(hbond[2], hbond[3], donor_loc,
                                                       hbond[6], hbond[7], acceptor_loc)
                if label in divies.keys():
                    if timesteps:
                        time=divies[label]
                        time.append(hbond[9])
                        divies[label]=time
                    else: divies[label]+=float(hbond[9])
                else:
                    if timesteps: divies[label]=[hbond[9]]
                    else: divies[label]=float(hbond[9])
            else:
                label='{0} {1} {2} {3} {4} {5}'.format(hbond[2], hbond[3], hbond[4],
                                                       hbond[6], hbond[7], hbond[8])
                if timesteps:
                    if label in hbonds.keys():
                        time=hbonds[label]['times']
                        time.append(hbond[9])
                        hbonds[label]['times']=time
                    else: hbonds[label]={'donor_restype':hbond[2], 'donor_residx':hbond[3],
                                             'donor_HA':hbond[4], 'acceptor_restype':hbond[6],
                                             'acceptor_residx':hbond[7], 'acceptor_HA':hbond[8],
                                             'times':[hbond[9]]}
                else:
                    if label in hbonds.keys():
                        hbonds[label]['occupancy']+=float(hbond[9])
                    else:
                        hbonds[label]={'donor_restype':hbond[2], 'donor_residx':hbond[3],
                                           'donor_HA':hbond[4], 'acceptor_restype':hbond[6],
                                           'acceptor_residx':hbond[7], 'acceptor_HA':hbond[8],
                                           'occupancy':float(hbond[9])}

    if divy: 
        for key,val in divies.items():
            infos=key.split(' ')
            if timesteps:
                hbonds[key]={'donor_restype':infos[0], 'donor_residx':infos[1],
                                 'donor_loc':infos[2], 'acceptor_restype':infos[3],
                                 'acceptor_residx':infos[4], 'acceptor_loc':infos[5],
                                 'times':val}
            else:
                hbonds[key]={'donor_restype':infos[0], 'donor_residx':infos[1],
                                 'donor_loc':infos[2], 'acceptor_restype':infos[3],
                                 'acceptor_residx':infos[4], 'acceptor_loc':infos[5],
                                 'occupancy':float(val)}
    for hb in hbonds.keys():
        if timesteps:
            hbonds[hb]['occupancy']=float('{:0.5f}'.format(len(
                hbonds[hb]['times'])/float(num_steps)))
        else:
            hbonds[hb]['occupancy']=float('{:0.5f}'.format(hbonds[hb]['occupancy']))
    return hbonds



