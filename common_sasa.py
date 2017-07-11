#!/usr/bin/python -i

import numpy as np
import re



def filter_sasas(data,work,SASA_keys,sasa_type='abs_sasa',
                 base_restype=None,comp_restype=None,res_list=None):

    """
    This function unpacks the data from a sasa file according to use specifications.
    """

    SASAs={}
    for sn in SASA_keys:
        raw_sasa=zip(data[sn]['data']['resname'],data[sn]['data']['resid'],
                       data[sn]['data'][sasa_type])
        SASAs[sn]={'name':' '.join(sn.split('_')), 'kcat':work.meta[sn]['kcat'],
                   'base_sasa':{}}
        if comp_restype: SASAs[sn]['comp_sasa']={}
        if type(res_list)!=list: 
            res_list=data[sn]['data']['resid']
            SASA_sums=np.sum([i[2] for i in raw_sasa if i[1] in res_list],axis=0)
        else:
            base_sasa=[];comp_sasa=[];SASA_sums=[]
            for elt in raw_sasa:
                resname=elt[0];resid=elt[1];sasa=elt[2]
                if resid in res_list:
                    if base_restype and resname in base_restype:
                        SASAs[sn]['base_sasa'][resid]={'resid':resid, 'resname':resname,
                                                         'sasa_vals': sasa}
                    if comp_restype and resname in comp_restype:
                        SASAs[sn]['comp_sasa'][resid]={'resid':resid, 'resname':resname,
                                                       'sasa_vals': sasa}
                    if not base_restype:                        
                        SASA_sums.append(sasa)
                        SASAs[sn]['base_sasa'][resid]={'resid':resid, 'resname':resname,
                                                         'sasa_vals': sasa}
    return SASAs

def combine_replicates(unpacked_SASAs):

    """
    This function will combine SASA values from replicate simulations. The names of the keys
    in the incoming unpacked_SASAs dict must contain the string 'active_wt', 'inactive_wt',
    or the mutation code in the format (from_AA)(residue_num)(to_AA). It will give a warning
    (but not an error!) if not all names have the same number of replicates.
    """

    sns=unpacked_SASAs.keys()
    unique_keys={}
    search_string='(\w\d{1,5}\w|active_wt|inactive_wt)'
    for sn in sns:
        name=re.search(search_string,sn).groups()[0]
        if name in unique_keys.keys(): unique_keys[name]+=1
        else:  unique_keys[name]=1
    combined_SASAs={}
    if not all([val==unique_keys.values()[0] for val in unique_keys.values()]):
        print "[warning] not all SASA replicates have same number of replicates"
    for name in unique_keys.keys():
        for key,val in unpacked_SASAs.items():
            if re.search("\b{0}".format(name),"\b{0}".format(key)):
                if name in combined_SASAs.keys():
                    for resid1 in val['base_sasa'].keys():
                        for resid2 in combined_SASAs[name]['base_sasa'].keys():
                            if resid1==resid2:
                                combined_SASAs[name]['base_sasa'][resid2][
                                    'sasa_vals']=np.append(val['base_sasa'][
                                        resid1]['sasa_vals'],combined_SASAs[
                                            name]['base_sasa'][resid2]['sasa_vals'])
                    if 'comp_sasa' in combined_SASAs[name].keys():
                        for resid1 in val['comp_sasa'].keys():
                            for resid2 in combined_SASAs[name]['comp_sasa'].keys():
                                if resid1==resid2:
                                    combined_SASAs[name]['comp_sasa'][resid2][
                                        'sasa_vals']=np.append(
                                            val['comp_sasa'][resid1]['sasa_vals'],
                                            combined_SASAs[name]['comp_sasa'][resid2]['sasa_vals'])
                else: 
                    combined_SASAs[name]=val
                    combined_SASAs[name]['name']=name
    return combined_SASAs


def sasa_stats(SASAs):
    SASA_stats={}
    for sn,info in SASAs.items():
        SASA_stats[sn]={'name':' '.join(info['name'].split('_')), 'kcat':info['kcat']}
        SASA_stats[sn]['means']=[np.mean(val['sasa_vals']) for val in info['base_sasa'].values()]
        SASA_stats[sn]['stds']=[np.std(val['sasa_vals']) for val in info['base_sasa'].values()]
        SASA_stats[sn]['median']=[np.median(val['sasa_vals']) for val in info['base_sasa'].values()]
        #!!!need to do something with comp_sasa
        SASA_stats[sn]['mean']=np.sum(SASA_stats[sn]['means'])
        SASA_stats[sn]['std']=np.sum(SASA_stats[sn]['stds'])
    means=[[k,v['mean']] for k,v in SASA_stats.items()]
    max_=max([SASA_stats[key]['mean'] for key in SASA_stats.keys()])
    min_=min([SASA_stats[key]['mean'] for key in SASA_stats.keys()])
    SASA_stats['max_sasa']=max_
    SASA_stats['min_sasa']=min_
    sorted_keys=[i[0] for i in sorted(means,key=lambda x: x[1],reverse=True)]
    return SASA_stats,sorted_keys
