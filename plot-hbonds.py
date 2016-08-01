#!/usr/bin/python -i
#-*- coding: utf-8 -*-

#---plot prep
execfile('./omni/base/header.py')
from plotter import *
from base.store import plotload
from common_plot import *
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.patches as mpatches

#---settings
plotname = 'hbonds'

#---colors
import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)

#---load the upstream data
data,calc = plotload(plotname,work)

def read_hbonds(my_data,hbond_keys,timesteps=False,donor_restype=None,acceptor_restype=None,
                           donor_reslist=None,acceptor_reslist=None,divy=False):
    hbonds={}
    main_donor='N'
    main_acceptor='O'
    for sn in hbond_keys:
        if timesteps: 
            raw_hbonds=my_data[sn]['data']['hbonds_timesteps'][0]
            num_steps=len(set([i[9] for i in raw_hbonds]))
        else: raw_hbonds=my_data[sn]['data']['hbonds_occupancies']
        hbonds[sn]={}
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
                        if label in hbonds[sn].keys():
                            time=hbonds[sn][label]['times']
                            time.append(hbond[9])
                            hbonds[sn][label]['times']=time
                        else: hbonds[sn][label]={'donor_restype':hbond[2], 'donor_residx':hbond[3],
                                                 'donor_HA':hbond[4], 'acceptor_restype':hbond[6],
                                                 'acceptor_residx':hbond[7], 'acceptor_HA':hbond[8],
                                                 'times':[hbond[9]]}
                    else:
                        if label in hbonds[sn].keys():
                            hbonds[sn][label]['occupancy']+=float(hbond[9])
                        else:
                            hbonds[sn][label]={'donor_restype':hbond[2], 'donor_residx':hbond[3],
                                               'donor_HA':hbond[4], 'acceptor_restype':hbond[6],
                                               'acceptor_residx':hbond[7], 'acceptor_HA':hbond[8],
                                               'occupancy':float(hbond[9])}

        if divy: 
            for key,val in divies.items():
                infos=key.split(' ')
                if timesteps:
                    hbonds[sn][key]={'donor_restype':infos[0], 'donor_residx':infos[1],
                                     'donor_loc':infos[2], 'acceptor_restype':infos[3],
                                     'acceptor_residx':infos[4], 'acceptor_loc':infos[5],
                                     'times':val}
                else:
                    hbonds[sn][key]={'donor_restype':infos[0], 'donor_residx':infos[1],
                                     'donor_loc':infos[2], 'acceptor_restype':infos[3],
                                     'acceptor_residx':infos[4], 'acceptor_loc':infos[5],
                                     'occupancy':float(val)}

        for hb in hbonds[sn].keys():
            if timesteps:
                hbonds[sn][hb]['occupancy']=float('{:0.5f}'.format(len(
                    hbonds[sn][hb]['times'])/float(num_steps)))
            else:
                hbonds[sn][hb]['occupancy']=float('{:0.5f}'.format(hbonds[sn][hb]['occupancy']))
        if timesteps:
            hbonds[sn]['steps']=num_steps
    for name in hbonds:
        hbonds[name]['active']=work.meta[name]['active']
    return hbonds

def combine_hbonds(unpacked_hbonds,sorted_keys,num_replicates=2,timesteps=False,divy=False):
    # this will only work with a sorted list of keys that contain num_replicate number of 
    # replicates sequentially
    if len(sorted_keys)%num_replicates!=0:
        print "[Error] number of keys is not a multiple of number of replicates; returning."
        return
    combined_hbonds={}
    for i in range(len(sorted_keys)/num_replicates):
        bonds={};key_list=[]
        for j in range(num_replicates):
            name=sorted_keys.pop(0)
            key_list.append(name)
            for label in unpacked_hbonds[name]:
                if label=='active': continue
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
        combined_hbonds[name]['active']=unpacked_hbonds[name]['active']
    return combined_hbonds
        
def occupancy_diff(hbonds,reference=None,threshold=0.4):
    altered_donors={}
    bond_list=list(set([label for name in hbonds for label in hbonds[name]]))
    for sn in hbonds:
        if sn==reference: continue
        altered_donors[sn]={'delta':0, 'name':hbonds[sn]['name'], 'active':hbonds[sn]['active']}
        for bond in bond_list:
            if bond=='name' or bond=='active': continue
            if bond in hbonds[sn].keys() or bond in hbonds[reference].keys():
                if bond not in hbonds[reference].keys(): hbonds[reference][bond]={'mean':0}
                if bond not in hbonds[sn].keys(): hbonds[sn][bond]={'mean':0}
                delta_wt=float(hbonds[sn][bond]['mean'])-float(hbonds[reference][bond]['mean'])
                if abs(delta_wt)>threshold:
                    altered_donors[sn]['delta']+=delta_wt
    means=[[k,v['delta']] for k,v in altered_donors.items()]
    sorted_keys=[i[0] for i in sorted(means,key=lambda x: x[1],reverse=True)]
    return altered_donors,sorted_keys


def histofusion(deltas,keys,mode='values',title=None, plot=True, out_file=None, y_limits=False,
                ylabel='H-bonds occupancy difference',meta=None):
    sorted_names=keys
    avgs=np.array([deltas[name]['delta'] for name in sorted_names])
    #maxs=np.array([data[name]['max'] for name in sorted_names])
    #mins=np.array([data[name]['min'] for name in sorted_names])
    labels=[deltas[name]['active'] for name in sorted_names]
    mutations=[deltas[name]['name'] for name in sorted_names]
    fig, ax = plt.subplots()
    x_ticks = np.arange(len(labels))
    width=0.8
    color_list=[color_dict[i] for i in labels]
    label_list=[label_dict[i] for i in labels]
    alpha=0.6
    bar = ax.bar(x_ticks-width/2, avgs, width, color=color_list, alpha=alpha)
    #ax.scatter(x_ticks, mins, 40, color='k', marker='*')
    #ax.scatter(x_ticks, maxs, 40, color='k', marker='*')
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(mutations, rotation='vertical', ha='center',size='large')
    ax.set_xlim(-1,len(avgs))
    plt.subplots_adjust(bottom=0.2)
    plt.ylabel(ylabel,size='x-large')
    if title!=None:
        plt.title(title,size='x-large')
    plt.grid(axis='y')
    
    patches={
        True:mpatches.Patch(color=color_dict[True], label='activating', alpha=alpha),
        False:mpatches.Patch(color=color_dict[False], label='non-activating', alpha=alpha),
        'wt':mpatches.Patch(color=color_dict['wt'], label='wild type', alpha=alpha),
        'maybe':mpatches.Patch(color=color_dict['maybe'], label='unknown', alpha=alpha)}    
    used_patch=[patches[label] for label in set(labels)]
    used_label=[label_dict[label] for label in set(labels)]
    ax.legend(used_patch,used_label)

    if plot:
        plt.show(block=False)
    else:
        picturesave('fig.histogram-%s'%(plotname),work.plotdir,backup=False,
                    version=True,meta=meta,dpi=200)



def confusion(data):
    names=[name for name in data if data[name]['active'] not in ['wt','maybe']]
    values=[data[name]['delta'] for name in names]
    labels=[data[name]['active'] for name in names]
    confusion={}
    for threshold in values:
        predictions=[]
        for value in values:
            if value>=threshold:
                predictions.append(1)
            else:
                predictions.append(0)

        tp=0;fp=0;tn=0;fn=0
        for i in range(len(predictions)):
            if predictions[i]==labels[i]:
                if labels[i]==1:
                    tp+=1
                elif labels[i]==0:
                    tn+=1
            elif predictions[i]!=labels[i]:
                if labels[i]==1:
                    fn+=1
                elif labels[i]==0:
                    fp+=1
        TPR=tp/float(tp+fn)
        FPR=fp/float(fp+tn)
        #roc=skm.roc_auc_score(np.array(labels,dtype=int),np.array(predictions))
        confusion[threshold]={'tpr':TPR, 'fpr':FPR, 'tp':tp, 'fp':fp, 'fn':fn, 'tn':tn}
    return confusion
    
def parameter_sweep1D(in_data, reference='inactive_wt', limits=[0,1,6], title=None,plot=True):
    sweep_limits=np.linspace(limits[0],limits[1],limits[2])
    params=[]
    for i in sweep_limits:
        deltas,keys=occupancy_diff(in_data, reference=reference, threshold=i)
        preds=confusion(deltas)
        rate=max([preds[pred]['tp']-preds[pred]['fp'] for pred in preds])
        params.append([i,rate])
    params=np.array(params)
    x=np.transpose(params)[0]
    y=np.transpose(params)[1]
    fig, ax = plt.subplots()
    ax.plot(x,y)
    if title!=None:
        plt.title(title,size='x-large')
    plt.ylabel('(True Positive) - (False Positive)',size='x-large')
    plt.xlabel('H-bond differential threshold',size='x-large')
    if plot:
        plt.show(block=False)
    else:
        picturesave('fig.parameter-sweep-%s'%(plotname),work.plotdir,backup=False,
                    version=True,meta={},dpi=200)
    best_rate=max([i[1] for i in params])
    best_threshold=min([i[0] for i in params if i[1]==best_rate])
    return best_threshold



#finalize the data and call the plotter
protein=work.c
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit


#dhbo=divy_hbonds_occ=read_hbonds(sort_keys,divy=True,timesteps=False,donor_reslist=None,acceptor_reslist=None)
#hbo=hbonds_occ=read_hbonds(sort_keys,timesteps=False)
#hbts=hbonds_ts=read_hbonds(sort_keys,timesteps=True)
#dhbts=_divy_hbonds_ts=read_hbonds(sort_keys,divy=True,timesteps=True)
#combos=combine_hbonds(divy_hbonds_occ,sort_keys,divy=True)
#best_thresh=parameter_sweep1D(combos, reference='inactive_wt', limits=[0,1.2,49],plot=False)
best_thresh=0.75
for key,val in domains.items():
    if type(val)==int: #need a list of residx
        continue
    sort_keys=sorted(data.keys())
    my_data=data
    dhbo=read_hbonds(my_data,sort_keys,divy=True,timesteps=False,donor_reslist=val,acceptor_reslist=val)
    combos=combine_hbonds(dhbo,sort_keys,divy=True)
    deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=best_thresh)
    histofusion(deltas,keys,title=u'Threshold = {0:1.3f}{1} Investigating: {2}'.format(best_thresh,'\n',key),plot=True,meta={key:val,'threshold':best_thresh})

