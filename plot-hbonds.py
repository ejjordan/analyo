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
        #if timesteps: hbonds[sn]['steps']=num_steps
    for name in hbonds: hbonds[name]['active']=work.meta[name]['active']
    return hbonds

def chunk_hbonds(unpacked_hbonds,hbond_keys,num_chunks=10):
    chunk_hbonds={}
    for key in hbond_keys:
        hbonds=unpacked_hbonds[key]
        timeseries=data[key]['slice']['timeseries']
        timeseries_len=len(timeseries)
        chunk_len=timeseries_len/num_chunks
        chunk_steps=range(0,timeseries_len,chunk_len)
        chunk_hbonds[key]={'active':hbonds['active']}
        for i in range(len(chunk_steps)-1):
            start_time=timeseries[chunk_steps[i]]
            stop_time=timeseries[chunk_steps[i+1]]
            for bond in hbonds:
                if bond=='active': continue
                if bond not in chunk_hbonds[key]:
                    chunk_hbonds[key][bond]={'donor_residx':hbonds[bond]['donor_residx'],
                                             'donor_restype':hbonds[bond]['donor_restype'],
                                             'donor_HA':hbonds[bond]['donor_HA'],
                                             'acceptor_residx':hbonds[bond]['acceptor_residx'],
                                             'acceptor_restype':hbonds[bond]['acceptor_restype'],
                                             'acceptor_HA':hbonds[bond]['acceptor_HA'],
                                             'times':{}}
                times=[time for time in hbonds[bond]['times'] if 
                       float(time)>float(start_time) and float(time)<float(stop_time)]
                if times:
                    chunk_hbonds[key][bond]['times'][i]={'start_time':start_time,
                                                         'stop_time':stop_time, 'times':times,
                                                         'occupancy':len(times)/float(chunk_len)}
    return chunk_hbonds


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
        altered_donors[sn]={'delta':0, 'name':hbonds[sn]['name'], 'active':hbonds[sn]['active'],
                            'bonds':[]}
        for bond in bond_list:
            if bond=='name' or bond=='active': continue
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

def occupancy_thresholder(hbonds,threshold=0.4):
    threshold_donors={}
    for sn in hbonds:
        threshold_donors[sn]={'name':' '.join(sn.split('_')), 'active':hbonds[sn]['active']}
        for bond in hbonds[sn]:
            if bond=='name' or bond=='active' or bond=='start_time' or bond=='stop_time': continue
            occupancy=hbonds[sn][bond]['occupancy']
            if occupancy>threshold:
                threshold_donors[sn][bond]={'occupancy':occupancy}
    return threshold_donors

def occupancy_stats_thresholder(hbonds,threshold_label='mean',threshold=0.4):
    threshold_donors={}
    for sn in hbonds:
        threshold_donors[sn]={'name':' '.join(sn.split('_')), 'active':hbonds[sn]['active']}
        for bond in hbonds[sn]:
            if bond=='name' or bond=='active' or bond=='start_time' or bond=='stop_time': continue
            value=hbonds[sn][bond][threshold_label]
            if value>threshold:
                threshold_donors[sn][bond]={'mean':hbonds[sn][bond]['mean'],
                                            'std':hbonds[sn][bond]['std'],
                                            'min':hbonds[sn][bond]['min'],
                                            'max':hbonds[sn][bond]['max'],
                                            'delta':hbonds[sn][bond]['delta']}
    return threshold_donors

def occupancy_variance(hbonds,sorted_keys,num_chunks=10):
    accumulator={}
    for sn in sorted_keys:
        accumulator[sn]={}
        for chunk in range(num_chunks):
            for bond in hbonds[sn]:
                if bond=='active' or bond=='start_time' or bond=='stop_time': continue
                if chunk not in hbonds[sn][bond]['times']: continue
                if bond in accumulator[sn]:
                    accumulator[sn][bond].append(float(
                        hbonds[sn][bond]['times'][chunk]['occupancy']))
                else:
                    accumulator[sn][bond]=[float(hbonds[sn][bond]['times'][chunk]['occupancy'])]
        for bond in accumulator[sn]:
            bond_frames=len(accumulator[sn][bond])
            non_bond_frames=num_chunks-bond_frames
            bond_frames_list=np.append(np.array(accumulator[sn][bond]),[0]*non_bond_frames)
            accumulator[sn][bond]={'mean':np.mean(bond_frames_list),
                                   'std':np.std(bond_frames_list),
                                   'min':min(bond_frames_list),
                                   'max':max(bond_frames_list)}
            accumulator[sn][bond]['delta']=accumulator[sn][bond]['max']-accumulator[sn][bond]['min']
        accumulator[sn]['active']=hbonds[sn]['active']
    return accumulator

def stat_plotter(var,threshold_label='std',threshold=0.2):
    stat=occupancy_stats_thresholder(var,threshold_label,threshold)
    keys=stat.keys()
    stats=[[key,stat[key]['active'],len(stat[key].keys())-2] for key in keys]
    labels=[];label_colors=[];values=[]
    for stat in stats:
        labels.append(' '.join(stat[0].split('_')))
        label_colors.append(stat[1])
        values.append(stat[2])
    fig, ax = plt.subplots()
    x_ticks = np.arange(len(labels))
    color_list=np.array([color_dict[i] for i in label_colors])
    bar = ax.bar(x_ticks, values, color=color_list)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(labels, rotation='vertical', ha='center',size='large')
    plt.subplots_adjust(bottom=0.4)
    plt.show(block=False)


def thresh_plotter(thresh,stats=False):
    labels=[];label_colors=[];values=[];means=[];stds=[]
    for sn in thresh:
        for bond in thresh[sn]:
            if bond=='name' or bond=='active': continue
            labels.append(thresh[sn]['name']+' '+bond)
            label_colors.append(thresh[sn]['active'])
            if stats:
                means.append(thresh[sn][bond]['mean'])
                stds.append(thresh[sn][bond]['std'])
            else:
                values.append(thresh[sn][bond]['occupancy'])
    fig, ax = plt.subplots()
    x_ticks = np.arange(len(labels))
    width=0.8
    color_list=np.array([color_dict[i] for i in label_colors])
    label_list=[label_dict[i] for i in label_colors]
    alpha=0.6
    if stats:
        for i,color in enumerate(color_list):
            ax.errorbar(x=x_ticks[i]-width/2, y=means[i], yerr=stds[i], color=color,
                        elinewidth=1, capthick=1, capsize=6, fmt='ko')
    else:
        bar = ax.bar(x_ticks-width/2, values, color=color_list, alpha=alpha)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(labels, rotation='vertical', ha='center',size='large')
    plt.subplots_adjust(bottom=0.4)
    plt.show(block=False)

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


sort_keys=sorted(data.keys())
my_data=data
hbts=read_hbonds(my_data,sort_keys,timesteps=False,donor_reslist=domains['$\\alpha$C helix, activation loop'],acceptor_reslist=domains['$\\alpha$C helix, activation loop'],divy=True)
#hbts={'active_wt_replicate':hbts['active_wt_replicate']}
#sort_keys=['active_wt_replicate']
best_thresh=0.8
combos=combine_hbonds(hbts,sort_keys,divy=True,num_replicates=2)
#chunks=chunk_hbonds(hbts,sort_keys)
#del(hbts)
#var=occupancy_variance(chunks,sort_keys)
#tstats=occupancy_stats_thresholder(var,threshold_label='std',threshold=0.35)
#thresh=occupancy_thresholder(chunks,threshold=best_thresh)
#thresh_plotter(tstats,stats=True)
deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=best_thresh)
histofusion(deltas,keys,title=u'Threshold = {0:1.3f}{1} Investigating: {2}'.format(best_thresh,'\n','all'),plot=True)

def plot_num_bond_scaling(num_bond_scaling=101,bonds=combos,plot=True,xlog=True):
    num_bonds=[]
    for thresh in np.linspace(0,1,num_bond_scaling):
        deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=thresh)
        bond_list=list(set(flatten(np.array([deltas[key]['bonds'] for key in deltas]))))
        num_bonds.append([thresh,len(bond_list)])
    fig, ax = plt.subplots()
    ax.plot([i[1] for i in num_bonds], [i[0] for i in num_bonds])
    if xlog:
        ax.set_xscale('log')
    ax.invert_yaxis()
    plt.show()
    if plot:
        plt.show(block=False)
    else:
        picturesave('fig.num_bond_scaling-%s'%(plotname),work.plotdir,backup=False,
                    version=True,dpi=200)
