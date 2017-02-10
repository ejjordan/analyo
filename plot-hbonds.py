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
import json

#---settings
plotname = 'hbonds'

#---colors
import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)

#---load the upstream data
data,calc = plotload(plotname,work)

def hbonds_timesteps(read_hbonds,hbond_keys,donor_restype=None,acceptor_restype=None,
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
        

def chunk_hbonds(unpacked_hbonds,hbond_keys,num_chunks=10,bond_list=None,
                 divy=False,deltas=False):
    chunk_hbonds={}
    for key in hbond_keys:
        hbonds=unpacked_hbonds[key]
        timeseries=data[key]['slice']['timeseries']
        timeseries_len=len(timeseries)
        chunk_len=timeseries_len/num_chunks
        chunk_steps=range(0,timeseries_len,chunk_len)
        chunk_hbonds[key]={'kcat':hbonds['kcat'],'name':' '.join(key.split('_'))}
        for i in range(len(chunk_steps)-1):
            start_time=timeseries[chunk_steps[i]]
            stop_time=timeseries[chunk_steps[i+1]]
            for bond in hbonds:
                if bond_list and bond not in bond_list: continue
                if bond=='kcat': continue
                if bond not in chunk_hbonds[key]:
                    chunk_hbonds[key][bond]={'donor_residx':hbonds[bond]['donor_residx'],
                                             'acceptor_residx':hbonds[bond]['acceptor_residx'],
                                             'times':{},'deltas':{}}
                    if divy:
                        chunk_hbonds[key][bond]['donor_loc']=hbonds[bond]['donor_loc']
                        chunk_hbonds[key][bond]['acceptor_loc']=hbonds[bond]['acceptor_loc']
                    else:
                        chunk_hbonds[key][bond]['donor_restype']=hbonds[bond]['donor_restype']
                        chunk_hbonds[key][bond]['donor_HA']=hbonds[bond]['donor_HA']
                        chunk_hbonds[key][bond]['acceptor_restype']=hbonds[bond]['acceptor_restype']
                        chunk_hbonds[key][bond]['acceptor_HA']=hbonds[bond]['acceptor_HA'],
                times=[time for time in hbonds[bond]['times'] if 
                       float(time)>float(start_time) and float(time)<float(stop_time)]
                if times:
                    chunk_hbonds[key][bond]['times'][i]={'start_time':start_time,
                                                         'stop_time':stop_time, 'times':times,
                                                         'occupancy':len(times)/float(chunk_len)}
                if deltas and i>0:
                    one=0;two=0
                    if i in chunk_hbonds[key][bond]['times']:
                        one=chunk_hbonds[key][bond]['times'][i]['occupancy']
                    if i-1 in chunk_hbonds[key][bond]['times']:
                        two=chunk_hbonds[key][bond]['times'][i-1]['occupancy']
                    chunk_hbonds[key][bond]['deltas'][i]=one-two

    return chunk_hbonds


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

def occupancy_thresholder(hbonds,threshold=0.4):
    threshold_donors={}
    for sn in hbonds:
        threshold_donors[sn]={'name':' '.join(sn.split('_')), 'kcat':hbonds[sn]['kcat']}
        for bond in hbonds[sn]:
            if bond=='name' or bond=='kcat' or bond=='start_time' or bond=='stop_time': continue
            occupancy=hbonds[sn][bond]['occupancy']
            if occupancy>threshold:
                threshold_donors[sn][bond]={'occupancy':occupancy}
    return threshold_donors

def occupancy_stats_thresholder(hbonds,threshold_label='mean',threshold=0.4):
    threshold_donors={}
    for sn in hbonds:
        threshold_donors[sn]={'name':' '.join(sn.split('_')), 'kcat':hbonds[sn]['kcat']}
        for bond in hbonds[sn]:
            if bond=='name' or bond=='kcat' or bond=='start_time' or bond=='stop_time': continue
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
                if bond=='kcat' or bond=='start_time' or bond=='stop_time' or bond=='name':
                    continue
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
        accumulator[sn]['kcat']=hbonds[sn]['kcat']
    return accumulator

def stat_plotter(var,threshold_label='std',threshold=0.2):
    stat=occupancy_stats_thresholder(var,threshold_label,threshold)
    keys=stat.keys()
    stats=[[key,stat[key]['kcat'],len(stat[key].keys())-2] for key in keys]
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


"""
def thresh_plotter(thresh,stats=False,deltas=False,chunks=10,plot=True,title=None,meta=None,plot_threshold=None,residue_to_highlight=None):
    labels=[];label_colors=[];label_marks=[];values=[];means=[];stds=[]
    for sn in thresh:
        for bond in thresh[sn]:
            if bond=='name' or bond=='kcat': continue
            if stats:
                labels.append(thresh[sn]['name']+' '+bond)
                means.append(thresh[sn][bond]['mean'])
                stds.append(thresh[sn][bond]['std'])
                label_colors.append(thresh[sn]['active'])
            elif deltas:
                for chunk in range(chunks):
                    chunk+=1
                    if bond in thresh[sn] and chunk in thresh[sn][bond]['times'] and abs(thresh[sn][bond]['deltas'][chunk])>plot_threshold:
                        values.append(thresh[sn][bond]['deltas'][chunk])
                        mark='black'
                        if residue_to_highlight:
                            bond_data=bond.split()
                            if str(residue_to_highlight) in bond_data:
                                mark='red'
                        labels.append(' '.join(sn.split('_'))+' '+bond+' '+str(chunk))
                        label_marks.append(mark)
                        label_colors.append(thresh[sn]['active'])
            else:
                for chunk in range(chunks):
                    if bond in thresh[sn] and chunk in thresh[sn][bond]['times']:
                        values.append(thresh[sn][bond]['times'][chunk]['occupancy'])
                    else: values.append(0)
                    labels.append(' '.join(sn.split('_'))+' '+bond+' '+str(chunk+1))
                    label_colors.append(thresh[sn]['active'])
    fig, ax = plt.subplots(figsize=(18,14))
    x_ticks = np.arange(len(labels))
    width=1
    color_list=np.array([color_dict[i] for i in label_colors])
    label_list=[label_dict[i] for i in label_colors]
    alpha=0.6
    if stats:
        for i,color in enumerate(color_list):
            ax.errorbar(x=x_ticks[i]-width/2, y=means[i], yerr=stds[i], color=color,
                        elinewidth=1, capthick=1, capsize=6, fmt='ko')
    else:
        bar = ax.bar(x_ticks+width, values, color=color_list, alpha=alpha)
    ax.set_xticks(x_ticks+1.5*width)
    ax.set_xticklabels(labels, rotation='vertical', ha='center',size='large')
    plt.subplots_adjust(bottom=0.4)
    if residue_to_highlight:
        [label_text.set_color(label_marks[i]) for i,label_text in enumerate(plt.gca().get_xticklabels())]
        print "number of bonds: {0}\nnumber of bonds highlighted: {1}".format(len(label_marks),len([i for i in label_marks if i=='red']))
    if title!=None:
        plt.title(title,size='x-large')
    patches={
        True:mpatches.Patch(color=color_dict[True], label='activating', alpha=alpha),
        False:mpatches.Patch(color=color_dict[False], label='non-activating', alpha=alpha),
        'wt':mpatches.Patch(color=color_dict['wt'], label='wild type', alpha=alpha),
        'maybe':mpatches.Patch(color=color_dict['maybe'], label='unknown', alpha=alpha)}    
    used_patch=[patches[label] for label in set(label_colors)]
    used_label=[label_dict[label] for label in set(label_colors)]
    #ax.legend(used_patch,used_label)
    ax.legend(used_patch,used_label, loc=8, ncol=3)
    if plot:
        plt.show(block=False)
    else:
        picturesave('fig.delta_bonds-%s'%(plotname),work.plotdir,backup=False,
                    version=True,meta=meta,dpi=200)
"""

def label_maker(deltas, kcat_cut=33, name_list=None):

    """
    This function takes a 'delta' object and a kcat cut-off and returns activation (or 
    non-activation) labels. If supplied a name_list the labels will be returned in the
    specified order.
    """
    
    if not name_list: name_list=[deltas[sn]['name'] for sn in deltas]
    kcats=[deltas[name]['kcat'] for name in name_list]
    labels=[]
    for kcat in kcats:
        if kcat=='X': labels.append('maybe')
        elif kcat=='WT': labels.append('wt')
        elif float(kcat)>kcat_cut:
            labels.append(True)
        else:
            labels.append(False)
    return labels

def histofusion(deltas,keys,mode='values',title=None, plot=True, out_file=None, y_limits=False,
                ylabel='H-bonds occupancy difference', meta=None, kcat_cut=20,zero=False):
    sorted_names=keys
    avgs=np.array([deltas[name]['delta'] for name in sorted_names])
    if zero: avgs=[avg if avg!=0 else 0.1 for avg in avgs]
    #maxs=np.array([data[name]['max'] for name in sorted_names])
    #mins=np.array([data[name]['min'] for name in sorted_names])
    labels=label_maker(deltas, kcat_cut, name_list=sorted_names)
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



def confusion(data, kcat_cut):
    names=[name for name in data if data[name]['kcat'] not in ['WT','X']]
    values=[data[name]['delta'] for name in names]
    labels=label_maker(data, kcat_cut,name_list=names)

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
        confusion[threshold]={'tpr':TPR, 'fpr':FPR, 'tp':tp, 'fp':fp, 'fn':fn, 'tn':tn}
    return confusion

    
def parameter_sweep1D(in_data, reference='inactive_wt', limits=[0,1,6], title=None, plot=True,
                      parameter='threshold', metric='ROC', meta=None, alt_param=None):
    sweep_limits=np.linspace(limits[0],limits[1],limits[2])
    fig, ax = plt.subplots()
    figdat=[]
    for i in sweep_limits:
        if parameter=='threshold':
            kcat=alt_param;thresh=i
        if parameter=='kcat':
            kcat=i;thresh=alt_param
        deltas,keys=occupancy_diff(in_data, reference=reference, threshold=thresh)
        preds=confusion(deltas,kcat_cut=kcat)
        #histofusion(deltas,keys,plot=True,kcat_cut=kcat,title="kcat {0} fp {1} tp {2} fn {3} tn {4}".format(int(kcat),preds['fp'],preds['tp'],preds['fn'],preds['tn']))
        if metric=='ROC':
            x,y=[],[]
            sort_preds_keys=sorted(preds.keys(),reverse=True)
            for pred in sort_preds_keys:
                x.append(preds[pred]['fpr'])
                y.append(preds[pred]['tpr'])
            figdat.append([x,y,i])
            plt.plot(x,y,label=i)
            plt.ylabel('True Positive Rate',size='x-large')
            plt.xlabel('False Positive Rate',size='x-large')
            ax.set_xlim(0,1.01)
            ax.set_ylim(0,1.01)
            ax.plot([0,1.1],[0,1.1],'--')
        if metric=='tp-fp':
            y.append(preds['tp']-preds['fp'])
            x.append(i)
            plt.ylabel('(True Positive) - (False Positive)',size='x-large')
            if parameter=='kcat':
                plt.xlabel('kcat',size='x-large')
            if parameter=='threshold':
                plt.xlabel('H-bond differential threshold',size='x-large')

    plt.legend(loc='best')
    if title!=None:
        plt.title(title,size='x-large')
    if plot:
        plt.show(block=False)
    else:
        picturesave('fig.parameter-sweep-%s'%(plotname),work.plotdir,backup=False,
                    version=True,meta=meta,dpi=200)
    #best_rate=max([i[1] for i in params])
    #best_threshold=min([i[0] for i in params if i[1]==best_rate])
    #return best_threshold


#finalize the data and call the plotter
protein=work.c
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit


sort_keys=sorted(data.keys())
hbts=hbonds_timesteps(data,sort_keys,donor_reslist=domains['$\\alpha$C helix, activation loop'],acceptor_reslist=domains['$\\alpha$C helix, activation loop'],divy=True)
#intradomain_bonds_discard(hbts,domains['$\\alpha$C helix'],print_bonds=False)
#intradomain_bonds_discard(hbts,domains['activation loop'],print_bonds=False)
#hbts=hbonds_timesteps(data,sort_keys,donor_reslist=domains['activation loop'],acceptor_reslist=domains['$\\alpha$C helix'],divy=True,reslist_XOR=True)
#timesteps_discard(hbts,start_time=51000,stop_time=101000)
threshold=0.75
combos=combine_hbonds(hbts,sort_keys,divy=True,num_replicates=2)
deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=threshold)
bond_list=list(set(flatten([deltas[key]['bonds'] for key in deltas])))
#chunks=chunk_hbonds(hbts,sort_keys,bond_list=bond_list,divy=True,num_chunks=2,deltas=True)
#var=occupancy_variance(chunks,sort_keys)
#tstats=occupancy_stats_thresholder(var,threshold_label='std',threshold=0.35)
#thresh=occupancy_thresholder(chunks,threshold=best_thresh)
#thresh_plotter(chunks,stats=False,chunks=1,deltas=True)
hili_res=1284
#thresh_plotter(chunks,stats=False,chunks=2,deltas=True,plot=True,plot_threshold=0.4,title='Significantly altered H-bonds',meta={'occupancy_diff threshold':threshold,'plot threshold':0.4,'bond_list':bond_list,'highlighted_residue':hili_res},residue_to_highlight=hili_res)

kcat=30;metric='ROC';param='kcat'
title=u'Threshold = {0:1.3f}\tkcat: {1}\nResidues: {2}'.format(threshold,kcat,'$\\alpha$C helix, activation loop')
#parameter_sweep1D(combos, limits=[10,40,4],title='{0} sweep\nthreshold={1}'.format(param,threshold),meta={'parameter':param,'metric':metric,'alt_param':{'threshold':threshold}},alt_param=threshold,parameter=param,plot=False)
histofusion(deltas,keys,title=title,plot=False,kcat_cut=kcat,meta={'occupancy_diff threshold':threshold,'donor_residues':'$\\alpha$C helix, activation loop','acceptor_residues':'$\\alpha$C helix, activation loop','kcat':kcat,'intradomain bonds excluded':'no'})
