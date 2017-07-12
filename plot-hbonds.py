#!/usr/bin/env python

#---plot prep
from common_plot import *
from common_hbonds import *
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



def thresh_plotter(thresh,stats=False,deltas=False,chunks=10,plot=True,title=None,meta=None,plot_threshold=None,residue_to_highlight=None):
    labels=[];label_colors=[];label_marks=[];values=[];means=[];stds=[]
    for sn in thresh:
        for bond in thresh[sn]:
            if bond=='name' or bond=='kcat': continue
            """
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
            """
            for chunk in range(chunks):
                chunk+=1
                if bond in thresh[sn] and chunk in thresh[sn][bond]['times'] and abs(thresh[sn][bond]['deltas'][chunk])>plot_threshold:
                    values.append(thresh[sn][bond]['deltas'][chunk])
                    mark='black'
                    if residue_to_highlight:
                        bond_data=bond.split()
                        if str(residue_to_highlight) in bond_data:
                            mark='red'
                    labels.append(' '.join(sn.split('_'))+', '+bond)
                    label_marks.append(mark)
                    label_colors.append(thresh[sn]['kcat'])

    label_colors=label_maker(label_colors,val_type='kcat')
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
    plt.ylabel('Change in H-bond occupancy',size='x-large')
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
    ax.legend(used_patch,used_label, loc=8, ncol=3)
    if plot:
        plt.show(block=False)
    else:
        picturesave('fig.delta_bonds-%s'%(plotname),work.plotdir,backup=False,
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
protein='alk'
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit


sort_keys=sorted(data.keys())
hbts=hbonds_timesteps(data,sort_keys,donor_reslist=domains['$\\alpha$C helix, activation loop'],acceptor_reslist=domains['$\\alpha$C helix, activation loop'],divy=True)
threshold=0.75
combos=combine_hbonds(hbts,sort_keys,divy=True,num_replicates=2)
deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=threshold)
bond_list=[item for sublist in [deltas[key]['bonds'] for key in deltas] for item in sublist]
bond_list=list(set(bond_list))
#chunks=chunk_hbonds(hbts,sort_keys,bond_list=bond_list,divy=True,num_chunks=2,deltas=True)
#var=occupancy_variance(chunks,sort_keys)
#tstats=occupancy_stats_thresholder(var,threshold_label='std',threshold=0.35)
#thresh=occupancy_thresholder(chunks,threshold=best_thresh)
#thresh_plotter(chunks,stats=False,chunks=1,deltas=True)
hili_res=868;plot_thresh=0.3
#thresh_plotter(chunks,stats=False,chunks=2,deltas=True,plot=False,plot_threshold=plot_thresh,title='$>${0} $\%$ occupancy change H-bonds'.format(plot_thresh*100),meta={'occupancy_diff threshold':threshold,'plot threshold':plot_thresh,'bond_list':bond_list,'highlighted_residue':hili_res},residue_to_highlight=None)

kcat=3;metric='ROC';param='kcat'
title=u'Threshold = {0:1.3f}\tkcat: {1}x\nResidues: {2}'.format(threshold,kcat,'$\\alpha$C helix, activation loop')
#parameter_sweep1D(combos, limits=[10,40,4],title='{0} sweep\nthreshold={1}'.format(param,threshold),meta={'parameter':param,'metric':metric,'alt_param':{'threshold':threshold}},alt_param=threshold,parameter=param,plot=False)
histofusion(deltas,keys,title=title,plot=False,kcat_cut=kcat,meta={'occupancy_diff threshold':threshold,'donor_residues':'$\\alpha$C helix, activation loop','acceptor_residues':'$\\alpha$C helix, activation loop','kcat':kcat,'intradomain bonds excluded':'no'},zero=True)
