#!/usr/bin/python -i

#---plot prep
execfile('./omni/base/header.py')
from plotter import *
from base.store import plotload
from common_plot import *
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.patches as mpatches
import re

#---settings
plotname = 'sasa'

#---colors
import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)

#---load the upstream data
data,calc = plotload(plotname,work)
sasa_type='abs_sasa'

protein=work.c
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit

def filter_sasas(SASA_keys,sasa_type=sasa_type,base_restype=None,comp_restype=None,res_list=None):

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

def chunk_sasa(unpacked_SASAs,chunk_length=500,sasa_length=5001):
    chunked_SASAs={}
    chunks=range(0,sasa_length,chunk_length)
    for name,info in unpacked_SASAs.items():
        for i in range(len(chunks)-1):
            chunk_name='{0}_{1}'.format(name,i)
            chunked_SASAs[chunk_name]={'active':info['active'],'base_sasa':{},'name':chunk_name}
            for resid,resinfo in info['base_sasa'].items():
                chunk=resinfo['sasa_vals'][chunks[i]:chunks[i+1]]
                chunked_SASAs[chunk_name]['base_sasa'][resid]={'sasa_vals':chunk,'resid':resid,
                                                               'resname':resinfo['resname']}
            if 'comp_sasa' in info.keys():
                chunked_SASAs[chunk_name]['comp_sasa']={}
                for resid,resinfo in info['comp_sasa'].items():
                    chunk=resinfo['sasa_vals'][chunks[i]:chunks[i+1]]
                    chunked_SASAs[chunk_name]['comp_sasa'][resid]={'sasa_vals':chunk,'resid':resid,
                                                                   'resname':resinfo['resname']}
    return chunked_SASAs


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
    #max_=max([SASA_stats[key]['mean'] for key in SASA_stats.keys()])
    #min_=min([SASA_stats[key]['mean'] for key in SASA_stats.keys()])
    #SASA_stats['max_sasa']=max_
    #SASA_stats['min_sasa']=min_
    means=[[k,v['mean']] for k,v in SASA_stats.items()]
    sorted_keys=[i[0] for i in sorted(means,key=lambda x: x[1],reverse=True)]
    return SASA_stats,sorted_keys

def combine_SASAs(unpacked_SASAs,sorted_keys,num_replicates=2):
    # this will only work with a sorted list of keys that contain num_replicate number of 
    # replicates sequentially
    if len(sorted_keys)%num_replicates!=0:
        print "[Error] number of keys is not a multiple of number of replicates; returning."
        return
    combined_SASAs={}
    for i in range(len(sorted_keys)/num_replicates):
        sums=[];key_list=[]
        for j in range(num_replicates):
            name=sorted_keys.pop(0)
            key_list.append(name)
            sums.append(unpacked_SASAs[name]['sums'])
        name=min(key_list,key=len)
        sums=np.array(sums)
        sums=np.concatenate(sums)
        active=unpacked_SASAs[name]['active']
        combined_SASAs[name]={'name':' '.join(name.split('_')), 'sums':sums, 'mean':np.mean(sums),
                              'std':np.std(sums), 'median':np.median(sums),
                              'min':min(sums), 'max':max(sums), 'active':active}
    means=[[k,v['mean']] for k,v in combined_SASAs.items()]
    new_sorted_keys=[i[0] for i in sorted(means,key=lambda x: x[1],reverse=True)]
    combined_SASAs['max_sasa']=unpacked_SASAs['max_sasa']
    combined_SASAs['min_sasa']=unpacked_SASAs['min_sasa']
    
    return combined_SASAs,new_sorted_keys
    

def each_SASA(data,sasa_type=sasa_type,base_restype=polar,comp_restype=hydrophobic,res_list=None):
    sort_keys=sorted(data.keys())
    sasas=unpack_sasas(sort_keys,sasa_type=sasa_type,base_restype=comp_restype,
                       comp_restype=comp_res_type,res_list=res_list)
    num_sims=len(sort_keys)
    base_size = 20.
    wide_factor = 1.5
    color_dict={True:'r', False:'g', 'maybe':'b', 'wt':'m'}
    ncols = int(np.ceil(np.sqrt(num_sims)))
    nrows = int(np.ceil(float(num_sims)/ncols))
    fig = plt.figure(figsize=(base_size,base_size*(float(nrows)/ncols)/wide_factor))
    gs = gridspec.GridSpec(nrows,ncols,hspace=0.65,wspace=0.8)
    axes = [plt.subplot(gs[plot_num/ncols,plot_num%ncols]) for plot_num in range(num_sims)]
    for plot_num,ax in enumerate(axes):
        SASA=sasas[sort_keys[plot_num]]
        name=SASA['name'];activity=SASA['active'];SASA_sums=SASA['sums']
        ts = np.array(range(len(SASA_sums)))
        ax.plot(ts,SASA_sums,color=color_dict[activity])
        ax.set_title(name)
        ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
        ax.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='on')
    max_SASA=sasas['max_sasa']
    min_SASA=sasas['min_sasa']
    for plot_num,ax in enumerate(axes):
        ax.set_ylim(min_SASA,max_SASA)
    plt.show()


def stack_SASA(data,sasa_type=sasa_type,base_restype=None,comp_restype=None,res_list=None):
    #---prepare an axis
    axes,fig = panelplot(
        layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},
	figsize=(12,8))

    color_dict={True:'r', False:'g', 'maybe':'b', 'wt':'m'}
    sort_keys=sorted(data.keys())
    sasas=unpack_sasas(sort_keys,sasa_type=sasa_type,base_restype=base_restype,
                       comp_restype=comp_restype,res_list=res_list)

    #---PLOT
    counter,xpos,xlim_left = 0,[],0
    max_SASA=sasas['max_sasa']
    min_SASA=sasas['min_sasa']
    SASA_bins = np.arange(min_SASA,max_SASA,(max_SASA-min_SASA)/100)
    for snum,sn in enumerate(sort_keys):
	#---unpack
        SASAs=sasas[sn]['sums'];name=sasas[sn]['name'];activity=sasas[sn]['active']
	color = color_dict[activity]
	ts = np.array(range(len(SASAs)))
	ts -= ts.min()
	#---histograms
	ax = axes[0]
	counts,bins = np.histogram(SASAs,bins=SASA_bins,normed=True)
	ax.fill_betweenx((bins[1:]+bins[:-1])/2.,counter,counter+counts,alpha=1.0,
                         color=color,lw=0)
	max_SASA = max([max_SASA,max(SASAs)])
	#---compute drift
	m,b = np.polyfit(ts[len(SASAs)/10:],SASAs[len(SASAs)/10:],1)
	drift = m*(ts[-1]-ts[len(SASAs)/10])
	drift_start = ts[len(SASAs)/10]*m+b
	#---drift arrows
	ax.arrow(counter,drift_start,0,drift,
                 head_width=max(counts)*0.1,head_length=max_SASA*0.05,fc=color,ec='w',lw=1.5,
                 path_effects=[path_effects.Stroke(linewidth=4,foreground=color),
                               path_effects.Normal()],zorder=3)
	xpos.append(counter)
	counter += max(counts)*1.1
	if snum == 0: xlim_left = -1*max(counts)*0.1
    ax.set_xticks(xpos)
    ax.set_xticklabels([sasas[sn]['name'] for sn in sort_keys],rotation=45,ha='right')
    ax.set_xlim(xlim_left,counter)
    for ax in axes: 
	ax.set_ylabel(r'relative SASA (A.U.)')
    fig.show()
    #picturesave('fig.%s'%plotname,work.plotdir,backup=False,version=True,meta={})


def error_SASA(data,sort_keys=None,sasa_type=sasa_type,plot=True,meta=None,title=None,kcat=30):
    #---prepare an axis
    axes,fig = panelplot(
        layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},
	figsize=(12,8))

    sasas=data
    if not sort_keys: sort_keys=data.keys()

    #---PLOT
    counter,xpos,xlim_left = 0,[],0
    #labels=[sasas[name]['active'] for name in sort_keys]
    labels=label_maker(sasas,kcat_cut=kcat,name_list=sort_keys)
    #import pdb;pdb.set_trace()
    for snum,sn in enumerate(sort_keys):
	#---unpack
        mean=sasas[sn]['mean'];std=sasas[sn]['std']
        name=sasas[sn]['name'];activity=labels[snum]
        kcat=sasas[sn]['kcat']
	color = color_dict[activity]
	#---boxes
	ax = axes[0]
        ax.errorbar(x=counter,y=mean,yerr=std,color=color,elinewidth=4,capthick=4,
                    capsize=6,fmt='ko')
	xpos.append(counter)
        counter+=1
    ax.set_xticks(xpos)
    ax.set_xticklabels([sasas[sn]['name'] for sn in sort_keys],rotation=45,ha='right')
    ax.set_xlim(xlim_left-1,counter)
    alpha=1
    patches={
        True:mpatches.Patch(color=color_dict[True], label='activating', alpha=alpha),
        False:mpatches.Patch(color=color_dict[False], label='non-activating', alpha=alpha),
        'wt':mpatches.Patch(color=color_dict['wt'], label='wild type', alpha=alpha),
        'maybe':mpatches.Patch(color=color_dict['maybe'], label='unknown', alpha=alpha)}    
    used_patch=[patches[label] for label in set(labels)]
    used_label=[label_dict[label] for label in set(labels)]
    if title!=None:
        plt.title(title,size='x-large')
    ax.legend(used_patch,used_label)
    for ax in axes: 
	ax.set_ylabel(r'SASA (A.U.)')
    if plot:
        fig.show()
    else: picturesave('fig.error-%s'%plotname,work.plotdir,backup=False,version=True,meta=meta)

def label_maker(sasas, kcat_cut=33, name_list=None):

    """
    This function takes a 'sasa' object and a kcat cut-off and returns activation (or 
    non-activation) labels. If supplied a name_list the labels will be returned in the
    specified order.
    """
    
    if not name_list: name_list=[sasas[sn]['name'] for sn in sasas]
    kcats=[sasas[name]['kcat'] for name in name_list]
    labels=[]
    for kcat in kcats:
        if kcat=='X': labels.append('maybe')
        elif kcat=='WT': labels.append('wt')
        elif float(kcat)>kcat_cut:
            labels.append(True)
        else:
            labels.append(False)
    return labels


sasas=filter_sasas(data.keys(),sasa_type=sasa_type,base_restype=hydrophobic,comp_restype=None,res_list=hydrophobic_core[protein])
combos=combine_replicates(sasas)
#keys=sorted(sasas.keys())
#combos,new_keys=combine_SASAs(sasas,keys)
#chunks=chunk_sasa(sasas)
stats,keys=sasa_stats(combos)
kcat=30
error_SASA(stats,sort_keys=keys,sasa_type=sasa_type,plot=False,title='kcat {0}'.format(kcat),meta={'kcat':kcat},kcat=kcat)
