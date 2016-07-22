#!/usr/bin/python -i

#---plot prep
execfile('./omni/base/header.py')
from plotter import *
from base.store import plotload
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.patches as mpatches

#---settings
plotname = 'sasa'
rmsd_bin_step = 0.1

#---colors
import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)
color_dict={True:'r', False:'g', 'maybe':'b', 'wt':'m'}

#---load the upstream data
data,calc = plotload(plotname,work)
sasa_type='abs_sasa'
label_dict={True:'activating', False:'non-activating', 'maybe':'unknown', 'wt':'wild type'}


def get_subdomains(protein,domainfile):
    with open(domainfile,'r') as fp: domain_lines=fp.readlines(); fp.close()
    for line in domain_lines:
        line.strip('\n')
        dom_info=line.split(' ')
        if dom_info[0].upper()==protein.upper():
            domains={'kd_start':int(dom_info[1]), 'kd_end':int(dom_info[2]),
                     'kd':range(int(dom_info[1]), int(dom_info[2])),
                     'ploop_start':int(dom_info[3]), 'ploop_end':int(dom_info[4]),
                     'ploop':range(int(dom_info[3]), int(dom_info[4])),
                     'alphac_start':int(dom_info[5]), 'alphac_end':int(dom_info[6]),
                     'alphac':range(int(dom_info[5]), int(dom_info[6])),
                     'catloop_start':int(dom_info[7]), 'catloop_end':int(dom_info[8]),
                     'cloop':range(int(dom_info[7]), int(dom_info[8])),
                     'activation_start':int(dom_info[9]), 'activation_end':int(dom_info[10]),
                     'aloop':range(int(dom_info[9]), int(dom_info[10]))}
            return domains
    return
domainfile='./calcs/kinase_subdomains'
protein=work.c
domains=get_subdomains(protein,domainfile)
if not domains: print "[ERROR] no subdomains found"; exit

hydrophobic=['PHE','TYR','ILE','LEU','VAL','TRP']
polar=['ARG','LYS','GLU','ASP','HIS','SER','THR','ASN','GLN']
hydrophobic_core = [773, 774, 777, 782, 835, 836, 841, 864]

def unpack_sasas(SASA_keys,sasa_type=sasa_type,base_restype=None,comp_restype=None,res_list=None):
    SASAs={}
    for sn in SASA_keys:
        raw_sasa=zip(data[sn]['data']['resname'],data[sn]['data']['resid'],
                       data[sn]['data'][sasa_type])
        SASAs[sn]={'name':' '.join(sn.split('_')), 'active':work.meta[sn]['active'],
                   'raw_sasa':raw_sasa}
        if type(res_list)!=list: res_list=data[sn]['data']['resid']
        if base_restype:
            base_sasa=[]
            for elt in raw_sasa:
                resname=elt[0];resid=elt[1];sasa=elt[2]
                for res in base_restype:
                    if resname==res: base_sasa.append(sasa)
            if not comp_restype and resid in res_list: SASA_sums=np.sum(base_sasa,axis=0)
        if base_restype and comp_restype:
            comp_sasa=[]
            for elt in raw_sasa:
                resname=elt[0];resid=elt[1];sasa=elt[2]
                for res in comp_restype:
                    if resname==res and resid in res_list: comp_sasa.append(sasa)
            SASA_sums=np.sum(base_sasa,axis=0)-np.sum(comp_sasa,axis=0)
        else: SASA_sums=np.sum([i[2] for i in raw_sasa if i[1] in res_list],axis=0)
        SASAs[sn]['sums']=np.array(SASA_sums)
        SASAs[sn]['mean']=np.mean(SASA_sums)
        SASAs[sn]['std']=np.std(SASA_sums)
        SASAs[sn]['median']=np.median(SASA_sums)
        SASAs[sn]['max']=max(SASA_sums)
        SASAs[sn]['min']=min(SASA_sums)
    max_=max([SASAs[sn]['max'] for sn in SASAs.keys()])
    min_=min([SASAs[sn]['min'] for sn in SASAs.keys()])
    SASAs['max_sasa']=max_
    SASAs['min_sasa']=min_
    return SASAs

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


def error_SASA(data,sasa_type=sasa_type,base_restype=None,comp_restype=None,res_list=None,
               combine=False,plot=True):
    #---prepare an axis
    axes,fig = panelplot(
        layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},
	figsize=(12,8))

    sort_keys=sorted(data.keys())
    sasas=unpack_sasas(sort_keys,sasa_type=sasa_type,base_restype=base_restype,
                       comp_restype=comp_restype,res_list=res_list)
    if combine: 
        sasas,sort_keys=combine_SASAs(sasas,sort_keys,num_replicates=combine)
        if not sasas: return
    #---PLOT
    counter,xpos,xlim_left = 0,[],0
    max_SASA=sasas['max_sasa']
    min_SASA=sasas['min_sasa']
    labels=[sasas[name]['active'] for name in sort_keys]
    for snum,sn in enumerate(sort_keys):
	#---unpack
        mean=sasas[sn]['mean'];std=sasas[sn]['std']
        name=sasas[sn]['name'];activity=sasas[sn]['active']
	color = color_dict[activity]
	#---boxes
	ax = axes[0]
        ax.errorbar(x=counter,y=mean,yerr=std,color=color,elinewidth=4,capthick=3,
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
    ax.legend(used_patch,used_label)
    for ax in axes: 
	ax.set_ylabel(r'relative SASA (A.U.)')
    if plot:
        fig.show()
    else: picturesave('fig.error-%s'%plotname,work.plotdir,backup=False,version=True,meta={})


error_SASA(data,sasa_type,base_restype=None,comp_restype=None,res_list=hydrophobic_core,combine=2,plot=False)

