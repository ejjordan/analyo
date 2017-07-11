#!/usr/bin/env python

#---plot prep
from common_plot import *
from common_sasa import *
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.patches as mpatches
import re,math

#---settings
plotname = 'sasa'

#---colors
import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)

#---load the upstream data
data,calc = plotload(plotname,work)


protein='alk'#work.c
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit


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

def each_SASA(sasas,sort_keys,kcat_cut=30,plot=True,meta=None):
    num_sims=len(sort_keys)
    labels=label_maker(sasas,kcat_cut=kcat_cut,name_list=sort_keys)
    base_size = 20.
    wide_factor = 1.5
    color_dict={True:'r', False:'g', 'maybe':'b', 'wt':'m'}
    ncols = int(np.ceil(np.sqrt(num_sims)))
    nrows = int(np.ceil(float(num_sims)/ncols))
    fig = plt.figure(figsize=(base_size,base_size*(float(nrows)/ncols)/wide_factor))
    gs = gridspec.GridSpec(nrows,ncols,hspace=0.65,wspace=0.8)
    axes = [plt.subplot(gs[plot_num/ncols,plot_num%ncols]) for plot_num in range(num_sims)]
    max_SASA=0;ts_scaling=0.02
    for plot_num,ax in enumerate(axes):
        SASA=sasas[sort_keys[plot_num]]
        ts_sasa=np.sum([SASA['base_sasa'][res]['sasa_vals'] for res in SASA['base_sasa']],axis=0)
        name=SASA['name'];activity=labels[plot_num]
        ts = np.array(range(len(ts_sasa)))*ts_scaling
        ax.plot(ts,ts_sasa,color=color_dict[activity])
        ax.set_title(name)
        ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
        ax.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='on')
        max_SASA=max(max_SASA,max(ts_sasa))

    min_SASA=0
    if meta:
        meta['kcat cut']=kcat_cut
        meta['max sasa']=max_sasa
        meta['ts scaling']=ts_scaling
    else: meta={'kcat cut':kcat_cut,'max sasa':max_SASA,'ts scaling':ts_scaling}
    for plot_num,ax in enumerate(axes):
        ax.set_ylim(min_SASA,max_SASA)
    if plot:
        plt.show(block=False)
    else: picturesave('fig.each-%s'%plotname,work.plotdir,backup=False,version=True,meta=meta)

def one_SASA(SASA,cut_thresh=25,meta=None,plot=True):
    ax=plt.subplot()
    for res in SASA['base_sasa']:
        res_vals=SASA['base_sasa'][res]['sasa_vals']
        res_type=SASA['base_sasa'][res]['resname']
        ts = np.array(range(len(res_vals)))
        max_sasa=max(res_vals)
        if max_sasa>cut_thresh:
            ax.plot(ts,res_vals,label="{0}-{1}".format(res_type,res))
    ax.legend()
    plt.title(SASA['name'],size='x-large')
    if meta:
        meta['cut thresh']=cut_thresh
        meta['name']=SASA['name']
    else: meta={'cut thresh':cut_thresh,'name':SASA['name']}
    if plot:
        plt.show(block=False)
    else: picturesave('fig.single-%s'%plotname,work.plotdir,backup=False,version=True,meta=meta)



def stack_SASA(data,sasa_type='abs_sasa',base_restype=None,comp_restype=None,res_list=None):
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
    ax.set_ylabel(r'relative SASA (A.U.)')
    fig.show()
    #picturesave('fig.%s'%plotname,work.plotdir,backup=False,version=True,meta={})


def error_SASA(data,sort_keys=None,sasa_type='abs_sasa',plot=True,meta=None,title=None,kcat=30):
    #---prepare an axis
    axes,fig = panelplot(
        layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},
    figsize=(12,8))

    sasas=data
    if not sort_keys: sort_keys=data.keys()

    #---PLOT
    counter,xpos,xlim_left = 0,[],0
    labels=label_maker(sasas,kcat_cut=kcat,name_list=sort_keys)
    for snum,sn in enumerate(sort_keys):
        mean=sasas[sn]['mean']
        std=sasas[sn]['std']
        name=sasas[sn]['name']
        activity=labels[snum]
        kcat=sasas[sn]['kcat']
        color = color_dict[activity]
        #---boxes
        ax = axes[0]
        ax.errorbar(x=counter,y=mean,yerr=std,ecolor=color,elinewidth=4,capthick=4,
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
    ax.set_ylabel('SASA (\AA$^2$)')
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
        elif float(kcat)>=kcat_cut:
            labels.append(True)
        else:
            labels.append(False)
    return labels

sasa_type='abs_sasa'
#res_list=[1169, 1179, 1185, 1226, 1239, 1366]
#[1147, 1162, 1233, 1246, 1291, 1293, 1301, 1366, 1376]
#[1155, 1171, 1245, 1278, 1291, 1293, 1295, 1301, 1313, 1325, 1358, 1366]
#[1170, 1240, 1265, 1277, 1278, 1295, 1320, 1327, 1376, 1385]
#[1187, 1190, 1233, 1268, 1282, 1305, 1306, 1322, 1338, 1376, 1383]
#[1246,1366,1224,1190,1383,1282,1152]
#[1118, 1185, 1187, 1193, 1224, 1278, 1283, 1295, 1325, 1327, 1339, 1378]
#[1278,1162,1378,1339,1323,1295,1118,1283]#Rspine[protein]
sasas=filter_sasas(data,work,data.keys(),base_restype=hydrophobic,comp_restype=None,res_list=res_list)
combined='no'
combos=combine_replicates(sasas);combined='yes'
keys=sorted(sasas.keys())
#chunks=chunk_sasa(sasas)
stats,keys=sasa_stats(combos)
kcat=30;title='kcat {0}x\n{1}'.format(kcat,'regulatory spine')
error_SASA(stats,sort_keys=keys,sasa_type=sasa_type,plot=False,title=title,meta={'kcat':kcat,'residues':res_list,'combined':'yes'},kcat=kcat)
#each_SASA(sasas,keys,plot=False)
#for key in keys:
#    one_SASA(sasas[key],plot=False)
