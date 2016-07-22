#!/usr/bin/python -i
#-*- coding: utf-8 -*-

#---plot prep
execfile('./omni/base/header.py')
from plotter import *
from base.store import plotload
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.patches as mpatches

#---settings
plotname = 'hbonds'
rmsd_bin_step = 0.1

#---colors
import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)
color_dict={True:'r', False:'g', 'maybe':'b', 'wt':'m'}

#---load the upstream data
data,calc = plotload(plotname,work)
hydrophobic=['PHE','TYR','ILE','LEU','VAL','TRP']
polar=['ARG','LYS','GLU','ASP','HIS','SER','THR','ASN','GLN']
charged=['ARG','LYS','GLU','ASP','HIS']
hydrophobic_core=[1096,1098,1170,1171,1174,1179,1239,1245,1271,1240]
label_dict={True:'activating', False:'non-activating', 'maybe':'unknown', 'wt':'wild type'}

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


"""
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

    color_dict={True:'r', False:'g', 'maybe':'b', 'wt':'m'}
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
    for ax in axes: 
	ax.set_ylabel(r'relative SASA (A.U.)')
    if plot:
        fig.show()
    else: picturesave('fig.%s'%plotname,work.plotdir,backup=False,version=True,meta={})
"""

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

def get_subdomains(protein,domainfile):
    with open(domainfile,'r') as fp: domain_lines=fp.readlines(); fp.close()
    for line in domain_lines:
        line.strip('\n')
        dom_info=line.split(' ')
        if dom_info[0].upper()==protein.upper():
            domains={'kd_start':int(dom_info[1]), 'kd_end':int(dom_info[2]),
                     'kinase domain':range(int(dom_info[1]), int(dom_info[2])),
                     'ploop_start':int(dom_info[3]), 'ploop_end':int(dom_info[4]),
                     'nucleotide binding loop':range(int(dom_info[3]), int(dom_info[4])),
                     'alphac_start':int(dom_info[5]), 'alphac_end':int(dom_info[6]),
                     u"$\\alpha$C helix":range(int(dom_info[5]), int(dom_info[6])),
                     'catloop_start':int(dom_info[7]), 'catloop_end':int(dom_info[8]),
                     'catalytic loop':range(int(dom_info[7]), int(dom_info[8])),
                     'activation_start':int(dom_info[9]), 'activation_end':int(dom_info[10]),
                     'activation loop':range(int(dom_info[9]), int(dom_info[10]))}
    domains[u"$\\alpha$C helix, activation loop"]=list(set(set(domains['activation loop'])|set(domains[u"$\\alpha$C helix"])))
    domains[u"$\\alpha$C helix, nucleotide binding loop"]=list(set(set(domains['nucleotide binding loop'])|set(domains[u"$\\alpha$C helix"])))
    domains[u"$\\alpha$C helix, catalytic loop"]=list(set(set(domains['catalytic loop'])|set(domains[u"$\\alpha$C helix"])))
    domains['nucleotide binding loop, activation loop']=list(set(set(domains['activation loop'])|set(domains['nucleotide binding loop'])))
    domains['nucleotide binding loop, catalytic loop']=list(set(set(domains['nucleotide binding loop'])|set(domains['catalytic loop'])))
    domains['activation loop, catalytic loop']=list(set(set(domains['catalytic loop'])|set(domains['activation loop'])))
    domains[u"$\\alpha$C helix, catalytic loop, activation loop"]=list(set(set(domains['activation loop'])|set(domains[u"$\\alpha$C helix"])|set(domains['catalytic loop'])))
    domains[u"$\\alpha$C helix, catalytic loop, nucleotide binding loop"]=list(set(set(domains['nucleotide binding loop'])|set(domains[u"$\\alpha$C helix"])|set(domains['catalytic loop'])))
    domains[u"$\\alpha$C helix, nucleotide binding loop, activation loop"]=list(set(set(domains['activation loop'])|set(domains[u"$\\alpha$C helix"])|set(domains['nucleotide binding loop'])))
    domains[u"nucleotide binding loop, $\\alpha$C helix, catalytic loop, activation loop"]=list(set(set(domains['activation loop'])|set(domains[u"$\\alpha$C helix"])|set(domains['catalytic loop'])|set(domains['nucleotide binding loop'])))
    domains[u"kinase domain - ($\\alpha$C helix, activation loop)"]=list(set(domains['kinase domain']) - set(domains[u"$\\alpha$C helix, activation loop"]))
    return domains


#finalize the data and call the plotter
domainfile='./calcs/kinase_subdomains'
protein=work.c
domains=get_subdomains(protein,domainfile)
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
    histofusion(deltas,keys,title=u'Threshold = {0:1.3f}{1} Investigating: {2}'.format(best_thresh,'\n',key),plot=False,meta={key:val,'threshold':best_thresh})

