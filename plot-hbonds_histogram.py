#!/usr/bin/env python

#---plot prep
from common_plot import *
from common_hbonds import *
import numpy as np
import matplotlib.patches as mpatches

#---settings
plotname = 'hbonds_histogram'

#---load the upstream data
data,calc = plotload(plotname,work)

def histofusion(deltas,keys,mode='values',title=None, plot=True, y_limits=False,
                meta=None,zero=False, max_inactive=2, min_active=4,
                ylabel='H-bonds occupancy difference'):
    sorted_names=keys
    avgs=np.array([deltas[name]['delta'] for name in sorted_names])
    if zero:
        pos_avgs=[avg if avg!=0 else 0.1 for avg in avgs]
        neg_avgs=[avg if avg!=0 else -0.1 for avg in avgs]
    labels=label_maker(data=deltas, name_list=sorted_names,
                       max_inactive=max_inactive, min_active=min_active)
    mutations=[deltas[name]['name'] for name in sorted_names]
    fig, ax = plt.subplots()
    x_ticks = np.arange(len(labels))
    width=0.8
    color_list=[color_dict[i] for i in labels]
    label_list=[label_dict[i] for i in labels]
    bar = ax.bar(x_ticks, avgs, width, color=color_list)
    if zero:
        bar = ax.bar(x_ticks, pos_avgs, width, color=color_list)
        bar = ax.bar(x_ticks, neg_avgs, width, color=color_list)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(mutations, rotation='vertical', ha='center',size='large')
    fig.subplots_adjust(bottom=0.2)
    plt.ylabel(ylabel,size='x-large')
    if title!=None:
        plt.title(title,size='x-large')
    ax.grid(axis='y')

    patches={
        True:mpatches.Patch(color=color_dict[True], label='activating'),
        False:mpatches.Patch(color=color_dict[False], label='non-activating'),
        'I':mpatches.Patch(color=color_dict['I'], label='mildly activating'),
        'wt':mpatches.Patch(color=color_dict['wt'], label='wild type'),
        'maybe':mpatches.Patch(color=color_dict['maybe'], label='unknown')}
    used_patch=[patches[label] for label in set(labels)]
    used_label=[label_dict[label] for label in set(labels)]
    ax.legend(used_patch,used_label)

    if plot:
        plt.show(block=False)
    else:
        picturesave('fig.%s'%(plotname),directory=work.plotdir,backup=False,
                    version=True,meta=meta)


#finalize the data and call the plotter
protein=work.meta['protein_name']
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit

sort_keys=sorted(data.keys())
min_active=work.plots[plotname]['specs']['min_active']
max_inactive=work.plots[plotname]['specs']['max_inactive']
threshold=work.plots[plotname]['specs']['threshold']
if protein=='alk':
    for name in sort_keys:
        if work.meta[name]['kcat'] not in ['WT','X','I']:
            work.meta[name]['kcat']=float(work.meta[name]['kcat'])/9

for domain_name,domain_residxs in domains.items():
	hbts=hbonds_timesteps(work,data,sort_keys,donor_reslist=domain_residxs,
						  acceptor_reslist=domain_residxs,divy=True)
	combos=combine_hbonds(hbts,sort_keys,divy=True,num_replicates=2)
	deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=threshold)
	bond_list=[item for sublist in [deltas[key]['bonds'] for key in deltas] for item in sublist]
	bond_list=list(set(bond_list))
	title=u'inactive kcat $<=$ {0}x, active kcat $>=$ {1}x, Threshold = {2:1.3f}\nResidues: {3}'.format(max_inactive,min_active,threshold,domain_name)
	print title
	histofusion(deltas, keys, title=title, zero=True, plot=False,
                max_inactive=max_inactive, min_active=min_active,
				meta={'occupancy_diff threshold':threshold, 'domain name':domain_name,
					  'donor_residues':domain_residxs, 'acceptor_residues':domain_residxs,
					  'max inactive':max_inactive, 'min active':min_active,
                      'intradomain bonds excluded':'no'})
