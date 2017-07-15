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
                ylabel='H-bonds occupancy difference', meta=None, kcat_cut=20,zero=False):
    sorted_names=keys
    avgs=np.array([deltas[name]['delta'] for name in sorted_names])
    if zero:
        pos_avgs=[avg if avg!=0 else 0.1 for avg in avgs]
        neg_avgs=[avg if avg!=0 else -0.1 for avg in avgs]
    labels=label_maker(deltas, kcat_cut, name_list=sorted_names)
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
kcat=work.plots[plotname]['specs']['kcat_cut']
threshold=work.plots[plotname]['specs']['threshold']

for domain_name,domain_residxs in domains.items():
	hbts=hbonds_timesteps(work,data,sort_keys,donor_reslist=domain_residxs,
						  acceptor_reslist=domain_residxs,divy=True)
	combos=combine_hbonds(hbts,sort_keys,divy=True,num_replicates=2)
	deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=threshold)
	bond_list=[item for sublist in [deltas[key]['bonds'] for key in deltas] for item in sublist]
	bond_list=list(set(bond_list))
	title=u'Threshold = {0:1.3f}\tkcat: {1}x\nResidues: {2}'.format(threshold,4,domain_name)
	print title
	histofusion(deltas,keys,title=title,kcat_cut=kcat,zero=True,plot=False,
				meta={'occupancy_diff threshold':threshold,'domain name':domain_name,
					  'donor_residues':domain_residxs,'acceptor_residues':domain_residxs,
					  'kcat cutoff':kcat,'intradomain bonds excluded':'no'})
