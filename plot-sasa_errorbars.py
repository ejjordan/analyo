#!/usr/bin/env python

#---plot prep
from common_plot import *
from common_sasa import *
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.patches as mpatches

#---settings
plotname = 'sasa_errorbars'

#---load the upstream data
data,calc = plotload(plotname,work)

def error_SASA(data,sort_keys=None,plot=True,meta=None,title=None,
			   max_inactive=2, min_active=4):

	#get the data ready
	sasas=data
	if not sort_keys: sort_keys=data.keys()
	labels=label_maker(data,name_list=sort_keys,max_inactive=max_inactive, min_active=min_active)
	mutations=[sasas[sn]['name'].upper() for sn in sort_keys]

	xaxis_size=max(len(mutations)/3.,8) #scale x axis width by number of mutations
	#---prepare an axis
	axes,fig = panelplot(
		layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},
		figsize=(xaxis_size,8))

	#---PLOT
	counter,x_ticks,lower_yvals = 0,[],[]
	for snum,sn in enumerate(sort_keys):
		mean=sasas[sn]['mean']
		std=sasas[sn]['std']
		name=sasas[sn]['name']
		activity=labels[snum]
		kcat=sasas[sn]['kcat']
		color = color_dict[activity]
		#---boxes
		ax = axes[0]
		ax.errorbar(x=counter,y=mean,yerr=std,ecolor=color,elinewidth=8,capthick=6,
					capsize=6,fmt='wo')
		x_ticks.append(counter)
		lower_yvals.append(mean-std)
		counter+=1

	ax.set_xticks(x_ticks)
	for ii,mut in enumerate(mutations):
		if mut=='ACTIVE WT':mutations[ii]='ACTIVE'
		if mut=='INACTIVE WT':mutations[ii]='INACTIVE'
	ax.set_xticklabels(mutations,rotation='vertical',ha='center',size='xx-large')
	ymin,ymax=ax.get_ylim()
	plot_size=ymax-ymin;buf_size=0.005*plot_size
	xtick_lines=[[[xval,xval],[ymin+buf_size,ymax-buf_size]] for xval,ymax in zip(
		x_ticks,lower_yvals)]
	for line in xtick_lines:
		ax.plot(line[0],line[1],color='k',linestyle='-.',zorder=1.5, linewidth=1)
	ax.set_ylim(bottom=ymin)
	ax.set_xlim(x_ticks[0]-1,x_ticks[-1]+1)
	fig.subplots_adjust(bottom=0.2)
	alpha=1
	patches={
		True:mpatches.Patch(color=color_dict[True], label='activating', alpha=alpha),
		False:mpatches.Patch(color=color_dict[False], label='non-activating', alpha=alpha),
		'I':mpatches.Patch(color=color_dict['I'], label='mildly activating'),
		'wt':mpatches.Patch(color=color_dict['wt'], label='wild type', alpha=alpha),
		'maybe':mpatches.Patch(color=color_dict['maybe'], label='unknown', alpha=alpha)}    
	used_patch=[patches[label] for label in set(labels)]
	used_label=[label_dict[label] for label in set(labels)]
	if title!=None:
		plt.title(title,size='xx-large')
	ax.legend(used_patch,used_label,fontsize=14)
	ax.set_ylabel('SASA (\AA$^2$)')
	if plot:
		fig.show()
	else: picturesave('fig.%s'%plotname,work.plotdir,backup=False,version=True,meta=meta)


#finalize the data and call the plotter
protein=work.meta['protein_name']
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit

sort_keys=sorted(data.keys())
min_active=work.plots[plotname]['specs']['min_active']
max_inactive=work.plots[plotname]['specs']['max_inactive']
if protein=='alk':
	for name in sort_keys:
		if work.meta[name]['kcat'] not in ['WT','X','I']:
			work.meta[name]['kcat']=float(work.meta[name]['kcat'])/9
res_list=Rspine[protein]
sasas=filter_sasas(data,work,data.keys(),base_restype=None,comp_restype=None,res_list=res_list)
combined='no'
combos=combine_replicates(sasas);combined='yes'
keys=sorted(sasas.keys())
ordering='std'
stats,keys=sasa_stats(combos,order=ordering)
title=u'inactive kcat $<=$ {0}x, active kcat $>=$ {1}x\nResidues: {2}'.format(max_inactive,min_active,"regulatory spine")
error_SASA(stats,sort_keys=keys,title=title,plot=False,
max_inactive=max_inactive, min_active=min_active,
		   meta={'max inactive':max_inactive, 'min active':min_active,'order':ordering,
				 'residues':res_list,'combined':'yes'})
