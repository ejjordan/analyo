#!/usr/bin/env python

from common_plot import *
from common_hbonds import *
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.patches as mpatches
import json

#---settings                                                                                                                                                  
plotname = 'hbonds_altered'
data,calc = plotload(plotname,work)

def construct_timeseries(work,sn):
	slices=work.slice_meta.slices[sn]
	slice_group=slices.keys()[0]
	slice_dict=slices[slice_group].flat()
	timeseries=np.arange(slice_dict['start'],slice_dict['end']+1,slice_dict['skip'])
	return timeseries

def chunk_hbonds(work,unpacked_hbonds,hbond_keys,num_chunks=10,bond_list=None,
				 divy=False,deltas=False):
	chunk_hbonds={}
	for sn in hbond_keys:
		hbonds=unpacked_hbonds[sn]
		timeseries=construct_timeseries(work,sn)
		timeseries_len=len(timeseries)
		chunk_len=timeseries_len/num_chunks
		chunk_steps=range(0,timeseries_len,chunk_len)
		chunk_hbonds[sn]={'kcat':hbonds['kcat'],'name':sn}
		for i in range(len(chunk_steps)-1):
			start_time=timeseries[chunk_steps[i]]
			stop_time=timeseries[chunk_steps[i+1]]
			for bond in hbonds:
				if bond_list and bond not in bond_list: continue
				if bond=='kcat': continue
				if bond not in chunk_hbonds[sn]:
					chunk_hbonds[sn][bond]={'donor_residx':hbonds[bond]['donor_residx'],
											'acceptor_residx':hbonds[bond]['acceptor_residx'],
											'times':{},'deltas':{}}
					if divy:
						chunk_hbonds[sn][bond]['donor_loc']=hbonds[bond]['donor_loc']
						chunk_hbonds[sn][bond]['acceptor_loc']=hbonds[bond]['acceptor_loc']
					else:
						chunk_hbonds[sn][bond]['donor_restype']=hbonds[bond]['donor_restype']
						chunk_hbonds[sn][bond]['donor_HA']=hbonds[bond]['donor_HA']
						chunk_hbonds[sn][bond]['acceptor_restype']=hbonds[bond]['acceptor_restype']
						chunk_hbonds[sn][bond]['acceptor_HA']=hbonds[bond]['acceptor_HA'],
				times=[time for time in hbonds[bond]['times'] if
					   float(time)>float(start_time) and float(time)<float(stop_time)]
				if times:
					chunk_hbonds[sn][bond]['times'][i]={'start_time':start_time,
														'stop_time':stop_time, 'times':times,
														'occupancy':len(times)/float(chunk_len)}
				if deltas and i>0:
					one=0;two=0
					if i in chunk_hbonds[sn][bond]['times']:
						one=chunk_hbonds[sn][bond]['times'][i]['occupancy']
					if i-1 in chunk_hbonds[sn][bond]['times']:
						two=chunk_hbonds[sn][bond]['times'][i-1]['occupancy']
					chunk_hbonds[sn][bond]['deltas'][i]=one-two

	return chunk_hbonds


def thresh_plotter(thresh,deltas=False,chunks=10,plot=True,title=None,meta=None,plot_threshold=None,residue_to_highlight=None):
	xtick_labels=[];xtick_labels_sns=[];label_colors=[];label_marks=[];values=[]
	for sn in thresh:
		for bond in thresh[sn]:
			if bond=='name' or bond=='kcat': continue
			for chunk in range(chunks):
				chunk+=1
				if bond in thresh[sn] and chunk in thresh[sn][bond]['times'] and abs(thresh[sn][bond]['deltas'][chunk])>plot_threshold:
					values.append(thresh[sn][bond]['deltas'][chunk])
					xtick_labels_sns.append(sn)
					mark='black'
					if residue_to_highlight:
						bond_data=bond.split()
						if str(residue_to_highlight) in bond_data:
							mark='red'
					xtick_labels.append(' '.join(sn.split('_'))+', '+bond)
					label_marks.append(mark)

	labels=label_maker(data=thresh, name_list=xtick_labels_sns,
					   max_inactive=max_inactive, min_active=min_active)
	fig, ax = plt.subplots(figsize=(45,20))
	x_ticks = np.arange(len(xtick_labels))
	ax.set_xlim([-1,len(xtick_labels)+1])
	width=0.75
	color_list=[color_dict[i] for i in labels]
	bar = ax.bar(x_ticks, values, width, color=color_list)
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(xtick_labels, rotation='vertical', ha='center',size='x-small')
	ax.xaxis.grid(True)
	[line.set_linestyle('-.') for line in ax.get_xgridlines()]
	ax.set_axisbelow(True)
	plt.ylabel('Change in H-bond occupancy',size='x-large')
	if residue_to_highlight:
		[label_text.set_color(label_marks[i]) for i,label_text in enumerate(plt.gca().get_xticklabels())]
		total_bonds=len(label_marks)
		hili_bonds=len([i for i in label_marks if i=='red'])
		print "number of bonds: {0}\nnumber of bonds highlighted: {1}".format(total_bonds,hili_bonds)
		title+="\nhighlighting residue %s-%d accounting for %d of %d hbonds"%('ARG',residue_to_highlight,
																			  hili_bonds,total_bonds)
		if title!=None:
			plt.title(title,size='x-large')
	alpha=1
	patches={
		True:mpatches.Patch(color=color_dict[True], label='activating', alpha=alpha),
		False:mpatches.Patch(color=color_dict[False], label='non-activating', alpha=alpha),
		'I':mpatches.Patch(color=color_dict['I'], label='mildly activating'),
		'wt':mpatches.Patch(color=color_dict['wt'], label='wild type', alpha=alpha),
		'maybe':mpatches.Patch(color=color_dict['maybe'], label='unknown', alpha=alpha)}
	used_patch=[patches[label] for label in set(label_colors)]
	used_label=[label_dict[label] for label in set(label_colors)]
	ax.legend(used_patch,used_label, loc='upper right')
	if plot:
		plt.show(block=False)
	else:
		picturesave('fig.%s'%(plotname),work.plotdir,backup=False,
					version=True,meta=meta,dpi=200)



#finalize the data and call the plotter
protein=work.meta['protein_name']
domains=get_subdomains(protein)
if not domains: print "[ERROR] no subdomains found"; exit

sort_keys=sorted(data.keys())
min_active=work.plots[plotname]['specs']['min_active']
max_inactive=work.plots[plotname]['specs']['max_inactive']
hbts=hbonds_timesteps(work,data,sort_keys,divy=True,
					  donor_reslist=domains['$\\alpha$C helix, activation loop'],
					  acceptor_reslist=domains['$\\alpha$C helix, activation loop'])


if protein=='alk':
	alk_kcat_normalize(work,hbts)
hili_res=None;plot_thresh=0.3
#!note that altered is the reverse of what it says it is
altered=False;altered_thresh=0.75
if altered:
	bond_list=[item for sublist in [hbts[sn].keys() for sn in hbts] for item in sublist]
else:
	combos=combine_hbonds(hbts,sort_keys,divy=True,num_replicates=2)
	deltas,keys=occupancy_diff(combos,reference='inactive_wt',threshold=altered_thresh)
	bond_list=[item for sublist in [deltas[key]['bonds'] for key in deltas] for item in sublist]
bond_list=list(set(bond_list))
chunks=chunk_hbonds(work,hbts,sort_keys,bond_list=bond_list,divy=True,num_chunks=2,deltas=True)
title='$>${0} $\%$ occupancy change H-bonds'.format(plot_thresh*100)

thresh_plotter(chunks,chunks=2,deltas=True,residue_to_highlight=hili_res,title=title,
			   plot=False,plot_threshold=plot_thresh,
			   meta={'plot threshold':plot_thresh,'bond_list':bond_list,
					 'altered theshold':altered_thresh,'highlighted_residue':hili_res,
					 'only show delta hbonds':altered})
