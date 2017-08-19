#!/usr/bin/env python

#---plot prep
from common_plot import *
import numpy as np

#---load data for plotting
plotname = 'rmsd_histogram'
data,calc = plotload(plotname,work)

#---GLOBALS
act_data=['CA_act','ach_act','aloop_act','ach_aloop_act']
inact_data=['CA_inact','ach_inact','aloop_inact','ach_aloop_inact']


def plot_rmsd(rmsd_stats,store=True):
	sns=sorted(rmsd_stats.keys())
	names=[' '.join(sn.split('_')).upper() for sn in sns]
	counter=0
	fig, ax_array = plt.subplots(2,2,figsize=(14,16))
	title='{0} C$\\alpha$ RMSDs'.format(protein)
	name_list=["all residues","$\\alpha$C helix","activation loop",
			   "$\\alpha$C helix & activation loop"]
	ind1=[0,0,1,1];ind2=[0,1,0,1];width=0.4
	for name,act,inact,i,j in zip(name_list,act_data,inact_data,ind1,ind2):
		act_ticks=np.arange(len(sns))
		inact_ticks=np.arange(len(sns))+0.5
		#counter+=1
		act_means=[];inact_means=[]
		act_stds=[];inact_stds=[]
		for sn in sns:
			act_means.append(rmsd_stats[sn][act]['mean'])
			inact_means.append(rmsd_stats[sn][inact]['mean'])
			act_stds.append(rmsd_stats[sn][act]['std'])
			inact_stds.append(rmsd_stats[sn][inact]['std'])
			#act_label="RMSD to active\nvariance {0:.2f}".format(act_variance)
			#inact_label="RMSD to inactive\nvariance {0:.2f}".format(inact_variance)
			
			#ax_array[i,j].set_xlim(act_time[0],act_time[-1])
			#ax_array[i,j].plot(act_time, act_rmsd, 'g', label=act_label)
			#ax_array[i,j].plot(inact_time, inact_rmsd, 'r', label=inact_label)
			#ax_array[i,j].legend(loc='best')
			#ax_array[i,j].set_xlabel("Time (ns)")
			#ax_array[i,j].set_ylabel("RMSD $\mathrm{(\\AA)}$")
			#plt.suptitle(title,size='large')
			#plt.tight_layout(rect=[0, 0, 1, 0.99])
		idx=range(len(sns))		
		lower_act=[act_means[k]-act_stds[k] for k in idx]
		upper_act=[act_means[k]+act_stds[k] for k in idx]
		lower_inact=[inact_means[k]+inact_stds[k] for k in idx]
		upper_inact=[inact_means[k]-inact_stds[k] for k in idx]
		
		ax_array[i,j].bar(left=act_ticks, height=upper_act, bottom=lower_act,
						  width=width, color='g')
		ax_array[i,j].bar(left=inact_ticks, height=upper_inact, bottom=lower_inact,
						  width=width, color='r')
		title_label="{0} RSMD".format(name)
		ax_array[i,j].set_title(title_label)
		ax_array[i,j].grid(axis='y')
		ax_array[i,j].grid(axis='x',linestyle='-.')
		import ipdb;ipdb.set_trace()
		ax_array[i,j].set_xticks(inact_ticks)
		ax_array[i,j].set_xticklabels(names, rotation='vertical', ha='left',size='small')

	if store:
		picturesave('fig.%s'%(plotname),work.plotdir,backup=False,
					version=True,meta={},dpi=200)
	else:
		plt.show(block=False)

def combine_replicates(sns,data,num_replicates=2):

	"""
	combine rmsd data from replicates
	"""
	combined={};count=0;unique_sns=len(sns)/num_replicates
	for i in range(unique_sns):
		name=sns[i*num_replicates]
		combined[name]={}
		for j in range(num_replicates):
			for act,inact in zip(act_data,inact_data):
				rmsd_data=data[sns[i*num_replicates+j]]['data']
				act_rmsd=rmsd_data[act][:,3]
				inact_rmsd=rmsd_data[inact][:,3]
				if act in combined[name]: 
					combined[name][act].append(act_rmsd)
				else: combined[name][act]=[act_rmsd]
				if inact in combined[name]: 
					combined[name][inact].append(inact_rmsd)
				else: combined[name][inact]=[inact_rmsd]
	for sn in combined:
		for act,inact in zip(act_data,inact_data):
			act_rmsd=combined[sn][act]
			inact_rmsd=combined[sn][inact]
			combined[sn][act]=[j for i in combined[sn][act] for j in i]
			combined[sn][inact]=[j for i in combined[sn][inact] for j in i]
	return combined


def rmsd_stats(rmsd_data):
	stats={}
	for sn in rmsd_data.keys():
		stats[sn]={}
		for act,inact in zip(act_data,inact_data):			
			act_rmsd=rmsd_data[sn][act]
			inact_rmsd=rmsd_data[sn][inact]
			act_mean=np.mean(act_rmsd);act_std=np.std(act_rmsd)
			inact_mean=np.mean(inact_rmsd);inact_std=np.std(inact_rmsd)
			cov_mat=np.cov(act_rmsd,inact_rmsd)
			act_variance=cov_mat[0][0]
			inact_variance=cov_mat[1][1]
			act_inact_covariance=cov_mat[0][1]
			stats[sn][act]={'mean':act_mean,'variance':act_variance,
							'std':act_std,'inact_covariance':act_inact_covariance}
			stats[sn][inact]={'mean':inact_mean,'variance':inact_variance,
							'std':inact_std,'act_covariance':act_inact_covariance}
			#import ipdb;ipdb.set_trace()
	return stats

protein=work.meta['protein_name']
sns=sorted(data.keys())
combos=combine_replicates(sns,data,num_replicates=2)
stats=rmsd_stats(combos)
plot_rmsd(stats,store=True)
