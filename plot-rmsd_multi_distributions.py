#!/usr/bin/python -i

#---plot prep
execfile('./omni/base/header.py')
from plotter import *
from base.store import plotload
import numpy as np
import re
from matplotlib.patches import Rectangle

#---load data for plotting
plotname = 'rmsd_multi_distributions'
data,calc = plotload(plotname,work)


import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)

#---MAIN

#---settings
#import pdb;pdb.set_trace()

def plot_rmsd(rmsd,sim_name,avg_window=20,store=True):
	save_name=sim_name
	#if "_replicate" in sim_name: sim_name=sim_name[:-10]
	sim_name=re.sub('_',' ',sim_name.upper())
	title='ALK {0} C$\\alpha$ RMSDs'.format(sim_name)
	fig, ax_array = plt.subplots(2,2,figsize=(12,10))
        name_list=["all residues","$\\alpha$C helix","activation loop",
                   "$\\alpha$C helix & activation loop"]
        act_data=['CA_act','ach_act','aloop_act','ach_aloop_act']
        inact_data=['CA_inact','ach_inact','aloop_inact','ach_aloop_inact']
        ind1=[0,0,1,1];ind2=[0,1,0,1]
        alpha=0.6
        act_rect = Rectangle((0, 0), 1, 1, fc="green", alpha=alpha)
        inact_rect = Rectangle((0, 0), 1, 1, fc="red", alpha=alpha)
        for name,act,inact,i,j in zip(name_list,act_data,inact_data,ind1,ind2):
            act_rmsd=rmsd_data[act][:,3]
            inact_rmsd=rmsd_data[inact][:,3]
            act_mean=np.mean(act_rmsd);act_sd=np.std(act_rmsd)
            inact_mean=np.mean(inact_rmsd);inact_sd=np.std(inact_rmsd)
            act_time=rmsd_data[act][:,1]/1000 #put in ns
            inact_time=rmsd_data[inact][:,1]/1000
            max_rmsd=max(max(act_rmsd),max(inact_rmsd))
            rmsd_bins = np.arange(0,max_rmsd*1.1,0.1)
            act_label="RMSD to active\n$\mu$={0:.2f}\t$\sigma$={1:.2f}".format(act_mean,act_sd)
            inact_label="RMSD to inactive\n$\mu$={0:.2f}\t$\sigma$={1:.2f}".format(
                inact_mean,inact_sd)
            title_label=name
            act_counts,act_bins = np.histogram(act_rmsd,bins=rmsd_bins,normed=True)
            inact_counts,inact_bins = np.histogram(inact_rmsd,bins=rmsd_bins,normed=True)
            ax_array[i,j].fill_between((act_bins[1:]+act_bins[:-1])/2.,0,act_counts,
                                       alpha=alpha,color='g',lw=0)
            ax_array[i,j].fill_between((inact_bins[1:]+inact_bins[:-1])/2.,0,inact_counts,
                                       alpha=alpha,color='r',lw=0)
            ax_array[i,j].set_xlim(0,max_rmsd+max_rmsd*0.1)
            ax_array[i,j].set_ylim(0,4)
            ax_array[i,j].legend([inact_rect,act_rect],[inact_label,act_label],loc='best')
            ax_array[i,j].grid(axis='y')
            ax_array[i,j].set_ylabel("(A. U.)")
            ax_array[i,j].set_xlabel("RMSD $\mathrm{(\\AA)}$")
            ax_array[i,j].set_title(title_label)

	plt.suptitle(title,size='large')
	plt.tight_layout(rect=[0, 0, 1, 0.99])
	if store:
		picturesave('fig.%s.%s'%(plotname,save_name),work.plotdir,backup=False,
                            version=True,meta={'name':sim_name},dpi=200)
	else:
		plt.show(block=False)


for key in data.keys():
    rmsd_data=data[key]['data']
    plot_rmsd(rmsd_data,key,store=True)
