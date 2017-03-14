#!/usr/bin/python -i

#---plot prep
#if 'plotload' not in globals(): execfile('/etc/pythonstart')
execfile('./omni/base/header.py')
from plotter import *
from base.store import plotload
#execfile('./calcs/specs/figures.py')
#execfile('./calcs/specs/colors.py')
import numpy as np
import re

#---load data for plotting
plotname = 'rmsd_multi'
data,calc = plotload(plotname,work)

#---FUNCTIONS
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N 


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
        for name,act,inact,i,j in zip(name_list,act_data,inact_data,ind1,ind2):
            act_rmsd=rmsd_data[act][:,3]
            inact_rmsd=rmsd_data[inact][:,3]
            cov_mat=np.cov(act_rmsd,inact_rmsd)
            act_variance=cov_mat[0][0]
            inact_variance=cov_mat[1][1]
            act_inact_covariance=cov_mat[0][1]
            act_time=rmsd_data[act][:,1]/1000 #put in ns
            inact_time=rmsd_data[inact][:,1]/1000
            act_time=act_time[avg_window-1:]
            inact_time=inact_time[avg_window-1:]
            act_rmsd=running_mean(act_rmsd,avg_window)
            inact_rmsd=running_mean(inact_rmsd,avg_window)
            act_label="RMSD to active\nvariance {0:.2f}".format(act_variance)
            inact_label="RMSD to inactive\nvariance {0:.2f}".format(inact_variance)
            title_label="{0}\ncovariance {1:.2f}".format(name,act_inact_covariance)
            ax_array[i,j].set_xlim(act_time[0],act_time[-1])
            ax_array[i,j].plot(act_time, act_rmsd, 'g', label=act_label)
            ax_array[i,j].plot(inact_time, inact_rmsd, 'r', label=inact_label)
            ax_array[i,j].legend(loc='best')
            ax_array[i,j].grid(axis='y')
            ax_array[i,j].set_xlabel("Time (ns)")
            ax_array[i,j].set_ylabel("RMSD $\mathrm{(\\AA)}$")
            ax_array[i,j].set_title(title_label)

	plt.suptitle(title,size='large')
	plt.tight_layout(rect=[0, 0, 1, 0.99])
	if store:
		picturesave('fig.%s.%s'%(plotname,save_name),work.plotdir,backup=False,
					version=True,meta={'avg_window':avg_window, 'name':sim_name},dpi=200)
	else:
		plt.show(block=False)


for key in data.keys():
    rmsd_data=data[key]['data']
    plot_rmsd(rmsd_data,key,store=False)

