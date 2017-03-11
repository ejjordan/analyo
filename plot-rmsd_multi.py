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
	act_rmsd=rmsd_data['CA_act'][:,3]
	inact_rmsd=rmsd_data['CA_inact'][:,3]
        cov_mat=np.cov(act_rmsd,inact_rmsd)
        act_variance=cov_mat[0][0]
        inact_variance=cov_mat[1][1]
        act_inact_covariance=cov_mat[0][1]
	act_time=rmsd_data['CA_act'][:,1]/1000 #put in ns
	inact_time=rmsd_data['CA_inact'][:,1]/1000
	ax_array[0,0].set_xlim(act_time[0],act_time[-1])
	act_time=act_time[avg_window-1:]
	inact_time=inact_time[avg_window-1:]
	act_rmsd=running_mean(act_rmsd,avg_window)
	inact_rmsd=running_mean(inact_rmsd,avg_window)
        act_label="RMSD to active\nvariance {0:.2f}".format(act_variance)
        inact_label="RMSD to inactive\nvariance {0:.2f}".format(inact_variance)
        title_label="all residues\ncovariance {0:.2f}".format(act_inact_covariance)
	ax_array[0,0].plot(act_time, act_rmsd, 'g', label=act_label)
	ax_array[0,0].plot(inact_time, inact_rmsd, 'r', label=inact_label)
	ax_array[0,0].legend(loc='best')
	ax_array[0,0].grid(axis='y')
	ax_array[0,0].set_xlabel("Time (ns)")
	ax_array[0,0].set_ylabel("RMSD $\mathrm{(\\AA)}$")
	ax_array[0,0].set_title(title_label)

	act_rmsd=rmsd_data['ach_act'][:,3]
	inact_rmsd=rmsd_data['ach_inact'][:,3]
        cov_mat=np.cov(act_rmsd,inact_rmsd)
        act_variance=cov_mat[0][0]
        inact_variance=cov_mat[1][1]
        act_inact_covariance=cov_mat[0][1]
	act_time=rmsd_data['ach_act'][:,1]/1000 #put in ns
	inact_time=rmsd_data['ach_inact'][:,1]/1000
	ax_array[0,1].set_xlim(act_time[0],act_time[-1])
	act_time=act_time[avg_window-1:]
	inact_time=inact_time[avg_window-1:]
	act_rmsd=running_mean(act_rmsd,avg_window)
	inact_rmsd=running_mean(inact_rmsd,avg_window)
        act_label="RMSD to active\nvariance {0:.2f}".format(act_variance)
        inact_label="RMSD to inactive\nvariance {0:.2f}".format(inact_variance)
        title_label="$\\alpha$C helix\ncovariance {0:.2f}".format(act_inact_covariance)
	ax_array[0,1].plot(act_time, act_rmsd, 'g', label=act_label)
	ax_array[0,1].plot(inact_time, inact_rmsd, 'r', label=inact_label)
	ax_array[0,1].legend(loc='best')
	ax_array[0,1].grid(axis='y')
	ax_array[0,1].set_xlabel("Time (ns)")
	ax_array[0,1].set_ylabel("RMSD $\mathrm{(\\AA)}$")
	ax_array[0,1].set_title(title_label)

	act_rmsd=rmsd_data['aloop_act'][:,3]
	inact_rmsd=rmsd_data['aloop_inact'][:,3]
        cov_mat=np.cov(act_rmsd,inact_rmsd)
        act_variance=cov_mat[0][0]
        inact_variance=cov_mat[1][1]
        act_inact_covariance=cov_mat[0][1]
	act_time=rmsd_data['aloop_act'][:,1]/1000 #put in ns
	inact_time=rmsd_data['aloop_inact'][:,1]/1000
	ax_array[1,0].set_xlim(act_time[0],act_time[-1])
	act_time=act_time[avg_window-1:]
	inact_time=inact_time[avg_window-1:]
	act_rmsd=running_mean(act_rmsd,avg_window)
	inact_rmsd=running_mean(inact_rmsd,avg_window)
        act_label="RMSD to active\nvariance {0:.2f}".format(act_variance)
        inact_label="RMSD to inactive\nvariance {0:.2f}".format(inact_variance)
        title_label="activation loop\ncovariance {0:.2f}".format(act_inact_covariance)
	ax_array[1,0].plot(act_time, act_rmsd, 'g', label=act_label)
	ax_array[1,0].plot(inact_time, inact_rmsd, 'r', label=inact_label)
	ax_array[1,0].legend(loc='best')
	ax_array[1,0].grid(axis='y')
	ax_array[1,0].set_xlabel("Time (ns)")
	ax_array[1,0].set_ylabel("RMSD $\mathrm{(\\AA)}$")
	ax_array[1,0].set_title(title_label)

	act_rmsd=rmsd_data['ach_aloop_act'][:,3]
	inact_rmsd=rmsd_data['ach_aloop_inact'][:,3]
        cov_mat=np.cov(act_rmsd,inact_rmsd)
        act_variance=cov_mat[0][0]
        inact_variance=cov_mat[1][1]
        act_inact_covariance=cov_mat[0][1]
	act_time=rmsd_data['ach_aloop_act'][:,1]/1000 #put in ns
	inact_time=rmsd_data['ach_aloop_inact'][:,1]/1000
	ax_array[1,1].set_xlim(act_time[0],act_time[-1])
	act_time=act_time[avg_window-1:]
	inact_time=inact_time[avg_window-1:]
	act_rmsd=running_mean(act_rmsd,avg_window)
	inact_rmsd=running_mean(inact_rmsd,avg_window)
        act_label="RMSD to active\nvariance {0:.2f}".format(act_variance)
        inact_label="RMSD to inactive\nvariance {0:.2f}".format(inact_variance)
        title_label="$\\alpha$C helix & activation loop\ncovariance {0:.2f}".format(
            act_inact_covariance)
	ax_array[1,1].plot(act_time, act_rmsd, 'g', label=act_label)
	ax_array[1,1].plot(inact_time, inact_rmsd, 'r', label=inact_label)
	ax_array[1,1].legend(loc='best')
	ax_array[1,1].grid(axis='y')
	ax_array[1,1].set_xlabel("Time (ns)")
	ax_array[1,1].set_ylabel("RMSD $\mathrm{(\\AA)}$")
	ax_array[1,1].set_title(title_label)


	plt.suptitle(title,size='large')
	plt.tight_layout(rect=[0, 0, 1, 0.99])
	if store:
		picturesave('fig.%s.%s'%(plotname,save_name),work.plotdir,backup=False,
					version=True,meta={'avg_window':avg_window, 'name':sim_name},dpi=200)
	else:
		plt.show(block=False)


for key in data.keys():
    rmsd_data=data[key]['data']
    plot_rmsd(rmsd_data,key,store=True)

    """
    A_ca=rmsd_data['CA_act'][:,3]
    I_ca=rmsd_data['CA_inact'][:,3]
    A_ac=rmsd_data['ach_act'][:,3]
    I_ac=rmsd_data['ach_inact'][:,3]
    A_al=rmsd_data['aloop_act'][:,3]
    I_al=rmsd_data['aloop_inact'][:,3]
    A_aa=rmsd_data['ach_aloop_act'][:,3]
    I_aa=rmsd_data['ach_aloop_inact'][:,3]

    for time in range(len(I_ac)):
        if I_ac[time]>A_ac[time]:
            print key, time
            break
    


    print key.split('_')[0],
    print "\tCA\tach\taloop\taa"
    print "act\t{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}".format(np.mean(A_ca),np.mean(A_ac),np.mean(A_al),np.mean(A_aa))
    print "inact\t{0:.2f}\t{1:.2f}\t{2:.2f}\t{3:.2f}".format(np.mean(I_ca),np.mean(I_ac),np.mean(I_al),np.mean(I_aa))

    for time in range(len(act_rmsd)):
        if act_rmsd[time]-inact_rmsd[time]<0.5:
            print key, time
            break
    """
