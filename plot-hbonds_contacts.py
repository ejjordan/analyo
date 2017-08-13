#!/usr/bin/env python

"""
ACTINLINK BONDS
see lab notebook for details
"""

import time
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.tools import status,framelooper
import common_plots

#---block: what to plot
routine = ['contact_map']
#relevant_
residue_types = {'ARG':'basic','HIS':'basic','LYS':'basic',
	'ASP':'acidic','GLU':'acidic','SER':'polar','THR':'polar','ASN':'polar','GLN':'polar',
	'ALA':'hydrophobic','VAL':'hydrophobic','ILE':'hydrophobic','LEU':'hydrophobic',
		'MET':'hydrophobic','PHE':'hydrophobic','TYR':'hydrophobic','TRP':'hydrophobic',
	'CYS':'special','SEC':'special','GLY':'special','PRO':'special'}

residue_type_colors = {'basic':'Blues','acidic':'Reds','hydrophobic':'Greens',
	'polar':'Purples','special':'Oranges'}
plotspec = {'fs_xlabel':14,'fs_ylabel':20,'fs_title':20,
	'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,'fs_legend_title':20,
	'binary_color_intensity':0.5}
residue_colors = dict([(name,residue_type_colors[residue_types[name]]) for name in residue_types])
ticks_font = mpl.font_manager.FontProperties(family='Latin Modern Mono',style='normal',
	size=14,weight='normal',stretch='normal')

#---block: mimic the coda in contacts.py
def hydrogen_bond_compactor():
	global data_hbonds,bonds,obs
	for sn in sns:
		#---custom mapping for collecting hydrogen bonds
		bonds,obs = [data_hbonds[sn]['data'][i] for i in ['bonds','observations']]
		#---! get resids for the protein and lipid_resnames from contact maps
		lipid_resnames = np.unique(
			data_contacts[sn]['data']['bonds'][:,rowspec.index('target_resname')])
		resids = data_contacts[sn]['data']['subject_residues_resids']
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		#---compute loop
		looper = [{'resid':resid,'resname_set':resname_set} 
			for resid in resids for resname_name,resname_set in resname_combos]
		compute_function = bond_counter
		incoming = basic_compute_loop(compute_function,looper,run_parallel=True)
		#---tacking on compacted data to mimic the form of the contact maps
		data_hbonds[sn]['data']['hbonds_compacted'] = np.array(incoming)
		data_hbonds[sn]['data']['pairs_resid_resname'] = np.array([(resid,resname_name) 
			for resid in resids for resname_name,resname_set in resname_combos]).astype(str)

#---block: post-post processing
def make_postdat():
	global sns,data_hbonds,data_contacts,bond_mappings
	#---easy lookup for multiple upstream data types
	pointers = {'hbonds':data_hbonds,'contacts':data_contacts}
	#---postdat is indexed by simulation, then bond type
	postdat = dict([(sn,{}) for sn in sns])
	#---loop over simulations
	for sn in sns:
		#---common features
		lipid_resnames = np.unique(
			data_contacts[sn]['data']['bonds'][:,rowspec.index('target_resname')])
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		lipid_resnames = list(zip(*resname_combos))[0]
		postdat[sn]['lipid_resnames'] = lipid_resnames
		postdat[sn]['resids'] = data_contacts[sn]['data']['subject_residues_resids']
		postdat[sn]['subject_residues_resnames'] = data_contacts[sn]['data']['subject_residues_resnames']
		postdat[sn]['times'] = data_contacts[sn]['data']['times']
		#---loop over bond mappings
		for bond in bond_mappings:
			this_data = pointers[bond['upstream']]
			#---raw data are timewise bond counts for resid,lipid_name pairs
			raw = dict([((int(i[0]),i[1]),j)
				for i,j in zip(
					this_data[sn]['data']['pairs_resid_resname'],
					#---this line specifies the two kinds of data to extract from contacts: either
					#---...either the explicit or reduced bonds
					this_data[sn]['data'][bond['post_key']])])
			#---post is the raw data organized by lipid_name into matrices suitable for plotting
			postdat[sn][bond['name']] = dict([(resname,
				np.array([raw[(r,resname)] 
				for r in postdat[sn]['resids']])) for resname in postdat[sn]['lipid_resnames']])
	return postdat

#---block: contact map code
def colorstreak_contact_map(sns,postdat,bond_name,plotspec,fn,**kwargs):
	"""
	Final resting place for the contact map code.
	"""
	do_crazy_colors = kwargs.get('do_crazy_colors',True)
	cmap = kwargs.get('cmap','Greys')

	#---same residue comparisons
	lipid_resnames = list(set([tuple(postdat[sn]['lipid_resnames']) for sn in sns]))
	if not len(lipid_resnames)==1: raise Exception('cannot use different resnames across simulations')
	else: lipid_resnames = lipid_resnames[0]
	if kwargs.get('lipid_resnames',None):
		lipid_resnames = [i for i in lipid_resnames if i in kwargs['lipid_resnames']]
	ceiling = float(max([max([i.max() for i in postdat[sn][bond_name].values()]) for sn in sns]))
	#---plot
	axes,fig = panelplot(
		layout={'out':{'hspace':0.05,'wspace':0.05,'grid':[1,len(lipid_resnames)]},
		'ins':[{'hspace':0.2,'wspace':0.05,'grid':[len(sns),1]} 
		for i in lipid_resnames]},figsize=plotspec['figsize'])
	for ss,sn in enumerate(sns):
		for rr,resname in enumerate(lipid_resnames):
			#only show residues that have at least one hydrogen bond during 10% of 
			#...the simulation
			hbond_residxs = [i for i in range(len(postdat[sn][bond_name][resname])) if
							 sum(postdat[sn][bond_name][resname][i])>len(postdat[sn]['times'])/10.]
			hbond_residxs.reverse()
			try:
				ax = axes[rr][ss]
			except:
				ax = axes[rr]
			resnames = postdat[sn]['subject_residues_resnames'][hbond_residxs]
			resids = postdat[sn]['resids'][hbond_residxs]
			resinfo = zip(resnames,resids)
			if rr==0: 
				ax.set_ylabel('residues (resname-resnum)',fontsize=plotspec['fs_ylabel'])
				if work.plots[plotname].get('settings',{}).get('show_residue_names',True):
					#---never set the ticklabels before the ticks
					ax.set_yticks(np.arange(len(resnames))+0.5)
					ax.set_yticklabels(['%s-%s'%(res,idx) for res,idx in resinfo],
						fontsize=plotspec['fs_yticks'])
					for label in ax.get_yticklabels():
						label.set_color(mpl.cm.__dict__[
							residue_colors.get(residue_codes[label._text[:3]],
							'Greys')](plotspec['label_color_strength'])[:-1])
			else: ax.set_yticks([])
			if ss==0:
				ax.set_title('%s\nhydrogen bonding persistence map'%' '.join(sn.split('_')),
					fontsize=plotspec['fs_title'])

			if do_crazy_colors:
				image_data = np.array([mpl.cm.__dict__[
					residue_colors.get(resnames[rnum],'Greys')](
						postdat[sn][bond_name][resname][hbond_residxs[rnum]].astype(float)/ceiling).T[:-1].T 
									   for rnum in range(len(hbond_residxs))])
			else: image_data = postdat[sn][bond_name][resname]
			duration = (postdat[sn]['times'].max()-postdat[sn]['times'].min())*10**-3
			xtick_interval = work.plots[plotname].get('settings',{}).get('xtick_interval',
				plotspec['time_tick_interval'])
			ax.set_xticks(np.arange(xtick_interval,duration,xtick_interval))
			ax.set_xlabel('time (ns)',fontsize=plotspec['fs_xlabel'])
			#---! removed division by the following, designed to make things lighter when reduced
			#---! ...{'explicit':1.0,'reduced':plotspec['binary_color_intensity']}[mode]
			im = ax.imshow(image_data,
				origin='lower',interpolation='nearest',extent=[0,duration,0,len(resids)],
				aspect=float(duration)/len(resids),
				vmax=ceiling,vmin=0,cmap=mpl.cm.__dict__[cmap])
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		if ss==0:
			#---color bar on the side
			axins = inset_axes(ax,width="8%",height="70%",loc=3,
				bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
			cbar_title = kwargs.get('cbar_title','score')
			cbar.ax.set_title(cbar_title)

	#---title on the legend
	title = kwargs.get('legend_title','contacts\n'+r'$\mathrm{r \leq %.1f \AA}$'%(cutoff))
	patches,labels = zip(*[(mpl.patches.Rectangle((0,0),1.0,1.0,
		fc=mpl.cm.__dict__[residue_type_colors[kind]](plotspec['legend_color_strength'])),kind) 
		for kind in residue_type_colors])
	legend = axes[0].legend(patches,labels,loc='upper left',fontsize=plotspec['fs_legend']-1,title='residue\ntype',
							bbox_to_anchor=(1.005,0.0,1.,0.97),shadow=True,fancybox=True)
	frame = legend.get_frame()
	frame.set_edgecolor('black')
	frame.set_facecolor('white')
	plt.setp(legend.get_title(),fontsize=plotspec['fs_legend_title']-1)
	picturesave(fn,work.plotdir,backup=False,version=True,meta={},extras=[])


#---block: counting hydrogen bonds
def bond_counter(resid,resname_set):
	"""
	Adapted from the explicit section of the contacts.py function called "count_reduced_contact".
	The hydrogen bonding currently uses the same rowspec as the contacts.
	"""
	global bonds,obs,rowspec
	#---filter the observations by the protein residue (subject_resid) and target resname
	#---...providing a result
	which = np.where(np.all((bonds[:,rowspec.index('subject_resid')].astype(int)==resid,
		np.in1d(bonds[:,rowspec.index('target_resname')],resname_set)),axis=0))
	result = obs.T[which].sum(axis=0)
	return result

#---block: stolen from contacts.py
def basic_compute_loop(compute_function,looper,run_parallel=True,debug=False):
	"""Canonical form of the basic compute loop."""
	start = time.time()
	if run_parallel:
		incoming = Parallel(n_jobs=8,verbose=10 if debug else 0)(
			delayed(compute_function,has_shareable_memory)(**looper[ll]) 
			for ll in framelooper(len(looper),start=start))
	else: 
		incoming = []
		for ll in framelooper(len(looper)):
			incoming.append(compute_function(**looper[ll]))
	return incoming

#---block: import the post-processed data	
if 'data' not in globals(): 
	sns_contacts,(data_contacts,calc_contacts) = work.sns(),plotload('contacts',work)
	sns_hbonds,(data_hbonds,calc_hbonds) = work.sns(),plotload('hydrogen_bonding',work)
	if sns_hbonds!=sns_contacts: 
		raise Exception('collections for hydrogen_bonding and contacts are not equal')
	else: sns = sns_hbonds
	
	#---set the cutoff in the yaml file to make this plot because we pull from multiple upstream sources
	cutoff = calc_contacts['calcs']['specs']['specs']['cutoff']
	#---map data type onto keys
	bond_mappings = [
		{'name':'hbonds','post_key':'hbonds_compacted','upstream':'hbonds'},]
	#---specify the upstream data and naming scheme
	rowspec = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom']
	lipid_label = lambda x: dict([(i,'$$\mathrm{{PIP}_{2}}$$') 
		for i in work.vars.get('selectors',{}).get('resnames_PIP2',{})]).get(x,x)
	sn_title = lambda sn: '%s%s'%(work.meta[sn].get('label',re.sub('_','-',sn)),
		'' if not work.meta[sn].get('cholesterol',False) else '\n(cholesterol)')

	#---post-post processing
	if 'postdat' not in globals(): 
		hydrogen_bond_compactor()
		postdat = make_postdat()

#---block: all chemical bonds on one contact map
plotspec = {'fs_xlabel':14,'fs_ylabel':20,'fs_title':20,
	'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,
	'label_color_strength':1.0,'fs_legend_title':20,
	'binary_color_intensity':0.5,'figsize':(14,14),'time_tick_interval':20,
	'fs_xticks':11,'fs_yticks':11}
plotspec_small = {'figsize':(8,8),'fs_legend_title':14,'fs_legend':12,'fs_xticks':12,'fs_yticks':12}
specs = {'all':{'sns':sns,'lipid_resnames':['all lipids','all lipids']}}
for sn in sns:
	for bond in bond_mappings:
		for specname,spec in specs.items():		
			this_plotspec = dict(plotspec)
			this_plotspec.update(**spec.get('plotspec',{}))
			kwargs = dict(postdat=postdat,bond_name=bond['name'],plotspec=this_plotspec)
			kwargs.update(sns=[sn])
			kwargs.update(lipid_resnames=spec.get('lipid_resnames',None))
			figname = 'hbonds-contacts'
			fn = 'fig.%s.%s'%(figname,sn)
			kwargs.update(legend_title='residue\ntypes',cbar_title='counts')
			colorstreak_contact_map(fn=fn,**kwargs)
