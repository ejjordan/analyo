#!/usr/bin/env python

import os

hydrophobic=['PHE','TYR','ILE','LEU','VAL','TRP']
polar=['ARG','LYS','GLU','ASP','HIS','SER','THR','ASN','GLN']
charged=['ARG','LYS','GLU','ASP','HIS']
hydrophobic_core={'alk':[1096,1098,1170,1171,1174,1179,1239,1245,1271,1240],
                  'egfr':[765,766,769,774,827,828,833,856],
                  'braf':[504,505,513,566,567,572,595],
                  'her2':[773,774,777,782,835,836,841,864]}
Rspine={'braf':[505,516,574,595],
        'alk':[1171,1182,1247,1271],
        'her2':[774,785,843,864]}
hydrophobic_core_noJM={'alk':[1170,1171,1174,1179,1239,1245,1271,1240]}
label_dict={True:'activating', False:'non-activating', 'I':'mildly activating',
			'maybe':'unknown', 'wt':'wild type'}
color_dict={True:'r', False:'g', 'I':'m', 'maybe':'b', 'wt':'k'}
pattern_dict = {'nucleotide binding loop':'//', '$\\alpha$C helix':'..',
				'catalytic loop':'\\', 'activation loop':'oo'}
#ipython notebooks are executed from calcs but make calls are a level up
domain_fh=os.path.abspath('calcs/kinase_subdomains')
if not os.path.isfile(domain_fh): domain_fh=os.path.abspath('kinase_subdomains')
def get_subdomains(protein,domainfile=domain_fh):
    with open(domainfile,'r') as fp: domain_lines=fp.readlines(); fp.close()
    for line in domain_lines:
        line.strip('\n')
        dom_info=line.split(' ')
        if dom_info[0].upper()==protein.upper():
            domains_borders={'kd_start':int(dom_info[1]), 'kd_end':int(dom_info[2]),
                     'ploop_start':int(dom_info[3]), 'ploop_end':int(dom_info[4]),
                     'alphac_start':int(dom_info[5]), 'alphac_end':int(dom_info[6]),
                     'catloop_start':int(dom_info[7]), 'catloop_end':int(dom_info[8]),
                     'activation_start':int(dom_info[9]), 'activation_end':int(dom_info[10])}
            domains={'kinase domain':range(int(dom_info[1]), int(dom_info[2])),
                     'nucleotide binding loop':range(int(dom_info[3]), int(dom_info[4])),
                     u"$\\alpha$C helix":range(int(dom_info[5]), int(dom_info[6])),
                     'catalytic loop':range(int(dom_info[7]), int(dom_info[8])),
                     'activation loop':range(int(dom_info[9]), int(dom_info[10]))}
    domains[u"$\\alpha$C helix, activation loop"]=list(set(set(domains['activation loop'])|
                                                           set(domains[u"$\\alpha$C helix"])))
    domains[u"$\\alpha$C helix, nucleotide binding loop"]=list(set(set(domains['nucleotide binding loop'])|
                                                                   set(domains[u"$\\alpha$C helix"])))
    domains[u"$\\alpha$C helix, catalytic loop"]=list(set(set(domains['catalytic loop'])|
                                                          set(domains[u"$\\alpha$C helix"])))
    domains['nucleotide binding loop, activation loop']=list(set(set(domains['activation loop'])|
                                                                 set(domains['nucleotide binding loop'])))
    domains['nucleotide binding loop, catalytic loop']=list(set(set(domains['nucleotide binding loop'])|
                                                                set(domains['catalytic loop'])))
    domains['activation loop, catalytic loop']=list(set(set(domains['catalytic loop'])|
                                                        set(domains['activation loop'])))
    domains[u"$\\alpha$C helix, catalytic loop, activation loop"]=list(set(set(domains['activation loop'])|
                                                                           set(domains[u"$\\alpha$C helix"])|
                                                                           set(domains['catalytic loop'])))
    domains[u"$\\alpha$C helix, catalytic loop, nucleotide binding loop"]=list(set(set(
        domains['nucleotide binding loop'])|set(domains[u"$\\alpha$C helix"])|set(domains['catalytic loop'])))
    domains[u"$\\alpha$C helix, nucleotide binding loop, activation loop"]=list(set(set(
        domains['activation loop'])|set(domains[u"$\\alpha$C helix"])|set(domains['nucleotide binding loop'])))
    domains[u"nucleotide binding loop, $\\alpha$C helix, catalytic loop, activation loop"]=list(set(set(
        domains['activation loop'])|set(domains[u"$\\alpha$C helix"])|set(
            domains['catalytic loop'])|set(domains['nucleotide binding loop'])))
    domains[u"kinase domain - ($\\alpha$C helix, activation loop)"]=list(set(
        domains['kinase domain']) - set(domains[u"$\\alpha$C helix, activation loop"]))
    return domains

def label_maker(data, name_list=None, max_inactive=4, min_active=4):

    """
    This function takes a data object (e.g. form sasa or hbonds) and a kcat cut-off and
	returns activation (or non-activation) labels. If supplied a name_list the labels
	will be returned in the specified order.
    """

    if not name_list: name_list=[data[sn]['name'] for sn in data]
    labels=[]
    for name in name_list:
        if 'kcat' in data[name]:
            kcat=data[name]['kcat']
            if kcat=='X': labels.append('maybe')
            elif kcat=='WT': labels.append('wt')
            elif float(kcat)>=min_active:
		        labels.append(True)
            elif float(kcat)<=max_inactive:
                labels.append(False)
            elif float(kcat)>max_inactive and float(kcat)<min_active:
		        labels.append('I') #for intermediate activity
            else: raise Exception('[ERROR]: something went wrong with the kcat data, '
                                  'see the meta file')            
        else: raise Exception('[ERROR]: no kcat information in "data" provided')
    return labels
