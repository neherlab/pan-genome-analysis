import os
import numpy as np
from Bio import Phylo
from collections import defaultdict
from sf_coreTree_json import Metadata
from sf_miscellaneous import write_pickle

class PresenceAbsenceAssociation(object):
    """docstring for Association"""
    def __init__(self, tree, meta_info):
        super(PresenceAbsenceAssociation, self).__init__()
        self.meta_info = meta_info
        if type(tree)==str and os.path.isfile(tree):
            self.tree = Phylo.load(tree_name, 'newick')
        else:
            self.tree = tree

    def set_gain_loss(self, gain_loss):
        clade = 0
        for n in self.tree.find_clades(order='preorder'):
            if n==self.tree.root:
                continue
            n.present = 'present' if gain_loss[clade] in [1, 3] else 'absent'
            n.event = (gain_loss[clade] in [1, 2])
            clade+=1

        rn = self.tree.root
        rn.present = 'present' if np.mean([c.present=='present' for c in rn])>=0.5 else 'absent'
        rn.event = True

    def calc_association_by_event(self,meta_column, transform=None, pc=3):
        '''
        calculate the mean value of the phenotype of leaves upstream and down stream
        of each branch in the tree.
        '''
        def mean_var(node, pc=3, total_var=0.0):
            nval = node.meta_count
            m = node.meta_value/nval
            var = ((node.meta_sq_value - nval*m**2) + pc*total_var)/(nval-1.0+pc)
            return (m,var, nval)

        if transform is None:
            transform = lambda x:x
        all_values = []
        values_by_state = {'present':[], 'absent':[]}
        for n in self.tree.find_clades(order='postorder'):
            if n.is_terminal():
                n.strain = n.name.split('|')[0]
                try:
                    n.meta_value = transform(self.meta_info[n.strain][meta_column])
                except:
                    print("WARNING: reading field %s for strain %s failed"%(meta_column, n.strain))
                    n.meta_value = np.nan

                if not np.isnan(n.meta_value):
                    all_values.append(n.meta_value)
                    values_by_state[n.present].append(n.meta_value)
                    n.meta_count = 1
                    n.meta_sq_value = n.meta_value*n.meta_value
                else:
                    n.meta_count = 0
                    n.meta_sq_value = np.nan

            else:
                n.meta_count = np.sum([c.meta_count for c in n if c.meta_count and not c.event])
                n.meta_value = np.sum([c.meta_value for c in n if c.meta_count and not c.event])
                n.meta_sq_value = np.sum([c.meta_sq_value for c in n if c.meta_count and not c.event])

        self.averages = {'present':[], 'absent':[]}
        rn = self.tree.root
        if rn.meta_count:
            self.averages[rn.present].append(mean_var(rn, pc=3, total_var = np.var(all_values)))
        for n in self.tree.find_clades(order='preorder'):
            if n.event and n.meta_count:
                self.averages[n.present].append(mean_var(n, pc=3, total_var = np.var(all_values)))

        from scipy.stats import ttest_ind
        self.averages = {k:np.array(v) for k,v in self.averages.iteritems()}
        p = self.averages['present']
        a = self.averages['absent']
        if len(a) and len(p) and len(values_by_state['present']) and len(values_by_state['absent']):
            #return (np.sum(p[:,0]/p[:,1])/np.sum(1./p[:,1]) - np.sum(a[:,0]/a[:,1])/np.sum(1./a[:,1]))*np.sqrt(1.0/(1.0/p.shape[0] + 1.0/a.shape[0]))
            #return (np.mean(p[:,0]) - np.mean(a[:,0]))*np.sqrt(1.0/(1.0/p.shape[0] + 1.0/a.shape[0]))
            #return (np.mean(values_by_state['present']) - np.mean(values_by_state['absent']))*np.sqrt(1.0/(1.0/len(values_by_state['present']) + 1.0/len(values_by_state['absent'])))
            return (np.mean(p[:,0]) - np.mean(a[:,0]))*np.sqrt(1.0/(1.0/p.shape[0] + 1.0/a.shape[0]))
        else:
            return np.nan

    def calc_association_simple(self,meta_column, transform=None, pc=3):
        '''
        calculate the mean value of the phenotype of leaves upstream and down stream
        of each branch in the tree.
        '''

        values_by_state = {'present':[], 'absent':[]}
        for n in self.tree.get_terminals():
            if n.is_terminal():
                n.strain = n.name.split('|')[0]
                n.meta_value = transform(self.meta_info[n.strain][meta_column])
                if not np.isnan(n.meta_value):
                    values_by_state[n.present].append(n.meta_value)


        n_events = len([n for n in self.tree.find_clades() if n.event])



        if len(values_by_state['present'])>1 and len(values_by_state['absent'])>1:
            # return (np.mean(values_by_state['present']) - np.mean(values_by_state['absent']))\
            #         *np.sqrt(1.0/(1.0/len(values_by_state['present']) + 1.0/len(values_by_state['absent'])))\
            #         /np.std(values_by_state['present']+values_by_state['absent'])
            #return (np.mean(values_by_state['present']) - np.mean(values_by_state['absent']))\
            #        /np.sqrt(np.var(values_by_state['present'])+np.var(values_by_state['present']))
            return (np.mean(values_by_state['present']) - np.mean(values_by_state['absent']))\
                   /np.std(values_by_state['present']+values_by_state['absent'])*np.sqrt(n_events)
        else:
            return np.nan


class BranchAssociation(object):
    """docstring for Association"""
    def __init__(self, tree, meta_info):
        super(BranchAssociation, self).__init__()
        self.meta_info = meta_info
        if type(tree)==str and os.path.isfile(tree):
            self.tree = Phylo.load(tree_name, 'newick')
        else:
            self.tree = tree


    def calc_up_down_averages(self,meta_column, transform=None, pc=3):
        '''
        calculate the mean value of the phenotype of leaves upstream and down stream
        of each branch in the tree.
        '''
        if transform is None:
            transform = lambda x:x
        for n in self.tree.find_clades(order='postorder'):
            for c in n: c.up=n # add up links for convenience

            if n.is_terminal():
                n.strain = n.name.split('|')[0]
                n.meta_value = transform(self.meta_info[n.strain][meta_column])
                if not np.isnan(n.meta_value):
                    n.meta_count = 1
                    n.meta_sq_value = n.meta_value*n.meta_value
                else:
                    n.meta_count = 0
                    n.meta_sq_value = np.nan

            else:
                n.meta_count = np.sum([c.meta_count for c in n if c.meta_count])
                n.meta_value = np.sum([c.meta_value for c in n if c.meta_count])
                n.meta_sq_value = np.sum([c.meta_sq_value for c in n if c.meta_count])


        root_node = self.tree.root
        n = root_node
        n.meta_derived_average = n.meta_value/n.meta_count
        n.meta_derived_var = n.meta_count/(n.meta_count-1.0)\
                        *(n.meta_sq_value/n.meta_count - n.meta_derived_average**2)
        n.meta_derived_SSEM = n.meta_derived_var/n.meta_count
        pseudo_var = self.tree.root.meta_derived_var

        for n in self.tree.find_clades(order='preorder'):
            if n==root_node:
                continue

            # calculate average and standard deviation of meta data of child nodes
            if n.meta_count==0:
                n.meta_derived_average = np.nan
                n.meta_derived_var = np.nan
                n.meta_derived_SSEM = np.inf
            else:
                n.meta_derived_average = n.meta_value/n.meta_count
                if n.meta_count==1:
                    n.meta_derived_var = np.nan
                    n.meta_derived_SSEM = np.inf
                else:
                    n.meta_derived_var = n.meta_count/(n.meta_count-1.0)\
                                *(n.meta_sq_value/n.meta_count - n.meta_derived_average**2)
                    n.meta_derived_SSEM = (n.meta_derived_var+pc*pseudo_var)/n.meta_count

            # calculate average and standard deviation of meta data of all non child nodes
            n_non_child = root_node.meta_count - n.meta_count
            n.meta_ancestral_average = (root_node.meta_value-n.meta_value)/n_non_child
            n.meta_ancestral_var = n_non_child/(n_non_child-1.0)\
                            *((root_node.meta_sq_value - n.meta_sq_value)/n_non_child
                               - n.meta_ancestral_average**2)
            n.meta_ancestral_SSEM = (n.meta_ancestral_var+pc*pseudo_var)/n_non_child


    def calc_significance(self):
        max_score = 0
        for n in self.tree.find_clades():
            if n==self.tree.root:
                n.z_score=np.nan
            else:
                n.z_score = np.abs(n.meta_derived_average - n.meta_ancestral_average)/\
                        np.sqrt(n.meta_ancestral_SSEM + n.meta_derived_SSEM)

                if (not np.isnan(n.z_score)) and n.z_score>max_score:
                    max_score=n.z_score

        return max_score

def infer_branch_associations(path, metainfo_fpath, meta_data_config,
    total_strains_count, strain_fraction_branch_association):
    from sf_geneCluster_align_makeTree import load_sorted_clusters
    from sf_coreTree_json import metadata_load
    data_description = meta_data_config
    association_dict = defaultdict(dict)
    metadata = Metadata(metainfo_fpath, data_description)
    metadata_dict = metadata.to_dict()

    sorted_genelist = load_sorted_clusters(path)
    ## sorted_genelist: [(clusterID, [ count_strains,[memb1,...],count_genes]),...]
    for clusterID, gene in sorted_genelist:
        if gene[-1]>=total_strains_count*strain_fraction_branch_association: # and clusterID=='GC00001136':
            print(clusterID)
            tree = Phylo.read("%s/geneCluster/%s.nwk"%(path, clusterID), 'newick')
            assoc = BranchAssociation(tree, metadata_dict)
            for col, d  in metadata.data_description.iterrows():
                if d['associate']=='yes':
                    if 'log_scale' in d and d['log_scale']=='yes':
                        t = lambda x:np.log(x)
                    else:
                        t = lambda x:x
                    assoc.calc_up_down_averages(d["meta_category"], transform = t)
                    max_assoc = assoc.calc_significance()
                    association_dict[clusterID][d["meta_category"]] = max_assoc

    write_pickle("%s/branch_association.cpk"%path, association_dict)


def load_gain_loss(path, clusterID):
    with open('%s/geneCluster/%s_patterns.json'%(path, clusterID), 'r') as ifile:
        tmp = ifile.readlines()[-1].strip().split(':')[-1].split('"')[-2]
    return map(int, list(tmp))


def infer_presence_absence_associations(path, metainfo_fpath, meta_data_config,
    total_strains_count, min_strain_fraction_association, max_strain_fraction_association):
    from sf_geneCluster_align_makeTree import load_sorted_clusters
    from sf_coreTree_json import metadata_load
    data_description = meta_data_config
    association_dict = defaultdict(dict)
    metadata = Metadata(metainfo_fpath, data_description)
    metadata_dict = metadata.to_dict()
    min_strains_association = total_strains_count*min_strain_fraction_association
    max_strains_association = total_strains_count*max_strain_fraction_association
    sorted_genelist = load_sorted_clusters(path)
    ## sorted_genelist: [(clusterID, [ count_strains,[memb1,...],count_genes]),...]
    # TODO fix vis
    tree = Phylo.read("%sgeneCluster/strain_tree.nwk"%(path), 'newick')
    assoc = PresenceAbsenceAssociation(tree, metadata_dict)
    for clusterID, gene in sorted_genelist:
        if gene[-1]>min_strains_association and gene[-1]<max_strains_association:
            print(clusterID)
            gl = load_gain_loss(path, clusterID)
            for col, d  in metadata.data_description.iterrows():
                if d['associate']=='yes':
                    if 'log_scale' in d and d['log_scale']=='yes':
                        t = lambda x:np.log(x)
                    else:
                        t = lambda x:x
                    assoc.set_gain_loss(gl)
                    score = assoc.calc_association_simple(d["meta_category"], transform = t)
                    if np.isinf(score):
                        association_dict[clusterID][d["meta_category"]] = 0.0
                    else:
                        association_dict[clusterID][d["meta_category"]] = np.abs(score)

    write_pickle("%s/presence_absence_association.cpk"%path, association_dict)

