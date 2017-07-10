from sf_geneCluster_align_makeTree import cluster_align_makeTree
from sf_core_diversity import estimate_core_gene_diversity
from sf_split_long_branch import postprocess_split_long_branch
from sf_split_paralogy import postprocess_paralogs_iterative
from sf_unclustered_genes import postprocess_unclustered_genes

class clusterCollector(object):
    """docstring for clusterCollector"""
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)

    def estimate_raw_core_diversity(self):
        """ computing raw core gene diversity which's refined as post-processing cutoff """
        if self.split_long_branch_cutoff==0.0:
            self.raw_core_diversity, self.split_long_branch_cutoff= estimate_core_gene_diversity(self.path,
                self.folders_dict, self.strain_list, self.threads, self.core_genome_threshold, self.factor_core_diversity, self.species)

    def make_geneCluster_alignment_and_tree(self):
        """ align genes in gene cluster and building gene tree """
        cluster_align_makeTree(self.path, self.folders_dict, self.threads, self.disable_cluster_postprocessing, self.simple_tree)

    def postprocessing_split_long_branch(self):
        """ postprocessing: split long branch via refined core gene diversity """
        postprocess_split_long_branch(self.threads, self.path, self.simple_tree, self.split_long_branch_cutoff)

    def postprocessing_split_paralogs(self):
        """ postprocessing: split paralogs"""
        if self.paralog_branch_cutoff==0.0:
            self.paralog_branch_cutoff=self.split_long_branch_cutoff
        postprocess_paralogs_iterative(self.threads, self.path, self.nstrains,
            self.simple_tree, self.paralog_branch_cutoff, self.disable_long_branch_splitting,
            self.paralog_frac_cutoff, self.explore_paralog_plot)

    def postprocess_merge_underclustered_genes(self):
        """ postprocessing: integrate under_clustered genes """
        postprocess_unclustered_genes(self.threads, self.path, self.nstrains, self.simple_tree,
            self.split_long_branch_cutoff, self.window_size_smoothed, self.strain_proportion, self.sigma_scale)
