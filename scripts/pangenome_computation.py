import os, glob
from collections import defaultdict
from cluster_collective_processing import clusterCollector
from sf_miscellaneous import times, write_in_fa, load_pickle, write_pickle
from sf_extract_sequences import extract_sequences
from sf_extract_metadata import extract_metadata
from sf_cluster_protein import clustering_protein
from sf_preclustering import preclustering_protein
from sf_cluster_protein_divide_conquer import clustering_divide_conquer
from sf_cluster_RNA import RNA_cluster
from sf_geneCluster_align_makeTree import write_final_cluster
from sf_RNAcluster_align_makeTree import RNAclusters_align_makeTree
from sf_core_SNP_matrix import create_core_SNP_matrix
from sf_core_tree_build import aln_to_Newick
from sf_gene_presence import make_genepresence_alignment
from sf_gain_loss import process_gain_loss
from sf_geneCluster_json import geneCluster_to_json
from sf_coreTree_json import json_parser

class pangenome:
    """
    pangenome analysis based on diamond and MCL
    This creates orthologous clusters from genes in a collection of reference
    genomes by first finding homologous genes using diamond, clustering those
    using MCL, and post-processing these clusters by building phylogenies using
    fasttree. Finally, alignments, meta data, and annotated phylogenies are
    exported fro visualization in a web browser.
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
            #self.params_dict[k]=v todos
        self.folders_dict=defaultdict( str,
            gbk_path='input_GenBank/',
            nucleotide_path='nucleotide_fna/',
            protein_path='protein_faa/',
            clustering_path='protein_faa/diamond_matches/',
            RNA_path='RNA_fna/',
            cluster_seq_path='geneCluster/',
            tmp_core_seq_path='tmp_core/',
            vis_json_path='vis/',
            vis_cluster_path='vis/geneCluster/',
            log_path='log/')

        # set up folder structure and files names
        self.organize_folders()
        self.specify_filepath()
        if os.path.exists(self.fpaths_dict['strain_cpk']):
            self.strain_list=load_pickle(self.fpaths_dict['strain_cpk'])
            self.nstrains=len(self.strain_list)

    def organize_folders(self):
        """ create folders for pangenome analysis """
        folders_dict=self.folders_dict
        command_mkdir='mkdir -p '
        for k,v in folders_dict.iteritems():
            folders_dict[k]='%s%s'%(self.path,v)
        for key, folder_path in folders_dict.iteritems():
            if not os.path.isdir(folder_path):
                os.system(''.join([command_mkdir,folder_path]))
        self.nucleotide_path=folders_dict['nucleotide_path']
        self.protein_path=folders_dict['protein_path']
        self.clustering_path=folders_dict['clustering_path']


    def specify_filepath(self):
        """
        create names for files written or accesses by different parts of the pipeline.
        Having a central store for all files names helps keeping the code in order.
        Note that self.organize_folders needs to be called first.
        """
        fpaths_dict=defaultdict( str,
            strain_cpk='%s%s'%(self.path,'strain_list.cpk'),
            geneID_to_geneSeqID_path='%s%s'%(self.path,'geneID_to_geneSeqID.cpk'),
            geneID_to_description_path='%s%s'%(self.path,'geneID_to_description.cpk'),
            RNAID_to_SeqID_path='%s%s'%(self.path,'RNAID_to_SeqID.cpk'),
            RNAID_to_description_path='%s%s'%(self.path,'RNAID_to_description.cpk'),
            protein_dict_path='%s%s'%(self.protein_path,'all_protein_seq.cpk'),
            nucleotide_dict_path='%s%s'%(self.nucleotide_path,'all_nucleotide_seq.cpk'),
            cluster_fpath='%s%s'%(self.clustering_path,'allclusters.tsv'),
            cluster_final_fpath='%s%s'%(self.clustering_path,'allclusters_final.tsv'),
            cluster_cpk_fpath='%s%s'%(self.clustering_path,'allclusters.cpk'),
            cluster_cpk_final_fpath='%s%s'%(self.clustering_path,'allclusters_postprocessed.cpk'))
        self.fpaths_dict=fpaths_dict


    def make_strain_list(self):
        """ make strainID list and harmonize input filenames"""
        path=self.path
        folders_dict=self.folders_dict
        ## load input strains from all gbk or fasta files in self.path
        if self.gbk_present==1:
            glob_item='.gbk'
            gbk_path=folders_dict['gbk_path']
            glob_list=glob.glob('%s*%s'%(path,glob_item))
            if len(glob_list)!=0:
                harmonize_filename(path,glob_list)
                strain_list= [i.split('/')[-1].split(glob_item)[0]
                              for i in glob.iglob('%s*%s'%(path,glob_item))]
                ## move gbk files in folder input_GenBank
                command_organize_gbk_input=' '.join(['mv', path+'*gbk',gbk_path])
                os.system(command_organize_gbk_input)
            else:
                gbk_glob=glob.iglob('%s*%s'%(gbk_path,glob_item))
                strain_list= [i.split('/')[-1].split(glob_item)[0] for i in gbk_glob]
        else:
            glob_item='.faa'
            glob_list=glob.glob('%s*%s'%(path,glob_item))
            if len(glob_list)!=0:
                harmonize_filename(path,glob_list)
                strain_list=[i.split('/')[-1].split(glob_item)[0]
                             for i in glob.iglob('%s*%s'%(path,glob_item))]
            else:
                protein_glob=glob.iglob('%s*%s'%(folders_dict['protein_path'],glob_item))
                strain_list= [i.split('/')[-1].split(glob_item)[0] for i in protein_glob]
            command_organize_aa_input= 'mv %s*.faa %s'%(path,folders_dict['protein_path'])
            command_organize_nuc_input='mv %s*.fna %s'%(path,folders_dict['nucleotide_path'])
            os.system(command_organize_nuc_input)
            os.system(command_organize_aa_input)
        # write the list of strains to a pickle file and store the list in self
        write_pickle(self.fpaths_dict['strain_cpk'], strain_list)
        self.strain_list=strain_list
        self.nstrains=len(strain_list)


    def extract_gbk_sequences(self):
        """ extract nucleotide and protein sequences from GenBank file """
        extract_sequences(self.path, self.strain_list, self.folders_dict, self.gbk_present, self.enable_RNA_clustering)
        #gene_aa_dict, gene_na_dict= extract_sequences()

    def extract_gbk_metadata(self):
        """ extract metainfo (collection date, country, etc.) from GenBank file """
        extract_metadata(self.path, self.strain_list, self.folders_dict, self.gbk_present, self.metainfo_reconcile)

    def clustering_protein_sequences(self):
        """ clustering protein sequences"""
        clustering_protein(self.path, self.folders_dict, self.threads, self.blast_fpath, self.roary_fpath,
            self.orthofinder_fpath, self.other_tool_fpath, self.diamond_evalue, self.diamond_max_target_seqs,
            self.diamond_identity, self.diamond_query_cover, self.diamond_subject_cover, self.diamond_path, self.mcl_inflation)

    def finding_gene_copies(self):
        """ running Diamond on each strain and find the duplicates """
        preclustering_protein(self.path, self.folders_dict, self.threads,
                self.diamond_evalue, self.diamond_max_target_seqs,
                self.diamond_identity_subproblem, self.diamond_query_cover_subproblem,
                self.diamond_subject_cover_subproblem)

    def clustering_protein_divide_conquer(self):
        """ applying divide and conquer algorithm to speed up all-against-all protein alignment """
        clustering_divide_conquer(self.path, self.folders_dict, self.threads,
             self.diamond_evalue, self.diamond_max_target_seqs,
            self.diamond_identity_subproblem, self.diamond_query_cover_subproblem,
            self.diamond_subject_cover_subproblem, self.mcl_inflation, self.diamond_path, self.diamond_dc_subset_size)

    def RNA_clustering(self):
        """clustering RNA sequences """
        RNA_cluster(self.path, self.threads, self.blastn_RNA_max_target_seqs, self.mcl_inflation)

    def process_clusters(self):
        """"""
        myClusterCollector= clusterCollector(
            path=self.path, folders_dict=self.folders_dict,
            strain_list=self.strain_list,
            nstrains=self.nstrains, threads=self.threads,
            core_genome_threshold=self.core_genome_threshold,
            factor_core_diversity=self.factor_core_diversity,
            disable_cluster_postprocessing=self.disable_cluster_postprocessing,
            disable_long_branch_splitting=self.disable_long_branch_splitting,
            simple_tree=self.simple_tree,
            split_long_branch_cutoff=self.split_long_branch_cutoff,
            paralog_branch_cutoff=self.paralog_branch_cutoff,
            paralog_frac_cutoff=self.paralog_frac_cutoff,
            explore_paralog_plot=self.explore_paralog_plot,
            window_size_smoothed=self.window_size_smoothed,
            strain_proportion=self.strain_proportion,
            sigma_scale=self.sigma_scale,
            species=self.species
            )

        myClusterCollector.estimate_raw_core_diversity()
        #self.clusterCollector= myClusterCollector.cluster_align_makeTree()
        myClusterCollector.make_geneCluster_alignment_and_tree()
        if not self.disable_cluster_postprocessing:
            if not self.disable_long_branch_splitting:
                myClusterCollector.postprocessing_split_long_branch()
                myClusterCollector.postprocess_merge_underclustered_genes()
            myClusterCollector.postprocessing_split_paralogs()
            write_final_cluster(self.path)

    def make_RNACluster_alignment_and_tree(self):
        """ aligning RNA clusters and building RNA tree """
        RNAclusters_align_makeTree(self.path, self.folders_dict, self.threads, self.simple_tree)

    def create_SNP_alignment(self):
        """ build pseudo-alignment based on core gene SNPs """
        create_core_SNP_matrix(self.path, self.core_genome_threshold, self.core_gene_strain_fpath)

    def build_core_tree(self):
        """ build core tree based on SNP alignment """
        aln_to_Newick(self.path, self.folders_dict, self.raxml_max_time, self.raxml_path, self.threads)

    def compute_gene_presence_pattern(self):
        """
        compute gene presence/absence pattern for each gene cluster
        (displayed on core tree)
        """
        make_genepresence_alignment(self.path, self.disable_gain_loss, self.merged_gain_loss_output)

    def infer_gene_gain_loss_pattern(self):
        """
        infer gene gain/loss pattern for each gene cluster
        (displayed on core tree)
        """
        process_gain_loss(self.path, self.merged_gain_loss_output)

    def inferAssociations(self):
        from sf_association import infer_branch_associations, infer_presence_absence_associations
        infer_branch_associations(self.path, self.metainfo_fpath, self.meta_data_config, self.nstrains, self.min_strain_fraction_branch_association)
        infer_presence_absence_associations(self.path, self.metainfo_fpath, self.meta_data_config, self.nstrains,
            self.min_strain_fraction_presence_association, self.max_strain_fraction_presence_association)
        # TODO: gain loss associations

    def export_geneCluster_json(self):
        """ export gene clusters as json file which is loaded in the cluster datatable """
        geneCluster_to_json(self.path, self.enable_RNA_clustering, self.store_locus_tag, self.raw_locus_tag, self.optional_table_column)

    def export_coreTree_json(self):
        """ export core tree as json file for core tree visualization"""
        json_parser(self.path, self.folders_dict, self.fpaths_dict, self.metainfo_fpath, self.meta_data_config, self.clean_temporary_files)


def harmonize_filename(path, glob_list):
    """ force '-' to be replaced as '_' in input filename """
    for fpath in glob_list:
        gbk_fname=fpath.split('/')[-1]
        if '-' in gbk_fname:
            ## force to replace '-' with '_' in GenBank filename
            gbk_fname= gbk_fname.replace('-','_')
            print ''.join(['filename harmonized: ',fpath,' -> ',gbk_fname])
            os.system(''.join(['mv ',fpath,' ',path,gbk_fname]))
