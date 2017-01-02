import argparse
import os, sys, time, glob
from collections import Counter
from SF00_miscellaneous import times,load_pickle, write_pickle, load_strains
from SF02_get_acc_single_serverDown import accessionID_single
from SF03_diamond_input import diamond_input
from SF04_gbk_metainfo import gbk_To_Metainfo
from SF05_cluster_protein import clustering_protein_sequences
from SF05_preclustering import preclustering_protein_sequences
from SF05_cluster_protein_divide_conquer import clustering_divide_conquer
from SF05_diamond_orthamcl import diamond_orthamcl_cluster
from SF05_2_blastRNA import RNA_cluster
from SF06_x_geneCluster_correl_stats import cluster_align_makeTree_correl_stats 
from SF06_geneCluster_align_makeTree import cluster_align_makeTree, postprocess_paralogs_iterative
from SF06_0_core_diversity import estimate_core_gene_diversity
from SF06_1_split_overclusters import postprocess_split_overclusters
from SF06_2_unclustered_genes import postprocess_unclustered_genes
from SF06_3_clusterRNA import RNAclusters_align_makeTree
from SF07_core_SNP_matrix import create_core_SNP_matrix
from SF08_core_tree_build import aln_to_Newick
from SF09_1_gene_presence import make_genepresence_alignment
from SF09_2_gain_loss import process_gain_loss
from SF10_geneCluster_export import geneCluster_to_json
from SF11_tree_metadata_export import json_parser
#command line example
#python run-pipeline.py -fn data/Pat3 -sl Pat3-RefSeq.txt -st 1 3 4 5 6 7 8 9 10 > Pat3.log 2>&1

#python run-pipeline.py -fn data/Pmnb -sl Pat3-RefSeq.txt -st 5 -bp Pm_clustering/Pm_allclusters.tab > Pmnb-5.log 2>&1

parser = argparse.ArgumentParser(description=\
    'Software for computing core and pan-genome from a set of genome sequences.',
    usage='%(prog)s')
parser.add_argument('-fn', '--folder_name', type = str, required=True,
    help='the absolute path for project folder ', metavar='')
parser.add_argument('-sl', '--strain_list', type = str, required=True,
    help='the file name for strain list (e.g.: Pa-RefSeq.txt)', metavar='')
parser.add_argument('-gbk', '--gbk_present', type = int, default = 1,
    help=' Using GenBank files by default. \
    Otherwise, nucleotide/amino_acid sequence fna/faa files should be provided.', metavar='')
parser.add_argument('-st', '--steps', nargs='+', type = int, default = ['all'],
    help='select specific steps to run or run all steps by default ', metavar='')
parser.add_argument('-rt', '--raxml_max_time', type = int, default = 30,
    help='RAxML tree optimization: maximal runing time (in minutes, default: 30 min)' , metavar='')
parser.add_argument('-t', '--threads', type = int, default = 1,
    help='number of threads', metavar='')

#/*==================================
#            clustering            
#==================================*/
parser.add_argument('-bp', '--blast_file_path', type = str, default = 'none',
    help='the absolute path for blast result (e.g.: /path/blast.out)' , metavar='')
parser.add_argument('-rp', '--roary_file_path', type = str, default = 'none',
    help='the absolute path for roary result (e.g.: /path/roary.out)' , metavar='')
parser.add_argument('-mi', '--meta_info_file_path', type = str, default = 'none',
    help='the absolute path for meta_information file (e.g.: /path/meta.out)' , metavar='')
parser.add_argument('-dme', '--diamond_evalue', type = str, default = '0.00001',
    help='default: e-value threshold below 0.001', metavar='')
parser.add_argument('-dmt', '--diamond_max_target_seqs', type = str, default = '600',
    help='Diamond: the maximum number of target sequences per query to keep alignments for.\
    Calculation: #strain * #max_duplication (40*15= 600)', metavar='')
parser.add_argument('-dmi', '--diamond_identity', type = str, default = '30',
    help='Diamond: sequence identity threshold to report an alignment. Default: empty.\
    When applied to species with low genetic diversity: 70 could be a decent starting point.\
    All alignments with identity below 0.7 of will not be reported, \
    thus also saving computational time. ', metavar='')
parser.add_argument('-dmqc', '--diamond_query_cover', type = str, default = '30',
    help='Diamond: sequence (query) coverage threshold to report an alignment.  Default: empty.\
    When applied to species with low genetic diversity: 70 could be a decent starting point.\
    All alignments with less than 0.7 of query coverage will not be reported, \
    thus also saving computational time. ', metavar='')
parser.add_argument('-dmsc', '--diamond_subject_cover', type = str, default = '30',
    help='Diamond: sequence (subject) coverage threshold to report an alignment.  Default: empty.\
    When applied to species with low genetic diversity: 70 could be a decent starting point.\
    All alignments with less than 0.7 of subject coverage will not be reported, \
    thus also saving computational time. ', metavar='')
parser.add_argument('-dmip', '--diamond_identity_precluster', type = str, default = '90',
    help='Diamond: sequence identity threshold to report an alignment.', metavar='')
parser.add_argument('-dmqcp', '--diamond_query_cover_precluster', type = str, default = '90',
    help='Diamond: sequence (query) coverage threshold to report an alignment.', metavar='')
parser.add_argument('-dmscp', '--diamond_subject_cover_precluster', type = str, default = '90',
    help='Diamond: sequence (subject) coverage threshold to report an alignment.', metavar='')
parser.add_argument('-dmdc', '--diamond_divide_conquer', type = int, default = 0,
    help='Default: running diamond alignment in divide-and-conquer(DC) method is not activated;\
    this option is designed for large dataset.', metavar='')
parser.add_argument('-dcs', '--subset_size', type = int, default = 50,
    help='Default:subset_size (number of strains in a subset) for divide-and-conquer(DC) method',\
    metavar='')
parser.add_argument('-imcl', '--mcl_inflation', type = float, default = 2.0,
    help='MCL: inflation parameter (varying this parameter affects granularity) ', metavar='')
parser.add_argument('-ortha', '--orthAgogue_used', type = int, default = 0,
    help='Default: orthAgogue not used', metavar='')
parser.add_argument('-bmt', '--blastn_RNA_max_target_seqs', type = str, default = '100',
    help='Blastn on RNAs: the maximum number of target sequences per query to keep alignments for.\
    Calculation: #strain * #max_duplication', metavar='')
parser.add_argument('-cst', '--enable_cluster_correl_stats', type = int, default = 0,
    help='default: calculate statistics on each cluster for correlation test', metavar='')

#/*=======================================
#            post-processing            
#=======================================*/
parser.add_argument('-np', '--disable_cluster_postprocessing', type = int, default = 0,
    help='default: run post-processing procedure for splitting overclustered genes and paralogs,\
    and clustering un-clustered genes', metavar='')
parser.add_argument('-nrna', '--disable_RNA_clustering', type = int, default = 1,
    help='default:  disabled, not cluster rRNAs and tRNAs', metavar='')
## split tree via breaking up long branches (resolving over-clustering)
parser.add_argument('-cbc', '--cut_branch_threshold_customized', type = str, default = False,
    help='[Resolving over-clustering]\
    default: use cut_tree branch_length threshold calculated from core gene diversity;\
    if this is set to True: use the threshold provided by user', metavar='')
parser.add_argument('-cb', '--cut_branch_threshold', type = float, default = 0.3,
    help='[Resolving over-clustering]\
    branch_length threshold to indicate whether to cut sub-trees:',metavar='')
## split paralogy
parser.add_argument('-pep', '--explore_paralog_plot', type = int, default = 0,
    help='default: not plot paralog statistics', metavar='')
parser.add_argument('-pc', '--paralog_cutoff', type = float, default = 0.3,
    help='Fraction of strains required for splitting paralogy. Default: 0.3', metavar='')
parser.add_argument('-pbc', '--paralog_branch_length_cutoff', type = float, default = 1,
    help='Branch_length cutoff used in paralogy splitting', metavar='')
## resolve peaks (unclustered records)
parser.add_argument('-ws', '--window_size_smoothed', type = int, default = 5,
    help='postprocess_unclustered_genes: window_size for smoothed cluster length distribution',
    metavar='')
parser.add_argument('-spr', '--strain_proportion', type = float, default = 0.3,
    help='postprocess_unclustered_genes: strain_proportion', metavar='')
parser.add_argument('-ss', '--sigma_scale', type = int, default = 3,
    help='postprocess_unclustered_genes: sigma_scale', metavar='')

## core genome cutoff
parser.add_argument('-cg', '--core_genome_threshold', type = float, default = 1.0,
    help='percentage of strains used to decide whether a gene is core.\
    Default: 1.0 for strictly core gene; customized instance: 0.9 for soft core genes',
    metavar='')

parser.add_argument('-kt', '--keep_temporary_file', type = str, default = True,
    help='default keep_temporary_file', metavar='')

params = parser.parse_args()
path = params.folder_name
strain_list_file= params.strain_list
## run all steps
if params.steps[0]=='all':
    params.steps=range(1,12)

if path[:-1]!='/':
    path=path+'/'
print path
species=strain_list_file.split('-RefSeq')[0]

if 1 in params.steps: #step 01:
    load_strains(path, params.gbk_present)
    print '======  step01: refSeq strain list successfully found.'

## load strain_list.cpk file and give the total number of strains
if os.path.isfile(path+'strain_list.cpk'):
    strain_list= load_pickle(path+'strain_list.cpk')
nstrains =len([ istrain for istrain in strain_list ])

## deactivated  (activation:'2' -> 2)
if '2' in params.steps:# step02:
    print '======  starting step02: download NCBI refseq GenBank file from strain list'
    start = time.time()
    accessionID_single(path, strain_list)
    print '======  time for step02: download NCBI refseq GenBank file from strain list'
    print times(start),'\n'

if 3 in params.steps:# step03:
    print '======  starting step03: create input file for Diamond from GenBank file (.gb)'
    start = time.time()
    diamond_input(path, strain_list,
        params.gbk_present, params.disable_RNA_clustering)
    print '======  time for step03: create input file for Diamond from GenBank file (.gb)'
    print times(start),'\n'

if 4 in params.steps:# step04:
    print '======  starting step04: extra meta_info (country,date, host...) from GenBank file'
    start = time.time()
    gbk_To_Metainfo(path, params.gbk_present)
    print '======  time for step04: extra meta_info (country,date, host...) from GenBank file'
    print times(start),'\n'

if 5 in params.steps:# step05:
    print '======  starting step05: cluster genes'
    start = time.time()
    if params.diamond_divide_conquer==1:
        ## clustering with divide_and_conquer method for large dataset
        clustering_divide_conquer(path, params.threads,
             params.diamond_evalue, params.diamond_max_target_seqs,
            params.diamond_identity, params.diamond_query_cover,
            params.diamond_subject_cover, params.mcl_inflation, params.subset_size)
    elif params.orthAgogue_used==0:
        ## pre-clustering
        if 0:
            preclustering_protein_sequences(path, params.threads,
                params.diamond_evalue, params.diamond_max_target_seqs,
                params.diamond_identity_precluster, params.diamond_query_cover_precluster, params.diamond_subject_cover_precluster)
        ## clustering
        if 1:
            clustering_protein_sequences(path, params.threads,
                params.blast_file_path, params.roary_file_path,
                params.diamond_evalue, params.diamond_max_target_seqs,
                params.diamond_identity, params.diamond_query_cover, params.diamond_subject_cover,
                params.mcl_inflation
                )
    else:
        ## test ortha+mcl
        diamond_orthamcl_cluster(path, params.threads,
            params.blast_file_path, params.roary_file_path,
            params.diamond_evalue, params.diamond_max_target_seqs,
            params.diamond_identity, params.diamond_query_cover, params.diamond_subject_cover,
            params.mcl_inflation)
    ## clustering RNA when option activated
    if params.disable_RNA_clustering==0:
        RNA_cluster(path, params.threads, params.blastn_RNA_max_target_seqs, params.mcl_inflation)
    print '======  time for step05: cluster genes'
    print times(start),'\n'

if 6 in params.steps:# step06:
    print '======  starting step06: align genes in geneCluster by mafft and build gene trees'
    start = time.time()

    ## core gene diveristy
    # tmp (should be uncommented)
    #cut_branch_threshold= params.cut_branch_threshold
    if 1: #or params.cut_branch_threshold_customized==False:
        cut_branch_threshold=\
            estimate_core_gene_diversity(path, params.threads, params.core_genome_threshold)
        #cut_branch_threshold= params.cut_branch_threshold
    #else:
    #    cut_branch_threshold= params.cut_branch_threshold

    ## align and make tree
    if params.enable_cluster_correl_stats==1:
        cluster_align_makeTree_correl_stats(path, params.threads, params.disable_cluster_postprocessing)
    else:
        cluster_align_makeTree(path, params.threads, params.disable_cluster_postprocessing)

    ## with/without post-processing
    if params.disable_cluster_postprocessing==0:
        postprocess_split_overclusters(params.threads, path, cut_branch_threshold)
        if 1:
            postprocess_paralogs_iterative(params.threads, path, nstrains,\
                params.paralog_cutoff, params.paralog_branch_length_cutoff,\
                params.explore_paralog_plot)
        postprocess_unclustered_genes(params.threads, path, nstrains,\
            params.window_size_smoothed, params.strain_proportion, params.sigma_scale )
    if params.disable_RNA_clustering==0:
        RNAclusters_align_makeTree(path, params.threads)
    print '======  time for step06: align genes in geneCluster by mafft and build gene trees'
    print times(start),'\n'

if 7 in params.steps:# step07:
    print '======  starting step07: call SNPs from core genes'
    start = time.time()
    create_core_SNP_matrix(path, params.core_genome_threshold)
    print '======  time for step07: call SNPs from core genes'
    print times(start),'\n'

if 8 in params.steps:# step08:
    print '======  starting step08: run fasttree and raxml for tree construction'
    start = time.time()
    aln_to_Newick(path, params.raxml_max_time, params.threads)
    print '======  time for step08: run fasttree and raxml for tree construction'
    print times(start),'\n'

if 9 in params.steps:# step09:
    print '======  starting step09: infer presence/absence and gain/loss patterns of all genes'
    start = time.time()
    make_genepresence_alignment(path)
    process_gain_loss(path)
    print '======  time for step09: infer presence/absence and gain/loss patterns of all genes'
    print times(start),'\n'

if 10 in params.steps:# step10:
    print '======  starting step10: create json file for geneDataTable visualization'
    start = time.time()
    geneCluster_to_json(path, params.disable_RNA_clustering)
    print '======  time for step10: create json file for geneDataTable visualization'
    print times(start),'\n'

if 11 in params.steps:# step11:
    print '======  starting step11: extract json files for tree and treeDataTable visualization,etc'
    start = time.time()
    json_parser(path, species, params.meta_info_file_path)
    print '======  time for step11: extract json files for tree and treeDataTable visualization,etc'
    print times(start),'\n'

if 12 in params.steps:# step12:
    from SF12_create_panGenome import ete_tree_split
    print '======  starting step12: create pan-genome'
    start = time.time()
    ete_tree_split(path, species)
    print '======  time for step12: create pan-genome'
    print times(start),'\n'
