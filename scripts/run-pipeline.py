import argparse
import os, sys, time
sys.path.append('./')
sys.path.append('./scripts/')
from pangenome_computation import pangenome
from sf_miscellaneous import times, load_pickle, write_pickle, organize_folder, load_strains
from sf_fetch_refseq import fetch_refseq
from sf_extract_sequences import extract_sequences
from sf_extract_metadata import extract_metadata
from sf_cluster_protein import clustering_protein
from sf_preclustering import preclustering_protein
from sf_cluster_protein_divide_conquer import clustering_divide_conquer
from sf_cluster_orthamcl import diamond_orthamcl_cluster
from sf_cluster_RNA import RNA_cluster
from sf_geneCluster_align_makeTree import cluster_align_makeTree
from sf_core_diversity import estimate_core_gene_diversity
from sf_split_paralogy import postprocess_paralogs_iterative
from sf_split_long_branch import postprocess_split_long_branch
from sf_unclustered_genes import postprocess_unclustered_genes
from sf_RNAcluster_align_makeTree import RNAclusters_align_makeTree
from sf_core_SNP_matrix import create_core_SNP_matrix
from sf_core_tree_build import aln_to_Newick
from sf_gene_presence import make_genepresence_alignment
from sf_gain_loss import process_gain_loss
from sf_geneCluster_json import geneCluster_to_json
from sf_coreTree_json import json_parser
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
parser.add_argument('-op', '--orthofinder_file_path', type = str, default = 'none',
    help='the absolute path for orthofinder result (e.g.: /path/orthofinder.out)' , metavar='')
parser.add_argument('-otp', '--other_tool_fpath', type = str, default = 'none',
    help='the absolute path for result from other orthology inference tool  (e.g.: /path/other_tool.out)' , metavar='')
parser.add_argument('-mi', '--meta_info_file_path', type = str, default = 'none',
    help='the absolute path for meta_information file (e.g.: /path/meta.out)' , metavar='')
parser.add_argument('-dme', '--diamond_evalue', type = str, default = '0.001',
    help='default: e-value threshold below 0.001', metavar='')
parser.add_argument('-dmt', '--diamond_max_target_seqs', type = str, default = '600',
    help='Diamond: the maximum number of target sequences per query to keep alignments for.\
    Calculation: #strain * #max_duplication (40*15= 600)', metavar='')
parser.add_argument('-dmi', '--diamond_identity', type = str, default = '0',
    help='Diamond: sequence identity threshold to report an alignment. Default: empty.\
    When applied to species with low genetic diversity: 70 could be a decent starting point.\
    All alignments with identity below 0.7 of will not be reported, \
    thus also saving computational time. ', metavar='')
parser.add_argument('-dmqc', '--diamond_query_cover', type = str, default = '0',
    help='Diamond: sequence (query) coverage threshold to report an alignment.  Default: empty.\
    When applied to species with low genetic diversity: 70 could be a decent starting point.\
    All alignments with less than 0.7 of query coverage will not be reported, \
    thus also saving computational time. ', metavar='')
parser.add_argument('-dmsc', '--diamond_subject_cover', type = str, default = '0',
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
parser.add_argument('-imcl', '--mcl_inflation', type = float, default = 1.5,
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
    help='default: disabled, not cluster rRNAs', metavar='')
## split tree via breaking up long branches (resolving over-clustering)
#parser.add_argument('-sf', '--split_long_branch_factor', type = float, default = 3.0,
#    help='use (0.1+3.0*core_diversity)/(1+3.0*core_diversity) to decide split_long_branch_cutoff',metavar='')
parser.add_argument('-cb', '--split_long_branch_cutoff', type = float, default = 0.0,
    help='split long branch cutoff provided by user (by default: 0.0 as not given):',metavar='')
## split paralogy
parser.add_argument('-pep', '--explore_paralog_plot', type = int, default = 0,
    help='default: not plot paralog statistics', metavar='')
parser.add_argument('-pfc', '--paralog_frac_cutoff', type = float, default = 0.3,
    help='fraction of strains required for splitting paralogy. Default: 0.3', metavar='')
parser.add_argument('-pbc', '--paralog_branch_cutoff', type = float, default = 0.0,
    help='branch_length cutoff used in paralogy splitting', metavar='')
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
## core gene strain constraint
parser.add_argument('-csf', '--core_gene_strain_fpath', type = str, default = '',
    help='file path for user-provided subset of strains( core genes should be present in all strains in this list)',
    metavar='')
## gene gain loss inference
parser.add_argument('-gl', '--enable_gain_loss', type = int, default = 1,
    help='default: not enable gene gain and loss inference', metavar='')

parser.add_argument('-sitr', '--simple_tree', type = int, default = 0,
    help='default: enable rich tree annotation; setting simple to 1 does not conduct ancestral inference via treetime.', metavar='')

parser.add_argument('-sp', '--species_name', type = str, default = '',
    help='default:', metavar='')
parser.add_argument('-lo', '--large_output', type = int, default = 0,
    help='default: not split gene presence/absence and gain/loss pattern into separate files for each cluster', metavar='')
parser.add_argument('-rl', '--raw_locus_tag', type = int, default = 0,
    help='default: use strain_ID + locus_tag as new locus_tag; if set to 0, use raw locus_tag from GenBank file', metavar='')

parser.add_argument('-kt', '--keep_temporary_file', type = int, default = 1,
    help='default: keep temporary files', metavar='')

params = parser.parse_args()
path = params.folder_name
if path[:-1]!='/':
    path='%s/'%path
if params.steps[0]=='all':
    ## run all steps
    params.steps=range(1,12)

print 'Running panX in main folder: %s'%path
species=params.strain_list.split('-RefSeq')[0]
#species=params.species_name TODO

myPangenome=pangenome(
    path=params.folder_name,
    gbk_present=params.gbk_present,
    threads=params.threads,
    raxml_limit=params.raxml_max_time,
    blast_fpath=params.blast_file_path,
    roary_fpath=params.roary_file_path,
    orthofinder_fpath=params.orthofinder_file_path,
    other_tool_fpath=params.other_tool_fpath,
    metainfo_fpath=params.meta_info_file_path,
    diamond_evalue=params.diamond_evalue,
    diamond_max_target=params.diamond_max_target_seqs,
    diamond_identity=params.diamond_identity,
    diamond_q_cover=params.diamond_query_cover,
    diamond_s_cover=params.diamond_subject_cover,
    diamond_identity_subprob=params.diamond_identity_precluster,
    diamond_q_cover_subprob=params.diamond_query_cover_precluster,
    diamond_s_cover_subprob=params.diamond_subject_cover_precluster,
    diamond_dc_used=params.diamond_divide_conquer,
    diamond_dc_subsize=params.subset_size,
    mcl_inflation=params.mcl_inflation,
    blastn_max_target=params.blastn_RNA_max_target_seqs,
    disable_postprocess=params.disable_cluster_postprocessing,
    disable_RNA_cluster=params.disable_RNA_clustering,
    split_long_branch_cutoff=params.split_long_branch_cutoff,
    paralog_frac_cutoff=params.paralog_frac_cutoff,
    paralog_branch_cutoff=params.paralog_branch_cutoff,
    windowsize_smoothed=params.window_size_smoothed,
    strain_proportion=params.strain_proportion,
    sigma_scale=params.sigma_scale,    
    core_genome_threshold=params.core_genome_threshold,
    core_gene_strain_fpath=params.core_gene_strain_fpath,
    enable_gain_loss=params.enable_gain_loss,
    simple_tree=params.simple_tree,
    large_output=params.large_output,
    raw_locus_tag=params.raw_locus_tag,
    keep_temporary_file=params.keep_temporary_file
    )
folders_dict=myPangenome.organize_folder()
if 1 in params.steps:#step 01:
    #load_strains(path, params.gbk_present,folders_dict)
    folders_dict=myPangenome.organize_folder()
    myPangenome.specify_filepath()
    myPangenome.make_strain_list()
    print '======  step01: strain list successfully loaded'
## load strain_list.cpk file and give the total number of strains
if os.path.isfile(path+'strain_list.cpk'):
    strain_list= load_pickle(path+'strain_list.cpk')
nstrains =len([ istrain for istrain in strain_list ])

## deactivated  (activation:'2' -> 2)
if '2' in params.steps:# step02:
    print '======  starting step02: download NCBI refseq GenBank file from strain list'
    start = time.time()
    fetch_refseq(path, strain_list)
    print '======  time for step02: download NCBI refseq GenBank file from strain list'
    print times(start),'\n'

if 3 in params.steps:# step03:
    print '======  extract sequences from GenBank file'
    start = time.time()
    gene_aa_dict, gene_na_dict= extract_sequences(path, strain_list, folders_dict, params.gbk_present,
        params.disable_RNA_clustering)
    print '======  extract sequences from GenBank file'
    print times(start),'\n'

if 4 in params.steps:# step04:
    print '======  extract metadata from GenBank file'
    start = time.time()
    extract_metadata(path, strain_list, folders_dict, params.gbk_present)
    print '======  extract metadata from GenBank file'
    print times(start),'\n'

if 5 in params.steps:# step05:
    print '======  cluster proteins'
    start = time.time()
    if params.diamond_divide_conquer==1:
        ## clustering with divide_and_conquer method for large dataset
        clustering_divide_conquer(path, folders_dict, params.threads,
             params.diamond_evalue, params.diamond_max_target_seqs,
            params.diamond_identity, params.diamond_query_cover,
            params.diamond_subject_cover, params.mcl_inflation, params.subset_size)
    elif params.orthAgogue_used==0:
        ## pre-clustering
        if 0:
            preclustering_protein(path, folders_dict, params.threads,
                params.diamond_evalue, params.diamond_max_target_seqs,
                params.diamond_identity_precluster, params.diamond_query_cover_precluster,
                params.diamond_subject_cover_precluster)
        ## clustering
        if 1:
            clustering_protein(path, folders_dict, params.threads,
                params.blast_file_path, params.roary_file_path,
                params.orthofinder_file_path,params.other_tool_fpath,
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
    print '======  cluster proteins'
    print times(start),'\n'

if 6 in params.steps:# step06:
    print '======  starting step06: align genes in geneCluster by mafft and build gene trees'
    start = time.time()

    ## core gene diveristy
    if 1:
        if params.split_long_branch_cutoff==0.0:
            params.split_long_branch_cutoff=\
            estimate_core_gene_diversity(path, folders_dict, strain_list, params.threads, params.core_genome_threshold)
    if 1:
        ## align and make tree
        if params.enable_cluster_correl_stats==1:
            cluster_align_makeTree_correl_stats(path, params.threads, params.disable_cluster_postprocessing)
        else:
            cluster_align_makeTree(path, folders_dict, params.threads, params.disable_cluster_postprocessing, params.simple_tree)

    ## with/without post-processing
    if 1 and params.disable_cluster_postprocessing==0:
        if 1:
            postprocess_split_long_branch(params.threads, path, params.simple_tree, params.split_long_branch_cutoff)
        if 1:
            if params.paralog_branch_cutoff==0.0: ## using default setting (core gene diverstiy)
                params.paralog_branch_cutoff=params.split_long_branch_cutoff
            postprocess_paralogs_iterative(params.threads, path, nstrains, params.simple_tree,\
                params.paralog_branch_cutoff, params.paralog_frac_cutoff, params.explore_paralog_plot)
        if 1:
            postprocess_unclustered_genes(params.threads, path, nstrains, params.simple_tree,\
            params.window_size_smoothed, params.strain_proportion, params.sigma_scale )
    if params.disable_RNA_clustering==0:
        RNAclusters_align_makeTree(path, folders_dict, params.threads)
    print '======  time for step06: align genes in geneCluster by mafft and build gene trees'
    print times(start),'\n'

if 7 in params.steps:# step07:
    print '======  starting step07: call SNPs from core genes'
    start = time.time()
    create_core_SNP_matrix(path, params.core_genome_threshold, params.core_gene_strain_fpath)
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
    make_genepresence_alignment(path, params.enable_gain_loss, params.large_output)
    if params.enable_gain_loss==1:
        process_gain_loss(path, params.large_output)
    print '======  time for step09: infer presence/absence and gain/loss patterns of all genes'
    print times(start),'\n'

if 10 in params.steps:# step10:
    print '======  starting step10: create json file for geneDataTable visualization'
    start = time.time()
    geneCluster_to_json(path, params.disable_RNA_clustering, params.large_output, params.raw_locus_tag)
    print '======  time for step10: create json file for geneDataTable visualization'
    print times(start),'\n'

if 11 in params.steps:# step11:
    print '======  starting step11: extract json files for tree and treeDataTable visualization,etc'
    start = time.time()
    json_parser(path, species, params.meta_info_file_path, params.large_output)
    print '======  time for step11: extract json files for tree and treeDataTable visualization,etc'
    print times(start),'\n'

if 12 in params.steps:# step12:
    from SF12_create_panGenome import ete_tree_split
    print '======  starting step12: create pan-genome'
    start = time.time()
    ete_tree_split(path, species)
    print '======  time for step12: create pan-genome'
    print times(start),'\n'
