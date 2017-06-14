#!/usr/bin/env python
import argparse
import os, sys, time
panX_script_dir= os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,panX_script_dir+'/scripts/')
from pangenome_computation import pangenome
from sf_miscellaneous import times
'''
panX clusters genes from a set of microbial genomes into orthologous clusters
and exports alignments, phylogenies, and meta data for interactive web visualization
and powerful pan-genomic data exploration.

USAGE EXAMPLE: ./panX -fn data/TestSet -sl TestSet-RefSeq.txt 1>TestSet.log 2>TestSet.err
For help, type:
    ./panX -h
'''
parser = argparse.ArgumentParser(description=\
    'panX: Software for computing core and pan-genome from a set of genome sequences.'
    ' The results will be exported as json files for visualization in the browser.',
    usage='./%(prog)s'+' -h (help)')
parser.add_argument('-fn', '--folder_name', type = str, required=True,
    help='the absolute path for project folder ', metavar='')
parser.add_argument('-sl', '--strain_list', type = str, required=True,
    help='the file name for strain list (e.g.: Pa-RefSeq.txt)', metavar='')
parser.add_argument('-ngbk', '--gbk_present', action='store_false',
    help='use nucleotide/amino acid sequence files (fna/faa) when no genBank files given (this option does not consider annotations)')
parser.add_argument('-st', '--steps', nargs='+', type = int, default = ['all'],
    help='run specific steps or run all steps by default', metavar='')
parser.add_argument('-mo', '--metainfo_organism', action='store_true',
    help='add organism information in metadata table.')
parser.add_argument('-rt', '--raxml_max_time', type = int, default = 30,
    help='RAxML tree optimization: maximal runing time (minutes, default:30min)' , metavar='')
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
parser.add_argument('-mi', '--metainfo_fpath', type = str, default = 'none',
    help='the absolute path for meta_information file (e.g.: /path/meta.out)' , metavar='')
parser.add_argument('-dmp', '--diamond_path', type = str, default = '',
    help='alternative diamond path provided by user', metavar='')
parser.add_argument('-dme', '--diamond_evalue', type = str, default = '0.001',
    help='default: e-value threshold below 0.001', metavar='')
parser.add_argument('-dmt', '--diamond_max_target_seqs', type = str, default = '600',
    help='Diamond: maximum number of target sequences per query\
    Calculation: #strain * #max_duplication (40*15= 600)', metavar='')
parser.add_argument('-dmi', '--diamond_identity', type = str, default = '0',
    help='Diamond: sequence identity threshold to report an alignment. Default: empty.', metavar='')
parser.add_argument('-dmqc', '--diamond_query_cover', type = str, default = '0',
    help='Diamond: sequence (query) coverage threshold to report an alignment.  Default: empty', metavar='')
parser.add_argument('-dmsc', '--diamond_subject_cover', type = str, default = '0',
    help='Diamond: sequence (subject) coverage threshold to report an alignment.  Default: empty', metavar='')
parser.add_argument('-dmip', '--diamond_identity_subproblem', type = str, default = '90',
    help='Diamond: sequence identity threshold to report an alignment.', metavar='')
parser.add_argument('-dmqcp', '--diamond_query_cover_subproblem', type = str, default = '90',
    help='Diamond: sequence (query) coverage threshold to report an alignment', metavar='')
parser.add_argument('-dmscp', '--diamond_subject_cover_subproblem', type = str, default = '90',
    help='Diamond: sequence (subject) coverage threshold to report an alignment', metavar='')
parser.add_argument('-dmdc', '--diamond_divide_conquer', action='store_true',
    help='running diamond alignment in divide-and-conquer(DC) algorithm for large dataset')
parser.add_argument('-dcs', '--subset_size', type = int, default = 50,
    help='subset_size (number of strains in a subset) for divide-and-conquer(DC) algorithm. Default:50',\
    metavar='')
parser.add_argument('-imcl', '--mcl_inflation', type = float, default = 1.5,
    help='MCL: inflation parameter (this parameter affects granularity) ', metavar='')
parser.add_argument('-bmt', '--blastn_RNA_max_target_seqs', type = str, default = '100',
    help='Blastn on RNAs: the maximum number of target sequences per query\
    Calculation: #strain * #max_duplication', metavar='')

#/*=======================================
#            post-processing
#=======================================*/
parser.add_argument('-np', '--disable_cluster_postprocessing', action='store_true',
    help='disable postprocessing (split overclustered genes and paralogs, and cluster unclustered genes)')
parser.add_argument('-rna', '--enable_RNA_clustering', action='store_true',
    help='cluster rRNAs')
## split tree via breaking up long branches (resolving over-clustering)
#parser.add_argument('-sf', '--split_long_branch_factor', type = float, default = 3.0,
#    help='use (0.1+3.0*core_diversity)/(1+3.0*core_diversity) to decide split_long_branch_cutoff',metavar='')
parser.add_argument('-fcd', '--factor_core_diversity', type = float, default = 2.0,
    help='default: factor used to refine raw core genome diversity, \
    apply (0.1+2.0*core_diversity)/(1+2.0*core_diversity) to decide split_long_branch_cutoff', metavar='')
parser.add_argument('-slb', '--split_long_branch_cutoff', type = float, default = 0.0,
    help='split long branch cutoff provided by user (by default: 0.0 as not given):',metavar='')
## split paralogy
parser.add_argument('-pep', '--explore_paralog_plot', action='store_true',
    help='default: not plot paralog statistics')
parser.add_argument('-pfc', '--paralog_frac_cutoff', type = float, default = 0.33,
    help='fraction of strains required for splitting paralogy. Default: 0.33', metavar='')
parser.add_argument('-pbc', '--paralog_branch_cutoff', type = float, default = 0.0,
    help='branch_length cutoff used in paralogy splitting', metavar='')
## resolve peaks (unclustered records)
parser.add_argument('-ws', '--window_size_smoothed', type = int, default = 5,
    help='postprocess_unclustered_genes: window size for smoothed cluster length distribution',
    metavar='')
parser.add_argument('-spr', '--strain_proportion', type = float, default = 0.3,
    help='postprocess_unclustered_genes: strain proportion', metavar='')
parser.add_argument('-ss', '--sigma_scale', type = int, default = 3,
    help='postprocess_unclustered_genes: sigma scale', metavar='')

## core genome cutoff
parser.add_argument('-cg', '--core_genome_threshold', type = float, default = 1.0,
    help='percentage of strains used to decide whether a gene is core.\
    Default: 1.0 for strictly core gene; < 1.0 for soft core genes',
    metavar='')
## core gene strain constraint
parser.add_argument('-csf', '--core_gene_strain_fpath', type = str, default = '',
    help='file path for user-provided subset of strains (core genes should be present in all strains in this list)',
    metavar='')
parser.add_argument('-sitr', '--simple_tree', action='store_true',
    help='simple tree: does not use treetime for ancestral inference')
## gene gain loss inference
parser.add_argument('-dgl', '--disable_gain_loss', action='store_true',
    help='disable enable gene gain and loss inference (not recommended)')
parser.add_argument('-mglo', '--merged_gain_loss_output', action='store_true',
    help='not split gene presence/absence and gain/loss pattern into separate files for each cluster')
## branch association inference
parser.add_argument('-iba', '--infer_branch_association', action='store_true',
    help='infer branch association')
## other options
parser.add_argument('-slt', '--store_locus_tag', action='store_true',
    help='store locus_tags in a separate file instead of saving locus_tags in gene cluster json for large dataset')
parser.add_argument('-rlt', '--raw_locus_tag', action='store_true',
    help='use raw locus_tag from GenBank instead of strain_ID + locus_tag')
parser.add_argument('-otc', '--optional_table_column', action='store_true',
    help='add customized column in gene cluster json file for visualization.')
parser.add_argument('-mtf', '--meta_tidy_fpath', type = str, default = '',
    help='file path for pre-defined metadata structure (discrete/continuous data type, etc.)', metavar='')
parser.add_argument('-rxm', '--raxml_path', type = str, default = '',
    help='absolute path of raxml', metavar='')
parser.add_argument('-ct', '--clean_temporary_files', action='store_true',
    help='default: keep temporary files')

params = parser.parse_args()
path = os.path.abspath(params.folder_name)+'/'
if params.steps[0]=='all':
    ## run all steps
    params.steps=range(1,12)

print 'Running panX in main folder: %s'%path
#species=params.species_name TODO

myPangenome=pangenome(
    path=path,
    species=params.strain_list.split('-RefSeq')[0],
    gbk_present=params.gbk_present,
    threads=params.threads,
    metainfo_organism=params.metainfo_organism,
    raxml_max_time=params.raxml_max_time,
    blast_fpath=params.blast_file_path,
    roary_fpath=params.roary_file_path,
    orthofinder_fpath=params.orthofinder_file_path,
    other_tool_fpath=params.other_tool_fpath,
    metainfo_fpath=params.metainfo_fpath,
    diamond_path=params.diamond_path,
    diamond_evalue=params.diamond_evalue,
    diamond_max_target_seqs=params.diamond_max_target_seqs,
    diamond_identity=params.diamond_identity,
    diamond_query_cover=params.diamond_query_cover,
    diamond_subject_cover=params.diamond_subject_cover,
    diamond_identity_subproblem=params.diamond_identity_subproblem,
    diamond_query_cover_subproblem=params.diamond_query_cover_subproblem,
    diamond_subject_cover_subproblem=params.diamond_subject_cover_subproblem,
    diamond_dc_used=params.diamond_divide_conquer,
    diamond_dc_subset_size=params.subset_size,
    mcl_inflation=params.mcl_inflation,
    blastn_RNA_max_target_seqs=params.blastn_RNA_max_target_seqs,
    disable_cluster_postprocessing=params.disable_cluster_postprocessing,
    enable_RNA_clustering=params.enable_RNA_clustering,
    factor_core_diversity=params.factor_core_diversity,
    split_long_branch_cutoff=params.split_long_branch_cutoff,
    explore_paralog_plot=params.explore_paralog_plot,
    paralog_frac_cutoff=params.paralog_frac_cutoff,
    paralog_branch_cutoff=params.paralog_branch_cutoff,
    window_size_smoothed=params.window_size_smoothed,
    strain_proportion=params.strain_proportion,
    sigma_scale=params.sigma_scale,
    core_genome_threshold=params.core_genome_threshold,
    core_gene_strain_fpath=params.core_gene_strain_fpath,
    disable_gain_loss=params.disable_gain_loss,
    simple_tree=params.simple_tree,
    merged_gain_loss_output=params.merged_gain_loss_output,
    store_locus_tag=params.store_locus_tag,
    raw_locus_tag=params.raw_locus_tag,
    meta_tidy_fpath=params.meta_tidy_fpath,
    raxml_path=params.raxml_path,
    optional_table_column=params.optional_table_column,
    clean_temporary_files=params.clean_temporary_files
    )

if 1 in params.steps:#step 01:
    myPangenome.make_strain_list()
    print '======  step01: strain list successfully loaded'

## deactivated  (activation:'2' -> 2)
if '2' in params.steps:# step02:
    print '======  starting step02: download NCBI RefSeq GenBank file from strain list'
    start = time.time()
    fetch_refseq(path, strain_list)
    print '======  time for step02: download NCBI RefSeq GenBank file from strain list'
    print times(start),'\n'

if 3 in params.steps:# step03:
    print '======  extract sequences from GenBank file'
    start = time.time()
    myPangenome.extract_gbk_sequences()
    print '======  extract sequences from GenBank file'
    print times(start),'\n'

if 4 in params.steps:# step04:
    print '======  extract metadata from GenBank file'
    start = time.time()
    myPangenome.extract_gbk_metadata()
    print '======  extract metadata from GenBank file'
    print times(start),'\n'

if 5 in params.steps:# step05:
    print '======  cluster proteins'
    start = time.time()
    if params.diamond_divide_conquer:
        ## clustering with divide_and_conquer method for large dataset
        myPangenome.clustering_protein_divide_conquer()
    else:
        if 0:## pre-clustering
            myPangenome.finding_gene_copies()
        ## clustering
        myPangenome.clustering_protein_sequences()
    ## clustering RNA when option activated
    if params.enable_RNA_clustering:
        myPangenome.RNA_clustering()
    print '======  cluster proteins'
    print times(start),'\n'

if 6 in params.steps:# step06:
    print '======  starting step06: align genes in geneCluster by mafft and build gene trees'
    start = time.time()
    myPangenome.process_clusters()

    if params.enable_RNA_clustering:
        myPangenome.make_RNACluster_alignment_and_tree()
    print '======  time for step06: align genes in geneCluster by mafft and build gene trees'
    print times(start),'\n'

if 7 in params.steps:# step07:
    print '======  starting step07: call SNPs from core genes'
    start = time.time()
    myPangenome.create_SNP_alignment()
    print '======  time for step07: call SNPs from core genes'
    print times(start),'\n'

if 8 in params.steps:# step08:
    print '======  starting step08: run fasttree and raxml for tree construction'
    start = time.time()
    myPangenome.build_core_tree()
    print '======  time for step08: run fasttree and raxml for tree construction'
    print times(start),'\n'

if 9 in params.steps:# step09:
    print '======  starting step09: infer presence/absence and gain/loss patterns of all genes'
    start = time.time()
    myPangenome.compute_gene_presence_pattern()
    if not params.disable_gain_loss:
        myPangenome.infer_gene_gain_loss_pattern()
    if params.infer_branch_association:
        myPangenome.inferAssociations()

    print '======  time for step09: infer presence/absence and gain/loss patterns of all genes'
    print times(start),'\n'

if 10 in params.steps:# step10:
    print '======  starting step10: create json file for geneDataTable visualization'
    start = time.time()
    myPangenome.export_geneCluster_json()
    print '======  time for step10: create json file for geneDataTable visualization'
    print times(start),'\n'

if 11 in params.steps:# step11:
    print '======  starting step11: extract json files for tree and treeDataTable visualization,etc'
    start = time.time()
    myPangenome.export_coreTree_json()
    print '======  time for step11: extract json files for tree and treeDataTable visualization,etc'
    print times(start),'\n'

if 12 in params.steps:# step12:
    from SF12_create_panGenome import ete_tree_split
    print '======  starting step12: create pan-genome'
    start = time.time()
    ete_tree_split(path, species)
    print '======  time for step12: create pan-genome'
    print times(start),'\n'
