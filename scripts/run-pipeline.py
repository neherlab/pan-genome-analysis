import argparse
import os, sys, time, glob
from collections import Counter
from SF00_miscellaneous import times,load_pickle, write_pickle
from SF02_get_acc_single_serverDown import accessionID_single
from SF03_diamond_input import diamond_input
from SF04_gbk_metainfo import gbk_To_Metainfo
from SF05_diamond_orthamcl import diamond_orthamcl_cluster
from SF05_2_blastRNA import RNA_cluster
from SF06_geneCluster_align_makeTree import cluster_align_makeTree, postprocess_paralogs_iterative
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

parser = argparse.ArgumentParser(description='Pipeline for identifying core and pan-genome from NCBI RefSeq genome sequences.')
parser.add_argument('-fn', '--folder_name', type = str, required=True, help='the absolute path for project folder ')
parser.add_argument('-sl', '--strain_list', type = str, required=True, help='the file name for strain list (e.g.: Pa-RefSeq.txt)')
parser.add_argument('-st', '--steps', nargs='+', type = int, default = ['all'], help='select specific steps to run or run all steps by default ')
parser.add_argument('-rt', '--raxml_max_time', type = int, default = 30, help='RAxML tree optimization: maximal runing time (in minutes, default: 30 min)')
parser.add_argument('-t', '--threads', type = int, default = 1, help='number of threads')
parser.add_argument('-bp', '--blast_file_path', type = str, default = 'none', help='the absolute path for blast result (e.g.: /path/blast.out)')
parser.add_argument('-rp', '--roary_file_path', type = str, default = 'none', help='the absolute path for roary result (e.g.: /path/roary.out)')
parser.add_argument('-mi', '--meta_info_file_path', type = str, default = 'none', help='the absolute path for meta_information file (e.g.: /path/meta.out)')
parser.add_argument('-np', '--disable_cluster_postprocessing', type = int, default = 0, help='default: run post-processing procedure for splitting paralogs and clustering un-clustered genes')
parser.add_argument('-nrna', '--disable_RNA_clustering', type = int, default = 1, help='default:  disabled, not cluster rRNAs and tRNAs')
parser.add_argument('-dmt', '--diamond_max_target_seqs', type = str, default = '600', help='Diamond: the maximum number of target sequences per query to keep alignments for. Defalut: #strain * #max_duplication= 40*15= 600 ')
parser.add_argument('-imcl', '--mcl_inflation', type = float, default = 2.0, help='MCL: inflation parameter (varying this parameter affects granularity) ')
parser.add_argument('-bmt', '--blastn_RNA_max_target_seqs', type = str, default = '600', help='Blastn on RNAs: the maximum number of target sequences per query to keep alignments for. Defalut: #strain * #max_duplication= 40*15= 600 ')
parser.add_argument('-ws', '--window_size_smoothed', type = int, default = 5, help='postprocess_unclustered_genes: window_size for smoothed cluster length distribution')
parser.add_argument('-spr', '--strain_proportion', type = float, default = 0.3, help='postprocess_unclustered_genes: strain_proportion')
parser.add_argument('-ss', '--sigma_scale', type = int, default = 3, help='postprocess_unclustered_genes: sigma_scale')
parser.add_argument('-kt', '--keep_temporary_file', type = str, default = True, help='default keep_temporary_file')
parser.add_argument('-cg', '--core_genome_threshold', type = float, default = 1.0, help='percentage of strains used to decide whether a gene is core. Default: 1.0 for strictly core gene; customized instance: 0.9 for soft core genes')

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

def load_strains():
    """ load input strains in strain_list """
    if os.path.isfile(path+strain_list_file):
        with open(path+strain_list_file,'rb') as infile:
            strain_list=[]
            for gbk_filepath in glob.glob(path+'*gbk'):
                gbk_filename=gbk_filepath.split('/')[-1]
                if '-' in gbk_filename:
                    new_gbk_filename=gbk_filename.replace('-','_')
                    #print'mv %s %s%s'%(gbk_filepath,path,new_gbk_filename)
                    print new_gbk_filename.split('.gbk')[0]
                    strain_list.append(new_gbk_filename.split('.gbk')[0])
                    os.system('mv %s %s%s'%(gbk_filepath,path,new_gbk_filename))
                else:
                    strain_list.append(gbk_filename.split('.gbk')[0])
            print strain_list
            write_pickle(path+'strain_list.cpk', strain_list )

if 1 in params.steps: #step 01:
    load_strains()
    print 'step01-refSeq strain list successfully found.'

## load strain_list.cpk file and give the total number of strains
if os.path.isfile(path+'strain_list.cpk'):
    strain_list= load_pickle(path+'strain_list.cpk')
nstrains =len([ istrain for istrain in strain_list ])

if 2 in params.steps:# step02:
    start = time.time()
    accessionID_single(path, strain_list)
    print 'step02-download NCBI refseq GenBank file from strain list:'
    print times(start)

if 3 in params.steps:# step03:
    start = time.time()
    diamond_input(path, strain_list, params.disable_RNA_clustering)
    print 'step03-create input file for Diamond from genBank file (.gb):'
    print times(start)

if 4 in params.steps:# step04:
    start = time.time()
    gbk_To_Metainfo(path)
    print 'step04-extra meta_info (country,date, host...) from genBank file:'
    print times(start)

if 5 in params.steps:# step05:
    start = time.time()
    diamond_orthamcl_cluster(path, params.threads, params.blast_file_path, params.roary_file_path, params.diamond_max_target_seqs, params.mcl_inflation )
    if params.disable_RNA_clustering==0:
        RNA_cluster( path, params.threads, params.blastn_RNA_max_target_seqs, params.mcl_inflation )
    print 'step05-run diamond and ortha-mcl to cluster genes: '
    print times(start)

if 6 in params.steps:# step06:
    start = time.time()
    cluster_align_makeTree(path, params.threads, params.disable_cluster_postprocessing)
    if params.disable_cluster_postprocessing==0:
        postprocess_paralogs_iterative(params.threads, path, nstrains)
        postprocess_unclustered_genes(params.threads, path, nstrains, params.window_size_smoothed, params.strain_proportion, params.sigma_scale )
    if params.disable_RNA_clustering==0:
        RNAclusters_align_makeTree(path, params.threads)
    print 'step06-align genes in geneCluster by mafft and build gene trees:'
    print times(start)

if 7 in params.steps:# step07:
    start = time.time()
    create_core_SNP_matrix(path, params.core_genome_threshold)
    print 'step07-call SNPs from core genes:'
    print times(start)

if 8 in params.steps:# step08:
    start = time.time()
    aln_to_Newick(path, params.raxml_max_time, params.threads)
    print 'step08-run fasttree and raxml for tree construction:'
    print times(start)

if 9 in params.steps:# step09:
    start = time.time()
    make_genepresence_alignment(path)
    process_gain_loss(path)
    print 'step09-infer presence/absence and gain/loss patterns of all genes:'
    print times(start)

if 10 in params.steps:# step10:
    start = time.time()
    geneCluster_to_json(path, params.disable_RNA_clustering)
    print 'step10-creates json file for geneDataTable visualization:'
    print times(start)

if 11 in params.steps:# step11:
    start = time.time()
    json_parser(path, species, params.meta_info_file_path)
    print 'step11-extract json files for tree&treeDataTable visualization,etc:'
    print times(start)

if 12 in params.steps:# step12:
    from SF12_create_panGenome import ete_tree_split
    start = time.time()
    ete_tree_split(path, species)
    print 'step12-create pan-genome:'
    print times(start)
