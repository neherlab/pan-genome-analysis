import os, sys, glob, time, shutil
import numpy as np
from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from SF00_miscellaneous import times, read_fasta, write_in_fa, load_pickle, write_pickle
from SF06_geneCluster_align_makeTree import multips, mpm_tree, prepare_cluster_seq

def calculate_diversity(file_path, files_list):
    """ calculate_diversity
    """
    diversity_dict=defaultdict()
    for input_filepath in files_list:
        try:
            myTree = mpm_tree(input_filepath)
            myTree.codon_align()
            myTree.diversity_statistics_nuc()
            diversity_dict[input_filepath.split('/')[-1]]=round(myTree.diversity_nuc,4)
        except:
            print('problem in sequence diversity calculation:', input_filepath)

    ## write records in gene_diversity file
    with open(file_path+'tmp_core_diversity.txt', 'a') as tmp_core_diversity_file:
        for clusterID, diversity in diversity_dict.iteritems():
            tmp_core_diversity_file.write('%s\t%f\n'%(clusterID, diversity))

def tmp_average_core_diversity(file_path):
    """ calculate mean diversity """
    with open(file_path+'tmp_core_diversity.txt', 'r') as tmp_core_diversity_file:
        diversity_lst=[float(iline.split('\t')[1]) for iline in tmp_core_diversity_file]
    return round(np.mean(diversity_lst),4)

def estimate_core_gene_diversity(path, parallel=32, core_cutoff=1.0):
    """
    estimate core gene diversity before gene cluster alignment
    and cluster post-processing
    """
    ## total number of strains
    strain_list=load_pickle(path+'strain_list.cpk')
    totalStrain= len(strain_list)

    ## load clusters
    cluster_path='%s%s'%(path,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(cluster_path+'allclusters.cpk')

    ## create core gene list
    core_geneCluster_dt= defaultdict()
    # geneCluster_dt: {clusterID:[ count_strains,[memb1,...],count_genes }
    for clusterID, cluster_stats in geneCluster_dt.iteritems():
        if core_cutoff==1.0:
            strain_core_cutoff=totalStrain
        else:
            strain_core_cutoff=int(totalStrain*core_cutoff)
        ## check whether #genes == #strains and it's a core/soft-core gene
        if cluster_stats[0]==cluster_stats[2] and cluster_stats[0]>=strain_core_cutoff:
            core_geneCluster_dt[clusterID]=cluster_stats

    ## remove the folder from previous run
    cluster_seqs_path= '%s%s'%(path,'geneCluster/')
    if os.path.exists(cluster_seqs_path):
        os.system(''.join(['rm -r ',cluster_seqs_path]))
    os.system('mkdir %s'%cluster_seqs_path)

    tmp_core_cluster_seqs_path= '%s%s'%(cluster_seqs_path,'tmp_core/');
    if os.path.exists(tmp_core_cluster_seqs_path):
        os.system(''.join(['rm -r ',tmp_core_cluster_seqs_path]))
    os.system('mkdir %s'%tmp_core_cluster_seqs_path)

    ## load geneID_to_geneSeqID geneSeqID cpk file
    geneID_to_geneSeqID_dict= load_pickle(path+'geneID_to_geneSeqID.cpk')

    ## create dict storing all genes' translation
    faa_path= path+'protein_faa/'
    gene_aa_dict= defaultdict(list)
    for ifaa in glob.glob(faa_path+"*.faa"):
        gene_aa_dict.update(read_fasta(ifaa))
    write_pickle('%s%s'%(cluster_seqs_path,'all_protein_seq.cpk'),\
        gene_aa_dict )

    ## create dict for all gene's nucleotide sequence
    strain_nuc_seq_cpk={}
    nucleotide_dict_path= '%s%s'%(path,'nucleotide_fna/')
    for istrain in strain_list:
        strain_nuc_seq_cpk[istrain]=load_pickle(nucleotide_dict_path+istrain+'_gene_nuc_dict.cpk')
    write_pickle('%s%s'%(cluster_seqs_path,'all_nucleotide_seq.cpk'),\
        strain_nuc_seq_cpk )

    ## write nucleotide and amino-acid sequences for each gene cluster 
    prepare_cluster_seq(tmp_core_cluster_seqs_path, core_geneCluster_dt, \
        strain_nuc_seq_cpk, geneID_to_geneSeqID_dict, gene_aa_dict)

    tmp_fa_files=glob.glob(tmp_core_cluster_seqs_path+"*.fna")
    multips(calculate_diversity, parallel, tmp_core_cluster_seqs_path, tmp_fa_files)

    calculated_core_diversity=tmp_average_core_diversity(tmp_core_cluster_seqs_path)
    #refined_core_diversity= max(0.3, 3*calculated_core_diversity)
    refined_core_diversity= \
        (0.1+3*calculated_core_diversity)/(1+3*calculated_core_diversity)
    print('estimated_core_diversity: ',refined_core_diversity, \
           'calculated_core_diversity: ',calculated_core_diversity)
    ## move folder tmp_core to the central data folder
    new_cluster_seqs_path= '%stmp_core'%path
    if os.path.exists(new_cluster_seqs_path):
        os.system(''.join(['rm -r ',new_cluster_seqs_path]))
    os.system('mv %s %s'%(tmp_core_cluster_seqs_path, path))
    return refined_core_diversity
