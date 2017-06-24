import os, sys, glob, time, shutil
import numpy as np
from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from sf_miscellaneous import times, read_fasta, write_in_fa, load_pickle, write_pickle, multips
from sf_geneCluster_align_makeTree import mpm_tree#, prepare_cluster_seq

def export_cluster_seq_tmp(cluster_seqs_path, geneCluster_dt,
    geneID_to_geneSeqID_dict, gene_na_dict, gene_aa_dict):
    """ write nuc/aa sequences for each cluster  """
    for clusterID, gene in geneCluster_dt.iteritems():
        ## geneCluster file name
        gene_cluster_nu_filename="%s%s"%(clusterID,'.fna')
        with open( cluster_seqs_path+gene_cluster_nu_filename, 'wb') as gene_cluster_nu_write:
            ## write nucleotide sequences into geneCluster files
            for gene_memb in gene[1]:
                ## gene_name format: strain_1|locusTag
                strain_name= gene_memb.split('|')[0]
                geneSeqID=geneID_to_geneSeqID_dict[gene_memb]
                write_in_fa(gene_cluster_nu_write, geneSeqID, gene_na_dict[strain_name][gene_memb] )

def calculate_diversity( files_list, file_path, species):
    """ calculate_diversity
    """
    diversity_dict=defaultdict()
    for input_filepath in files_list:
        try:
            myTree = mpm_tree(input_filepath, speciesID=species)
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

def estimate_core_gene_diversity(path, folders_dict, strain_list, parallel, core_cutoff, factor_core_diversity, species):
    """
    estimate core gene diversity before gene cluster alignment
    and cluster post-processing
    """
    totalStrain= len(strain_list)

    ## load clusters
    clustering_path= folders_dict['clustering_path']
    geneCluster_dt= load_pickle(clustering_path+'allclusters.cpk')
    protein_path= folders_dict['protein_path']
    nucleotide_path= folders_dict['nucleotide_path']
    protein_dict_path= '%s%s'%(protein_path,'all_protein_seq.cpk')
    nucleotide_dict_path= '%s%s'%(nucleotide_path,'all_nucleotide_seq.cpk')
    tmp_core_seq_path= '%s%s'%(clustering_path,'tmp_core/')
    ## load geneID_to_geneSeqID geneSeqID cpk file
    geneID_to_geneSeqID_dict= load_pickle(path+'geneID_to_geneSeqID.cpk')

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
    if os.path.exists(tmp_core_seq_path):
        os.system(''.join(['rm -rf ',tmp_core_seq_path]))
    os.system('mkdir %s'%tmp_core_seq_path)

    ## create dict storing all genes' translation
    if 0:
        gene_aa_dict= defaultdict(dict)
        for accession_id in strain_list:
            gene_aa_dict[accession_id]= read_fasta(''.join([protein_path,accession_id,'.faa']))
        write_pickle(protein_dict_path, gene_aa_dict)

        ## create dict for all gene's nucleotide sequence
        gene_na_dict= defaultdict(dict)
        for accession_id in strain_list:
            gene_na_dict[accession_id]=read_fasta(''.join([nucleotide_path,accession_id,'.fna']))
        write_pickle(nucleotide_dict_path, gene_na_dict)

    gene_aa_dict= load_pickle(protein_dict_path)
    gene_na_dict= load_pickle(nucleotide_dict_path)

    ## write nucleotide and amino-acid sequences for each gene cluster
    export_cluster_seq_tmp(tmp_core_seq_path, core_geneCluster_dt, geneID_to_geneSeqID_dict,
        gene_na_dict, gene_aa_dict)

    tmp_fa_files=glob.glob(tmp_core_seq_path+"*.fna")
    multips(calculate_diversity, parallel, tmp_fa_files, tmp_core_seq_path, species)

    calculated_core_diversity=tmp_average_core_diversity(tmp_core_seq_path)
    refined_core_diversity= round((0.1+factor_core_diversity*calculated_core_diversity)/(1+factor_core_diversity*calculated_core_diversity),4)
    print('factor used: '+str(factor_core_diversity))
    print('average core genome diversity: '+str(calculated_core_diversity))
    print('defined core genome diversity cutoff for splitting long branches: '+str(refined_core_diversity))

    ## move folder tmp_core to the central data folder
    new_clustering_path= '%stmp_core'%path
    if os.path.exists(new_clustering_path):
        os.system(''.join(['rm -r ',new_clustering_path]))
    os.system('mv %s %s'%(tmp_core_seq_path, path))
    return calculated_core_diversity, refined_core_diversity
