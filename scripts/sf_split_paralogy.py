import os, glob, time
import numpy as np
from collections import Counter
from Bio import Phylo
from sf_miscellaneous import read_fasta, write_in_fa, load_pickle, multips
from sf_geneCluster_align_makeTree import align_and_makeTree, find_best_split, update_diversity_cpk, load_sorted_clusters, update_geneCluster_cpk, mem_check

def split_cluster(tree, nstrains, max_branch_length, max_paralogs):
    '''
    determine which clusters to split
    return tree/false depending on whether the cluster should be split or not
    '''
    # determine the optimal split
    best_split = find_best_split(tree)
    # explore linear discriminator
    #return best_split.branch_length/max_branch_length + float(len(best_split.para_nodes))/max_paralogs > 1.0 and len(best_split.para_nodes) > 1
    core_genome_diversity=max_branch_length
    condition1= True if len(best_split.para_nodes)>=nstrains and best_split.split_bl> core_genome_diversity else False
    condition2= True if len(best_split.para_nodes)>max_paralogs and (len(best_split.leaf)==nstrains or len(best_split.not_leafs)==nstrains) and best_split.split_bl>core_genome_diversity else False
    to_split= True if condition1 or condition2 else False
    #return best_split.branch_length/max_branch_length > 1.0 and float(len(best_split.para_nodes))/max_paralogs >= 1.0
    return to_split

def explore_paralogs(path, nstrains, paralog_branch_cutoff, paralog_frac_cutoff=0.3, plot=0):
    '''
    gather paralog statistics for all trees and plot if desired
    parameters:
    paralog_branch_cutoff -- cutoff used to determined whether or not to split cluster.
    paralog_frac_cutoff  -- cutoff for paralog splitting as fraction of total strains.
                            (default 0.3 -- that is 30%)
    '''
    cluster_seqs_path=path+'geneCluster/'
    fname_list =glob.iglob(cluster_seqs_path+'*nwk')
    paralog_stat = []
    for fi,fname in enumerate(fname_list):
        try:
            tree = Phylo.read(fname, 'newick')
        except:
            print 'abc ', fname
        best_split = find_best_split(tree)
        if best_split is not None:
            paralog_stat.append([fname, best_split.branch_length, len(best_split.para_nodes)])
            #paralog_stat.append([fname, best_split.branch_length/mean_branch_length, len(best_split.para_nodes)])
    with open(cluster_seqs_path+'paralogy_statistics.txt','wb') as paralogy_statistics:
        for x,y,z in paralog_stat:
            paralogy_statistics.write('%s\t%s\t%s\n'%(x.split('/')[-1],y,z))

    def plot_paralogs(path):
        '''
        plot branch length against # of paralogs across trees
        '''
        import matplotlib.pyplot as plt
        plt.ion()
        plt.figure()
        plt.scatter([x[1] for x in paralog_stat], [x[2] for x in paralog_stat])
        plt.ylabel('# paralogs')
        plt.xlabel('branch length')
        plt.savefig(path+'explore_paralogs.pdf')

    if plot: plot_paralogs(cluster_seqs_path)
    #return paralog_split_list

def create_split_cluster_files(file_path, fname,
    gene_list1, gene_list2, geneCluster_dt):
    """
    delete the old cluster and create two new clusters
    params:
        new_fa_files: list to which new file names are appeneded
        gene_list1/2: lists containing the genes in the new split clusters
        geneCluster_dt: cluster dictionary to be updated
    """
    orgin_nwk_name = fname.split('/')[-1]
    clusterID = orgin_nwk_name.replace('.nwk','')
    origin_cluster_nu_fa = orgin_nwk_name.replace('nwk','fna')
    origin_cluster_aa_fa = orgin_nwk_name.replace('nwk','faa')

    split_fa_files_set=set()
    try:
        print('deleting:',orgin_nwk_name)
        ##debug:
        ##print('deleting:',orgin_nwk_name,gene_list1,gene_list2, clusterID)
        del geneCluster_dt[clusterID]
    except:
        print("can't delete",orgin_nwk_name)
        ##debug:
        ##print("can't delete",orgin_nwk_name,gene_list1,gene_list2, clusterID)

    ## write new cluster fa files
    origin_nu_fa_dt = read_fasta(file_path+origin_cluster_nu_fa)
    origin_aa_fa_dt = read_fasta(file_path+origin_cluster_aa_fa)
    sgs_index=0

    ## split_gene_list has geneSeqID instead of geneID
    for split_gene_list in (list(gene_list1), list(gene_list2)):
        sgs_index+=1
        newClusterId="%s_%s"%(clusterID,sgs_index)
        gene_cluster_nu_filename="%s%s"%(newClusterId,'.fna')
        gene_cluster_aa_filename="%s%s"%(newClusterId,'.faa')
        gene_cluster_nu_write=open( file_path+gene_cluster_nu_filename, 'wb')
        gene_cluster_aa_write=open( file_path+gene_cluster_aa_filename, 'wb')

        split_fa_files_set |=  set([file_path+gene_cluster_nu_filename])

        ## write new split cluster files
        for gene_memb in split_gene_list:
            if "\\'" in gene_memb:
                gene_memb=gene_memb.replace("\\'","'")
            try:
                write_in_fa(gene_cluster_nu_write, gene_memb, origin_nu_fa_dt[gene_memb])
                write_in_fa(gene_cluster_aa_write, gene_memb, origin_aa_fa_dt[gene_memb])
            except:
                print 'xxxx', fname, gene_memb, gene_list1, gene_list2

        gene_cluster_nu_write.close(); gene_cluster_aa_write.close();

        geneCluster_dt[ newClusterId ] = [0,[],0]
        ## num_stains
        geneCluster_dt[ newClusterId ][0]=len(dict(Counter([ ig.split('|')[0] for ig in split_gene_list])).keys())
        ## num_genes
        geneCluster_dt[ newClusterId ][2]=len(dict(Counter([ ig for ig in split_gene_list])).keys())
        ## gene members
        geneCluster_dt[ newClusterId ][1]=[ ig.split('-')[0] for ig in split_gene_list ]
    return split_fa_files_set

def postprocess_paralogs_iterative(parallel, path, nstrains, simple_tree,
   	paralog_branch_cutoff, paralog_frac_cutoff=0.3, plot=0):

    cluster_path= path+'protein_faa/diamond_matches/'
    geneCluster_dt=load_pickle(cluster_path+'allclusters_postprocessed.cpk')

    split_result= postprocess_paralogs( parallel, path, nstrains, simple_tree,
                                            geneCluster_dt, set(),
                                            paralog_branch_cutoff=paralog_branch_cutoff,
                                            paralog_frac_cutoff=paralog_frac_cutoff, plot=plot)
    n_split_clusters, new_fa_files_set = split_result
    iteration=0
    while(n_split_clusters):
        print('---- split a total of ',n_split_clusters, 'in iteration', iteration)
        split_result= postprocess_paralogs( parallel, path, nstrains, simple_tree,
                                                geneCluster_dt, new_fa_files_set,
                                                paralog_branch_cutoff=paralog_branch_cutoff,
                                                paralog_frac_cutoff=paralog_frac_cutoff, plot=plot)
        n_split_clusters, new_fa_files_set = split_result
        iteration+=1

    ## write gene_diversity_Dt cpk file
    update_diversity_cpk(path)

    ## remove old gene cluster and create new split cluster
    update_geneCluster_cpk(path, geneCluster_dt)

def postprocess_paralogs(parallel, path, nstrains, simple_tree, geneCluster_dt,
    new_fa_files_set,  paralog_branch_cutoff, paralog_frac_cutoff=0.3, plot=0):
    """
    splitting paralogs, discarding old gene clusters and creating new clusters of split paralogs
    params:
        parallel: number of threads to use
        nstrains: total number of strains
        paralog_branch_cutoff: branch length to split (E.g.: core gene diversity as cutoff)
        paralog_frac_cutoff:  fraction of nstrains required for splitting
        plot:      save figure with paralog statistics
    """

    ## exploring paralogs, default: False (not explore and plot), otherwise figure with statistics will saved
    if plot==1:
        explore_paralogs(path, nstrains, paralog_branch_cutoff=paralog_branch_cutoff,
                         paralog_frac_cutoff=paralog_frac_cutoff, plot=plot)

    file_path = path+'geneCluster/'

    if len(new_fa_files_set)==0:
        fname_list =glob.iglob(file_path+'*nwk')
    else:
        fname_list = [ new_fa.replace('.fna','.nwk') for new_fa in new_fa_files_set ]
        print fname_list

    new_fa_files_set= set()
    n_split_clusters = 0

    for fname in fname_list:
        try:
            tree = Phylo.read(fname, 'newick')
        except:
            print 'debug: ',fname, ' ', os.getcwd()

        best_split = find_best_split(tree)

        if best_split is not None:
            do_split = split_cluster(tree, nstrains,
                                     max_branch_length = paralog_branch_cutoff,
                                     #max_branch_length = paralog_branch_cutoff*mean_branch_length,
                                     max_paralogs = paralog_frac_cutoff*nstrains)
            if do_split:
                print('will split:', fname,
                    '#leaves:', tree.count_terminals(),
                    '#best_split.para_nodes:',len(best_split.para_nodes),
                    '#best_split.branch_length:', best_split.branch_length)

                all_genes = set([n.name for n in tree.get_terminals()])
                gene_list1 = set([n.name for n in best_split.get_terminals()])
                gene_list2 = all_genes.difference(gene_list1)
                #print all_genes, gene_list1, gene_list2

                new_fa_files = create_split_cluster_files(file_path, fname, gene_list1, gene_list2, geneCluster_dt)
                new_fa_files_set |= new_fa_files
                n_split_clusters+=1

    print 'new_split_fasta_files', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), new_fa_files_set
    ## make new aln and tree
    #mem_check('multips(align_and_')
    multips(align_and_makeTree, parallel, list(new_fa_files_set), file_path, simple_tree)
    return n_split_clusters, new_fa_files_set
