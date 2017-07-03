import os, sys, glob, time, shutil
import numpy as np
from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from sf_miscellaneous import times, read_fasta, write_in_fa, load_pickle, write_pickle, check_dependency, multips
from sf_geneCluster_align_makeTree import mpm_tree, \
    align_and_makeTree, write_final_cluster, update_geneCluster_cpk, update_diversity_cpk, mem_check

def update_geneCluster_dt(path,geneCluster_dt):
    """
    add new cluster statistics in folder update_long_branch_splits
    geneCluster_dt: geneCluster dict to be updated
    """
    update_long_branch_splits=''.join([path,'geneCluster/update_long_branch_splits/'])
    for ifile in glob.iglob(update_long_branch_splits+'*.cpk'):
        for k,v in load_pickle(ifile).iteritems():
            #print('adding newly split clusters %s'%k)
            geneCluster_dt[k] = v

def cut_tree_gather_clades(tree, cut_branch_threshold):
    """
    cut tree via breaking up long branches and gather genes in cutted clades,
    the rest ones are collected in leftover set.
    """
    # gene list containing the genes in the each new split cluster
    gene_list=[]
    # gene list containing the left-over not included in split gene list (gene_list)
    rest_genes=[]
    ## save children(leaves) for each node
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            node.leafs = set([node.name])
        else:
            node.leafs = set.union(*[c.leafs for c in node])

    leaves = set([leaf.name for leaf in tree.get_terminals()])

    ## find long branches and their children
    # loop through tree in post-order, cut at long branch,
    # store leaves of node intersected with leaves not yet dealt with
    try:
        for node in tree.find_clades(order='postorder'):
            if node.branch_length > cut_branch_threshold:
                gene_list.append(set.intersection(node.leafs, leaves))
                leaves=leaves-node.leafs
            elif False and node==tree.root: # check the sum of children.branch_length
            #Notes: not encourage to use based on test using simulation data
                root_children_bl=[child.branch_length for child in node]
                if max(root_children_bl) < cut_branch_threshold and sum(root_children_bl) > cut_branch_threshold:
                    gene_list.append(set.intersection(node[0].leafs, leaves))
                    leaves=leaves-node[0].leafs

        ## gather the rest unsplit genes in one cluster, filter for empty clusters
    except:
        #import ipdb; ipdb.set_trace();
        print("find_clades problem")

    ## gather leftover genes in
    rest_genes= leaves

    num_rest_genes= len(rest_genes)
    if num_rest_genes==1:
        ## singleton appended to the gene_list
        gene_list.append(leaves)
        rest_genes=[]
    elif num_rest_genes>1:
        ## keep the data structure consistent for the loop in output_cutted_clusters()
        rest_genes= [rest_genes]

    ## ignore empty sets
    gene_list = [x for x in gene_list if len(x)]

    return gene_list, rest_genes

def delete_original_clusters(file_path, geneCluster_dt):
    """
    Delete records in old_clusters_longSplit.txt to update geneCluster_dt
    Param:
    geneCluster_dt:
        a dict storing cluster statistics
    return:
        updated geneCluster_dt
    """
    cwd = os.getcwd()
    os.chdir(file_path)
    with open('old_clusters_longSplit.txt', 'rb') as delete_cluster_file:
        uncluster_filename_list= [ uncluster_filename.split('.fna')[0] for  uncluster_filename in delete_cluster_file]
        #geneCluster_dt.keys()[0]
        suffix_list=['_aa_aln.fa','_na_aln.fa','.fna','.faa','.nwk','_tree.json']
        for uncluster_filename in uncluster_filename_list:
            if uncluster_filename in geneCluster_dt:
                del geneCluster_dt[uncluster_filename]
                tmp_files=' '.join([ uncluster_filename+suffix for suffix in suffix_list])
                command_move_deleted_clusters=' '.join(['mv', tmp_files, './deleted_clusters_longSplit/'])
                os.system(command_move_deleted_clusters)
    os.chdir(cwd)
    return geneCluster_dt

def output_cutted_clusters(file_path, uncluster_filename, gene_list, cut_branch_threshold, treefile_used=None, cut_leftover=None):
    """
    delete the unclustered file and create new clusters
    params:
        gene_list: lists containing the genes in the new split clusters
        geneCluster_dt: cluster dictionary to be updated
        cut_leftover: flag to indicate whether there are the leftover nodes
            after cutting long branches. Default: empty.
    """
    clusterID = uncluster_filename.replace('.fna','')
    origin_uncluster_nu_fa = uncluster_filename
    origin_uncluster_aa_fa = uncluster_filename.replace('fna','faa')

    new_fa_files=set()

    ## load origin cluster fa files
    origin_nu_fa_dt = read_fasta(file_path+origin_uncluster_nu_fa)
    origin_aa_fa_dt = read_fasta(file_path+origin_uncluster_aa_fa)

    ## split_gene_list has geneSeqID instead of geneID
    for sgs_index,split_gene_list in enumerate(gene_list,1):
        if cut_leftover==True:
            ## newClusterId for the rest genes (_r as identifier)
            newClusterId="%s_r%s"%(clusterID,sgs_index)
        else:
            newClusterId="%s_%s"%(clusterID,sgs_index)

        #=============================================
        ## write new divided/split cluster files
        gene_cluster_nu_filename="%s%s"%(newClusterId,'.fna')
        gene_cluster_nu_filepath= file_path+gene_cluster_nu_filename
        gene_cluster_nu_write=open(gene_cluster_nu_filepath , 'wb')

        gene_cluster_aa_filename="%s%s"%(newClusterId,'.faa')
        gene_cluster_aa_filepath= file_path+gene_cluster_aa_filename
        gene_cluster_aa_write=open( file_path+gene_cluster_aa_filename, 'wb')

        for gene_memb in split_gene_list:
            if "\\'" in gene_memb: # Replace '\' in node name:
                ## NC_018495|CM9_RS01675-1-guanosine-3',5'-... in fasta ID
                ## 'NC_018495|CM9_RS01675-1-guanosine-3\',5\'-...' in nwk node name
                ## Use origin_nu_fa_dt[gene_memb] will throw the KeyError:
                ## "NC_018495|CM9_RS01675-1-guanosine-3\\',5\\'"
                gene_memb=gene_memb.replace("\\'","'")

            write_in_fa(gene_cluster_nu_write, gene_memb, origin_nu_fa_dt[gene_memb])
            write_in_fa(gene_cluster_aa_write, gene_memb, origin_aa_fa_dt[gene_memb])
        gene_cluster_nu_write.close(); gene_cluster_aa_write.close();
        #=============================================

        if cut_leftover==True:
            ## align the rest genes, build tree, cut long branches till nothing can be cutted.
            cutTree_outputCluster([gene_cluster_nu_filepath],file_path, cut_branch_threshold, treefile_used)
        else:
            ## record the misclusters to be deleted (already addressed in cutTree_outputCluster )
            ## it will output the same cluster several times
            #with open(file_path+'old_clusters_longSplit.txt', 'a') as delete_cluster_file:
            #    delete_cluster_file.write('%s\n'%uncluster_filename)

            ## add record in new_clusters_longSplit.txt, which is used for align new clusters
            new_fa_files.add(gene_cluster_nu_filepath)

            ## write cluster statistics in folder update_long_branch_splits
            addin_geneCluster_dt=defaultdict(list)
            addin_geneCluster_dt[ newClusterId ] = [0,[],0]
            ## num_stains
            addin_geneCluster_dt[ newClusterId ][0]=len(dict(Counter([ ig.split('|')[0] for ig in split_gene_list])).keys())
            ## num_genes
            addin_geneCluster_dt[ newClusterId ][2]=len(dict(Counter([ ig for ig in split_gene_list])).keys())
            ## gene members
            addin_geneCluster_dt[ newClusterId ][1]=[ ig.split('-')[0] for ig in split_gene_list ]
            ## cPickle new cluster statistics
            write_pickle(''.join([file_path,'update_long_branch_splits/', newClusterId,'.cpk']),addin_geneCluster_dt)

    ## write records in gene_diversity file
    with open(file_path+'new_clusters_longSplit.txt', 'a') as refined_cluster_file:
        for i in new_fa_files:
            refined_cluster_file.write('%s\n'%i)

def quick_align_makeTree(input_cluster_filepath,fasttree_name):
    """
    make a quick alignment and tree-building on post-processed clusters,
    cut the tree and create new clusters for the cutted clades
    """
    ## make a temporary tree for the purpose of checking long branches
    ## (different from creating a concrete tree via function align_and_makeTree)

    try:
        myTree = mpm_tree(input_cluster_filepath)
        myTree.codon_align()
        myTree.translate()
        myTree.build(raxml=False,fasttree_program=fasttree_name)
    except:
        print('tree problem in', input_cluster_filepath)
    return myTree.tree

def cutTree_outputCluster( file_list, file_path, cut_branch_threshold, treefile_used):
    """
    process flow for parallelization to cut the tree and output the clades in new clusters
    """
    new_fa_files=set()
    fasttree_name= 'fasttree' if check_dependency('fasttree') else 'FastTree'
    for input_filepath in file_list:
        if treefile_used==True:
            ## read tree
            input_cluster_filename=input_filepath.split('/')[-1].replace('.nwk','.fna')
            try:
                tree= Phylo.read(input_filepath, 'newick')
            except:
                print 'reading tree failed: ',input_filepath
        else:
            ## make tree
            input_cluster_filename=input_filepath.split('/')[-1]
            tree= quick_align_makeTree(input_filepath,fasttree_name)

        ## attempt to cut the tree
        gene_list, rest_genes = cut_tree_gather_clades(tree,cut_branch_threshold)

        ## add to-be-deleted cluster records
        if len(gene_list)!=0 and '_r' not in input_cluster_filename:
            ## 1st check: original cluster has been split
            ## 2nd check: it's not a "further-split" cluster
            ##            from an already split cluster
            with open(file_path+'old_clusters_longSplit.txt', 'a') as delete_cluster_file:
                #print 'delete clusters that have been split: ',input_cluster_filename
                delete_cluster_file.write('%s\n'%input_cluster_filename)

        ## output cutted clusters
        if len(gene_list)==0:
            ## nothing can be further cutted,
            ## cutting process for current tree will stop.
            if '_r' not in input_cluster_filename:
                ## a tree does not need to be split, skip the following
                pass#continue
            else:
                ## this's a list of rest genes which cannot be split.
                ## fill gene_list with genes in rest_genes
                gene_list=rest_genes
                ## set the rest_genes to empty list
                rest_genes=[]
        else:
            ## further process on left-over genes
            if len(rest_genes)!=0:
                output_cutted_clusters(file_path, input_cluster_filename,
                                    rest_genes, cut_branch_threshold,
                                    treefile_used=False, cut_leftover=True)

        ## write clades in gene_list into clusters
        output_cutted_clusters(file_path, input_cluster_filename,
                            gene_list, cut_branch_threshold,
                            treefile_used=False, cut_leftover=False)

def postprocess_split_long_branch(parallel, path, simple_tree, cut_branch_threshold=0.3):
    """
    Split tree via breaking up long branches.
    Remote homology leads to over-clustering. This yields tree with long branches.
    """

    file_path = ''.join([path,'geneCluster/'])
    new_split_folder= ''.join([file_path,'update_long_branch_splits/'])
    if os.path.exists(new_split_folder):
        ## remove the folder from previous run
        os.system(''.join(['rm -r ',new_split_folder]))
    os.system(''.join(['mkdir ',new_split_folder]))
    deleted_clusters_folder=''.join([file_path,'deleted_clusters_longSplit/'])
    if os.path.exists(deleted_clusters_folder):
        os.system(''.join(['rm -r ',deleted_clusters_folder]))
    os.system(''.join(['mkdir ',deleted_clusters_folder]))

    ## load clusters
    cluster_path='%s%s'%(path,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(cluster_path+'allclusters.cpk')

    ## gather all trees generated before postprocessing
    tree_path = file_path
    tree_fname_list =glob.glob(tree_path+'*nwk')

    ## ensure that writing to new_clusters_longSplit starts at the beginning (for re-running)
    if os.path.exists(''.join([file_path,'new_clusters_longSplit.txt'])):
        os.system(''.join(['rm ',file_path,'new_clusters_longSplit.txt']))
    if os.path.exists(''.join([file_path,'old_clusters_longSplit.txt'])):
        os.system(''.join(['rm ',file_path,'old_clusters_longSplit.txt']))

    # =============================================
    # parallelization:
    # "post-clustering workflow for splitting trees on over-clustered records"
    treefile_used=True
    multips(cutTree_outputCluster, parallel, tree_fname_list, file_path, cut_branch_threshold, treefile_used)

    ## If new_clusters_longSplit.txt (over_split records) exists,
    ## then gather new clusters from new_clusters_longSplit.txt
    if os.path.exists(''.join([file_path,'new_clusters_longSplit.txt'])):
        with open(file_path+'new_clusters_longSplit.txt', 'rb') as new_clusters_longSplit:
            new_fa_files_list=[ clus.rstrip() for clus in new_clusters_longSplit ]
            print '#times of splitting long branches:',len(new_fa_files_list)-1
        with open(file_path+'old_clusters_longSplit.txt', 'rb') as delete_cluster_file:
            deleted_file_count=len([ clus for clus in delete_cluster_file ])
            print '#clusters split during the checking of long branches:',deleted_file_count

        ## parallelization of "align and make tree on new cluster"
        multips(align_and_makeTree, parallel, new_fa_files_list, file_path, simple_tree)
        # =============================================

        ## delete original clusters which are split
        delete_original_clusters(file_path, geneCluster_dt)
        ## add newly split clusters
        update_geneCluster_dt(path,geneCluster_dt)
        ## write updated gene clusters in cpk file
        update_geneCluster_cpk(path, geneCluster_dt)
        ## write gene_diversity_Dt cpk file
        update_diversity_cpk(path)

        os.system(' '.join(['mv ',file_path+'new_clusters_longSplit.txt' ,file_path+'added_clusters_split_long.txt' ]))
        os.system(' '.join(['mv ',file_path+'old_clusters_longSplit.txt', file_path+'deleted_clusters_split_long.txt']))
    else: # no clusters postprocessed
        os.system(' '.join(['cp',cluster_path+'allclusters.cpk',cluster_path+'allclusters_postprocessed.cpk']))



