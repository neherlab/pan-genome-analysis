import os,glob,sys,time,shutil
import numpy as np
from itertools import izip
from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from SF00_miscellaneous import times, read_fasta, load_pickle, write_pickle, write_in_fa, write_json
from SF06_geneCluster_align_makeTree import multips, align_and_makeTree, load_sorted_clusters, update_geneCluster_cpk, update_diversity_cpk, mpm_tree

sys.setrecursionlimit(2000)

def update_geneCluster_dt(path,geneCluster_dt):
    """
    add new cluster statistics in folder update_uncluster_splits
    """
    update_uncluster_splits=''.join([path,'geneCluster/update_uncluster_splits/'])
    for ifile in glob.glob(update_uncluster_splits+'*.cpk'):
        for k,v in load_pickle(ifile).iteritems():
            geneCluster_dt[k] = v

def concatenate_cluster_files(clusterID_list, index, file_path):
    """
    concatenate unclustered genes
    params:
        clusterID_list: list of all clusterIDs in detected cluster peak
        index: index of concatenated cluster (based on the order of peaks)
    """
    ## gather the clusterID needed to be deleted
    cluster_needed_deletion=[i_file for i_file in clusterID_list]
    ## for nucleotide
    filenames_str=' '.join(['%s%s.fna'%(file_path,i_file) for i_file in clusterID_list])
    merged_cluster_filename='GC_unclust%03d.fna'%index
    os.system('cat %s > %s%s'%(filenames_str, file_path, merged_cluster_filename))
    ## for amino_acid
    faa_filenames_str=' '.join(['%s%s.faa'%(file_path,i_file) for i_file in clusterID_list])
    faa_clusters_in_peak_filename='GC_unclust%03d.faa'%index
    os.system('cat %s > %s%s'%(faa_filenames_str, file_path, faa_clusters_in_peak_filename))
    return merged_cluster_filename, cluster_needed_deletion

def find_and_merge_unclustered_genes( path, nstrains, window_size=5, strain_proportion=0.3 , sigma_scale=3):
    """
    detect the unclustered genes and concatenate them
    params:
        nstrains: total number of strains
        window_size:
        strain_proportion:
        sigma_scale:
    """
    file_path='%s%s'%(path,'geneCluster/')
    gene_clusters = load_sorted_clusters(path)
    length_to_cluster = defaultdict(list)
    length_list = []
    ## calculate cluster length distribution, link clusterIDs with their clusterLength
    for gid, (clusterID, gene) in enumerate(gene_clusters):
        # average length of the cluster in amino acids
        clusterLength= int(np.mean([len(igene) for igene in read_fasta(file_path+'%s%s'%(clusterID,'.fna')).values()])/3.0)
        length_to_cluster[clusterLength].append(clusterID)
        length_list.append(clusterLength)
    cluster_length_distribution = np.bincount(length_list)

    ## calculate smoothed cluster length distribution
    window_size=5
    window = np.ones(window_size, dtype=float)/window_size
    smoothed_length_distribution = np.convolve(cluster_length_distribution, window, mode='same')

    ## detect peaks
    peaks = (cluster_length_distribution - smoothed_length_distribution)> np.maximum(strain_proportion*nstrains, sigma_scale*np.sqrt(smoothed_length_distribution))
    position_peaks =np.where(peaks)[0]; #cluster_len_peaks= position_peaks*3
    ## concatenate clusters with the same aver. length, return dict of these clusters
    merged_clusters_dict=defaultdict(dict)
    for index, i_peak in enumerate(position_peaks,1):
        merged_cluster_filename, cluster_needed_deletion=concatenate_cluster_files(length_to_cluster[i_peak], index,file_path)
        merged_clusters_dict[merged_cluster_filename]=cluster_needed_deletion
    return merged_clusters_dict

def cut_tree_gather_clades(tree,filepath):
    """
    cut tree via breaking up long branches and gather genes in cutted clades,
    the rest ones are collected in leftover set.
    """
    gene_list=[]
    rest_genes=set()
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
            if node.branch_length > 0.5:
                gene_list.append(set.intersection(node.leafs, leaves))
                leaves=leaves-node.leafs
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
    elif num_rest_genes!=0:
    ## keep the data structure consistent for the loop in output_cutted_clusters()
        rest_genes= [list(rest_genes)]
        
    ## ignore empty sets    
    gene_list = [x for x in gene_list if len(x)]
    if len(gene_list)==1:
        gene_list= [gene_list]
    
    return gene_list, rest_genes

def delete_old_clusters(file_path, geneCluster_dt, merged_clusters_dict):
    """
    delete mis-clustered records in geneCluster_dt
    """
    
    with open(file_path+'delete_misclusters.txt', 'rb') as delete_cluster_file:
        uncluster_filename_list= [ uncluster_filename for  uncluster_filename in delete_cluster_file]
        try:
            for uncluster_filename in uncluster_filename_list:
                uncluster_filename=uncluster_filename.rstrip()
                for cluster_needed_deletion in merged_clusters_dict[uncluster_filename]:
                    if cluster_needed_deletion in geneCluster_dt:
                        del geneCluster_dt[cluster_needed_deletion]
        except:
            print("can't delete"," mis-clusterd genes gathered in ",uncluster_filename)
    return geneCluster_dt

def output_cutted_clusters(file_path, uncluster_filename, gene_list, geneCluster_dt,cut_leftover=False):
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
        ## align the rest genes, build tree, cut long branches till nothing can be cutted.
            
            ## write new divided/split cluster files
            newClusterId="%s_res_%s"%(clusterID,sgs_index)
            gene_cluster_nu_filename="%s%s"%(newClusterId,'.fna')
            gene_cluster_nu_filepath= file_path+gene_cluster_nu_filename
            gene_cluster_nu_write=open(gene_cluster_nu_filepath , 'wb')
            
            gene_cluster_aa_filename="%s%s"%(newClusterId,'.faa')
            gene_cluster_aa_filepath= file_path+gene_cluster_aa_filename
            gene_cluster_aa_write=open( file_path+gene_cluster_aa_filename, 'wb')

            for gene_memb in split_gene_list:
                write_in_fa(gene_cluster_nu_write, gene_memb, origin_nu_fa_dt[gene_memb])
                write_in_fa(gene_cluster_aa_write, gene_memb, origin_aa_fa_dt[gene_memb])

            gene_cluster_nu_write.close(); gene_cluster_aa_write.close();

            ## align rest genes, build tree
            makeTree_cutTree_outputCluster(file_path,[gene_cluster_nu_filepath],geneCluster_dt)

        else:
            ## record the misclusters to be deleted
            with open(file_path+'delete_misclusters.txt', 'a') as delete_cluster_file:
                delete_cluster_file.write('%s\n'%uncluster_filename)
            
            ## define new cluster filenames for nucl/prot sequences
            newClusterId="%s_%s"%(clusterID,sgs_index)
            gene_cluster_nu_filename="%s%s"%(newClusterId,'.fna')
            gene_cluster_aa_filename="%s%s"%(newClusterId,'.faa')
            gene_cluster_nu_write=open( file_path+gene_cluster_nu_filename, 'wb')
            gene_cluster_aa_write=open( file_path+gene_cluster_aa_filename, 'wb')

            ## add record in refined_clusters.txt, which is used for align new clusters
            new_fa_files.add(file_path+gene_cluster_nu_filename)

            ## write new split cluster files
            for gene_memb in split_gene_list:
                write_in_fa(gene_cluster_nu_write, gene_memb, origin_nu_fa_dt[gene_memb])
                write_in_fa(gene_cluster_aa_write, gene_memb, origin_aa_fa_dt[gene_memb])
            gene_cluster_nu_write.close(); gene_cluster_aa_write.close();

            ## write cluster statistics in folder update_uncluster_splits
            addin_geneCluster_dt=defaultdict(list)
            addin_geneCluster_dt[ newClusterId ] = [0,[],0]
            ## num_stains
            addin_geneCluster_dt[ newClusterId ][0]=len(dict(Counter([ ig.split('|')[0] for ig in split_gene_list])).keys())
            ## num_genes
            addin_geneCluster_dt[ newClusterId ][2]=len(dict(Counter([ ig for ig in split_gene_list])).keys())
            ## gene members
            addin_geneCluster_dt[ newClusterId ][1]=[ ig.split('-')[0] for ig in split_gene_list ]
            ## cPickle new cluster statistics
            write_pickle(''.join([file_path,'update_uncluster_splits/', newClusterId,'.cpk']),addin_geneCluster_dt)

    ## write records in gene_diversity file
    with open(file_path+'refined_clusters.txt', 'a') as refined_cluster_file:
        for i in new_fa_files:
            refined_cluster_file.write('%s\n'%i)

def quick_align_makeTree(input_cluster_filepath):
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
        myTree.build(raxml=False)
        # myTree.ancestral(translate_tree=True)
        ###Phylo.write(myTree.tree, input_cluster_filepath+'.s.nwk', 'newick')
        ###print input_cluster_filepath, input_cluster_filepath+'.s.nwk'
        ###print(myTree.tree)
    except:
        print 'tree problem in', input_cluster_filepath
    return myTree.tree

def makeTree_cutTree_outputCluster(file_path, files_list, geneCluster_dt):
    """
    process flow for parallelization to make a tree based on merged cluster,
    cut the tree and output the clades in new clusters
    """
    new_fa_files=set()
    for input_cluster_filepath in files_list:
        input_cluster_filename=input_cluster_filepath.split('/')[-1]
        
        ## make tree
        tree= quick_align_makeTree(input_cluster_filepath)

        ## attempt to cut the tree
        gene_list, rest_genes = cut_tree_gather_clades(tree,input_cluster_filename)
        
        ## output cutted clusters
        if len(gene_list)==0:
            ## nothing can be further cutted, cutting procedure for current tree will stop.
            ## fill gene_list with genes in rest_genes
            gene_list=rest_genes
            ## set the rest_genes to empty set
            rest_genes=set()
        else:
            ## further process left-over genes
            if len(rest_genes)!=0:
                output_cutted_clusters(file_path,
                                    input_cluster_filename,
                                    rest_genes,  geneCluster_dt,
                                    cut_leftover=True)
                                
        ## write clades in gene_list into clusters
        output_cutted_clusters(file_path,
                            input_cluster_filename,
                            gene_list,  geneCluster_dt,
                            cut_leftover=False)

def cut_all_trees_from_merged_clusters(parallel, path, geneCluster_dt):
    """
    split tree from the unclustered genes and create new cluster files
    params:
        gene_list: lists containing the genes in the new split clusters
        geneCluster_dt: cluster dictionary to be updated
    """
    file_path='%s%s'%(path,'geneCluster/')
    merged_cluster_filelist=glob.glob(file_path+'GC_unclust*.fna')
    ## parallelization of "post-clustering workflow for merged unclustered records"
    multips(makeTree_cutTree_outputCluster, parallel, file_path, merged_cluster_filelist,geneCluster_dt)

    ## gather new clusters from refined_clusters.txt
    with open(file_path+'refined_clusters.txt', 'rb') as refined_clusters:
        new_fa_files_list=[ clus.rstrip() for clus in refined_clusters ]

    ## parallelization of "align and make tree on new cluster"
    multips(align_and_makeTree, parallel, file_path, new_fa_files_list)

def postprocess_unclustered_genes(parallel, path, nstrains, window_size_smoothed=5, strain_proportion=0.3 , sigma_scale=3):
    """
        1) detect suspicious peaks in the distribution of average length of genes in gene cluster (aa count)
           np.bincount([1,2,3,34,3]) -> how often each entry is found
           np.convolve([1,1,1,1,1], gene_length_count)
          ->  unclustered genes will contribute many small clusters (size 1)
              that result in peaks in the distribution
        2) for each peak detected, align the sequences of all genes in clusters in peak
        3) to cluster aligned genes, build tree. However, to ensure long branches
          ->  between unaligned sub-alignment, fill gaps with random sequence (skipped, not tested)
              importantly, this random sequence needs to be the same in different columns of the alignment.
                  - random sequence a la rseq = [alpha[ii] for ii in np.random.randint(len(alpha), size=aln.get_aligmentlength())]
                  - for seq in aln: seq[seq=='-'] = rseq[seq=='-']
        4) make and split tree at branches >.5
        5) for each subtree (ideally only one big tree), define new gene cluster and run
           maketree_align from standard step 6
    """
    geneCluster_fasta_path = ''.join([path,'geneCluster/'])
    os.system(''.join(['mkdir ',geneCluster_fasta_path,'update_uncluster_splits']))

    ## load clusters
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(geneClusterPath+'orthamcl-allclusters_final.cpk')

    ## merge unclustered genes
    merged_clusters_dict=defaultdict(list)
    merged_clusters_dict=find_and_merge_unclustered_genes(path, nstrains, window_size_smoothed, strain_proportion , sigma_scale)

    ## ensure that writing to refined_clusters starts at the beginning (for re-running)
    if os.path.exists(''.join([geneCluster_fasta_path,'refined_clusters.txt'])):
        os.system(''.join(['rm ',geneCluster_fasta_path,'refined_clusters.txt']))
    if os.path.exists(''.join([geneCluster_fasta_path,'delete_misclusters.txt'])):
        os.system(''.join(['rm ',geneCluster_fasta_path,'delete_misclusters.txt']))

    ## cut tree and make new clusters
    cut_all_trees_from_merged_clusters(parallel,path,geneCluster_dt)
    ## write new clusters in orthamcl-allclusters_final.cpk
    os.system('cp %sorthamcl-allclusters_final.cpk %s/orthamcl-allclusters_final.cpk.bk '%(geneClusterPath,geneClusterPath))

    delete_old_clusters(geneCluster_fasta_path, geneCluster_dt, merged_clusters_dict)
    update_geneCluster_dt(path,geneCluster_dt)
    update_geneCluster_cpk(path, geneCluster_dt )
    ## write gene_diversity_Dt cpk file
    update_diversity_cpk(path)
