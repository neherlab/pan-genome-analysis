import os,glob,sys,time,shutil
import numpy as np
from itertools import izip
from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from SF00_miscellaneous import times, read_fasta, load_pickle, write_pickle, write_in_fa, write_json
from SF06_geneCluster_align_makeTree import multips, align_and_makeTree, load_sorted_clusters, update_gene_cluster, update_diversity_cpk_file
sys.path.append('./scripts/')
sys.setrecursionlimit(2000)

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

def find_and_merge_unclustered_genes(n_threads, path, nstrains, window_size=5, strain_proportion=0.3 , sigma_scale=3):
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

def create_split_unclustered_files(file_path, fname, 
                               gene_list,  diamond_geneCluster_dt, merged_clusters_dict):
    """
    delete the unclustered file and create new clusters
    params:
        gene_list: lists containing the genes in the new split clusters
        diamond_geneCluster_dt: cluster dictionary to be updated
        merged_clusters_dict: merged_clusters_dict: dictionary of merged clusters (key) with original clusterIDs (value),
            which is used to delete the old unclustered items in diamond_geneCluster_dt
    """
    origin_uncluster_nwk_name = fname.split('/')[-1]
    clusterID = origin_uncluster_nwk_name.replace('.fna','')
    origin_uncluster_nu_fa = origin_uncluster_nwk_name
    origin_uncluster_aa_fa = origin_uncluster_nwk_name.replace('fna','faa')

    split_fa_files_set=set()
    try:
        ## delete under-clustered clusters
        for cluster_needed_deletion in merged_clusters_dict[origin_uncluster_nwk_name]:
            if cluster_needed_deletion in diamond_geneCluster_dt: 
                del diamond_geneCluster_dt[cluster_needed_deletion]
                print('deleting:', cluster_needed_deletion, ' gathered in ', origin_uncluster_nwk_name)
    except:
        #print("can't delete",origin_uncluster_nwk_name,gene_list, clusterID)
        print("can't delete"," under_clusterd genes gathered in ",origin_uncluster_nwk_name)

    ## write new cluster fa files
    origin_nu_fa_dt = read_fasta(file_path+origin_uncluster_nu_fa)
    origin_aa_fa_dt = read_fasta(file_path+origin_uncluster_aa_fa)
    #print gene_list
    ## split_gene_list has geneSeqID instead of geneID
    for sgs_index,split_gene_list in enumerate(gene_list,1):
        newClusterId="%s_%s"%(clusterID,sgs_index)
        gene_cluster_nu_filename="%s%s"%(newClusterId,'.fna')
        gene_cluster_aa_filename="%s%s"%(newClusterId,'.faa')
        gene_cluster_nu_write=open( file_path+gene_cluster_nu_filename, 'wb')
        gene_cluster_aa_write=open( file_path+gene_cluster_aa_filename, 'wb')

        split_fa_files_set |=  set([file_path+gene_cluster_nu_filename])

        ## write new split cluster files
        for gene_memb in split_gene_list:
            write_in_fa(gene_cluster_nu_write, gene_memb, origin_nu_fa_dt[gene_memb])
            write_in_fa(gene_cluster_aa_write, gene_memb, origin_aa_fa_dt[gene_memb])
        gene_cluster_nu_write.close(); gene_cluster_aa_write.close();

        diamond_geneCluster_dt[ newClusterId ] = [0,[],0]
        ## num_stains
        diamond_geneCluster_dt[ newClusterId ][0]=len(dict(Counter([ ig.split('|')[0] for ig in split_gene_list])).keys())
        ## num_genes
        diamond_geneCluster_dt[ newClusterId ][2]=len(dict(Counter([ ig for ig in split_gene_list])).keys())
        ## gene members
        diamond_geneCluster_dt[ newClusterId ][1]=[ ig.split('-')[0] for ig in split_gene_list ]
    return split_fa_files_set

def cut_tree_from_merged_clusters(parallel, path, diamond_geneCluster_dt, merged_clusters_dict):
    """
    split tree from the unclustered genes into new clusters
    params:
        gene_list: lists containing the genes in the new split clusters
        diamond_geneCluster_dt: cluster dictionary to be updated
        merged_clusters_dict: dictionary of merged clusters (key) with original clusterIDs (value),
            which is used to delete the old unclustered items in diamond_geneCluster_dt
    """
    file_path='%s%s'%(path,'geneCluster/')
    from SF06_geneCluster_align_makeTree import mpm_tree
    new_fa_files_set= set()
    for merged_cluster_file in glob.glob(file_path+'GC_unclust*.fna'):
        myTree = mpm_tree(merged_cluster_file)
        myTree.align()
        myTree.build(raxml=False,treetime_used=True)
        #Phylo.write(myTree.tree, file_path+myTree.fname_prefix+'.nwk', 'newick')

        gene_list=[]
        #uncluster_filename= merged_cluster_file.split('/')[-1]
        ## save children(leaves) for each node
        for node in myTree.tree.find_clades(order='postorder'):
            if node.is_terminal():
                node.leafs = set([node.name])
            else:
                node.leafs = set.union(*[c.leafs for c in node])

        leaves = set([leaf.name for leaf in myTree.tree.get_terminals()])

        ## find long branches and their children
        for node in myTree.tree.get_nonterminals('preorder'):
            if node.branch_length> 0.5:
                gene_list.append(node.leafs)
                leaves=leaves-node.leafs
        ## gather the rest unsplit genes in one cluster
        gene_list.append(leaves)
        ## split clades into clusters
        new_fa_files = create_split_unclustered_files(file_path, merged_cluster_file,
                               gene_list,  diamond_geneCluster_dt, merged_clusters_dict)
        new_fa_files_set |= new_fa_files
    ## align and make tree on new cluster
    multips(align_and_makeTree, file_path, parallel, list(new_fa_files_set))

def postprocess_unclustered_genes(parallel, n_threads, path, nstrains, window_size=5, strain_proportion=0.3 , sigma_scale=3):

    # 1) detect suspicious peaks in the distribution of average length of genes in gene cluster (aa count)
    #    np.bincount([1,2,3,34,3]) -> how often each entry is found
    #    np.convolve([1,1,1,1,1], gene_length_count)
    #   ->  unclustered genes will contribute many small clusters (size 1)
    #       that result in peaks in the distribution
    # 2) for each peak detected, align the sequences of all genes in clusters in peak
    # 3) to cluster aligned genes, build tree. However, to ensure long branches
    #   ->  between unaligned sub-alignment, fill gaps with random sequence (skipped, not tested)
    #       importantly, this random sequence needs to be the same in different columns of the alignment.
    #           - random sequence a la rseq = [alpha[ii] for ii in np.random.randint(len(alpha), size=aln.get_aligmentlength())]
    #           - for seq in aln: seq[seq=='-'] = rseq[seq=='-']
    # 4) make and split tree at branches >.5
    # 5) for each subtree (ideally only one big tree), define new gene cluster and run
    #    maketree_align from standard step 6

    ## load clusters
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    diamond_geneCluster_dt=load_pickle(geneClusterPath+'orthamcl-allclusters_final.cpk')
    ## merge unclustered genes
    merged_clusters_dict=defaultdict(list)
    merged_clusters_dict=find_and_merge_unclustered_genes(n_threads, path, nstrains, window_size, strain_proportion , sigma_scale)
    ## cut tree and make new clusters
    cut_tree_from_merged_clusters(parallel,path,diamond_geneCluster_dt,merged_clusters_dict)
    ## write new clusters in orthamcl-allclusters_final.cpk
    os.system('cp %sorthamcl-allclusters_final.cpk %s/orthamcl-allclusters_final.cpk.bk '%(geneClusterPath,geneClusterPath))
    update_gene_cluster(path, diamond_geneCluster_dt )
    ## write gene_diversity_Dt cpk file
    update_diversity_cpk_file(path)