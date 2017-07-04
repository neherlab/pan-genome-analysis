import os,glob,sys,time,shutil
import numpy as np
from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from sf_miscellaneous import times, read_fasta, write_in_fa, load_pickle, write_pickle,\
 write_json, multips
from sf_geneCluster_align_makeTree import align_and_makeTree, update_geneCluster_cpk,\
 update_diversity_cpk, mpm_tree, load_sorted_clusters#, write_final_cluster
from sf_split_long_branch import update_geneCluster_dt, cut_tree_gather_clades, \
 output_cutted_clusters, quick_align_makeTree, cutTree_outputCluster

sys.setrecursionlimit(50000)

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
    merged_cluster_filename='GC_un%03d.fna'%index
    os.system('cat %s > %s%s'%(filenames_str, file_path, merged_cluster_filename))
    ## for amino_acid
    faa_filenames_str=' '.join(['%s%s.faa'%(file_path,i_file) for i_file in clusterID_list])
    faa_clusters_in_peak_filename='GC_un%03d.faa'%index
    os.system('cat %s > %s%s'%(faa_filenames_str, file_path, faa_clusters_in_peak_filename))
    return merged_cluster_filename, cluster_needed_deletion

def find_and_merge_unclustered_genes( path, nstrains, window_size=5, strain_proportion=0.3 , sigma_scale=3):
    """
    detect the unclustered genes and concatenate them
    params:
        nstrains: total number of strains
        window_size
        strain_proportion
        sigma_scale
    return:
        a dict with key of the merged cluster and value of
        a list of related unclustered cluster-name for deletion

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

def delete_old_merged_clusters(file_path, geneCluster_dt, merged_clusters_dict):
    """
    Using given merged_clusters_dict to delete
    mis-clustered records in geneCluster_dt
    Param:
    geneCluster_dt:
        a dict storing cluster statistics
    merged_clusters_dict:
        a dict with key of the merged cluster and value of
        a list of related unclustered cluster-name for deletion
    return:
        updated geneCluster_dt
    """
    cwd = os.getcwd()
    os.chdir(file_path)
    with open('old_clusters_longSplit.txt', 'rb') as delete_cluster_file:
        uncluster_filename_list= [ uncluster_filename.rstrip() for  uncluster_filename in delete_cluster_file]
        try:
            for uncluster_filename in uncluster_filename_list:
                for cluster_needed_deletion in merged_clusters_dict[uncluster_filename]:
                    if cluster_needed_deletion in geneCluster_dt:
                        del geneCluster_dt[cluster_needed_deletion]

                        if os.path.exists(cluster_needed_deletion+'.nwk'):
                            suffix_list=['_aa_aln.fa','_na_aln.fa','.fna','.faa','.nwk','_tree.json']
                        else:
                            suffix_list=['_aa_aln.fa','_na_aln.fa','.fna','.faa']
                        tmp_files=' '.join([ cluster_needed_deletion+suffix for suffix in suffix_list])
                        command_move_deleted_clusters=' '.join(['mv', tmp_files, './deleted_clusters_peaks_splits/'])
                        os.system(command_move_deleted_clusters)
        except:
            print("can't delete"," mis-clusterd genes gathered in ",uncluster_filename)
    os.chdir(cwd)
    return geneCluster_dt

def cut_all_trees_from_merged_clusters(parallel, path, cut_branch_threshold, simple_tree):
    """
    split tree from the unclustered genes and create new cluster files
    params:
        gene_list: lists containing the genes in the new split clusters
    """
    geneCluster_fasta_path='%s%s'%(path,'geneCluster/')
    merged_cluster_filelist=glob.glob(geneCluster_fasta_path+'GC_un*.fna')
    ## parallelization of "post-clustering workflow for merged unclustered records"
    multips(cutTree_outputCluster, parallel, merged_cluster_filelist, geneCluster_fasta_path,
        cut_branch_threshold, treefile_used=False)

    ## gather new clusters from new_clusters_longSplit.txt
    #if os.path.exists(''.join([geneCluster_fasta_path,'new_clusters_longSplit.txt'])):
    with open(''.join([geneCluster_fasta_path,'new_clusters_longSplit.txt']),'rb') as new_clusters_longSplit:
        new_fa_files_list=[ clus.rstrip() for clus in new_clusters_longSplit ]

    ## parallelization of "align and make tree on new cluster"
    multips(align_and_makeTree, parallel, new_fa_files_list, geneCluster_fasta_path, simple_tree)

def postprocess_unclustered_genes(parallel, path, nstrains, simple_tree, split_long_branch_cutoff,
    window_size_smoothed=5, strain_proportion=0.3 , sigma_scale=3):
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
        4) make and split tree at branches > 0.5
        5) for each subtree (ideally only one big tree), define new gene cluster and run
           maketree_align from standard step 6
    """

    geneCluster_fasta_path = ''.join([path,'geneCluster/'])
    new_split_folder= ''.join([geneCluster_fasta_path,'update_long_branch_splits/'])
    if os.path.exists(new_split_folder):
        ## remove the folder from previous run
        os.system(''.join(['rm -r ',new_split_folder]))
    os.system(''.join(['mkdir ',new_split_folder]))

    deleted_clusters_folder=''.join([geneCluster_fasta_path,'deleted_clusters_peaks_splits/'])
    if os.path.exists(deleted_clusters_folder):
        os.system(''.join(['rm -r ',deleted_clusters_folder]))
    os.system(''.join(['mkdir ',deleted_clusters_folder]))

    ## load clusters
    ClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(ClusterPath+'allclusters_postprocessed.cpk')

    ## merge unclustered genes
    merged_clusters_dict=defaultdict(list)
    merged_clusters_dict=find_and_merge_unclustered_genes(path, nstrains, window_size_smoothed, strain_proportion , sigma_scale)

    if len(merged_clusters_dict)!=0:
        ## there are merged clusters corresponding to the cluster peaks

        ## ensure that writing to new_clusters_longSplit starts at the beginning (for re-running)
        if os.path.exists(''.join([geneCluster_fasta_path,'new_clusters_longSplit.txt'])):
            os.system(''.join(['rm ',geneCluster_fasta_path,'new_clusters_longSplit.txt']))
        if os.path.exists(''.join([geneCluster_fasta_path,'old_clusters_longSplit.txt'])):
            os.system(''.join(['rm ',geneCluster_fasta_path,'old_clusters_longSplit.txt']))

        cut_branch_threshold=split_long_branch_cutoff#0.3
        ## cut tree and make new clusters
        cut_all_trees_from_merged_clusters(parallel, path, cut_branch_threshold, simple_tree)

        ## update clusters in allclusters_final.cpk
        #os.system('cp %sallclusters_final.cpk %s/allclusters_final.cpk.bk '%(ClusterPath,ClusterPath))

        ## delete old clusters
        delete_old_merged_clusters(geneCluster_fasta_path, geneCluster_dt, merged_clusters_dict)
        ## add newly split clusters
        update_geneCluster_dt(path,geneCluster_dt)
        ## write updated gene clusters in cpk file
        update_geneCluster_cpk(path,geneCluster_dt)
        ## write gene_diversity_Dt cpk file
        update_diversity_cpk(path)
    #write_final_cluster(path)
