import os,glob,sys,time,shutil
import numpy as np
from itertools import izip
from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from SF00_miscellaneous import times, read_fasta, load_pickle, write_pickle, write_in_fa, write_json
from SF06_geneCluster_align_makeTree import load_sorted_clusters
sys.path.append('./scripts/')
sys.setrecursionlimit(2000)


def postprocess_unclustered_genes(n_threads, path, nstrains):
    pass
    # 1) detect suspicious peaks in the distribution of average length of genes in gene cluster (aa count)
    #    np.bincount([1,2,3,34,3]) -> how often each entry is found
    #    np.convolve([1,1,1,1,1], gene_length_count)
    #   ->  unclustered genes will contribute many small clusters (size 1)
    #       that result in peaks in the distribution
    # 2) for each peak detected, align the sequences of all genes in clusters in peak
    # 3) to cluster aligned genes, build tree. However, to ensure long branches
    #   ->  between unaligned sub-alignment, fill gaps with random sequence
    #       importantly, this random sequence needs to be the same in different
    #       columns of the alignment.
    #           - random sequence a la rseq = [alpha[ii] for ii in np.random.randint(len(alpha), size=aln.get_aligmentlength())]
    #           - for seq in aln: seq[seq=='-'] = rseq[seq=='-']
    # 4) make and split tree at branches >.5
    # 5) for each subtree (ideally only one big tree), define new gene cluster and run
    #    maketree_align from standard step 6

    output_path='%s%s'%(path,'geneCluster/')
    gene_clusters = load_sorted_clusters(path)
    length_to_cluster = defaultdict(list)
    length_list = []
    for gid, (clusterID, gene) in enumerate(gene_clusters):
        # average length of the cluster in amino acids
        clusterLength= int(np.mean([len(igene) for igene in read_fasta(output_path+'%s%s'%(clusterID,'.fna')).values()])/3.0)
        length_to_cluster[clusterLength].append(clusterID)
        length_list.append(clusterLength)

    cluster_length_distribution = np.bincount(length_list)
    ws=5
    window = np.ones(ws, dtype=float)/ws
    smoothed_length_distribution = np.convolve(cluster_length_distribution, window, mode='same')
    print 3*np.sqrt(smoothed_length_distribution)
    #print max(0.3*nstrains, 3*np.sqrt(smoothed_length_distribution))
    peaks = (cluster_length_distribution - smoothed_length_distribution)> np.maximum(0.3*nstrains, 3*np.sqrt(smoothed_length_distribution))
    print peaks

    def plot_gene_cluster_length_distr():
        '''
        plot branch length against # of gene_cluster_length_distr across trees
        '''
        import matplotlib.pyplot as plt
        plt.ion()
        plt.figure()
        plt.hist( cluster_length_distribution ,normed=True)
        plt.xlabel('gene_cluster_length_distr')
        plt.savefig('./explore_gene_cluster_length_distr.pdf')

    #plot=True
    if True: plot_gene_cluster_length_distr()
