import os, sys, glob, time, shutil
import numpy as np
from sf_miscellaneous import times, read_fasta, load_pickle, write_pickle, write_in_fa, write_json, check_dependency, multips
from sf_geneCluster_align_makeTree import mpm_tree, align_and_makeTree, load_sorted_clusters, update_diversity_cpk

sys.setrecursionlimit(50000)

def update_gene_cluster_with_RNA(path, diamond_RNACluster_dt, diamond_geneCluster_dt ):
    ## update gene cluster pickled file
    cluster_path = path+'protein_faa/diamond_matches/'

    diamond_geneCluster_dt.update(diamond_RNACluster_dt)

    write_pickle(cluster_path+'allclusters_postprocessed.cpk',diamond_geneCluster_dt)

def create_RNACluster_fa(path,folders_dict):
    """
        input: '.fna', '_RNA_nuc_dict.cpk', 'allclusters.cpk'
        output: '.aln', 'tree.json', etc
    """
    RNA_path= folders_dict['RNA_path']
    RNA_dict= load_pickle('%s%s'%(RNA_path,'all_RNA_seq.cpk'))

    ## load RNA cluster cpk file
    diamond_RNACluster_dt=load_pickle(RNA_path+'allclusters.cpk')

    ## load RNAID_to_RNASeqID RNASeqID cpk file
    RNAID_to_RNASeqID_dict=load_pickle(path+'RNAID_to_SeqID.cpk')

    ## create cluster-RNAs fasta files (By default: put RNAs in geneCluster folder)
    fasta_path=path+'geneCluster/';
    ## diamond_RNACluster_dt: {clusterID:[ count_strains,[memb1,...],count_RNAs }
    for clusterID, RNA in diamond_RNACluster_dt.iteritems():
        ## RNACluster file name
        RNA_cluster_nu_filename="%s%s"%(clusterID,'.fna')
        RNA_cluster_nu_write=open( fasta_path+RNA_cluster_nu_filename, 'wb')
        ## write nucleotide/amino_acid sequences into RNACluster files
        for RNA_memb in RNA[1]:
            ## RNA_name format: strain_1|locusTag
            strain_name= RNA_memb.split('|')[0]
            RNA_memb_seq=str(RNA_dict[strain_name][RNA_memb])
            RNASeqID=RNAID_to_RNASeqID_dict[RNA_memb]
            write_in_fa(RNA_cluster_nu_write, RNASeqID, RNA_memb_seq )
        RNA_cluster_nu_write.close()
    return diamond_RNACluster_dt

def single_RNACluster_align_and_makeTree(fa_files_list, alignFile_path, simple_tree):
    fasttree_name= 'fasttree' if check_dependency('fasttree') else 'FastTree'
    for RNA_cluster_nu_filename in fa_files_list:
        try:
            # extract GC_RNA002 from path/GC_RNA002.aln
            clusterID = RNA_cluster_nu_filename.split('/')[-1].split('.')[0]
            geneDiversity_file = open(alignFile_path+'gene_diversity.txt', 'a')
            if len( read_fasta(RNA_cluster_nu_filename) )==1: # nothing to do for singletons
                ## na.aln
                RNA_cluster_nu_aln_filename= RNA_cluster_nu_filename.replace('.fna','_na.aln')
                ## RNA SeqID separator '|' is replaced by '-' for msa viewer compatibility
                with open(RNA_cluster_nu_aln_filename,'wb') as write_file:
                    for SeqID, Sequence in read_fasta(RNA_cluster_nu_filename).iteritems():
                        write_in_fa(write_file, SeqID.replace('|','-'), Sequence)
                geneDiversity_file.write('%s\t%s\n'%(clusterID,'0.0'))
            else: # align and build tree
                print RNA_cluster_nu_filename
                myTree = mpm_tree(RNA_cluster_nu_filename)
                myTree.align()

                if simple_tree==False:
                    myTree.build(raxml=False,fasttree_program=fasttree_name,treetime_used=True)
                    myTree.ancestral(translate_tree=True)
                    myTree.refine(CDS=False)
                else:
                    myTree.build(raxml=False,fasttree_program=fasttree_name,treetime_used=False)
                myTree.diversity_statistics_nuc()
                myTree.export(path=alignFile_path, RNA_specific=True)

                RNA_diversity_values='{0:.3f}'.format(myTree.diversity_nuc)
                geneDiversity_file.write('%s\t%s\n'%(clusterID,RNA_diversity_values))
                print clusterID,RNA_diversity_values
        except:
            print("Aligning and tree building of RNA %s failed"%RNA_cluster_nu_filename)

def RNAclusters_align_makeTree( path, folders_dict, parallel, simple_tree ):
    """
    create RNA clusters as nucleotide fasta files
    and build individual RNA trees based on fna files
    """

    diamond_RNACluster_dt=create_RNACluster_fa(path,folders_dict)

    ## align, build_tree, make_RNATree_json
    fasta_path = path+'geneCluster/'
    fa_files=glob.glob(fasta_path+"*RC*.fna")
    multips(single_RNACluster_align_and_makeTree, parallel, fa_files, fasta_path, simple_tree)
    ## add RNA cluster in diamond_geneCluster_dt
    ### load gene cluster
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    os.system('cp %sallclusters_postprocessed.cpk %s/allclusters_postprocessed.cpk.bk '%(geneClusterPath,geneClusterPath))
    diamond_geneCluster_dt=load_pickle(geneClusterPath+'allclusters_postprocessed.cpk')
    ### update gene cluster with RNA cluster
    update_gene_cluster_with_RNA(path, diamond_RNACluster_dt, diamond_geneCluster_dt)
    ### update diversity file
    update_diversity_cpk(path)
