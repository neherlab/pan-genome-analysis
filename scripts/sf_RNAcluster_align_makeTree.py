import os, sys, glob, time, shutil
import numpy as np
from SF00_miscellaneous import times, read_fasta, load_pickle, write_pickle, write_in_fa, write_json
from SF06_geneCluster_align_makeTree import mpm_tree, multips, align_and_makeTree, load_sorted_clusters, update_diversity_cpk

#sys.path.append('./scripts/')
sys.setrecursionlimit(2000)

def update_gene_cluster_with_RNA(path, diamond_RNACluster_dt, diamond_geneCluster_dt ):
    ## update gene cluster pickled file
    cluster_path = path+'protein_faa/diamond_matches/'

    diamond_geneCluster_dt.update(diamond_RNACluster_dt)

    write_pickle(cluster_path+'allclusters_final.cpk',diamond_geneCluster_dt)

def create_RNACluster_fa(path):
    """
        input: '.fna', '_RNA_nuc_dict.cpk', 'allclusters.cpk'
        output: '.aln', 'tree.json', etc
    """
    if 0:
        ## make sure the RNACluster folder is empty
        if os.path.isdir(path+'RNACluster/')==True:
            print 'remove previous folder: ',path+'RNACluster/'
            os.system('rm -rf %s'%(path+'RNACluster/'))

    ## dict storing nucleotide Id/Seq from '_RNA_nuc_dict.cpk' files
    istrain_cpk={}; strain_list= load_pickle(path+'strain_list.cpk');
    nucleotide_dict_path= '%s%s'%(path,'nucleotide_fna/')
    for istrain in strain_list:
        istrain_cpk[istrain]=load_pickle(nucleotide_dict_path+istrain+'_RNA_nuc_dict.cpk')

    ## load RNA cluster cpk file
    RNACluster_path=path+'RNA_fna/'
    diamond_RNACluster_dt=load_pickle(RNACluster_path+'allclusters.cpk')

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
            RNA_memb_seq=str(istrain_cpk[strain_name][RNA_memb])
            RNASeqID=RNAID_to_RNASeqID_dict[RNA_memb]
            write_in_fa(RNA_cluster_nu_write, RNASeqID, RNA_memb_seq )
        RNA_cluster_nu_write.close()
    return diamond_RNACluster_dt

def single_RNACluster_align_and_makeTree(thread, alignFile_path, fa_files_list):
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
                myTree.build(raxml=False)
                myTree.ancestral(translate_tree=True)
                myTree.refine()
                myTree.export(path=alignFile_path, RNA_specific=True)
                myTree.diversity_statistics()
                diversity=myTree.diversity
                RNA_diversity_values='{0:.3f}'.format(diversity)
                geneDiversity_file.write('%s\t%s\n'%(clusterID,RNA_diversity_values))
        except:
            print("Aligning and tree building of %s failed"%RNA_cluster_nu_filename)
            print(myTree.tree)

def RNAclusters_align_makeTree( path, parallel ):
    """
    create RNA clusters as nucleotide fasta files
    and build individual RNA trees based on fna files
    """

    diamond_RNACluster_dt=create_RNACluster_fa(path)

    ## align, build_tree, make_RNATree_json
    fasta_path = path+'geneCluster/'
    fa_files=glob.glob(fasta_path+"*RNA*.fna")
    multips(single_RNACluster_align_and_makeTree, fasta_path, parallel, fa_files)
    ## add RNA cluster in diamond_geneCluster_dt
    ### load gene cluster
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    os.system('cp %sallclusters_final.cpk %s/allclusters_final.cpk.bk '%(geneClusterPath,geneClusterPath))
    diamond_geneCluster_dt=load_pickle(geneClusterPath+'allclusters_final.cpk')
    ### update gene cluster with RNA cluster
    update_gene_cluster_with_RNA(path, diamond_RNACluster_dt, diamond_geneCluster_dt)
    ### update diversity file
    update_diversity_cpk(path)
