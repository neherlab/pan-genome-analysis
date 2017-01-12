from collections import defaultdict
from sf_geneCluster_align_makeTree import load_sorted_clusters
from sf_miscellaneous import load_pickle, write_pickle, write_in_fa, write_json

def create_genePresence(dt_strainGene, totalStrain, set_totalStrain, all_gene_names):
    """
    append to dictionary for gene presence/absence string {strain: '010100010010...',}
    function is called for every gene cluster
    params:
        - totalStrain:      total number of strains
        - set_totalStrain:  set of all strain names
        - all_gene_names:    list of gene names (incl strain name) of all genes in cluster
    """
    # determine set of strains that are represented in gene cluster (first element of gene name)
    set_sharedStrain=set([ igl.split('|')[0] for igl in all_gene_names])
    for ist in set_sharedStrain:
        # append '1' to strains that have the gene
        dt_strainGene[ist]='%s1'%dt_strainGene[ist]
    if totalStrain!=len(set_sharedStrain):
        # append '0' to the remaining strains
        for ist0 in set_totalStrain-set_sharedStrain:
            dt_strainGene[ist0]='%s0'%dt_strainGene[ist0]

def make_genepresence_alignment(path, enable_gain_loss):
    '''
    loop over all gene clusters and append 0/1 to strain specific
    string used as pseudo alignment of gene presence absence
    '''
    geneClusterPath='%s%s'%(path,'protein_fna/diamond_matches/')
    output_path='%s%s'%(path,'geneCluster/');

    ## load strain list and prepare for gene presence/absence
    strain_list= load_pickle('%s%s'%(path,'strain_list.cpk'))
    set_totalStrain=set([ istrain for istrain in strain_list ])
    totalStrain=len(set_totalStrain)
    dt_strainGene= defaultdict(str)

    sorted_genelist = load_sorted_clusters(path)
    ## sorted_genelist: [(clusterID, [ count_strains,[memb1,...],count_genes]),...]
    for clusterID, gene in sorted_genelist:
        ## append 0/1 to each strain
        create_genePresence(dt_strainGene, totalStrain, set_totalStrain, gene[1])

    with open('%s%s'%(output_path,'genePresence.aln'),'wb') as presence_outfile:
        for istkey in dt_strainGene:
            dt_strainGene[istkey]=''.join(dt_strainGene[istkey] )
            write_in_fa( presence_outfile, istkey, dt_strainGene[istkey] )
    write_pickle('%s%s'%(output_path,'dt_genePresence.cpk'), dt_strainGene)

    if enable_gain_loss==0:
        gene_loss_fname='%s%s'%(output_path,'geneGainLossEvent.json')
        write_json(dt_strainGene, gene_loss_fname, indent=1)
        geneEvents_dt={ i:0 for i in range(len(sorted_genelist)) }
        write_pickle('%s%s'%(output_path,'dt_geneEvents.cpk'), geneEvents_dt)
