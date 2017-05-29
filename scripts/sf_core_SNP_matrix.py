from sf_geneCluster_align_makeTree import load_sorted_clusters

def create_core_SNP_matrix(path, core_cutoff=1.0, core_gene_strain_fpath=''):#1.0
    """ create SNP matrix using core gene SNPs
        input: strain_list.cpk, core_geneList.cpk
        output: SNP_whole_matrix.aln
        core_cutoff: percentage of strains used to decide whether a gene is core
            default: 1.0 (strictly core gene, which is present in all strains)
            customized: 0.9 ( soft core, considered as core if present in 90% of strains)
    """
    import os,sys,operator
    import numpy as np
    import numpy.ma as ma
    from collections import defaultdict
    from sf_miscellaneous import read_fasta, write_pickle, load_pickle, write_in_fa

    alnFilePath='%s%s'%(path,'geneCluster/')
    output_path= alnFilePath

    ## create core gene list
    corelist=[]
    strain_list=load_pickle(path+'strain_list.cpk')
    totalStrain= len(strain_list)
    sorted_geneList = load_sorted_clusters(path)
    if core_gene_strain_fpath!='':
        with open(core_gene_strain_fpath,'rb') as core_gene_strain_file:
            core_strain_set= set([i.rstrip().replace('-','_') for i in core_gene_strain_file])
    with open(output_path+'core_geneList.txt','wb') as outfile:
        for clusterID, vg in sorted_geneList:
            if core_cutoff==1.0:
                strain_core_cutoff=totalStrain
            else:
                strain_core_cutoff=int(totalStrain*core_cutoff)
            if vg[0]==vg[2] and vg[0]>=strain_core_cutoff:
                coreGeneName='%s%s'%(clusterID,'_na_aln.fa')
                ## sequences might be discarded because of premature stops
                coreGeneName_path= alnFilePath+coreGeneName
                if os.path.exists(coreGeneName_path) and len(read_fasta(coreGeneName_path)) >= strain_core_cutoff:
                    if core_gene_strain_fpath!='' and len(core_strain_set-set([i.split('|')[0] for i in vg[1]]))!=0:
                        continue
                    outfile.write(coreGeneName+'\n')
                    corelist.append(coreGeneName)
                else:
                    #print '%s%s%s'%('warning: ',coreGeneName_path,' is not a core gene')
                    pass

        write_pickle(output_path+'core_geneList.cpk',corelist)

    refSeqList=load_pickle(path+'strain_list.cpk');refSeqList.sort()

    snp_fre_lst=[]; snp_wh_matrix_flag=0
    snp_pos_dt=defaultdict(list); snp_whole_matrix=np.array([])
    snps_by_gene=[]
    for align_file in corelist:## core genes
        nuc_array=np.array([]) # array to store nucleotides for each gene
        gene_seq_dt=read_fasta(alnFilePath+align_file)
        if core_cutoff!=1.0:
            # set sequences for missing gene (space*gene_length)
            missing_gene_seq=' '*len(gene_seq_dt.values()[0])
            totalStrain_sorted_lst=sorted(strain_list)
        # build strain_seq_dt from gene_seq_dt
        strain_seq_dt=defaultdict()
        for gene, seq in gene_seq_dt.iteritems():
            strain_seq_dt[gene.split('-')[0]]=seq # strain-locus_tag-...
        strain_seq_sorted_lst=sorted(strain_seq_dt.items(), key=lambda x: x[0])

        start_flag=0
        if core_cutoff==1.0:
            for ka, va in strain_seq_sorted_lst:
                if start_flag==0:
                    nuc_array=np.array(np.fromstring(va, dtype='S1'))
                    start_flag=1
                else:
                    nuc_array=np.vstack((nuc_array,np.fromstring(va, dtype='S1')))
            ## find SNP positions
            position_polymorphic = np.any(nuc_array != nuc_array[0, :], axis = 0)
            position_has_gap = np.any(nuc_array=='-', axis=0)
            position_SNP = position_polymorphic&(~position_has_gap)
            snp_columns = nuc_array[:,position_SNP]
            snp_pos_dt[align_file]=np.where(position_SNP)[0]
        else:
        ## add '-' for missing genes when dealing with soft core genes
            core_gene_strain=[ gene for gene in strain_seq_dt.keys()]
            for strain in totalStrain_sorted_lst:
                if start_flag==0:
                    if strain in core_gene_strain:
                        nuc_array=np.array(np.fromstring(strain_seq_dt[strain], dtype='S1'))
                    else:
                        print 'Soft core gene: gene absent in strain %s on cluster %s'%(strain,align_file)
                        nuc_array=np.array(np.fromstring(missing_gene_seq, dtype='S1'))
                    start_flag=1
                else:
                    if strain in core_gene_strain:
                        nuc_array=np.vstack((nuc_array,np.fromstring(strain_seq_dt[strain], dtype='S1')))
                    else:
                        print 'Soft core gene: gene absent in strain %s on cluster %s'%(strain,align_file)
                        nuc_array=np.vstack((nuc_array,np.fromstring(missing_gene_seq, dtype='S1')))
            ## find SNP positions
            ## mask missing genes -- determine rows that have ' ' in every column
            is_missing = np.all(nuc_array==' ',axis=1)
            masked_non_missing_array= np.ma.masked_array(nuc_array, nuc_array==' ')
            position_polymorphic = np.any(masked_non_missing_array!= masked_non_missing_array[0, :],axis = 0)
            position_has_gap = np.any(masked_non_missing_array=='-',axis=0)
            position_SNP = position_polymorphic&(~position_has_gap)
            # the below seems duplicated from 5 lines above??
            if is_missing.sum()>0: # with missing genes
                nuc_array[is_missing]='-'
            snp_columns = nuc_array[:,position_SNP]
            snp_pos_dt[align_file]=np.where(position_SNP)[0]
            #print snp_columns

        if snp_wh_matrix_flag==0:
            snp_whole_matrix=snp_columns;
            snp_wh_matrix_flag=1
        else:
            snp_whole_matrix=np.hstack((snp_whole_matrix, snp_columns))
    write_pickle(output_path+'snp_pos.cpk',snp_pos_dt)

    with open(output_path+'SNP_whole_matrix.aln','wb') as outfile:
        for ind, isw in enumerate(snp_whole_matrix):
            write_in_fa( outfile, refSeqList[ind], isw.tostring() )

