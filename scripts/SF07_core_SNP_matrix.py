from SF06_geneCluster_align_makeTree import load_sorted_clusters

def create_core_SNP_matrix(path):
    """ create SNP matrix using core gene SNPs
        input: strain_list.cpk, core_geneList.cpk
        output: SNP_whole_matrix.aln
    """
    import sys,operator
    import numpy as np
    from collections import defaultdict
    from SF00_miscellaneous import read_fasta, write_pickle, load_pickle, write_in_fa

    alnFilePath='%s%s'%(path,'geneCluster/')
    output_path= alnFilePath

    ## create core gene list
    corelist=[];
    totalStrain= len(load_pickle(path+'strain_list.cpk'))
    sorted_geneList = load_sorted_clusters(path)
    with open(output_path+'core_geneList.txt','wb') as outfile:
        for clusterID, vg in sorted_geneList:
            if vg[0]==totalStrain and vg[2]==totalStrain:
                coreGeneName='%s%s'%(clusterID,'_na.aln')
                outfile.write(coreGeneName+'\n')
                corelist.append(coreGeneName)
        write_pickle(output_path+'core_geneList.cpk',corelist)

    refSeqList=load_pickle(path+'strain_list.cpk');refSeqList.sort()

    snp_fre_lst=[]; snp_wh_matrix_flag=0
    snp_pos_dt=defaultdict(list); snp_whole_matrix=np.array([])

    snps_by_gene=[]
    for align_file in corelist:## all core genes
        fa_dt=read_fasta(alnFilePath+align_file)
        if len(fa_dt) != totalStrain:
            print '%s%s%s'%('warning: problem in ',align_file,'; not a core gene')
            continue
        else:
            fa_sorted_lst=sorted(fa_dt.items(), key=lambda x: x[0].split('|')[0])
            nuc_array=np.array([])
            flag=0
            for ka, va in enumerate(fa_sorted_lst):
                if flag==0:
                    flag=1
                    nuc_array=np.array(np.fromstring(va[1], dtype='S1'))
                else:
                    nuc_array=np.vstack((nuc_array,np.fromstring(va[1], dtype='S1')))

            position_polymorphic = np.where(np.all(nuc_array== nuc_array[0, :],axis = 0)==False)[0]
            position_has_gap = np.where(np.any(nuc_array=='-',axis=0))[0]
            position_SNP = np.setdiff1d(position_polymorphic, position_has_gap)
            snp_columns = nuc_array[:,position_SNP]
            snp_pos_dt[align_file]=snp_columns

            if snp_wh_matrix_flag==0:
                snp_whole_matrix=snp_columns;
                snp_wh_matrix_flag=1
            else:
                snp_whole_matrix=np.hstack((snp_whole_matrix, snp_columns))

    write_pickle(output_path+'snp_pos.cpk',snp_pos_dt)

    with open(output_path+'SNP_whole_matrix.aln','wb') as outfile:
        for ind, isw in enumerate(snp_whole_matrix):
            write_in_fa( outfile, refSeqList[ind], isw.tostring() )

