import sys, glob
from SF00miscellaneous import load_pickle, write_pickle, read_fasta, write_in_fa
from operator import itemgetter
from collections import defaultdict, Counter

def consolidate_annotation(all_gene_name):
    """ determine annotation majority """
    annotations=dict(Counter( [ igi.split('|')[1].split('-',2)[2] for igi in all_gene_name ]) )
    annotations_sorted=sorted(annotations.iteritems(), key=itemgetter(1),reverse=True)
    if  len(annotations)==0 : 
        majority='hypothetical_protein'
        allAnn  ='hypothetical_protein'
    else:
        majority = annotations_sorted[0][0]
        if majority=='hypothetical_protein' and len(annotations)!=1: 
            majority=annotations_sorted[1][0]
        # "#" to delimit key/value ; "@" to seperate various annotations
        allAnn=''.join(['%s#%s@'%(i_ann,j_ann) for i_ann,j_ann in annotations_sorted ])[:-1]
    return allAnn,majority

def create_genePresence(dt_strainGene, totalStrain, set_totalStrain, all_gene_name):
    """ create dt for gene presence/absence string """
    set_sharedStrain=set([ igl.split('|')[0] for igl in all_gene_name])
    for ist in set_sharedStrain: 
        dt_strainGene[ist].append('1') 
    if totalStrain!=len(set_sharedStrain):
        for ist0 in set_totalStrain-set_sharedStrain:
            dt_strainGene[ist0].append('0')
    return dt_strainGene

def geneCluster_to_json(path, species,present_cf='0'):
    """
    create json file for gene cluster table visualzition
    input: pickled strain list (strain_list.cpk) ,pickled gene cluster file (species-orthamcl-allclusters.cpk), pickled gene diversity (gene_diversity.cpk) etc. 
    output: species-geneCluster.json
    """
    geneClusterPath='%s%s'%(path,'protein_fna/diamond_matches/')
    alnFilePath='%s%s'%(path,'geneCluster/'); output_path=alnFilePath  
    gene_diversity_Dt=load_pickle(alnFilePath+'gene_diversity.cpk')
    find_clusterFaName_Dt=load_pickle(alnFilePath+'find_clusterFaName_final.cpk')
    write_file_lst_json=open(output_path+species+'-geneCluster.json', 'wb')

    ## load strain list and prepare for gene presence/absence
    strain_list= load_pickle(path+'strain_list.cpk')
    set_totalStrain=set([ istrain for istrain in strain_list ])
    totalStrain=len(set_totalStrain)
    dt_strainGene= defaultdict(list)
    dt_geneID_to_geneName= defaultdict(list)
    dt_geneID_diversity=defaultdict(list)

    ## create core gene list
    corelist=[];
    diamond_geneCluster_dt=load_pickle(geneClusterPath+species+'-orthamcl-allclusters_final.cpk')
    with open(output_path+'core_geneList.txt','wb') as outfile:
        for kg, vg in diamond_geneCluster_dt.iteritems():
            if vg[0]==totalStrain and vg[2]==totalStrain:
                coreGeneName=find_clusterFaName_Dt[kg].replace('.fa','.aln')
                outfile.write(coreGeneName+'\n')
                corelist.append(coreGeneName)
        write_pickle(output_path+'core_geneList.cpk',corelist)

    ## sort by strain_counts
    #sorted_genelist=sorted(diamond_geneCluster_dt, key=lambda x: int(x[1][0]), reverse=True)
    from operator import itemgetter
    sorted_genelist=sorted(diamond_geneCluster_dt.iteritems(), key=lambda (k,v): itemgetter(0)(v), reverse=True) 

    ## prepare geneId_Dt_to_locusTag   
    locusTag_to_geneId_Dt=load_pickle(path+'locusTag_to_geneId.cpk')
    geneId_Dt_to_locusTag=defaultdict(list)
    geneId_Dt_to_locusTag={v:k for k,v in locusTag_to_geneId_Dt.items()}

    write_file_lst_json.write('['); begin=0
    for gid, gene in enumerate(sorted_genelist): 
    # data structure: [ ref,[ count,[ann_1,ann_2,...ann_n] ]
        if begin==0: 
            begin=1
        else:
            write_file_lst_json.write(',\n')
        ref_gene=gene[0]
        ref_strainName,full_ann=ref_gene.split('|')
        st, en=full_ann.split('-')[1].split(':')
        st, en=int(st), int(en)

        ## annotation majority
        allAnn,majority = consolidate_annotation(gene[1][1])

        ## average length
        geneLength_list= [ len(igene) for igene in read_fasta(output_path+find_clusterFaName_Dt[gene[0]]).values() ]
        geneClusterLength = sum(geneLength_list) / len(geneLength_list)
        #print geneLength_list,geneClusterLength

        ## msa
        geneCluster_aln=find_clusterFaName_Dt[gene[0]].replace('.nu.fa','.aa.aln')
        
        ## gene presence 
        dt_strainGene=create_genePresence(dt_strainGene, totalStrain, set_totalStrain, gene[1][1])

        ## geneID_to_geneName
        dt_geneID_to_geneName[gid+1]=geneCluster_aln.replace('.aa.aln','-core-event-tr.json')

        ## dt_geneID_diversity (Franz's plot)
        # dt_geneID_diversity[gid+1]=gene_diversity_Dt[ref_gene] 

        ## duplicate
        if len(gene[1])==2:
            geneCount=len(gene[1][1]) # no geneCount record stored
        else: 
            geneCount=gene[1][2]
        if geneCount>gene[1][0]: 
            duplicated_state='yes';
            dup_list=[ ig.split('|')[0] for ig in gene[1][1] ]
            # "#" to delimit (gene/geneCount)key/value ; "@" to seperate genes
            # Counter({'g1': 2, 'g2': 1})
            dup_detail=''.join(['%s#%s@'%(kd,vd) for kd, vd in dict(Counter(dup_list)).items() if vd>1 ])[:-1]
        else:
            duplicated_state='no';dup_detail=''

        ## locus_tag
        locus_tag_strain=' '.join([ geneId_Dt_to_locusTag[igl] for igl in gene[1][1] ])
        #locus_tag_strain=' '.join([ '%s_%s'%(igl.split('|')[0],geneId_Dt_to_locusTag[igl]) for igl in gene[1][1] ])

        ## write json 
        newline='{"geneId":%d,"geneLen":%d,"count": %d,"dupli":"%s","dup_detail": "%s","ann":"%s","msa":"%s","divers":"%s","allAnn":"%s","locus":"%s"}'
        write_file_lst_json.write(newline%(gid+1, geneClusterLength, gene[1][0], duplicated_state,dup_detail,majority, geneCluster_aln, gene_diversity_Dt[ref_gene],allAnn,locus_tag_strain))
    write_file_lst_json.write(']')
    write_file_lst_json.close()

    with open(output_path+'genePresence.aln','wb') as presence_outfile:
        for istkey in dt_strainGene: 
            dt_strainGene[istkey]=''.join(dt_strainGene[istkey] )
            write_in_fa( presence_outfile, istkey, dt_strainGene[istkey] )

    write_pickle(output_path+'dt_genePresence.cpk', dt_strainGene)
    write_pickle(output_path+'dt_geneEvent.cpk', dt_geneID_to_geneName)
    write_pickle(output_path+'dt_geneID_diversity.cpk', dt_geneID_diversity)