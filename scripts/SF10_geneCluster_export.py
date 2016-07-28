import sys, glob
from SF00miscellaneous import load_pickle, write_pickle, read_fasta, write_in_fa
from SF06geneCluster_align_makeTree import load_sorted_clusters
from operator import itemgetter
from collections import defaultdict, Counter

def consolidate_annotation(all_gene_names):
    """
    determine a consensus annotation if annotation of reference sequences
    is conflicting. hypothetical protein annotations are avoided is possible
    """
    # count annotation and sort by decreasing occurence
    annotations=dict(Counter( [ igi.split('|')[1].split('-',2)[2] for igi in all_gene_names ]) )
    annotations_sorted=sorted(annotations.iteritems(), key=itemgetter(1),reverse=True)
    if  len(annotations)==0 :
        majority='hypothetical_protein'
        allAnn  ='hypothetical_protein'
    else:
        # take majority, unless majority is hypothetical protein and other annotations exist
        majority = annotations_sorted[0][0]
        if majority=='hypothetical_protein' and len(annotations)!=1:
            majority=annotations_sorted[1][0]
        # "#" to delimit key/count ; "@" to seperate various annotations
        allAnn=''.join(['%s#%s@'%(i_ann,j_ann) for i_ann,j_ann in annotations_sorted ])[:-1]
    return allAnn,majority


def geneCluster_to_json(path, species):
    """
    create json file for gene cluster table visualzition
    input:  path to genecluster output
    output: species-geneCluster.json
    """
    output_path='%s%s'%(path,'geneCluster/');
    gene_diversity_Dt=load_pickle(output_path+'gene_diversity.cpk')
    write_file_lst_json=open(output_path+species+'-geneCluster.json', 'wb')

    ## sorted clusters
    sorted_genelist= load_sorted_clusters(path, species)

    ## prepare geneId_Dt_to_locusTag
    locusTag_to_geneId_Dt=load_pickle(path+'locusTag_to_geneId.cpk')
    geneId_Dt_to_locusTag=defaultdict(list)
    geneId_Dt_to_locusTag={v:k for k,v in locusTag_to_geneId_Dt.items()}

    ## load gain/loss event count dictionary
    dt_geneEvents= load_pickle(output_path+'dt_geneEvents.cpk')

    write_file_lst_json.write('['); begin=0
    ## sorted_genelist: [(clusterID, [ count_strains,[memb1,...],count_genes]),...]
    for gid, (clusterID, gene) in enumerate(sorted_genelist):
        strain_count, gene_list, gene_count = gene
        if begin==0:
            begin=1
        else:
            write_file_lst_json.write(',\n')

        ## annotation majority
        allAnn,majority = consolidate_annotation(gene_list)
        ## extract gain/loss event count
        gene_event= dt_geneEvents[gid]

        ## average length
        geneLength_list= [ len(igene) for igene in read_fasta(output_path+'%s%s'%(clusterID,'.nu.fa')).values() ]
        geneClusterLength = sum(geneLength_list) / len(geneLength_list)
        #print geneLength_list,geneClusterLength

        ## msa
        geneCluster_aln='%s%s'%(clusterID,'.nu.fa')

        ## check for duplicates
        if gene_count>strain_count:
            duplicated_state='yes';
            dup_list=[ ig.split('|')[0] for ig in gene_list]
            # "#" to delimit (gene/gene_count)key/value ; "@" to seperate genes
            # Counter({'g1': 2, 'g2': 1})
            dup_detail=''.join(['%s#%s@'%(kd,vd) for kd, vd in dict(Counter(dup_list)).items() if vd>1 ])[:-1]
        else:
            duplicated_state='no';dup_detail=''

        ## locus_tag
        locus_tag_strain=' '.join([ geneId_Dt_to_locusTag[igl] for igl in gene_list ])
        #locus_tag_strain=' '.join([ '%s_%s'%(igl.split('|')[0],geneId_Dt_to_locusTag[igl]) for igl in gene[1][1] ])

        ## write json
        newline='{"geneId":%d,"geneLen":%d,"count": %d,"dupli":"%s","dup_detail": "%s","ann":"%s","msa":"%s","divers":"%s","event":"%s","allAnn":"%s","locus":"%s"}'
        write_file_lst_json.write(newline%(gid+1, geneClusterLength, strain_count, duplicated_state,
                                           dup_detail,majority, geneCluster_aln,
                                           gene_diversity_Dt[clusterID],gene_event, allAnn,locus_tag_strain))
    write_file_lst_json.write(']')
    write_file_lst_json.close()

