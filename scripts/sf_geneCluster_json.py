import numpy as np
import os, sys,glob, time
from sf_miscellaneous import load_pickle, write_pickle, read_fasta, write_in_fa ,times
from sf_geneCluster_align_makeTree import load_sorted_clusters
from operator import itemgetter
from collections import defaultdict, Counter

def consolidate_annotation(path,all_gene_names, geneID_to_description_dict):
    """
    determine a consensus annotation if annotation of reference sequences
    is conflicting. hypothetical protein annotations are avoided is possible
    """
    # count annotation and sort by decreasing occurence

    #geneSeqID based on geneID-contig-geneName:Annotation
    annotations=dict(Counter( [ geneID_to_description_dict[igi]['annotation'] for igi in all_gene_names]) )
    ## create split cluster generating new ID conflict

    annotations_sorted=sorted(annotations.iteritems(), key=itemgetter(1),reverse=True)
    # take majority, unless majority is hypothetical protein and other annotations exist
    majority = annotations_sorted[0][0]
    if majority=='hypothetical_protein' and len(annotations)!=1:
        majority=annotations_sorted[1][0]
    majority=majority.replace('\\','\\\\').replace('_',' ')
    # "#" to delimit key/count ; "@" to seperate various annotations
    allAnn=''.join(['%s#%s@'%(i_ann,j_ann) for i_ann,j_ann in annotations_sorted ])[:-1]
    allAnn=allAnn.replace('\\','\\\\').replace('_',' ')
    return allAnn,majority

def consolidate_geneName(path,all_gene_names, geneID_to_description_dict):
    """
    determine a consensus geneName if geneName of reference sequences
    is conflicting.
    """
    # count geneName and sort by decreasing occurence

    #geneSeqID based on geneID-conitg-geneName:geneName
    geneNames=dict(Counter( [ geneID_to_description_dict[igi]['geneName'] for igi in all_gene_names ]) )
    ## create split cluster generating new ID conflict

    geneNames_sorted=sorted(geneNames.iteritems(), key=itemgetter(1),reverse=True)
    majority = geneNames_sorted[0][0]
    if majority=='':
        if len(geneNames)!=1:
            majority=geneNames_sorted[1][0]
        else:
            majority='None'

    # "#" to delimit key/count ; "@" to seperate various geneNames
    all_geneName=''.join(['%s#%s@'%(i_ann,j_ann) if i_ann!='' else 'None#%s@'%j_ann  for i_ann,j_ann in geneNames_sorted ])[:-1]
    #print all_geneName, ' ?', majority
    return all_geneName, majority

def optional_geneCluster_properties(gene_list, optional_table_column):
    strain_to_locustag = dict([igl.split('|')[:2] for igl in gene_list])
    ## add optional table column
    if optional_table_column:
        new_column_attr='PAO1'
        if 'NC_002516' in strain_to_locustag:
            new_column_data=strain_to_locustag['NC_002516']
        else: new_column_data=''

    if optional_table_column:
        return ['"'+new_column_attr+'":"'+new_column_data+'"']
    else:
        []


def geneCluster_associations(associations, suffix='BA'):
    return ['"%s %s":%1.2f'%(k.split()[0], suffix , np.abs(v)) for k,v in associations.iteritems() if not np.isnan(v)]


def geneCluster_to_json(path, enable_RNA_clustering, store_locus_tag,
                        raw_locus_tag, optional_table_column):
    """
    create json file for gene cluster table visualzition
    input:  path to genecluster output
    output: geneCluster.json
    """
    # define path and make output directory
    geneCluster_path='%s%s'%(path,'geneCluster/')
    output_path='%s%s'%(path,'vis/')

    # open files
    geneClusterJSON_outfile=open(output_path+'geneCluster.json', 'wb')
    ##store locus_tags in a separate file for large dataset
    if store_locus_tag:
        locus_tag_outfile=open(path+'search_locus_tag.tsv', 'wb')


    ### load precomputed annotations, diversity, associations etc
    # load geneID_to_descriptions
    geneID_to_descriptions=load_pickle(path+'geneID_to_description.cpk')

    if enable_RNA_clustering:
        # load RNAID_to_description_file
        geneID_to_descriptions.update(load_pickle(path+'RNAID_to_description.cpk'))

    gene_diversity_Dt = load_pickle(geneCluster_path+'gene_diversity.cpk')
    ## load gain/loss event count dictionary
    dt_geneEvents = load_pickle(geneCluster_path+'dt_geneEvents.cpk')
    ## load association
    branch_associations_path = path+'branch_association.cpk'
    if os.path.isfile(branch_associations_path):
        branch_associations = load_pickle(branch_associations_path)
    else:
        branch_associations={}
    presence_absence_associations_path = path+'presence_absence_association.cpk'
    if os.path.isfile(presence_absence_associations_path):
        presence_absence_associations = load_pickle(presence_absence_associations_path)
    else:
        presence_absence_associations={}

    ## load list of clustered sorted by strain count
    sorted_genelist = load_sorted_clusters(path)

    geneClusterJSON_outfile.write('[')
    ## sorted_genelist: [(clusterID, [ count_strains,[memb1,...],count_genes]),...]
    for gid, (clusterID, gene) in enumerate(sorted_genelist):
        strain_count, gene_list, gene_count = gene
        # #print strain_count, gene_count
        if gid!=0: ## begin
            geneClusterJSON_outfile.write(',\n')

        ## annotation majority
        allAnn, majority_annotation = consolidate_annotation(path, gene_list, geneID_to_descriptions)

        ## geneName majority
        all_geneName, majority_geneName =  consolidate_geneName(path, gene_list, geneID_to_descriptions)

        ## extract gain/loss event count
        gene_event= dt_geneEvents[gid]

        ## average length
        seqs = read_fasta(geneCluster_path+'%s%s'%(clusterID,'.fna')).values()
        geneClusterLength = int(np.mean([ len(igene) for igene in seqs]))

        ## msa
        #geneCluster_aln='%s%s'%(clusterID,'_aa.aln')
        geneCluster_aln=clusterID

        ## check for duplicates
        if gene_count>strain_count:
            duplicated_state='yes'
            dup_list=[ ig.split('|')[0] for ig in gene_list]
            # "#" to delimit (gene/gene_count)key/value ; "@" to seperate genes
            # Counter({'g1': 2, 'g2': 1})
            dup_detail=''.join(['%s#%s@'%(kd,vd) for kd, vd in Counter(dup_list).iteritems() if vd>1 ])[:-1]
        else:
            duplicated_state='no';dup_detail=''

        ## locus_tag
        if raw_locus_tag: # make a string of all locus tags [1] in igl.split('|')
            all_locus_tags=' '.join([ igl.split('|')[1] for igl in gene_list ])
        else: # in addition to locus tag, keep strain name (but replace '|')
            all_locus_tags=' '.join([ igl.replace('|','_') for igl in gene_list ])

        ## optionally store locus tags to file, remove from geneClusterJSON
        if store_locus_tag:
            locus_tag_outfile.write('%s\t%s\n'%(clusterID,all_locus_tags))
            all_locus_tags=''

        ## default cluster json fields
        cluster_json_line=['"geneId":'+str(gid+1),
                            '"geneLen":'+str(geneClusterLength),
                            '"count":'+str(strain_count),
                            '"dupli":"'+duplicated_state+'"',
                            '"dup_detail":"'+dup_detail+'"',
                            '"ann":"'+majority_annotation+'"',
                            '"msa":"'+geneCluster_aln+'"',
                            '"divers":"'+gene_diversity_Dt[clusterID]+'"',
                            '"event":"'+str(gene_event)+'"',
                            '"allAnn":"'+allAnn+'"',
                            '"GName":"'+majority_geneName+'"',
                            '"allGName":"'+all_geneName+'"',
                            '"locus":"'+all_locus_tags+'"'
                            ]

        if optional_table_column:
            cluster_json_line.extend(optional_geneCluster_properties(gene_list,optional_table_column))
        if clusterID in branch_associations:
            cluster_json_line.extend(geneCluster_associations(branch_associations[clusterID], suffix='BA'))
        if clusterID in presence_absence_associations:
            cluster_json_line.extend(geneCluster_associations(presence_absence_associations[clusterID], suffix='PA'))

        #write file
        cluster_json_line=','.join(cluster_json_line)
        geneClusterJSON_outfile.write('{'+cluster_json_line+'}')

    # close files
    geneClusterJSON_outfile.write(']')
    geneClusterJSON_outfile.close()
    if store_locus_tag: locus_tag_outfile.close()