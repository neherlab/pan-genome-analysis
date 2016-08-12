import os; from collections import defaultdict
from Bio import SeqIO
from SF00_miscellaneous import load_pickle, write_in_fa, write_pickle

def gbk_translation(each_gbk_path,nucleotide_dict_path, gb_file, output_filename, geneID_to_geneSeqID_dict, geneID_to_description_dict):
    '''
    extract sequences and meta informations of all genes in one reference genbank file
    params:
        - each_gbk_path:    path to the set of reference sequences used to construct
                            the core genome
        - nucleotide_dict_path:
                            path to the cPickled dicts of all nucleotide sequences 
                            for each genome
        - gb_file:          name of the reference to be analyzed
        - output_filename:  file into which all amino acid sequences are written
                            in fasta format. needed as input for diamond
        - geneID_to_geneSeqID_dict: dictionary linking geneID to gene sequence ID
                            modified in place (key: geneID; value: geneSeqID )
        - geneID_to_description_dict: dictionary linking geneID to description info
                            modified in place (key: geneID; value: a dict including
                            information on contig_index, annotation or more)
    '''
    
    reference_gb = '%s%s'%(each_gbk_path, gb_file)
    strainName=gb_file.split('.')[0]
    gene_nuc_seq_dict= '%s%s_gene_nuc_dict.cpk'%(nucleotide_dict_path, strainName)
    gene_nucleotide_sequences=defaultdict()
    aa_sequence_file=open(output_filename, 'wb')
    contig_index=0
    for contig in SeqIO.parse(reference_gb,'genbank'):
        contig_index+=1
        for feature in contig.features:
            if feature.type=='CDS':
                if 'product' in feature.qualifiers and 'translation' in feature.qualifiers :
                    if 'gene' in feature.qualifiers :
                        geneName='%s'%(feature.qualifiers['gene'][0]).replace(' ','_')
                    else: geneName='None'
                    product=feature.qualifiers['product'][0]
                    annotation= '_'.join(product.split(' '))
                    trans_seq=feature.qualifiers['translation'][0];
                    locus_tag=feature.qualifiers['locus_tag'][0];
                    if "PROKKA" in locus_tag:
                        locus_tag=locus_tag.replace('PROKKA_','')
                    if '%s_'%strainName in locus_tag:
                        locus_tag=locus_tag.split('%s_'%strainName)[1]
                    ## geneID is composed of strain_name and locus_tag
                    ## Keeping '|' separator is important, which is used later in orthAgogue.
                    geneID= '%s|%s'%(strainName,locus_tag)
                    write_in_fa(aa_sequence_file, geneID, trans_seq)
                    # give tag 'gname:' to genes which have gene name and separate it from annotation 
                    geneID_to_description_dict[geneID]={'geneName': geneName,
                                                        'contig': contig_index,
                                                        'annotation': annotation}                    
                    if geneName!='None':
                        geneName='%s_'%geneName
                    geneID_to_geneSeqID_dict[geneID]='%s|%s-%d-%s%s'%(strainName,
                                                    locus_tag, contig_index,
                                                    geneName, annotation)

                    gene_nucleotide_sequences[geneID] = feature.extract(contig.seq)
    write_pickle(gene_nuc_seq_dict, gene_nucleotide_sequences)
    aa_sequence_file.close()

def diamond_input(path, strain_lst):
    '''
        go through all GenBank files and extract sequences and metadata for each one
    ''' 
    each_gbk_path='%s%s'%(path,'input_GenBank/')
    os.system('mkdir %s;mv %s*gbk %s'%(each_gbk_path,path,each_gbk_path))
    protein_folder='%s%s'%(path,'protein_faa/')
    os.system('mkdir %s'%protein_folder)
    nucleotide_dict_path= '%s%s'%(path,'nucleotide_fna/')
    os.system('mkdir %s'%nucleotide_dict_path)
    geneID_to_geneSeqID_file= path+'geneID_to_geneSeqID.cpk'
    geneID_to_geneSeqID_dict= defaultdict()
    geneID_to_description_file= path+'geneID_to_description.cpk'
    geneID_to_description_dict= defaultdict()
    for istrain in strain_lst:
        strain_name='%s'%(istrain.split('.')[0])
        diamond_input_fname=protein_folder+'%s%s'%(strain_name,'.faa');
        gbk_translation(each_gbk_path,nucleotide_dict_path,'%s%s'%(strain_name,'.gbk'), diamond_input_fname, geneID_to_geneSeqID_dict, geneID_to_description_dict )
    write_pickle(geneID_to_geneSeqID_file, geneID_to_geneSeqID_dict)
    write_pickle(geneID_to_description_file, geneID_to_description_dict)
