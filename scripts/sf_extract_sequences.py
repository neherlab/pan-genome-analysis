import os
from collections import defaultdict
from Bio import SeqIO
from SF00_miscellaneous import read_fasta, write_in_fa, load_pickle, write_pickle

def gbk_translation(gbk_path,nucleotide_dict_path, gb_file,
    output_filename, output_filename2, geneID_to_geneSeqID_dict,geneID_to_description_dict,
    RNAID_to_SeqID_dict, RNAID_to_description_dict, disable_RNA_clustering):
    '''
    extract sequences and meta informations of all genes in one reference genbank file
    params:
        - gbk_path:    path to the set of reference sequences used to construct
                            the core genome
        - nucleotide_dict_path:
                            path to the cPickled dicts of all nucleotide sequences 
                            for each genome
        - gb_file:          name of the reference to be analyzed
        - output_filename:  file into which all amino acid sequences are written
                            in fasta format. needed as input for diamond
        - output_filename2: RNA nucleotide_sequences are written in fasta format.
                            Needed as RNA_blast_input
        - geneID_to_geneSeqID_dict: dictionary linking geneID to gene sequence ID
                            modified in place (key: geneID; value: geneSeqID )
        - geneID_to_description_dict: dictionary linking geneID to description info
                            modified in place (key: geneID; value: a dict including
                            information on contig_index, annotation or more)
        - RNAID_to_SeqID_dict: dictionary linking RNAID to RNA sequence ID
                            modified in place (key: RNAID; value: SeqID )
        - RNAID_to_description_dict: dictionary linking RNAID to description info
                            modified in place (key: RNAID; value: a dict including
                            information on contig_index, annotation or more)
        - disable_RNA_clustering: not cluster rRNA and tRNA (default: 0 -> cluster RNAs)
    '''
    
    reference_gb = '%s%s'%(gbk_path, gb_file)
    strainName=gb_file.split('.gbk')[0]
    gene_nuc_seq_dict= '%s%s_gene_nuc_dict.cpk'%(nucleotide_dict_path, strainName)
    gene_nucleotide_sequences=defaultdict()
    aa_sequence_file=open(output_filename, 'wb')

    if disable_RNA_clustering==0:
        RNA_nuc_seq_dict= '%s%s_RNA_nuc_dict.cpk'%(nucleotide_dict_path, strainName)
        RNA_nucleotide_sequences=defaultdict()
        RNA_sequence_file=open(output_filename2, 'wb')

    contig_index=0
    for contig in SeqIO.parse(reference_gb,'genbank'):
        contig_index+=1
        for feature in contig.features:
            if feature.type=='CDS':
                if 'product' in feature.qualifiers and 'translation' in feature.qualifiers :
                    if 'gene' in feature.qualifiers :
                        geneName='%s'%(feature.qualifiers['gene'][0]).replace(' ','_')
                    else: geneName=''
                    product=feature.qualifiers['product'][0]
                    annotation= '_'.join(product.split(' '))
                    trans_seq=feature.qualifiers['translation'][0];
                    locus_tag=feature.qualifiers['locus_tag'][0];
                    ## force to replace '-' with '_' in locus_tag
                    if '-' in locus_tag:
                        locus_tag=locus_tag.replace('-','_')
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
                    if geneName!='':
                        geneName='%s_'%geneName
                    geneID_to_geneSeqID_dict[geneID]='%s|%s-%d-%s%s'%(strainName,
                                                    locus_tag, contig_index,
                                                    geneName, annotation)

                    gene_nucleotide_sequences[geneID] = feature.extract(contig.seq)
            elif not disable_RNA_clustering and (feature.type=='rRNA' or feature.type=='tRNA'):
                if 'product' in feature.qualifiers:
                    geneName=''
                    product=feature.qualifiers['product'][0]
                    annotation= '_'.join(product.split(' '))
                    locus_tag=feature.qualifiers['locus_tag'][0];
                    if "PROKKA" in locus_tag:
                        locus_tag=locus_tag.replace('PROKKA_','')
                    if '%s_'%strainName in locus_tag:
                        locus_tag=locus_tag.split('%s_'%strainName)[1]
                    ## RNA is composed of strain_name and locus_tag
                    ## Keeping '|' separator is important, which is used later in orthAgogue.
                    RNAID= '%s|%s'%(strainName,locus_tag)
                    RNA_seq= str(feature.extract(contig.seq))
                    write_in_fa(RNA_sequence_file, RNAID, RNA_seq)
                    # give tag 'gname:' to genes which have gene name and separate it from annotation
                    RNAID_to_description_dict[RNAID]={
                                                    'geneName': '',
                                                    'contig': contig_index,
                                                    'annotation': annotation}
                    RNAID_to_SeqID_dict[RNAID]='%s|%s-%d-%s%s'%(strainName,
                                                    locus_tag, contig_index,
                                                    geneName, annotation)
                    RNA_nucleotide_sequences[RNAID] = RNA_seq

    write_pickle(gene_nuc_seq_dict, gene_nucleotide_sequences)
    if disable_RNA_clustering==0:
        write_pickle(RNA_nuc_seq_dict, RNA_nucleotide_sequences)
    aa_sequence_file.close()

def diamond_input(path, strain_lst, gbk_present, disable_RNA_clustering):
    '''
        go through all GenBank files and extract sequences and metadata for each one
    '''
    gbk_path='%s%s'%(path,'input_GenBank/')
    protein_folder='%s%s'%(path,'protein_faa/')
    os.system('mkdir %s'%protein_folder)
    nuc_seq_folder='%s%s'%(path,'nucleotide_fna/')
    nucleotide_dict_path=nuc_seq_folder
    os.system('mkdir %s'%nucleotide_dict_path)
    RNA_folder='%s%s'%(path,'RNA_fna/')
    os.system('mkdir %s'%RNA_folder)    
    ## CDS
    geneID_to_geneSeqID_file= path+'geneID_to_geneSeqID.cpk'
    geneID_to_geneSeqID_dict= defaultdict()
    geneID_to_description_file= path+'geneID_to_description.cpk'
    geneID_to_description_dict= defaultdict()
    ## RNA
    RNAID_to_SeqID_file= path+'RNAID_to_SeqID.cpk'
    RNAID_to_SeqID_dict= defaultdict()
    RNAID_to_description_file= path+'RNAID_to_description.cpk'
    RNAID_to_description_dict= defaultdict()
    if gbk_present==1:
        ## process gbk file
        for strain_name in strain_lst:
            diamond_input_fname=protein_folder+'%s%s'%(strain_name,'.faa')
            RNA_blast_input_fname=RNA_folder+'%s%s'%(strain_name,'.fna')
            gbk_translation(gbk_path,nucleotide_dict_path,'%s%s'%(strain_name,'.gbk'),\
                    diamond_input_fname, RNA_blast_input_fname,\
                    geneID_to_geneSeqID_dict,geneID_to_description_dict,\
                    RNAID_to_SeqID_dict, RNAID_to_description_dict, disable_RNA_clustering)
    else:
        ## process fna/faa files if gbk files are not given.
        command_organize_aa_input= 'mv %s*.faa %s'%(path,protein_folder)
        command_organize_nuc_input='mv %s*.fna %s'%(path,nuc_seq_folder)
        os.system(command_organize_nuc_input)
        os.system(command_organize_aa_input)
        for strain_name in strain_lst:
            ## amino acid sequences
            diamond_input_fname=protein_folder+'%s%s'%(strain_name,'.faa')
            aa_sequence_dt=read_fasta(diamond_input_fname)
            ## nucleotide sequences
            nuc_sequence_fname='%s%s%s'%(nucleotide_dict_path,strain_name,'.fna')
            nu_sequence_dt=read_fasta(nuc_sequence_fname)
            gene_nuc_seq_dict= '%s%s_gene_nuc_dict.cpk'%(nucleotide_dict_path, strain_name)
            write_pickle(gene_nuc_seq_dict, nu_sequence_dt)
            ## prepare geneSeqID and description
            for geneID in aa_sequence_dt.keys():
                geneName, annotation= '',''
                geneID_to_geneSeqID_dict[geneID]=geneID
                geneID_to_description_dict[geneID]={'geneName': geneName,
                                                    'annotation': annotation}
    write_pickle(geneID_to_geneSeqID_file, geneID_to_geneSeqID_dict)
    write_pickle(geneID_to_description_file, geneID_to_description_dict)
    ## option: process RNA sequences for RNA_clustering
    if disable_RNA_clustering==0:
        write_pickle(RNAID_to_SeqID_file, RNAID_to_SeqID_dict)
        write_pickle(RNAID_to_description_file, RNAID_to_description_dict)
