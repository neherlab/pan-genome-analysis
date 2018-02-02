import os
from collections import defaultdict
from Bio import SeqIO
from sf_miscellaneous import read_fasta, write_in_fa, load_pickle, write_pickle
import random,time

def gbk_translation(strainID, gbk_fname, protein_fname, nucleotide_fname, RNA_fname,
    geneID_to_geneSeqID_dict,geneID_to_description_dict,
    RNAID_to_SeqID_dict, RNAID_to_description_dict,
    gene_aa_dict, gene_na_dict, RNA_dict, enable_RNA_clustering):
    '''
    extract sequences and meta informations of all genes in one reference genbank file
    params:
        - gbk_fname:        Genbank filename
        - protein_fname:  file into which all amino acid sequences are written
                            in fasta format. needed as input for diamond
        - nucleotide_fname: file into which all nucleotide sequences are written
                            in fasta format. needed for cluster sequences
        - RNA_fname: RNA nucleotide_sequences are written in fasta format.
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
        - enable_RNA_clustering: not cluster rRNA
    '''

    aa_sequence_file=open(protein_fname, 'wb')
    nu_sequence_file=open(nucleotide_fname, 'wb')

    if enable_RNA_clustering:
        RNA_sequence_file=open(RNA_fname, 'wb')

    contig_index=0
    check_CDS_passed=0
    for contig in SeqIO.parse(gbk_fname,'genbank'):
        contig_index+=1
        for feature in contig.features:
            if feature.type=='CDS':
                if not check_CDS_passed:
                    check_CDS_passed=1
                if 'product' in feature.qualifiers and 'translation' in feature.qualifiers :
                    if 'gene' in feature.qualifiers :
                        geneName='%s'%(feature.qualifiers['gene'][0]).replace(' ','_')
                    else: geneName=''
                    product=feature.qualifiers['product'][0]
                    annotation= '_'.join(product.split(' '))
                    trans_seq=feature.qualifiers['translation'][0]
                    if 'locus_tag' in feature.qualifiers:
                        locus_tag=feature.qualifiers['locus_tag'][0]
                    else:
                        locus_tag=feature.qualifiers['db_xref'][0].split(':')[1]
                    ## force to replace '-' with '_' in locus_tag
                    if '-' in locus_tag:
                        locus_tag=locus_tag.replace('-','_')
                    if "PROKKA" in locus_tag:
                        locus_tag=locus_tag.replace('PROKKA_','')
                    if '%s_'%strainID in locus_tag:
                        locus_tag=locus_tag.split('%s_'%strainID)[1]
                    ## geneID is composed of strain_name and locus_tag
                    ## Keeping '|' separator is important, which is used later in orthAgogue.
                    geneID= '%s|%s'%(strainID,locus_tag)
                    na_seq=str(feature.extract(contig.seq))
                    write_in_fa(aa_sequence_file, geneID, trans_seq)
                    write_in_fa(nu_sequence_file, geneID, na_seq)
                    gene_aa_dict[strainID][geneID]=trans_seq
                    gene_na_dict[strainID][geneID]=na_seq
                    # give tag 'gname:' to genes which have gene name and separate it from annotation
                    geneID_to_description_dict[geneID]={'geneName': geneName,
                                                        'contig': contig_index,
                                                        'annotation': annotation}
                    if geneName!='':
                        geneName='%s_'%geneName
                    geneID_to_geneSeqID_dict[geneID]='%s|%s-%d-%s%s'%(strainID,
                                                    locus_tag, contig_index,
                                                    geneName, annotation)
            elif enable_RNA_clustering and (feature.type=='rRNA'):
            #elif not enable_RNA_clustering and (feature.type=='rRNA' or feature.type=='tRNA'):
                if 'product' in feature.qualifiers:
                    geneName=''
                    product=feature.qualifiers['product'][0]
                    annotation= '_'.join(product.split(' '))
                    try:
                        locus_tag=feature.qualifiers['locus_tag'][0]
                    except: # make a random string when locus_tag absent
                        locus_tag=time.strftime('%S',time.gmtime())+str(random.randint(0,10000000))
                    if "PROKKA" in locus_tag:
                        locus_tag=locus_tag.replace('PROKKA_','')
                    if '%s_'%strainID in locus_tag:
                        locus_tag=locus_tag.split('%s_'%strainID)[1]
                    ## RNA is composed of strain_name and locus_tag
                    ## Keeping '|' separator is important, which is used later in orthAgogue.
                    RNAID= '%s|%s'%(strainID,locus_tag)
                    RNA_seq= str(feature.extract(contig.seq))
                    write_in_fa(RNA_sequence_file, RNAID, RNA_seq)
                    RNA_dict[strainID][RNAID]=RNA_seq
                    # give tag 'gname:' to genes which have gene name and separate it from annotation
                    RNAID_to_description_dict[RNAID]={
                                                    'geneName': '',
                                                    'contig': contig_index,
                                                    'annotation': annotation}
                    RNAID_to_SeqID_dict[RNAID]='%s|%s-%d-%s%s'%(strainID,
                                                    locus_tag, contig_index,
                                                    geneName, annotation)

    aa_sequence_file.close(); nu_sequence_file.close()

    return check_CDS_passed

def extract_sequences(path, strain_list, folders_dict, gbk_present, enable_RNA_clustering):
    '''
        go through all GenBank files and extract sequences and metadata for each one
    '''
    gbk_path= folders_dict['gbk_path']
    protein_path= folders_dict['protein_path']
    nucleotide_path= folders_dict['nucleotide_path']
    RNA_path= folders_dict['RNA_path']

    geneID_to_geneSeqID_file= '%sgeneID_to_geneSeqID.cpk'%path
    geneID_to_description_file= '%sgeneID_to_description.cpk'%path
    RNAID_to_SeqID_file= '%sRNAID_to_SeqID.cpk'%path
    RNAID_to_description_file= '%sRNAID_to_description.cpk'%path

    protein_dict_path= '%s%s'%(protein_path,'all_protein_seq.cpk')
    nucleotide_dict_path= '%s%s'%(nucleotide_path,'all_nucleotide_seq.cpk')
    RNA_dict_path= '%s%s'%(RNA_path,'all_RNA_seq.cpk')

    geneID_to_geneSeqID_dict= defaultdict()
    geneID_to_description_dict= defaultdict()
    RNAID_to_SeqID_dict= defaultdict()
    RNAID_to_description_dict= defaultdict()
    gene_aa_dict= defaultdict(dict)
    gene_na_dict= defaultdict(dict)
    RNA_dict= defaultdict(dict)

    if gbk_present:
        ## clean up folder when data from previous run exist.
        os.system('rm -rf '+protein_path+'*.faa')
        os.system('rm -rf '+nucleotide_path+'*.fna')
        missing_CDS_list=[] ## a list containing strains which have no CDS (if any)
        ## process gbk file
        for strainID in strain_list:
            gbk_fname= ''.join([gbk_path,strainID,'.gbk'])
            protein_fname= ''.join([protein_path,strainID,'.faa'])
            nucleotide_fname= ''.join([nucleotide_path,strainID,'.fna'])
            RNA_fname= ''.join([RNA_path,strainID,'.fna'])
            check_CDS_passed= gbk_translation(strainID, gbk_fname, protein_fname, nucleotide_fname, RNA_fname,
                geneID_to_geneSeqID_dict,geneID_to_description_dict,
                RNAID_to_SeqID_dict, RNAID_to_description_dict,
                gene_aa_dict, gene_na_dict, RNA_dict, enable_RNA_clustering)
            if not check_CDS_passed:
                missing_CDS_list.append(strainID)
        if len(missing_CDS_list)!=0:
            print 'Warning: no CDS found in the following genome/genomes, please double-check\n', missing_CDS_list
            exit()
    else:
        ## process fna/faa files if gbk files are not given.
        for strainID in strain_list:
            ## amino acid sequences
            protein_fname=''.join([protein_path,strainID,'.faa'])
            nucleotide_fname=''.join([nucleotide_path,strainID,'.fna'])
            aa_sequence_dt=read_fasta(protein_fname)
            na_sequence_dt=read_fasta(nucleotide_fname)
            ## prepare geneSeqID and description
            for geneID in aa_sequence_dt.keys():
                geneName, annotation= '',''
                geneID_to_geneSeqID_dict[geneID]=geneID
                geneID_to_description_dict[geneID]={'geneName': geneName,
                                                    'annotation': annotation}
                gene_aa_dict[strainID][geneID]=aa_sequence_dt[geneID]
                gene_na_dict[strainID][geneID]=na_sequence_dt[geneID]
    write_pickle(geneID_to_geneSeqID_file, geneID_to_geneSeqID_dict)
    write_pickle(geneID_to_description_file, geneID_to_description_dict)
    write_pickle(protein_dict_path,gene_aa_dict)
    write_pickle(nucleotide_dict_path,gene_na_dict)
    ## option: process RNA sequences for RNA_clustering
    if enable_RNA_clustering:
        write_pickle(RNA_dict_path,RNA_dict)
        write_pickle(RNAID_to_SeqID_file, RNAID_to_SeqID_dict)
        write_pickle(RNAID_to_description_file, RNAID_to_description_dict)
    return gene_aa_dict, gene_na_dict
