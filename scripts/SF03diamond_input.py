import os; from collections import defaultdict
from Bio import SeqIO
from SF00miscellaneous import load_pickle, write_in_fa, write_pickle

def gbk_translation(each_gb_path, gb_file, output_filename, locusTag_to_geneId_Dt):
    '''
    extract meta informations of all genes in one reference genbank file
    params:
        - each_gb_path:     path to the set of reference sequences
                            used to construct the core genome
        - gb_file:          name of the reference to be analyzed
        - output_filename:  file into which all amino acid sequences are written
                            in fasta format. needed as input for diamond
        - locusTag_to_geneID: dictionary linking locusTag to geneID
                             modified in place.
    '''
    reference_gb = each_gb_path + gb_file
    strainName=gb_file.split('.')[0]
    geneFasta= each_gb_path + strainName +'.fastaCpk';
    gene_nucleotide_sequences=defaultdict()
    aa_sequence_file=open(output_filename, 'wb')
    flag=0; contig_index=0
    for contig in SeqIO.parse(reference_gb,'genbank'):
        contig_index+=1
        for feature in contig.features:
            loc=str(feature.location) # [543:1905](+): (1)0-base; (2)-: complement;
            # TODO: make consistent gene IDs

            if feature.type=='gene' and 'order' not in loc and 'pseudo' not in feature.qualifiers:
                rawLocation= str(feature.location);
                # complement(join(7099896..7099963,1..67))
                # TODO: refactor into function
                if 'join' in rawLocation:
                    flag_join=0; joinLen=0
                    for ira in rawLocation.split(','): # if 'join' includes more than 2 parts
                        lt_start=int(ira.split(':')[0].split('[')[1]);
                        lt_end=int(ira.split(':')[1].split(']')[0]);
                        joinLen+=lt_end-lt_start;
                        if flag_join==0:
                            lpos_start=lt_start; flag_join=1
                    # case1: join{[4395554:4396119](+), [0:11](+)}
                    # case2: join{[0:67](-), [7099895:7099963](-)} <in NZ_CP011369>
                    lpos_end=lpos_start+joinLen
                    print rawLocation, lpos_start, lpos_end;
                else:
                    lpos=str(feature.location).lstrip('[').split(']')[0].split(':');
                    if '<' in lpos[0]: #e.g.:[<363062:363251](-)
                        lpos_start=int(feature.location.nofuzzy_start)
                    else:
                        lpos_start=int(lpos[0])
                    if '>' in lpos[1]: #e.g.:[36550:>37009](-)
                        lpos_end=int(feature.location.nofuzzy_end)
                    else:
                        lpos_end=int(lpos[1]);
                geneID = strainName+'-'+str(contig_index)+'-'+str(lpos_start)
                gene_nucleotide_sequences[strainName+'-'+str(contig_index)+'-'+str(lpos_start)+'-'+str(lpos_end)] = feature.extract(genome.seq)
            elif feature.type=='CDS':
                if 'product' in feature.qualifiers and 'translation' in feature.qualifiers :
                    if 'gene' in feature.qualifiers :
                        geneName='%s_'%(feature.qualifiers['gene'][0]).replace(' ','_')
                    else: geneName=''
                    product=feature.qualifiers['product'][0];
                    trans_seq=feature.qualifiers['translation'][0];
                    #geneId='%s|%d-%d:%d-%s'%(gb_file.split('.gbk')[0], contig_index, lpos_start+1, lpos_end, '_'.join(product.split(' ')) )
                    geneId='%s|%d-%d:%d-%s%s'%(gb_file.split('.gbk')[0], contig_index, lpos_start+1, lpos_end, geneName, '_'.join(product.split(' ')) )
                    write_in_fa(aa_sequence_file, geneId, trans_seq)
                    ##if 'locus_tag' in feature.qualifiers :
                    locus_tag=feature.qualifiers['locus_tag'][0];
                    if "PROKKA" in locus_tag:
                        locus_tag=locus_tag.replace('PROKKA_','')
                        locusTag_to_geneId_Dt['%s-%s'%(strainName,locus_tag)]=geneId
                    else:
                        locusTag_to_geneId_Dt[locus_tag]=geneId
                flag=0;
    write_pickle(geneFasta, gene_nucleotide_sequences)
    aa_sequence_file.close()

def diamond_input(each_gb_path, strain_lst):
    # TODO: add comments
    folder=each_gb_path+ 'protein_fna/'
    os.system('mkdir '+folder);
    locusTag_to_geneId_file= each_gb_path +'locusTag_to_geneId.cpk'
    locusTag_to_geneId_Dt=defaultdict()
    for istrain in strain_lst:
        geach_name='%s'%(istrain.split('.')[0]);print geach_name;
        for_dmd_input=folder+'%s%s'%(geach_name,'.fna' );
        gbk_translation(each_gb_path,'%s%s'%(geach_name,'.gbk'), for_dmd_input, locusTag_to_geneId_Dt )
    write_pickle(locusTag_to_geneId_file, locusTag_to_geneId_Dt)


