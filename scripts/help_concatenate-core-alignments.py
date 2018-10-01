import os, json, gzip
import argparse
from Bio import SeqIO
from collections import defaultdict
from sf_miscellaneous import write_in_fa

parser = argparse.ArgumentParser(description='run this script in the ./scripts/ folder: concatenate all core gene alignments based on the core_geneList.txt file in the folder ./geneCluster/ ',\
         usage=' python %(prog)s'+'  -in ../data/TestSet/ -out /yourPath/concatenated_core_gene_alignments.fa')
parser.add_argument('-in',  '--input_filepath',  type = str, required=True, help='')
parser.add_argument('-out', '--output_filepath', type = str, required=True, help='')
params = parser.parse_args()
input_filepath=params.input_filepath
output_filepath=params.output_filepath

def concatenate_core_gene_alignments(input_path, output_path):
    core_genes_dt=defaultdict(str)
    with open(input_path+'/geneCluster/core_geneList.txt') as core_list:
    # all core gene alignments in FASTA files
        for gene in core_list:
            gene_path= input_path+'/vis/geneCluster/'+gene.rstrip()+'.gz'
            with gzip.open(gene_path, 'rb') as zip_file:
                for record in SeqIO.parse(zip_file, "fasta"):
                    #NC_018495-CM9_RS00390-1-hypothetical_protein
                    accession=record.id.split('-')[0]
                    core_genes_dt[accession]= '%s%s'%(core_genes_dt[accession], record.seq)

    with open(output_path,'wb') as output_file:
        for gene_id, gene_seq in core_genes_dt.iteritems():
            write_in_fa(output_file, gene_id, gene_seq)

concatenate_core_gene_alignments(input_filepath, output_filepath)
