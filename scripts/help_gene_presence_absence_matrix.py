import os
import argparse
from sf_geneCluster_align_makeTree import load_sorted_clusters
from sf_miscellaneous import read_fasta

parser = argparse.ArgumentParser(description='create the gene presence absence matrix',\
         usage=' python %(prog)s'+'  -in /yourPath/pan-genome-analysis/data/Escherichia_coli/ -out /yourPath/gene_presence_absence_matrix.csv')
parser.add_argument('-in', '--input_filepath',   type = str, required=True, help='/path/pan-genome-analysis/data/')
parser.add_argument('-out', '--output_filepath',   type = str, required=True, help='/outputPath/')
params = parser.parse_args()
input_filepath=params.input_filepath
output_filepath=params.output_filepath
if input_filepath[-1]!='/':
    input_filepath+='/'

def make_gene_presence_absence_matrix(input_filepath):
    os.chdir(input_filepath)
    gene_order= ','.join([gene.rstrip() for gene, content in load_sorted_clusters('./')])

    with open('./geneCluster/genePresence.aln') as inputf,\
         open(output_filepath,'wb') as outputf:
        outputf.write('accession,%s\n'%gene_order)
        for strain, genes in read_fasta(inputf).iteritems():
            outputf.write('%s,%s\n'%(strain,','.join(genes)))

make_gene_presence_absence_matrix(input_filepath)