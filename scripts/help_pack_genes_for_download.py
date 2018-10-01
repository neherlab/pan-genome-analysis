import os
import argparse
from sf_geneCluster_align_makeTree import load_sorted_clusters

parser = argparse.ArgumentParser(description='packing genes for downloading',\
         usage=' python %(prog)s'+'  -sn Escherichia_coli -path /path/pan-genome-analysis/data/')
parser.add_argument('-sn',  '--species_name', type = str, required=True, help='')
parser.add_argument('-path', '--analysis_folder',   type = str, required=True, help='/path/pan-genome-analysis/data/')
params = parser.parse_args()
species_name=params.species_name
analysis_folder=params.analysis_folder
if analysis_folder[-1]!='/':
    analysis_folder+='/'

def make_core_all_targz(species_name, analysis_folder):
    species_folder= analysis_folder+species_name
    cwd=os.getcwd()
    os.chdir(species_folder)

    ## packing core genes
    os.system('mkdir -p ./core_gene_alignments;')
    with open('./geneCluster/core_geneList.txt') as core_list:
     # all core gene alignments in FASTA files
        for gene in core_list:
            os.system('cp ./vis/geneCluster/'+gene.rstrip()+'.gz ./core_gene_alignments')
            os.system('cp ./vis/geneCluster/'+gene.rstrip().replace('_na','_aa')+'.gz ./core_gene_alignments')
        #os.system('gunzip ./core_gene_alignments/*')
        os.system('tar -zcf core_gene_alignments.tar.gz core_gene_alignments; rm -r ./core_gene_alignments; mv core_gene_alignments.tar.gz vis')

    ## packing all genes in pan-genome
    os.system('mkdir -p all_gene_alignments')
    for gene, content in load_sorted_clusters('./'):
        os.system('cp ./vis/geneCluster/'+gene.rstrip()+'_na_aln.fa.gz ./all_gene_alignments')
        os.system('cp ./vis/geneCluster/'+gene.rstrip()+'_aa_aln.fa.gz ./all_gene_alignments')
    #os.system('gunzip ./all_gene_alignments/*')
    os.system('tar -zcf all_gene_alignments.tar.gz all_gene_alignments; rm -r all_gene_alignments; mv all_gene_alignments.tar.gz ./vis')

make_core_all_targz(species_name, analysis_folder)