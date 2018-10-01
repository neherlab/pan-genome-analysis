import os, json
import argparse

parser = argparse.ArgumentParser(description='convert geneCluster.json into a csv file ',\
         usage=' python %(prog)s'+'  -in /path/geneCluster.json -out /path/geneCluster_stats.csv')
parser.add_argument('-in',  '--input_filepath',  type = str, required=True, help='')
parser.add_argument('-out', '--output_filepath', type = str, required=True, help='')
params = parser.parse_args()
input_filepath=params.input_filepath
output_filepath=params.output_filepath

def convert_cluster_statistics_csv(input_path, output_path):
    with open(input_path) as inputf,\
         open(output_path, 'w') as outputf:
        cluster_stats_list= json.load(inputf)
        header_list=['count','ann','dupli','GName','geneLen','event','msa','divers','geneId','allAnn','dup_detail','allGName','locus']
        header_description=['#strain count','major annotation','whether duplicated','major gene name','gene length','number of gene events','clusterID','diversity','gene ID','all annotations','duplication details','all gene names','all locus tags']
        outputf.write(','.join([ header for header in header_description ]))
        outputf.write('\n')
        for cluster_stats in cluster_stats_list:
            outputf.write(','.join([ str(cluster_stats[header]) for header in header_list ]))
            outputf.write('\n')

convert_cluster_statistics_csv(input_filepath, output_filepath)