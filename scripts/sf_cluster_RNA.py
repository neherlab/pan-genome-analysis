import os, sys
from collections import defaultdict, Counter
from sf_miscellaneous import read_fasta, write_pickle
from sf_cluster_protein import filter_hits_single, mcl_run

def parse_RNACluster(path,inputfile):
    """ store clusters as dictionary in cpk file """
    inputfile="%s%s"%(path,inputfile)
    with open(inputfile, 'rb') as infile:
        RNACluster_dt=defaultdict(list)
        for gid, iline in enumerate(infile,1): ##format: NC_022226|1-1956082:1956435
            col=iline.rstrip().split('\t')
            clusterID="RC%05d"%gid
            num_strains=len(dict(Counter([ ivg.split('|')[0] for ivg in col])).keys())
            num_RNAs=len(dict(Counter([ ivg for ivg in col])).keys())
            RNA_mem=[ icol for icol in col ]
            RNACluster_dt[clusterID]=[num_strains,RNA_mem,num_RNAs]
    write_pickle(path+'allclusters.cpk',RNACluster_dt)

def RNA_cluster(path, threads,blastn_RNA_max_target_seqs, mcl_inflation ):
    '''
    use the output of all-to-all blast comparison and MCL
    params:
        path:                    path to directory including data and output
        threads:                 number of parallel threads used to run blastn

    '''
    threads=str(threads)
    RNA_path='%sRNA_fna/'%path
    ## prepare query & reference file (all against all): merge all fna files
    all_RNA_filename='all_RNAs'
    os.system(''.join(['cat ',RNA_path,'*fna > ',RNA_path,all_RNA_filename,'.fna']))
    ## blastn on all RNA nucleotides
    ### make database
    run_makeblastdb=''.join(['makeblastdb -in ',RNA_path,all_RNA_filename,'.fna -dbtype nucl -out ',RNA_path,all_RNA_filename])
    os.system(run_makeblastdb)
    print 'rRNA command line record (makeblastdb): ', run_makeblastdb
    ### run blastn
    run_blastn=''.join(['blastn -db ',RNA_path,all_RNA_filename,' -query ',RNA_path,all_RNA_filename,'.fna -out ',RNA_path,'query_matches.m8 -evalue 0.001 -outfmt 6 -max_target_seqs ',blastn_RNA_max_target_seqs,' -num_threads ',threads,' > ',RNA_path,'blastn-output.log'])
    os.system(run_blastn)
    print 'rRNA command line record (blastn): ', run_blastn
    ## filtering hits via BS score
    filter_hits_single(RNA_path, threads)
    ## running mcl
    mcl_run(RNA_path, threads, mcl_inflation)
    ## run orthagogue and MCL
    all_cluster_file='allclusters.tsv'
    parse_RNACluster(RNA_path,all_cluster_file)
    ## remove database file
    os.system(''.join(['rm ',RNA_path,all_RNA_filename,'*']))
    os.system(''.join(['rm ',RNA_path,'*fna']))
    os.system(''.join(['rm ',RNA_path,'*.m8']))
    os.system(''.join(['rm ',RNA_path,'*.abc']))
