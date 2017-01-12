import os, sys
from collections import defaultdict, Counter
from sf_miscellaneous import write_pickle
from sf_cluster_orthamcl import ortha_mcl_run, orthagogue_singletons

def parse_RNACluster(path,inputfile):
    """ store clusters as dictionary in cpk file """
    from operator import itemgetter
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
    make all-against-all comparison using diamond
    THEN generate RNA clusters followed by orthoMCL/orthagogue+MCL
    OR use the output of all-to-all blast comparison and orthoMCL/orthagogue+MCL
    OR use the output of roary
    params:
        path:                    path to directory including data and output
        threads:                 number of parallel threads used to run blastn

    '''
    threads=str(threads)
    input_path=path+'RNA_fna/'
    output_path=input_path
    ## prepare query & reference file (all against all): merge all fna files
    all_RNA_filename='all_RNAs'
    os.system(''.join(['cat ',input_path,'*fna > ',input_path,all_RNA_filename,'.fna']))
    ## blastn on all RNA nucleotides
    ### make database
    run_makeblastdb=''.join(['makeblastdb -in ',input_path,all_RNA_filename,'.fna -dbtype nucl -out ',input_path,all_RNA_filename])
    os.system(run_makeblastdb)
    ### run blastn
    run_blastn=''.join(['blastn -db ',input_path,all_RNA_filename,' -query ',input_path,all_RNA_filename,'.fna -out ',input_path,'query_matches.m8 -evalue 0.001 -outfmt 6 -max_target_seqs ',blastn_RNA_max_target_seqs,' -num_threads ',threads,' > ',input_path,'blastn-output.log'])
    os.system(run_blastn)
    
    ## run orthagogue and MCL
    ortha_mcl_run(output_path, threads, mcl_inflation)
    ## save singeltons
    origin_cluster_file='cluster.output'
    orthagogue_singletons(output_path,origin_cluster_file,'%s.fna'%(all_RNA_filename))
    all_cluster_file='allclusters.tsv'
    parse_RNACluster(output_path,all_cluster_file)
    ## remove database file
    os.system(''.join(['rm ',input_path,all_RNA_filename,'*']))
    os.system(''.join(['rm ',input_path,'*fna']))
    os.system(''.join(['rm ',input_path,'all.abc']))
