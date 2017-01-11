import operator
import os, sys, time, glob
from collections import defaultdict, Counter
from SF00_miscellaneous import times, read_fasta, write_pickle

def diamond_run(output_path, dmd_ref_file, threads,
    diamond_evalue, diamond_max_target_seqs,diamond_identity,
    diamond_query_cover, diamond_subject_cover, diamond_no_self_hits=0):
    """ run diamond using sensitive alignment mode """
    diam='./tools/diamond'
    print 'diamond inputfile:', dmd_ref_file
    input_prefix= dmd_ref_file.split('.faa')[0]
    if input_prefix=='reference':
        output_m8_filename='query_matches.m8'
    else:
        output_m8_filename= '%s%s'%(input_prefix,'.m8')
    
    makedb_command= ''.join([diam,' makedb -p ',threads,
                        ' --in ', output_path, dmd_ref_file,
                        ' -d ',output_path,'nr_',input_prefix,
                        '> ',output_path,'diamond_makedb_',input_prefix,'.log 2>&1'
                        ])
    start = time.time()
    os.system(makedb_command)
    print 'diamond build index (makedb):', times(start)
    print 'command line record:', makedb_command
    ## option to enable --no-self-hits
    if diamond_no_self_hits==0:
        option_no_self_hits=''
    else:
        option_no_self_hits=' --no-self-hits'

    blastp_command= ''.join([diam,' blastp --sensitive -p ',threads,
                        ' -e ', diamond_evalue,
                        ' --id ', diamond_identity,
                        ' --query-cover ', diamond_query_cover,
                        ' --subject-cover ', diamond_subject_cover,
                        option_no_self_hits,
                        ' -k ', diamond_max_target_seqs,
                        ' -d ',output_path,'nr_', input_prefix,
                        ' -q ',output_path,dmd_ref_file,
                        ' -o ',output_path,output_m8_filename,
                        ' -t ./  > ',output_path,'diamond_blastp_',input_prefix,'.log  2>&1'
                        ])
    start = time.time()
    os.system(blastp_command)
    print 'diamond alignment (blastp):', times(start)
    print 'diamond_max_target_seqs used: %s'%diamond_max_target_seqs
    print 'command line record:', blastp_command

    ## remove diamond binary database (dmnd) file
    os.system(''.join(['rm ',output_path,'nr_', input_prefix,'.dmnd']))

def gather_seq_length(faa_path):
    """ """
    seq_length_dt=defaultdict()
    for faa_file in glob.glob(''.join([faa_path,'*faa'])):
        for gene_tag, seq in read_fasta(faa_file).iteritems():
            seq_length_dt[gene_tag]=len(seq)
    return seq_length_dt

def filter_hits_single(output_path, threads,
    input_prefix='', gather_seq_length_flag=0):
    """  """
    start = time.time()
    if input_prefix=='':
        ## default empty for all-vs.-all m8 file
        m8_filename='query_matches.m8'
        filtered_hits_filename='filtered_hits.abc'
    else:
        ## add input prefix for sub_problem abc files in divide-and-conquer method
        m8_filename='%s%s'%(input_prefix,'.m8')
        filtered_hits_filename='%s%s'%(input_prefix,'_filtered_hits.abc')

    with open(output_path+m8_filename,'rb') as m8_file,\
        open(output_path+filtered_hits_filename,'wb') as abc_file:
        print('filter hits:')
        for iline in m8_file:
            cols_list= iline.rstrip().split('\t')
            #query, ref, bit_score from column (0,1,-1)
            abc_file.write('%s\n'%'\t'.join([cols_list[ind] for ind in (0,1,-1)]))

    if input_prefix=='':
        print 'filter hits runtime:', times(start),'\n'
    else:
        print 'filter hits runtime for ',input_prefix,':', times(start),'\n'

def mcl_run(output_path, threads, mcl_inflation, input_prefix=''):
    """ """
    start = time.time()
    if input_prefix=='':
        ## default empty for abc file from all-vs.-all m8 file
        filtered_hits_filename='filtered_hits.abc'
    else:
        filtered_hits_filename='%s%s'%(input_prefix,'_filtered_hits.abc')

    command_mcl=''.join(['mcl ',output_path,filtered_hits_filename,' --abc ',\
                        '-o ',output_path,'allclusters.tsv -I ',str(mcl_inflation),\
                        ' -te ',str(threads),' > ',output_path,'mcl.log 2>&1'])
    print 'command line mcl:', command_mcl
    print 'mcl runtime:', times(start),'\n'
    os.system(command_mcl)

def cleanup_clustering(clustering_path):
    cwd = os.getcwd()
    os.chdir(clustering_path)
    os.system('rm -rf *.faa *.m8 *.abc *.output *dict.cpk ./subproblem_cluster_seqs/')
    os.chdir(cwd)

def parse_geneCluster(input_fpath, output_fpath, cluster_log=False):
    """ store clusters as dictionary in cpk file """
    with open(input_fpath, 'rb') as infile:
        geneCluster_dt=defaultdict(list)
        for gid, iline in enumerate(infile): ##format: NC_022226|1-1956082:1956435
            col=iline.rstrip().split('\t')
            clusterID="gc%08d"%gid
            num_stains=len(dict(Counter([ ivg.split('|')[0] for ivg in col])).keys())
            num_genes=len(dict(Counter([ ivg for ivg in col])).keys())
            gene_mem=[ icol for icol in col ]
            geneCluster_dt[clusterID]=[num_stains,gene_mem,num_genes]
    write_pickle(output_fpath,geneCluster_dt)
    return geneCluster_dt

def clustering_protein(path, folders_dict, threads,
    blast_fpath, roary_fpath,
    diamond_evalue, diamond_max_target_seqs, diamond_identity,
    diamond_query_cover, diamond_subject_cover, mcl_inflation):
    '''
    Procedure: all-against-all protein comparison + hits filtering + mcl clustering
    By default: DIAMOND -> BSAL -> MCL
    Alternatives: 
    1. Blastp (user-provided) -> BSAL -> MCL
    2. Roary
    params:
        path:                    path to directory including data and output
        threads:                 number of parallel threads used to run diamond
        blast_fpath: gene clusters by all-vs-all blast 
                                 comparison and other clusterings methods
        roary_fpath: gene clusters by roary
        diamond_max_target_seqs: Diamond setting: the maximum number of target sequences 
                                  per query to keep alignments for. Defalut: 
                                  #strain * #max_duplication= 40*15= 600 
    '''
    threads=str(threads)
    protein_path= folders_dict['protein_path']
    clustering_path= folders_dict['clustering_path']
    cluster_fpath= '%s%s'%(clustering_path,'allclusters.tsv')
    cluster_dt_cpk_fpath='%s%s'%(clustering_path,'allclusters.cpk')
    
    ## using standard pipeline (roary_fpath=='none')
    if blast_fpath=='none' and roary_fpath=='none':
        dmd_ref_file='reference.faa'
        ## prepare dmd_query_file (dmd_query_file is dmd_ref_file)
        os.system(''.join(['cat ',protein_path,'*faa > ',clustering_path,dmd_ref_file]))
        ## run diamond
        diamond_run(clustering_path, dmd_ref_file, threads, diamond_evalue,
            diamond_max_target_seqs, diamond_identity, diamond_query_cover, diamond_subject_cover)
        ## filtering hits via BS score
        filter_hits_single(clustering_path, threads)
        ## running mcl
        mcl_run(clustering_path, threads, mcl_inflation)
        ## clean up diamond_query_file
        os.system(''.join(['rm ',clustering_path,'*faa']))
    elif blast_file_path!='none': ## using user-given cluster file based on blast
            os.system(''.join(['cp ',blast_fpath,' ',clustering_path,'blastp.m8']))
            ## filtering hits via BS score
            filter_hits_single(clustering_path, threads, input_prefix='blastp')
            ## running mcl
            mcl_run(clustering_path, threads, mcl_inflation, input_prefix='blastp')
    elif roary_file_path!='none': ## using cluster files from roary
        os.system('ln -sf %s %sclustered_proteins'%(roary_fpath, clustering_path))
        with open(roary_fpath, 'rb') as cluster_external_file:
            with open(cluster_fpath, 'wb') as cluster_final_file:
                for cluster_line in cluster_external_file:
                     cluster_final_file.write( '%s\n'%'\t'.join([ gene_tag.replace('_','|') if '|' not in gene_tag else gene_tag for gene_tag in cluster_line.rstrip().split(': ')[1].split('\t')]) )
    cleanup_clustering(clustering_path)
    return parse_geneCluster(cluster_fpath, cluster_dt_cpk_fpath)
