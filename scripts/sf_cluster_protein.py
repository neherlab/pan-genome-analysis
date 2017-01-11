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
                        '-o ',output_path,'cluster.output -I ',str(mcl_inflation),\
                        ' -te ',str(threads),' > ',output_path,'mcl.log 2>&1'])
    print 'command line mcl:', command_mcl
    print 'mcl runtime:', times(start),'\n'
    os.system(command_mcl)

def parse_geneCluster(path,inputfile, cluster_log=False):
    """ store clusters as dictionary in cpk file """
    from operator import itemgetter
    inputfile="%s%s"%(path,inputfile)
    with open(inputfile, 'rb') as infile:
        geneCluster_dt=defaultdict(list)
        for gid, iline in enumerate(infile): ##format: NC_022226|1-1956082:1956435
            col=iline.rstrip().split('\t')
            clusterID="GC_%08d"%gid
            geneCluster_dt[clusterID]=[0,[],0]
            ## num_stains
            geneCluster_dt[clusterID][0]=len(dict(Counter([ ivg.split('|')[0] for ivg in col])).keys())
            ## num_genes
            geneCluster_dt[clusterID][2]=len(dict(Counter([ ivg for ivg in col])).keys())
            ## gene members
            geneCluster_dt[clusterID][1]=[ icol for icol in col ]
    write_pickle(path+'allclusters.cpk',geneCluster_dt)

    if cluster_log==True:
        with open(path+'clusters.log', 'wb') as write_fn_lst:
            geneCount_lst=sorted( geneCluster_dt.iteritems(), key=itemgetter(1), reverse=True);
            for kd, vd in geneCount_lst:
                write_fn_lst.write('%s%s\n'%(kd, vd));

def clustering_protein_sequences(path, threads,
    blast_cluster_file_path, roary_cluster_file_path,
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
        blast_cluster_file_path: gene clusters by all-vs-all blast 
                                 comparison and other clusterings methods
        roary_cluster_file_path: gene clusters by roary
        diamond_max_target_seqs: Diamond setting: the maximum number of target sequences 
                                  per query to keep alignments for. Defalut: 
                                  #strain * #max_duplication= 40*15= 600 
    '''
    input_path=path+'protein_faa/';
    output_path=input_path+'diamond_matches/';
    threads=str(threads)
    ## using standard pipeline (roary_cluster_file_path=='none')
    if roary_cluster_file_path=='none':
        if blast_cluster_file_path=='none':
            dmd_ref_file='reference.faa'#dmd_query_file='query.faa'
            ## prepare dmd_query_file
            os.system('mkdir '+output_path)
            os.system('cat '+input_path+'*faa > '+output_path+dmd_ref_file)
            ## dmd_query_file is dmd_ref_file
            #os.system('cp '+output_path+dmd_ref_filedmd_query_file+' '+output_path+dmd_query_file)
            ## run diamond
            diamond_run(output_path, dmd_ref_file, threads, diamond_evalue,
                diamond_max_target_seqs, diamond_identity, diamond_query_cover, diamond_subject_cover)
            ## filtering hits via BS score
            filter_hits_single(output_path, threads)
            ## running mcl
            mcl_run(output_path, threads, mcl_inflation)
            cluster_file='allclusters.tsv'
            os.system(''.join(['mv ',output_path,'cluster.output',\
                            ' ',output_path,'allclusters.tsv']))
            ## clean up diamond_query_file
            os.system(''.join(['rm ',output_path,'*faa']))
            parse_geneCluster(output_path,cluster_file)
        else: ## using user-given cluster file based on blast
            os.system('mkdir %s'%output_path)
            os.system('cp %s %sblastp.m8'%(blast_cluster_file_path, output_path))
            ## filtering hits via BS score
            filter_hits_single(output_path, threads, input_prefix='blastp')
            ## running mcl
            mcl_run(output_path, threads, mcl_inflation, input_prefix='blastp')
            cluster_file='allclusters.tsv'
            os.system(''.join(['mv ',output_path,'cluster.output',\
                            ' ',output_path,'allclusters.tsv']))
            ## clean up diamond_query_file
            os.system(''.join(['rm ',output_path,'*faa']))
            parse_geneCluster(output_path,cluster_file)
    else: ## using cluster files from roary
        os.system('mkdir %s'%output_path)
        os.system('ln -sf %s %sclustered_proteins'%(roary_cluster_file_path, output_path))
        with open(roary_cluster_file_path, 'rb') as cluster_external_file:
            with open(output_path+'allclusters.tsv', 'wb') as cluster_final_file:
                for cluster_line in cluster_external_file:
                     cluster_final_file.write( '%s\n'%'\t'.join([ gene_tag.replace('_','|') if '|' not in gene_tag else gene_tag for gene_tag in cluster_line.rstrip().split(': ')[1].split('\t')]) )
        all_cluster_file='allclusters.tsv';
        parse_geneCluster(output_path,all_cluster_file)


#path='/ebio/ag-neher/share/users/wding/pan-genome-analysis/data/M_geni/'
#control(path)