import os, sys, time, glob
from distutils.spawn import find_executable
from collections import defaultdict, Counter
from sf_miscellaneous import times, read_fasta, load_pickle, write_pickle

def diamond_run(output_path, dmd_ref_file, threads,
    diamond_evalue, diamond_max_target_seqs, diamond_identity,
    diamond_query_cover, diamond_subject_cover, diamond_path, diamond_no_self_hits=0):
    """ run diamond using sensitive alignment mode """
    if diamond_path=='':
        diamond_path=find_executable('diamond')
    diam=diamond_path
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
                        ' -f 6 qseqid sseqid bitscore',
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
    for faa_file in glob.iglob(''.join([faa_path,'*faa'])):
        for gene_tag, seq in read_fasta(faa_file).iteritems():
            seq_length_dt[gene_tag]=len(seq)
    return seq_length_dt

def filter_hits_single(output_path, threads,
    input_prefix='', gather_seq_length_flag=0, no_filtering=True):
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

    if no_filtering==True:
        m8_fpath=os.path.abspath(output_path+m8_filename)
        filtered_hits_fpath=os.path.abspath(output_path+filtered_hits_filename)
        os.system(' '.join(['ln -s',m8_fpath,filtered_hits_fpath]))
    else:
        with open(output_path+m8_filename,'rb') as m8_file,\
            open(output_path+filtered_hits_filename,'wb') as abc_file:
            print('filter hits:')
            #using_BS=1; using_BAL=0
            for iline in m8_file:
                cols_list= iline.rstrip().split('\t')
                #query, ref, bit_score from column (0,1,-1)
                abc_file.write('%s\n'%'\t'.join([cols_list[ind] for ind in (0,1,-1)]))
                ## BAL scoring option not used
                # if using_BS:## using bscore
                #     abc_file.write('%s\n'%'\t'.join([cols_list[ind] for ind in (0,1,-1)]))
                # elif using_BAL:## using bscore/aln_len
                #     abc_file.write('%s\n'%'\t'.join([cols_list[0], cols_list[1], \
                #         str(round( float(cols_list[-1])/float(cols_list[3]),5 ))]))
        if input_prefix=='':
            print 'filter hits runtime:', times(start),'\n'
        else:
            print 'filter hits runtime for ',input_prefix,':', times(start),'\n'

def mcl_run(output_path, threads, mcl_inflation, input_prefix=''):
    """ """
    start = time.time()
    if input_prefix=='': ## default empty for abc file from all-vs.-all m8 file
        filtered_hits_filename='filtered_hits.abc'
    else:
        filtered_hits_filename='%s%s'%(input_prefix,'_filtered_hits.abc')

    command_mcl=''.join(['mcl ',output_path,filtered_hits_filename,' --abc ',\
                        '-o ',output_path,'allclusters.tsv -I ',str(mcl_inflation),\
                        ' -te ',str(threads),' > ',output_path,'mcl.log 2>&1'])
    os.system(command_mcl)
    print 'command line mcl:', command_mcl
    print 'mcl runtime:', times(start),'\n'


def cleanup_clustering(clustering_path):
    cwd = os.getcwd()
    os.chdir(clustering_path)
    os.system('rm -rf *.faa *.m8 *.abc *.output *dict*.cpk ./subproblem_cluster_seqs/')
    os.chdir(cwd)

def parse_geneCluster(input_fpath, output_fpath, cluster_log=False):
    """ store clusters as dictionary in cpk file """
    with open(input_fpath, 'rb') as infile:
        geneCluster_dt=defaultdict(list)
        for gid, iline in enumerate(infile,1): ##format: NC_022226|1-1956082:1956435
            col=iline.rstrip().split('\t')
            clusterID="GC%08d"%gid
            num_strains=len(dict(Counter([ ivg.split('|')[0] for ivg in col])).keys())
            num_genes=len(dict(Counter([ ivg for ivg in col])).keys())
            gene_mem=[ icol for icol in col ]
            geneCluster_dt[clusterID]=[num_strains,gene_mem,num_genes]
    write_pickle(output_fpath,geneCluster_dt)
    return geneCluster_dt

def roary_cluster_process(locus_tag_to_geneID_dict, input_fpath, output_fpath):
    """ """
    print input_fpath,output_fpath
    with open(input_fpath,'rb') as in_clusters:
        with open(output_fpath,'wb') as out_clusters:
            for cluster in in_clusters:
                cluster_line_lst=[]
                #out_clusters.write('%s\n'%'\t'.join([ locus_tag_to_geneID_dict[gene.split('.')[0]] for gene in cluster.rstrip().split(': ')[1].split('\t') ]) )
                for gene in cluster.rstrip().split(': ')[1].split('\t'):
                    if gene.split('.')[0] not in locus_tag_to_geneID_dict:
                        ## Roary also clusters non-CDS genes (skip and record them)
                        #print gene.split('.')[0], 'not in locus_tag_to_geneID_dict'
                        print gene.split('.')[0], ' non-coding gene from Roary file (skipped) '
                    else:
                        cluster_line_lst.append(locus_tag_to_geneID_dict[gene.split('.')[0]])
                if len(cluster_line_lst)!=0: # skip clusters on single non-CDS gene
                    out_clusters.write('%s\n'%'\t'.join(cluster_line_lst))

def process_orthofinder(input_fpath,output_fpath):
    """
    Orthogroups.txt: inputfile_orthogroups
    Orthogroups_UnassignedGenes.csv: inputfile_singleton
    """
    #outfilename='orthofinder_cluster.tsv'
    with open(input_fpath,'rb') as orthogroups,\
        open(output_fpath,'wb') as outfile:
        for iline in orthogroups:
            cluster_line='\t'.join([gene for gene in iline.rstrip().split(': ')[1].split(' ')])
            outfile.write('%s\n'%cluster_line)

def clustering_protein(path, folders_dict, threads,
    blast_fpath, roary_fpath, orthofinder_fpath, other_tool_fpath,
    diamond_evalue, diamond_max_target_seqs, diamond_identity,
    diamond_query_cover, diamond_subject_cover, diamond_path, mcl_inflation):
    '''
    Procedure: all-against-all protein comparison + hits filtering + mcl clustering
    By default: DIAMOND -> BS -> MCL
    Alternatives:
    1. Blastp output (user-provided) -> BS -> MCL
    2. Roary
    3. OrthoFinder
    4. Other tools.
    params:
        path:                    path to directory including data and output
        threads:                 number of parallel threads used to run diamond
        blast_fpath: gene clusters by all-vs-all blast
                                 comparison and other clusterings methods
        roary_fpath: gene clusters by roary
        diamond_max_target_seqs: Diamond setting: the maximum number of target sequences
                                  per query to keep alignments for. Defalut:600
                                  #strain * #max_duplication= 50*10= 500
    '''
    threads=str(threads)
    protein_path= folders_dict['protein_path']
    clustering_path= folders_dict['clustering_path']
    cluster_fpath= '%s%s'%(clustering_path,'allclusters.tsv')
    cluster_dt_cpk_fpath='%s%s'%(clustering_path,'allclusters.cpk')
    if any( i!='none' for i in [blast_fpath, roary_fpath, orthofinder_fpath, other_tool_fpath]):
        geneID_to_geneSeqID_dict= load_pickle('%sgeneID_to_geneSeqID.cpk'%path)
        locus_tag_to_geneID_dict= defaultdict(list)
        for geneID in geneID_to_geneSeqID_dict.keys():
            locus_tag=geneID.split('|')[1]
            locus_tag_to_geneID_dict[locus_tag]=geneID
    ## using standard pipeline (roary_fpath=='none')
    if all( i=='none' for i in [blast_fpath, roary_fpath, orthofinder_fpath, other_tool_fpath]):
        dmd_ref_file='reference.faa'
        ## prepare dmd_query_file (dmd_query_file is dmd_ref_file)
        os.system(''.join(['cat ',protein_path,'*faa > ',clustering_path,dmd_ref_file]))
        ## run diamond
        diamond_run(clustering_path, dmd_ref_file, threads, diamond_evalue,
            diamond_max_target_seqs, diamond_identity, diamond_query_cover, diamond_subject_cover, diamond_path )
        ## filtering hits via BS score
        filter_hits_single(clustering_path, threads)
        ## running mcl
        mcl_run(clustering_path, threads, mcl_inflation)
        ## clean up diamond_query_file
        os.system(''.join(['rm ',clustering_path,'*faa']))
    elif blast_fpath!='none': ## using user-given cluster file based on blast
            os.system(''.join(['cp ',blast_fpath,' ',clustering_path,'blastp.m8']))
            ## filtering hits via BS score
            filter_hits_single(clustering_path, threads, input_prefix='blastp')
            ## running mcl
            mcl_run(clustering_path, threads, mcl_inflation, input_prefix='blastp')
    elif roary_fpath!='none': ## using cluster files from roary
        roary_cluster_process(locus_tag_to_geneID_dict, roary_fpath, cluster_fpath)
        # with open(roary_fpath, 'rb') as cluster_external_file:
        #     with open(cluster_fpath, 'wb') as cluster_final_file:
        #         for cluster_line in cluster_external_file:
        #              cluster_final_file.write( '%s\n'%'\t'.join([ gene_tag.replace('_','|') if '|' not in gene_tag else gene_tag for gene_tag in cluster_line.rstrip().split(': ')[1].split('\t')]) )
    elif orthofinder_fpath!='none':
        process_orthofinder(orthofinder_fpath,cluster_fpath)
    elif other_tool_fpath!='none':
        os.system('cp %s %s'%(other_tool_fpath,cluster_fpath))
    cleanup_clustering(clustering_path)
    return parse_geneCluster(cluster_fpath, cluster_dt_cpk_fpath)
