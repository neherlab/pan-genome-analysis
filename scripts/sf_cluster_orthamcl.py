import os, sys, time
from collections import defaultdict, Counter
from sf_miscellaneous import times, load_pickle, read_fasta, write_pickle
from sf_cluster_protein import diamond_run

def ortha_mcl_run(output_path, threads, mcl_inflation):
    """ run orthAgogue and MCL """
    #orthAgogue_path='./orthagogue/'
    #os.system(''.join([orthAgogue_path,"orthAgogue -c ",threads," -i ",output_path,"query_matches.m8 -s '|' -O ",output_path,"ortha > ",output_path,"orthAgogu.log  2>&1"]))
    current_path=os.getcwd()
    os.chdir(output_path)
    print 'Running orthAgogue and mcl under ',current_path
    command_orthagogue=''.join(["orthAgogue \
        -c ",threads," -i query_matches.m8 -s '|' -u -O ./ortha \
        > ./orthAgogue.log  2>&1"])
    print 'command line orthAgogue: ', command_orthagogue
    os.system(command_orthagogue)
    #os.system(''.join(['mv ortha/all.abc ./']))
    command_mcl=''.join(['mcl ortha/all.abc --abc -o cluster.output -I ',str(mcl_inflation),' > mcl.log 2>&1'])
    print 'command line MCL:', command_mcl
    os.system(command_mcl)
    os.chdir(current_path)
    #os.system(''.join(['rm ',output_path,'all.abc']))

def orthagogue_singletons(path,origin_cluster_file,all_faa_file):
    """ add singletons from original MCL output """
    #from operator import or_
    all_faa_file="%s%s"%(path,all_faa_file)
    origin_cluster_file="%s%s"%(path,origin_cluster_file)
    all_cluster_file="%s%s"%(path,'allclusters.tsv')

    # loop over cluster_file, each line is one cluster tab delimited geneIDs (strain-locusTag)
    # generate union of all genes in all clusters excluding singletons
    with open(origin_cluster_file, 'rb') as infile:
        orthagogue_set= set.union(*[ set(iline.rstrip().split('\t')) for iline in infile])
        #orthagogue_set=reduce(or_, [ set(iline.rstrip().split('\t')) for iline in infile ])

    # read all geneIDs from all genes from all strains, determine singletons as set difference
    all_faa_set=set( read_fasta(all_faa_file).keys() )
    singletons=all_faa_set-orthagogue_set
    print ''.join(['#all genes: ', str(len(all_faa_set)), \
        '; #genes in orthagogue: ', str(len(orthagogue_set)), '; #singletons: ', str(len(singletons))])

    # append singleton clusters to a copy of the original file
    os.system('cp '+origin_cluster_file+' '+all_cluster_file)
    with open(all_cluster_file, 'a') as outputfile:
        for isi in singletons:
            outputfile.write(isi+'\n')

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
        with open(path+'allclusters.log', 'wb') as write_fn_lst:
            orthagogue_geneCount_lst=sorted( geneCluster_dt.iteritems(), key=itemgetter(1), reverse=True);
            for kd, vd in orthagogue_geneCount_lst:
                write_fn_lst.write('%s%s\n'%(kd, vd));

def diamond_orthamcl_cluster(
    path, threads, blast_cluster_file_path, roary_cluster_file_path,
    diamond_evalue, diamond_max_target_seqs, diamond_identity,
    diamond_query_cover, diamond_subject_cover, mcl_inflation):
    ''' 
    make all-against-all comparison using diamond
    THEN generate gene clusters followed by orthoMCL/orthagogue+MCL
    OR use the output of all-to-all blast comparison and orthoMCL/orthagogue+MCL
    OR use the output of roary
    params:
        path:                    path to directory including data and output
        threads:                 number of parallel threads used to run diamond
        blast_cluster_file_path: gene clusters by all-to-all blast 
                                 comparison and orthoMCL/orthagogue+MCL
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
            #os.system('cp '+output_path+dmd_ref_file+' '+output_path+dmd_query_file)
            ## run diamond
            diamond_run(output_path, dmd_ref_file, threads, diamond_evalue,
                diamond_max_target_seqs, diamond_identity, diamond_query_cover, diamond_subject_cover)
            ortha_mcl_run(output_path, threads, mcl_inflation)
            ## save singletons
            origin_cluster_file='cluster.output';
            orthagogue_singletons(output_path,origin_cluster_file,dmd_query_file)
            ## clean up diamond_query_file
            os.system(''.join(['rm ',output_path,'*faa']))
            all_cluster_file='allclusters.tsv';
            os.system('cp %sallclusters.tsv %sorthamcl-allclusters.tsv'%(output_path,output_path))
            parse_geneCluster(output_path,all_cluster_file)
        else: ## using user-given cluster file based on blast
            os.system('mkdir %s'%output_path)
            os.system('ln -sf %s %sclustered_proteins'%(blast_cluster_file_path, output_path))
            from operator import itemgetter
            ## create gene cluster from blast output
            with open(blast_cluster_file_path, 'rb') as infile:
                geneCluster_dt=defaultdict(list)
                for gid, iline in enumerate(infile):
                    column=[  ico.replace('_','|') for ico in iline.rstrip().split('\t') ]
                    clusterID="GC_%08d"%gid
                    gene_list=[ ico for ico in column ]
                    geneCluster_dt[clusterID]=[0,[],0]
                    num_stains=len( dict(Counter([ ivg.split('|')[0] for ivg in gene_list ])) )
                    num_gene=len(dict(Counter([ ivg for ivg in column])))
                    geneCluster_dt[ clusterID ][0]=num_stains
                    geneCluster_dt[ clusterID ][2]=num_gene
                    geneCluster_dt[ clusterID ][1]=gene_list
            write_pickle(output_path+'allclusters.cpk', geneCluster_dt)

            orthagogue_geneCount_lst=sorted( geneCluster_dt.iteritems(), key=itemgetter(1), reverse=True)
            with open(output_path+'allclusters.log', 'wb') as write_fn_lst:
                for kd, vd in orthagogue_geneCount_lst:
                    write_fn_lst.write('%s%s\n'%(kd, vd))
    else: ## using cluster files from roary
        os.system('mkdir %s'%output_path)
        os.system('ln -sf %s %sclustered_proteins'%(roary_cluster_file_path, output_path))
        with open(roary_cluster_file_path, 'rb') as cluster_external_file:
            with open(output_path+'allclusters.tsv', 'wb') as cluster_final_file:
                for cluster_line in cluster_external_file:
                     cluster_final_file.write( '%s\n'%'\t'.join([ gene_tag.replace('_','|') if '|' not in gene_tag else gene_tag for gene_tag in cluster_line.rstrip().split(': ')[1].split('\t')]) )
        all_cluster_file='allclusters.tsv';
        parse_geneCluster(output_path,all_cluster_file)
