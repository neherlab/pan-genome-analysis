import os, sys, glob
from collections import defaultdict
from sf_cluster_protein import diamond_run

def gather_precluster(protein_path, m8_filename):
    """  """
    pre_cluster_dict= defaultdict(set)
    precluster_filename=m8_filename.replace('.m8','.preclust')
    with open('%s%s'%(protein_path,m8_filename),'rb') as m8_file,\
        open('%s%s'%(protein_path,precluster_filename),'wb') as precluster_file:
        for iline in m8_file:
            query, ref = iline.split('\t')[0:2]
            pre_cluster_dict[query].add(ref)
        for k in pre_cluster_dict.keys():
            if k in pre_cluster_dict:
                v=pre_cluster_dict[k]
                for iv in v:## set update
                    pre_cluster_dict[k]= pre_cluster_dict[k]|pre_cluster_dict[iv]
                for iv in v:
                    del pre_cluster_dict[iv]
        for k, v in pre_cluster_dict.iteritems():
            ## k is the first gene, v contains other genes in that cluster
            precluster_file.write('%s\t%s\n'%(k,'\t'.join([ iv for iv in v])))

def cleanup_preclustering(preclustering_path):
    cwd = os.getcwd()
    os.chdir(preclustering_path)
    os.system('rm -rf *.m8 diamond*.log')
    empty_list= [preclust for preclust in glob.iglob('*.preclust') if os.stat(preclust).st_size==0]
    os.system('rm -rf %s'%(' '.join(empty_list)))
    os.chdir(cwd)

def preclustering_protein(path, folders_dict, threads,
    diamond_evalue, diamond_max_target_seqs, diamond_identity,
    diamond_query_cover, diamond_subject_cover):
    ''' preclustering for finding highly similar genes (duplicates) '''
    threads=str(threads)
    protein_path= folders_dict['protein_path']

    for faa_file in glob.iglob('%s%s'%(protein_path,'*.faa')):
        dmd_ref_file=faa_file.split('/')[-1]
        ## run diamond
        diamond_run(protein_path, dmd_ref_file, threads, diamond_evalue, diamond_max_target_seqs,
            diamond_identity, diamond_query_cover, diamond_subject_cover, diamond_no_self_hits=1)
        m8_filename='%s%s'%(dmd_ref_file.split('.faa')[0],'.m8')
        print 'Highly similar genes blastp record: ', m8_filename,'\n'
        gather_precluster(protein_path, m8_filename)

    ## clean-up
    cleanup_preclustering(protein_path)
