import os, glob, time, cPickle
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord  
from itertools import izip_longest

def times(start): 
    
    t = int(time.time()-start)
    time_record=' %d minutes %d seconds (%d s)'%(t/60, t%60, t)
    return time_record

def read_fasta(filename):
    fa_dt={}
    for record in SeqIO.parse(filename, "fasta"):
        fa_dt[record.id]=str(record.seq)
    return fa_dt

def gbk_seq(each_gb_path, gbk_name):  
    # get sequence from gbk file
    gb = each_gb_path + gbk_name
    genome = SeqIO.read(gb,'genbank')
    allseq = genome.seq
    return allseq
    
def write_in_fa(write_file,Id, seq):
    #write_file=open(filename, 'wb');
    sequences = SeqRecord(Seq(seq), Id)
    SeqIO.write(sequences, write_file, "fasta")
    #write_file.close()

def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)

def write_pickle(filename, data_out):
    write_file=open(filename, 'wb')
    cPickle.dump(data_out,write_file,protocol=2)
    write_file.close()

def write_json(data, file_name, indent=1):
    import json
    try:
        handle = open(file_name, 'w')
    except IOError:
        pass
    else:
        json.dump(data, handle, indent=indent)
        handle.close()

def organize_folder(main_path):
    """ create folders for pangenome analysis """
    ## Genbank folder
    command_mkdir='mkdir -p '
    folders_dict=defaultdict( str,
        gbk_path='input_GenBank/', RNA_path='RNA_fna/',
        protein_path='protein_faa/', nucleotide_path='nucleotide_fna/',
        clustering_path='protein_faa/diamond_matches/',
        cluster_seq_path='geneCluster/', tmp_core_seq_path='tmp_core/'
        )
    for k,v in folders_dict.iteritems():
        folders_dict[k]='%s%s'%(main_path,v)
    for key, folder_path in folders_dict.iteritems():
        os.system(''.join([command_mkdir,folder_path]))
    return folders_dict

def harmonize_filename(path,glob_list):
    """ """
    for fpath in glob_list:
        gbk_fname=fpath.split('/')[-1]
        if '-' in gbk_fname:
            ## force to replace '-' with '_' in GenBank filename
            gbk_fname= gbk_fname.replace('-','_')
            print ''.join(['filename harmonized: ',fpath,' -> ',gbk_fname]) 
            os.system(''.join(['mv ',fpath,' ',path,gbk_fname]))

def load_strains(path,gbk_present,folders_dict):
    """ load input strains in strain_list """
    if gbk_present==1:
        glob_item='.gbk'
        gbk_path=folders_dict['gbk_path']
        glob_list=glob.glob('%s*%s'%(path,glob_item))
        if len(glob_list)!=0:
            harmonize_filename(path,glob_list)
            strain_list= [i.split('/')[-1].split(glob_item)[0] for i in glob.glob('%s*%s'%(path,glob_item))]
            ## move gbk files in folder input_GenBank
            command_organize_gbk_input=''.join(['mv ',path,'*gbk ',gbk_path])
            os.system(command_organize_gbk_input)
        else:
            glob_list=glob.glob('%s*%s'%(gbk_path,glob_item))
            strain_list= [i.split('/')[-1].split(glob_item)[0] for i in glob_list]
    else:
        glob_item='.faa'
        glob_list=glob.glob('%s*%s'%(path,glob_item))
        if len(glob_list)!=0:
            harmonize_filename(path,glob_list)
            strain_list=[i.split(glob_item)[0] for i in glob.glob('%s*%s'%(path,glob_item))]
        else:
            glob_list=glob.glob('%s*%s'%(folders_dict['protein_path'],glob_item))
            strain_list= [i.split('/')[-1].split(glob_item)[0] for i in glob_list]
        command_organize_aa_input= 'mv %s*.faa %s'%(path,folders_dict['protein_path'])
        command_organize_nuc_input='mv %s*.fna %s'%(path,folders_dict['nucleotide_path'])
        os.system(command_organize_nuc_input)
        os.system(command_organize_aa_input)
    write_pickle(path+'strain_list.cpk', strain_list)

def build_sublist_multithread(threads, full_list, pad_val=None):
    """ divide a list into sub_list for multi-threading """
    return izip_longest(*[iter(full_list)]*threads, fillvalue=pad_val)

def multips(function_in_use, threads, full_list, *args, **kwargs):
    """ running multiple threads """
    from multiprocessing import Process
    procs = []
    kwargs_dict= kwargs.copy()
    managed_dict_return=0
    index_used=0
    if 'manager_needed_dicts' in kwargs_dict:
        from multiprocessing import Manager
        managed_dicts=tuple(Manager().dict() for i in range(len(kwargs_dict['manager_needed_dicts'])))
        managed_dict_return=1
        del kwargs_dict['manager_needed_dicts']
    if 'index_needed' in kwargs_dict:
        index_used=1
        del kwargs_dict['index_needed']

    for ind, sublist in enumerate(build_sublist_multithread(threads, full_list)):
        if sublist[-1] is None:
            sublist=(i for i in sublist if i is not None)
        args_content= (sublist,)+ args
        if index_used==1:
            kwargs_dict['index']= ind*threads
        if managed_dict_return==1:
            args_content+= managed_dicts
        #print 'args_content: ', args_content, 'kwargs_dict: ',kwargs_dict
        p= Process(target=function_in_use, args=args_content, kwargs=kwargs_dict)
        p.Daemon = True; p.start(); procs.append(p)

    for p in procs:
        p.join()

    if managed_dict_return==1:
        return managed_dicts
