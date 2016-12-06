import os, glob
import numpy as np, gzip, cPickle
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord  

def times(start): 
    import time
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

def load_strains(path):
    """ load input strains in strain_list """
    gbk_path='%s%s'%(path,'input_GenBank/')
    ## move gbk files in folder input_GenBank
    os.system('mkdir %s;mv %s*gbk %s'%(gbk_path,path,gbk_path))
    strain_list=[]
    for gbk_filepath in glob.glob(gbk_path+'*gbk'):
        gbk_filename=gbk_filepath.split('/')[-1]
        ## harmonize GenBank file name 
        ## force '-' to be replaced as '_' in GenBank filename
        if '-' in gbk_filename:
            new_gbk_filename=gbk_filename.replace('-','_')
            print('Filename harmonized: ',\
                gbk_filename,' -> ', new_gbk_filename.split('.gbk')[0]) 
            strain_list.append(new_gbk_filename.split('.gbk')[0])
            os.system('mv %s %s%s'%(gbk_filepath,gbk_path,new_gbk_filename))
        else:
            strain_list.append(gbk_filename.split('.gbk')[0])
    write_pickle(path+'strain_list.cpk', strain_list )
