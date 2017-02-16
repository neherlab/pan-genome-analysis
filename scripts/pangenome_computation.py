import os, glob
from collections import defaultdict
from sf_miscellaneous import times, write_in_fa, load_pickle, write_pickle

class pangenome:
    """organize and streamline pangenome analysis """
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
        self.path= '%s/'%self.path

    def organize_folder(self):
        """ create folders for pangenome analysis """
        command_mkdir='mkdir -p '
        folders_dict=defaultdict( str,
            gbk_path='input_GenBank/',
            nucleotide_path='nucleotide_fna/',
            protein_path='protein_faa/',
            clustering_path='protein_faa/diamond_matches/',
            RNA_path='RNA_fna/',
            cluster_seq_path='geneCluster/',
            tmp_core_seq_path='tmp_core/')
        for k,v in folders_dict.iteritems():
            folders_dict[k]='%s%s'%(self.path,v)
        for key, folder_path in folders_dict.iteritems():
            os.system(''.join([command_mkdir,folder_path]))
        self.nucleotide_path=folders_dict['nucleotide_path']
        self.protein_path=folders_dict['protein_path']
        self.clustering_path=folders_dict['clustering_path']
        self.folders_dict=folders_dict
        return folders_dict

    def specify_filepath(self):
        """ organize file paths """
        fpaths_dict=defaultdict( str,
            strain_cpk='%s%s'%(self.path,'strain_list.cpk'),
            geneID_to_geneSeqID_path='%s%s'%(self.path,'geneID_to_geneSeqID.cpk'),
            geneID_to_description_path='%s%s'%(self.path,'geneID_to_description.cpk'),
            RNAID_to_SeqID_path='%s%s'%(self.path,'RNAID_to_SeqID.cpk'),
            RNAID_to_description_path='%s%s'%(self.path,'RNAID_to_description.cpk'),
            protein_dict_path='%s%s'%(self.protein_path,'all_protein_seq.cpk'),
            nucleotide_dict_path='%s%s'%(self.nucleotide_path,'all_nucleotide_seq.cpk'),
            cluster_fpath='%s%s'%(self.clustering_path,'allclusters.tsv'),
            cluster_final_fpath='%s%s'%(self.clustering_path,'allclusters_final.tsv'),
            cluster_cpk_fpath='%s%s'%(self.clustering_path,'allclusters.cpk'),
            cluster_cpk_final_fpath='%s%s'%(self.clustering_path,'allclusters_postprocessed.cpk'))
        self.fpaths_dict=fpaths_dict

    def make_strain_list(self):
        """ make strainID list and harmonize input filename """
        ## load input strains in strain_list 
        path=self.path
        folders_dict=self.folders_dict
        if self.gbk_present==1:
            glob_item='.gbk'
            gbk_path=folders_dict['gbk_path']
            glob_list=glob.glob('%s*%s'%(path,glob_item))
            if len(glob_list)!=0:
                harmonize_filename(path,glob_list)
                strain_list= [i.split('/')[-1].split(glob_item)[0] for i in glob.iglob('%s*%s'%(path,glob_item))]
                ## move gbk files in folder input_GenBank
                command_organize_gbk_input=''.join(['mv ',path,'*gbk ',gbk_path])
                os.system(command_organize_gbk_input)
            else:
                gbk_glob=glob.iglob('%s*%s'%(gbk_path,glob_item))
                strain_list= [i.split('/')[-1].split(glob_item)[0] for i in gbk_glob]
        else:
            glob_item='.faa'
            glob_list=glob.glob('%s*%s'%(path,glob_item))
            if len(glob_list)!=0:
                harmonize_filename(path,glob_list)
                strain_list=[i.split('/')[-1].split(glob_item)[0] for i in glob.iglob('%s*%s'%(path,glob_item))]
            else:
                protein_glob=glob.iglob('%s*%s'%(folders_dict['protein_path'],glob_item))
                strain_list= [i.split('/')[-1].split(glob_item)[0] for i in protein_glob]
            command_organize_aa_input= 'mv %s*.faa %s'%(path,folders_dict['protein_path'])
            command_organize_nuc_input='mv %s*.fna %s'%(path,folders_dict['nucleotide_path'])
            os.system(command_organize_nuc_input)
            os.system(command_organize_aa_input)
        write_pickle('%s%s'%(path,'strain_list.cpk'), strain_list)
        self.strain_list=strain_list
        self.strainCollector= defaultdict()

def harmonize_filename(path,glob_list):
    """ force '-' to be replaced as '_' in input filename """
    for fpath in glob_list:
        gbk_fname=fpath.split('/')[-1]
        if '-' in gbk_fname:
            ## force to replace '-' with '_' in GenBank filename
            gbk_fname= gbk_fname.replace('-','_')
            print ''.join(['filename harmonized: ',fpath,' -> ',gbk_fname]) 
            os.system(''.join(['mv ',fpath,' ',path,gbk_fname]))
