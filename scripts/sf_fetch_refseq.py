def fetch_refseq(path, strain_lst, species_to_search='Mycoplasma genitalium'):
    """ download NCBI refseq GenBank file from strain list """
    import os, sys, time, glob, csv
    from Bio import GenBank
    from sf_miscellaneous import write_pickle
    #species_to_search
    ## fetch the newest refseq assembly_summary file
    os.system('wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt > %sassembly_summary.txt'%path)
    with open('assembly_summary.txt','rb') as csvfile:
        outfile='downloadlink.txt'
        with open(path+outfile,'wb') as output:
            csv_reader = csv.reader(csvfile, delimiter='\t')
            headers = csv_reader.next()
            for icsv_line in csv_reader:
                # species name and complete
                if species_to_search in icsv_line[7] and 'Complete' in icsv_line[11]:
                    #os.system('wget -c %s/%s'%(icsv_line[19],'*_genomic.gbff.gz -P ./Refseq/Mt'))
                    output.write('%s/%s\n'%(icsv_line[19],'*_genomic.gbff.gz'))

    gbk_path='%sinput_GenBank/'%path
    command_download='wget -c --input %sdownloadlink.txt -P %s'%(path,gbk_path)
    os.system(command_download)
    command_gunzip='gunzip %s*.gz'%gbk_path
    os.system(command_gunzip)
    for each_gbk_path in glob.iglob('%s*gbff*'%gbk_path):
        with open(each_gbk_path) as gbk_file:
            for record in GenBank.parse(gbk_file):
                print(each_gbk_path,record.accession[0])
                break
            os.system('mv %s %s%s.gbk'%(each_gbk_path, gbk_path, record.accession[0]))

    if 0:
        os.chdir(path)
        species=glob.glob('*txt')[0].split('_list.')[0]
        os.system('rm *txt; rm *sh')
        os.system('gunzip *')
        while len(glob.glob('*.gz'))!=0:
            time.sleep(5)
        # rename gbk file
        for each_gbk_path in glob.iglob('*gbff*'):
            with open(each_gbk_path) as handle:
                print handle
                for record in GenBank.parse(handle):
                    print(each_gbk_path,record.accession[0])
                    break
            os.system('mv %s %s'%(each_gbk_path, record.accession[0]))
        for each_gbk_path in glob.iglob('*'):
            os.system('mv %s %s.gbk'%(each_gbk_path, each_gbk_path))
            #os.system('mv %s %s'%(each_gbk_path, each_gbk_path.split('.')[0]))
        os.system('ls *gbk > %s-RefSeq.txt; sed -i -- "s/.gbk//g" *txt'%species)
        os.system('wc -l *txt ; ls *gbk |wc -l')
        path='../../pan-genome-analysis/'
        os.system('cp %srun-TestSet-template.sh %srun-%s.sh; sed -i -- "s/TestSet/%s/g" %srun-%s.sh'%(path,path,species,species,path,species))
        os.system('mv ../%s/ %sdata/'%(species,path))
