def accessionID_single(path, strain_lst):
    """ download NCBI refseq GenBank file from strain list """
    import os, sys
    from Bio import Entrez
    from SF00miscellaneous import write_pickle 
            
    Entrez.email = "yourname@mail.com"       
    for gi in strain_lst:
        print gi
        if os.path.exists(path + gi + '.gbk'):
            continue
        handle = Entrez.efetch(db="nucleotide", id=gi,
                               rettype='gbwithparts', retmode="text")
        data = handle.read()
        handle.close()
        write_gi = open(path + gi + '.gbk', 'wb')
        write_gi.write(data)
        write_gi.close()