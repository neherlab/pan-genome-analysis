def accessionID_single(path, strain_lst, email=None):
    """ download NCBI refseq GenBank file from strain list """
    import os, sys
    from Bio import Entrez
    from SF00miscellaneous import write_pickle

    if email is None:
        Entrez.email = "yourname@mail.com"
    else:
        Entrez.email = email

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