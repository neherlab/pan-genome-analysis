def gbk_To_Metainfo(path):
    """
    extract metainfo (date/country) from genBank file
    This step is not necessary if the user provides a tab-delimited
    meta-information table as path/"metainfo_curated.tsv"
    Input: genBank file
    Output: metainfo csv file
    """
    import os, sys
    from Bio import SeqIO
    from SF00miscellaneous import load_pickle
    strainList=load_pickle(path+ 'strain_list.cpk')
    writeseq= open(path+'metainfo.txt', 'wb')
    # write the headers:
    # default: accName, strainName, antiBiotics, dateInfo, country, host
    writeseq.write( "%s\n"%('\t'.join(['accName' , 'collection_date', 'country', 'host'])) )
    # check each genBank file to get meta-type
    for eachstrain in strainList:
        for index, record in enumerate(SeqIO.parse(open(path+ eachstrain+'.gbk'), "genbank")):
            for i,feature in enumerate(record.features):
                if feature.type=='source':
                    host,datacolct,country ='', '', ''
                    if 'host' in feature.qualifiers:
                        host= feature.qualifiers['host'][0]
                    if 'collection_date' in feature.qualifiers:
                        datacolct= feature.qualifiers['collection_date'][0]
                    if 'country' in feature.qualifiers:
                        country= feature.qualifiers['country'][0]
                        country= country.split(':')[0] #USA: New...

                    # date processing
                    import re, calendar
                    datacolct= ''.join(datacolct.split('-'))
                    dates=re.findall('\d+', datacolct);
                    # two versions of date: 15-Seq-2011/2014-03-14
                    if sum([str.isalpha(ic) for ic in datacolct])!=0:
                        month_abbr=re.findall('[a-zA-Z]+', datacolct)[0]
                        month=str(list(calendar.month_abbr).index(month_abbr))
                        if len(datacolct)==9:
                            if len(month)==1: month='0'+month
                            datacolct=dates[1]+'-'+month+'-'+dates[0]
                        else:
                            if len(month)==1: month='0'+month
                            datacolct=dates[0]+'-'+month+'-01'#artificial day 01
                    elif datacolct!='':
                        if  len(datacolct)==8:
                            datacolct='%s-%s-%s'%(dates[0][:4], dates[0][4:6], dates[0][6:])
                        else: datacolct=dates[0]+'-01-01'

                    if datacolct=='': datacolct='unknown'
                    if host=='': host='unknown'
                    if country=='': country='unknown'
                    # just get the year
                    datacolct = datacolct.split('-')[0]
                    # antibiotic default: unknown
                    # antibio='unknown'
                    break
            #writeseq.write( "%s\n"%('\t'.join([eachstrain, antibio, datacolct, country, host])) )
            writeseq.write( "%s\n"%('\t'.join([eachstrain, datacolct, country, host])) )
    writeseq.close()
    # TODO change to .tsv -- also in later scripts
    os.system('cp '+path+'metainfo.txt'+' '+path+'metainfo_curated.txt')
