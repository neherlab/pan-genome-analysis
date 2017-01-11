def gbk_To_Metainfo(path, gbk_present):
    """
    extract metainfo (date/country) from genBank file
    This step is not necessary if the user provides a tab-delimited
    meta-information table as path/"metainfo_curated.tsv"
    path: species folder
    gbk_present: flag to indicate whether genBank is provided
    Input: genBank file
    Output: metainfo csv file
    """
    import os, sys
    from Bio import SeqIO
    from SF00_miscellaneous import load_pickle
    with open(path+'metainfo.tsv', 'wb') as writeseq:
        # write the headers:
        # (default: accName, strainName, dateInfo, country, host)
        writeseq.write( "%s\n"%('\t'.join(['accName' , 'strainName', 'collection_date', 'country', 'host'])) )
        strainList=load_pickle(path+ 'strain_list.cpk')
        if gbk_present==1:
            each_gbk_path='%s%s'%(path,'input_GenBank/')
            # check each genBank file to get meta-type
            for eachstrain in strainList:
                for index, record in enumerate(SeqIO.parse(open(each_gbk_path+ eachstrain+'.gbk'), "genbank")):
                    for i,feature in enumerate(record.features):
                        host, datacolct, country, strainName ='unknown', 'unknown', 'unknown', 'unknown'
                        if feature.type=='source':
                            if 'strain' in feature.qualifiers:
                                strainName= feature.qualifiers['strain'][0]
                            if 'host' in feature.qualifiers:
                                host= feature.qualifiers['host'][0]
                            if 'collection_date' in feature.qualifiers:
                                datacolct= feature.qualifiers['collection_date'][0]
                            if 'country' in feature.qualifiers:
                                country= feature.qualifiers['country'][0]
                                country= country.split(':')[0] #USA: New...

                            # date processing
                            if datacolct!='unknown':
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
                                    elif len(datacolct)==6: #'2010-05'
                                        datacolct='%s-%s-01'%(dates[0][:4], dates[0][4:6])
                                    else: 
                                        datacolct=dates[0]+'-01-01'

                            # just get the year
                            datacolct = datacolct.split('-')[0]
                            break
                    writeseq.write( "%s\n"%('\t'.join([eachstrain, strainName, datacolct, country, host])) )

        else: #gbk files are not provided
            strainName = datacolct = country = host='unknown'
            for eachstrain in strainList:
                writeseq.write( "%s\n"%('\t'.join([eachstrain, strainName, datacolct, country, host])) )
    
    os.system('cp %smetainfo.tsv %smetainfo_curated.tsv'%(path,path))    
    ## meta-info tsv table can be manually modified by user.