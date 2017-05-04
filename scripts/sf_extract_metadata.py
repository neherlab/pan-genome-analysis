def extract_metadata(path, strain_list, folders_dict, gbk_present, metainfo_organism):
    """
    extract metainfo (date/country) from genBank file
    This step is not necessary if the user provides a tab-delimited
    meta-information table as path/"metainfo_curated.tsv"
    path: species folder
    gbk_present: flag to indicate whether genBank is provided
    Input: genBank file
    Output: metainfo csv file
    """
    from Bio import SeqIO
    with open('%s%s'%(path,'metainfo.tsv'), 'wb') as writeseq:
        #headers: accession, strainName, dateInfo, country, host
        header=['accession' , 'strain', 'collection_date', 'country', 'host']
        if metainfo_organism:
            header.insert(1, 'organism')
        writeseq.write( "%s\n"%('\t'.join(header)) )
        if gbk_present==1:
            gbk_path=folders_dict['gbk_path']
            for strainID in strain_list:
                gbk_fpath=''.join([gbk_path,strainID,'.gbk'])
                for record in SeqIO.parse(gbk_fpath, "genbank"):
                    for feature in record.features:
                        host, datacolct, country, strainName ='unknown', 'unknown', 'unknown', 'unknown'
                        if feature.type=='source':
                            if metainfo_organism:
                                if 'organism' in feature.qualifiers:
                                    organismName= feature.qualifiers['organism'][0]
                                else:
                                    organismName= 'unknown'
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
                    metadata_row=[strainID, strainName, datacolct, country, host]
                    if metainfo_organism:
                        metadata_row.insert(1, organismName)
                    writeseq.write( "%s\n"%('\t'.join(metadata_row)) )
                    break
        else: #gbk files are not provided
            strainName = datacolct = country = host='unknown'
            for strainID in strain_list:
                writeseq.write( "%s\n"%('\t'.join([strainID, strainName, datacolct, country, host])) )