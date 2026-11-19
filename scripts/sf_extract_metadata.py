def extract_metadata(path, strain_list, folders_dict, gbk_present, metainfo_reconcile):
    """
    extract metainfo (date/country) from genBank file
    This step is not necessary if the user provides a tab-delimited
    meta-information table as path/"metainfo_curated.tsv"
    path: species folder
    gbk_present: flag to indicate whether genBank is provided
    Input: genBank file
    Output: metainfo csv file
    """
    def organism_format(raw_data):
        """"""
        #Pseudomonas putida ND6
        #[Pseudomonas syringae] pv. tomato str. DC3000
        new_tag=raw_data
        if ' ' in raw_data: #Binomial nomenclature
            new_tag=raw_data.replace('[','').replace(']','').split(' ')[:2]
            new_tag=new_tag[0][:1]+'.'+new_tag[1]
        return new_tag
    from Bio import SeqIO
    with open('%s%s'%(path,'metainfo.tsv'), 'wb') as writeseq:
        #headers: accession, strainName, dateInfo, country, host
        header=['accession' , 'strain', 'collection_date', 'country', 'host', 'organism']
        if metainfo_reconcile:
            header.insert(4,'region')
        writeseq.write( "%s\n"%('\t'.join(header)) )
        if gbk_present==True:
            if metainfo_reconcile:
                host_synonym_file= path+'../../metadata/host_synonym.tsv'
                country_to_region_file=path+'../../metadata/country2region.tsv'
                host_synonym_dict={}
                country_to_region_dict={}
                with open(host_synonym_file) as host_synonym_items:
                    host_synonym_items.next()
                    for host_synonym in host_synonym_items:
                        host_original, host_unified = host_synonym.rstrip().split('\t')
                        host_synonym_dict[host_original] = host_unified
                with open(country_to_region_file) as country_to_region_input:
                    country_to_region_input.next()
                    for country_to_region in country_to_region_input:
                        country, region= country_to_region.rstrip().split('\t')
                        country_to_region_dict[country]=region
            gbk_path=folders_dict['gbk_path']
            for strainID in strain_list:
                gbk_fpath=''.join([gbk_path,strainID,'.gbk'])
                for record in SeqIO.parse(gbk_fpath, "genbank"):
                    for feature in record.features:
                        host, datacolct, country, region, strainName, organismName ='unknown','unknown','unknown','unknown','unknown','unknown'
                        if feature.type=='source':
                            if 'organism' in feature.qualifiers:
                                raw_data=feature.qualifiers['organism'][0]
                                organismName= organism_format(raw_data)
                            else:
                                organismName= 'unknown'
                            if 'strain' in feature.qualifiers:
                                strainName= feature.qualifiers['strain'][0]
                            if 'host' in feature.qualifiers:
                                host= feature.qualifiers['host'][0]
                                #capitalize host string to harmonize GenBank meta-data
                                if host!='unknown':
                                    host= host.capitalize()
                                if metainfo_reconcile:
                                    if host in host_synonym_dict:
                                        host=host_synonym_dict[host]
                            if 'collection_date' in feature.qualifiers:
                                datacolct= feature.qualifiers['collection_date'][0]
                            if 'country' in feature.qualifiers:
                                country= feature.qualifiers['country'][0]
                                country= country.split(':')[0] #USA: New...
                                if metainfo_reconcile:
                                    if country in country_to_region_dict:
                                        region=country_to_region_dict[country]
                                    else:
                                        print 'country name %s not found in country2region.tsv'%country
                            # date processing
                            if datacolct!='unknown':
                                import re, calendar
                                # ignore time
                                match_time = re.search(r'\dT\d\d:\d\d:\d\dZ', datacolct) # for date-time format
                                if match_time:
                                    start_index_time = datacolct.index(match_time.group()) + 1
                                    datacolct = datacolct[:start_index_time]
                                match_time = re.search(r'\dT\d\dZ', datacolct) # for date-time format but time deleted manually
                                if match_time:
                                    start_index_time = datacolct.index(match_time.group()) + 1
                                    datacolct = datacolct[:start_index_time]
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
                    metadata_row=[strainID, strainName, datacolct, country, host, organismName]
                    if metainfo_reconcile:
                        metadata_row.insert(4,region)
                    writeseq.write( "%s\n"%('\t'.join(metadata_row)) )
                    break
        else: #gbk files are not provided
            strainName = datacolct = country = host= organismName ='unknown'
            for strainID in strain_list:
                writeseq.write( "%s\n"%('\t'.join([strainID, strainName, datacolct, country, host, organismName])) )
