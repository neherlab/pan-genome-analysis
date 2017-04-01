import os, sys, csv, json
from collections import defaultdict, Counter, OrderedDict
from ete2 import Tree
from itertools import izip
from sf_miscellaneous import write_json

def metadata_process(species, path, infile):
    """ extract meta-info from meta-info file
        input:  metainfo.tsv
        return: strain_meta_dict as dictionary storing meta-info for each strain
    """
    strain_meta_dict=defaultdict()
    metajson_dict = defaultdict(list)
    # for metatable
    metatable_strains_json_dict={"data":[]}
    with open(infile) as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t')
        headers = csv_reader.next()

        #store metadata for each accession
        for icsv_line in csv_reader:
            #assign unknown to empty meta-item
            icsv_line= ["unknown" if len(i)==0 else i for i in icsv_line ]
            #capitalize host string to harmonize GenBank meta-data
            host_rawdata= icsv_line[4]
            if host_rawdata!='unknown':
                icsv_line[4]= host_rawdata[0].upper() + host_rawdata[1:]
            accession=icsv_line[0]
            strain_meta_dict[ accession ] = icsv_line

            #store metadata belonging to each metadata_type (metadata_type as key)
            #headers[1:]: not consider accession as metadata
            metatable_strain_dt={}
            for header, meta_element in izip(headers, icsv_line):
                if len(meta_element)==0:
                    meta_element="unknown"
                metatable_strain_dt[header]=meta_element
                if header!='accession' and header!='strain':
                    metajson_dict[header].append(meta_element)
            #append metatable_strain_dt for each strain
            metatable_strains_json_dict["data"].append(metatable_strain_dt)

        #filter the redundant metadata for each metadata_type
        for k,v in metajson_dict.iteritems():
            metajson_dict[k]=sorted(dict(Counter(v)).keys())
        output_path= ''.join([path,'geneCluster/'])
        with open(output_path+'strainMetainfo.json', 'wb') as strain_metadata_json:
            strain_metadata_json.write(json.dumps(metatable_strains_json_dict))
    return strain_meta_dict, headers, metajson_dict

def core_tree_to_json( species, node, path, metadata_process_result, strain_list):
    ## tree json with meta-info labels
    strain_meta_dict, headers= metadata_process_result
    node.name = node.name.replace("'", '')
    core_tree_dict = OrderedDict()
    core_tree_dict["name"]=node.name
    core_tree_dict["clade"]=strain_list.index(node.name)

    if node.children: ##branchset
        core_tree_dict["children"] = []
        for child in node.children: ## recursively
            core_tree_dict["children"].append(core_tree_to_json(species, child, path, metadata_process_result, strain_list))
    else:
        ## accession, strainName, antiBio, dateInfo, country, host
        accession=core_tree_dict["name"]
        core_tree_dict['attr']={}
        for index_head, head in enumerate(headers):
            if index_head!=0: # skip the accName head
                core_tree_dict['attr'][head]=strain_meta_dict[accession][index_head]
    return core_tree_dict

def json_parser( path, species, meta_info_file_path, large_output ):
    """ create json files for web-visualiaztion
        input: tree_result.newick, metainfo.tsv
        output: json files for core gene SNP tree and strain metadata table
    """
    metaFile= '%s%s'%(path,'metainfo.tsv')
    if meta_info_file_path !='none':
        ## create a link of user-provided meta_info_file
        os.system('cp %s %s'%(meta_info_file_path, metaFile))

    output_path= ''.join([path,'geneCluster/'])

    tree = Tree('%s%s'%(output_path,'tree_result.newick'),format=1)
    strain_list=[node.name for node in tree.traverse("preorder")]

    ## create strain tree json file
    strain_meta_dict, headers, metajson_dict = metadata_process(species, path, metaFile)
    meta_display_dict = { h:h for h in headers }
    metadata_process_result=(strain_meta_dict, headers)

    coreTree_jsonString=json.dumps(core_tree_to_json(species, tree, path, metadata_process_result, strain_list))

    with open(output_path+'coreGenomeTree.json', 'wb') as core_tree_json:
        core_tree_json.write(coreTree_jsonString)

    ## write meta-dict.js in metajson
    for metatype in metajson_dict:
        if metatype=='country':

        elif metatype=='collection_date':

        elif metatype=='host':

        else:

    with open(''.join([path,'meta-dict-',species,'.js']),'wb') as meta_js_out:
        meta_js_out.write('var meta_set = ')
        meta_js_out.write('%s;'%json.dumps(metajson_dict))
        meta_js_out.write('var meta_display_set = ')
        meta_js_out.write('%s;'%json.dumps(meta_display_dict))

    ## move all *.cpk file to ./data/YourSpecies/ folder
    ##      coreGenomeTree.json and strainMetainfo.json file to ./data/YourSpecies/vis/ folder
    ##      GC*json file to ./data/YourSpecies/vis/geneCluster/ folder
    #os.system('ln -sf %s/*.cpk %s/../'%(output_path,output_path))
    os.chdir(output_path)
    #visualzition_path='%s%s'%(path,'vis/')
    os.system('mv coreGenomeTree.json strainMetainfo.json ../vis/')
    os.system('mv *C*.aln *C*_tree.json ../vis/geneCluster/')
    os.system('mv *C*.nwk ../vis/geneCluster/')
    os.system('cp tree_result.newick ../vis/strain_tree.nwk')
    if large_output==1:
        os.system('mv *C*patterns.json ../vis/geneCluster/')

    keep_temporary_file=0
    if keep_temporary_file:
        strain_protein_fa='./protein_faa/'
        strain_nucleotide_fa='./nucleotide_fna/'

        os.system('rm %s'%strain_protein_fa)
        os.system('rm %s'%strain_nucleotide_fa)
        # os.system('rm ')

    print('Pan-genome analysis is finished, your data can be transfered to the local server for data visualization and exploration via link-to-server.py in the main folder.')
