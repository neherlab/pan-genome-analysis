import os, sys, csv, json
from collections import defaultdict, Counter
from itertools import izip
from sf_miscellaneous import load_pickle, write_json

def metadata_process(species, path, infile):
    """ extract meta-info from meta-info file
        input:  metainfo.tsv
        return: metadata_dict as dictionary storing meta-info
    """
    metadata_dict=defaultdict()
    with open(infile) as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t')
        headers = csv_reader.next()
        meta_dict = defaultdict(list)
        meta_display_dict = { h:h for h in headers }
        for icsv_line in csv_reader:
            #capitalize host string to consolidate GenBank meta-data
            host_raw= icsv_line[4]
            icsv_line= ["unknown" if len(i)==0 else i for i in icsv_line ]
            if host_raw!='unknown':
                icsv_line[4]='%s%s'%(host_raw[0].upper(), host_raw[1:])
            metadata_dict[ icsv_line[0] ] = icsv_line
            for h, v in izip(headers, icsv_line):
                if len(v)==0:
                    v= "unknown"
                meta_dict[h].append(v)
        for k,v in meta_dict.iteritems():
            meta_item_list=sorted(dict(Counter(v)).keys())
            meta_dict[k]=meta_item_list

        del meta_dict[headers[0]]
        del meta_display_dict[headers[0]]

        ## write meta-dict.js
        with open(''.join([path,'meta-dict-',species,'.js']),'wb') as meta_js_out:
            meta_js_out.write('var meta_set = ')
            meta_js_out.write('%s;'%json.dumps(meta_dict))
            meta_js_out.write('var meta_display_set = ')
            meta_js_out.write('%s;'%json.dumps(meta_display_dict))
    return metadata_dict, headers

def core_tree_to_json( species, dt_genePattern, node, no_pres_flag, table_flag, path, metadata_process_result, strain_list=None):
    ## tree json with meta-info labels
    from collections import OrderedDict
    metadata_dict, headers= metadata_process_result
    node.name = node.name.replace("'", '')
    json0 = OrderedDict()
    if table_flag==0:
        json0["name"]=node.name
        json0["clade"]=strain_list.index(node.name)
    elif table_flag==1:
        if node.name.startswith('NODE_'):
            json0["name"]=""
        else:
            json0["name"]=node.name

    if node.children: ##branchset
        json0["children"] = []
        for ch in node.children: ## recursively
            json0["children"].append(core_tree_to_json(species, dt_genePattern, ch, no_pres_flag, table_flag, path, metadata_process_result, strain_list))
    if table_flag==0:
        json0["branch_length"]=float(node.dist)

    if json0["name"].startswith('NODE_') or json0["name"]=="":
        if no_pres_flag==0:
            json0["genePattern"]=dt_genePattern[json0["name"]]
    else:
        # accName, strainName, antiBio, dateInfo, country, host
        accName_each=json0["name"].split('-')[0]
        json0['attr']=defaultdict(dict)
        for index_head, head in enumerate(headers):
            if index_head!=0: # skip the accName head
                json0['attr'][head]=metadata_dict[accName_each][index_head]
        # load genePattern info
        if no_pres_flag==0:
            json0["genePattern"]=dt_genePattern[accName_each]
    return json0

def json_tnt_parser():
    read_json=open('coreGenomeTree-noBranch.json', 'rb')
    jsonString=''
    for ir in read_json:
        jsonString=ir.split('\n')[0]
    jsonString1=jsonString.replace('{"name": "", "children": [', '').replace(']}', '');
    jsonString_out1='[%s]'%jsonString1.replace('"name"','"accession"')
    tmp_meta_dict=json.loads(jsonString_out1)
    tmp_meta_list= [item_dict for item_dict in tmp_meta_dict]
    jsonString_out1= json.dumps(tmp_meta_list)
    jsonString_out1='{ "data": %s }'%jsonString_out1
    tmp_meta_dict=json.loads(jsonString_out1)
    with open('strainMetainfo.json', 'wb') as strain_metadata_json:
        strain_metadata_json.write(jsonString_out1)
    os.system('rm coreGenomeTree-noBranch.json')

def json_parser( path, species, meta_info_file_path, large_output ):
    """ create json file for web-visualiaztion
        input: tree_result.newick, *metainfo.tsv
        output: json files for gene cluster table and core gene SNP tree
    """
    from ete2 import Tree
    metaFile= '%s%s'%(path,'metainfo.tsv')
    if meta_info_file_path !='none':
        ## create a link of user-provided meta_info_file
        os.system('cp %s %s'%(meta_info_file_path, metaFile))

    output_path= ''.join([path,'geneCluster/'])
    if large_output==0:
        gene_pattern_dict_path= ''.join([output_path,'dt_genePattern.cpk'])
        dt_genePattern=load_pickle(gene_pattern_dict_path)
    else:
        dt_genePattern=None

    #visualzition_path='%s%s'%(path,'vis/')
    tree = Tree('%s%s'%(output_path,'tree_result.newick'),format=1)
    strain_list=[node.name for node in tree.traverse("preorder")]

    ## create tree json files
    no_presence = 0 if large_output==0 else 1
    metadata_process_result= metadata_process(species, path, metaFile);
    coreTree_jsonString=json.dumps(core_tree_to_json(species, dt_genePattern, tree, no_presence, 0, path, metadata_process_result, strain_list))
    coreTree_json_noBranch_string=json.dumps(core_tree_to_json(species, dt_genePattern, tree, 1, 1, path, metadata_process_result, None))
    os.chdir(output_path)
    with open('coreGenomeTree.json', 'wb') as core_tree_json:
        core_tree_json.write(coreTree_jsonString)
    with open('coreGenomeTree-noBranch.json', 'wb') as write_json1:
        write_json1.write(coreTree_json_noBranch_string)

    ## create tnt-nodeAttri-dataTable.json and tnt-nodeAttri.json for tree tables
    json_tnt_parser()

    ## move all *.cpk file to ./data/YourSpecies/ folder
    ##      coreGenomeTree.json and strainMetainfo.json file to ./data/YourSpecies/vis/ folder
    ##      GC*json file to ./data/YourSpecies/vis/geneCluster/ folder
    #os.system('ln -sf %s/*.cpk %s/../'%(output_path,output_path))
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
        # os.system('rm ')
        # os.system('rm ')

    print('Pan-genome analysis is finished, your data can be transfered to the local server for data visualization and exploration via link-to-server.py in the main folder.')
