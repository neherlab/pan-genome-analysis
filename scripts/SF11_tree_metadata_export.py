import os, sys, csv, json
from SF00_miscellaneous import load_pickle
def id_separation(species, path, infile):
    """ extract meta-info from meta-info file 
        input:  metainfo_curated.tsv
        return: id_dict as dictionary storing meta-info 
    """
    from collections import Counter
    id_dict={}; gene_dict={}
    with open(infile) as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t')
        headers = csv_reader.next()
        meta_dict = {}; meta_display_dict = {}
        for h in headers:
            meta_dict[h] = []
            meta_display_dict[h]=h

        for icsv_line in csv_reader:
            id_dict[ icsv_line[0] ] = icsv_line
            for h, v in zip(headers, icsv_line):
                meta_dict[h].append(v)
        for k,v in meta_dict.iteritems():
            meta_dict[k]=dict(Counter(v)).keys()
        del meta_dict[headers[0]]
        del meta_display_dict[headers[0]]

        ## write meta-dict.js
        with open(path+'meta-dict-'+species+'.js','wb') as meta_js_out:
            meta_js_out.write('var meta_set = ')
            meta_js_out.write('%s;'%json.dumps(meta_dict))
            meta_js_out.write('var meta_display_set = ')
            meta_js_out.write('%s;'%json.dumps(meta_display_dict))
    return id_dict, headers

def create_json_simple(node):
    ## tree json without meta-info
    from collections import OrderedDict
    node.name = node.name.replace("'", '')
    json0 = OrderedDict()
    json0["name"]=node.name.split('-')[0]
    if node.children: ##branchset 
        json0["children"] = []
        for ch in node.children: ## recursively 
            json0["children"].append(create_json_simple(ch))
    json0["branch_length"]=float(node.dist)
    return json0

def create_json_addLabel( species, node, flag, path, metaFile):
    ## tree json with meta-info labels
    from collections import OrderedDict
    id_dict,headers=id_separation(species, path, metaFile); 
    node.name = node.name.replace("'", '')
    dt_strainGene=load_pickle(path+'geneCluster/dt_genePresence.cpk')
    json0 = OrderedDict()
    if flag==0:
        json0["name"]=node.name
    elif flag==1:
        if node.name.startswith('NODE_'): 
            json0["name"]=""
        else:
            json0["name"]=node.name

    if node.children: ##branchset 
        json0["children"] = []
        for ch in node.children: ## recursively 
            json0["children"].append(create_json_addLabel(species, ch, flag, path, metaFile))
    if flag==0:
        json0["branch_length"]=float(node.dist)

    if json0["name"].startswith('NODE_') or json0["name"]=="": pass
    else:
        # accName, strainName, antiBio, dateInfo, country, host
        accName_each=json0["name"].split('-')[0]
        for index_head, head in enumerate(headers):
            if index_head!=0: # skip the accName head
                json0[head]=id_dict[accName_each][index_head]
        # load genePresence info
        if flag==0:
            json0["genePresence"]=dt_strainGene[accName_each]
    return json0

def json_tnt_parser():
    read_json=open('coreGenomeTree-noBranch.json', 'rb')
    jsonString='';
    for ir in read_json: 
        jsonString=ir.split('\n')[0]
    jsonString1=jsonString.replace('{"name": "", "children": [', '').replace('name', 'accession').replace(']}', '');
    jsonString_out1='%s%s%s'%('{ "data":[', jsonString1,']}')

    dt_tnt_nodeAttri=dict();
    for i_json in json.loads('[%s]'%jsonString1):
        # jsonString1 [{'acc':1,date:'1981'},{'acc':2,date:'1942'}]
        acc_key=i_json['accession']
        del i_json['accession']
        dt_tnt_nodeAttri[acc_key]=i_json
    jsonString_out2=json.dumps(dt_tnt_nodeAttri)

    with open('strainMetainfo.json', 'wb') as write_json:
        write_json.write(jsonString_out1)

def json_parser( path, species, meta_info_file_path ):
    """ create json file for web-visualiaztion
        input: tree_result.newick, *metainfo_curated.tsv
        output: json files for gene cluster table and core gene SNP tree
    """
    from ete2 import Tree;
    metaFile= path+'metainfo_curated.tsv'
    if meta_info_file_path =='none':
        metaFile= path+'metainfo_curated.tsv'
    else: ## create a link of user meta_info_file
        os.system('pwd')
        os.system('cp %s %s'%(meta_info_file_path, metaFile))

    output_path='%s%s'%(path,'geneCluster/')
    visualzition_path='%s%s'%(path,'Vis/')
    tree = Tree(output_path+'tree_result.newick',format=1)
    ## create tree json files
    jsonString=json.dumps(create_json_addLabel(species, tree, 0, path, metaFile))
    jsonString1=json.dumps(create_json_addLabel(species, tree, 1, path, metaFile))
    os.chdir(output_path)
    with open('coreGenomeTree.json', 'wb') as write_json:
        write_json.write(jsonString)
    with open('coreGenomeTree-noBranch.json', 'wb') as write_json1:
        write_json1.write(jsonString1)

    ## create tnt-nodeAttri-dataTable.json and tnt-nodeAttri.json for tree tables
    json_tnt_parser()

    ## move all *.cpk file to ./data/YourSpecies/ folder
    ##      coreGenomeTree.json and strainMetainfo.json file to ./data/YourSpecies/vis/ folder
    ##      GC*json file to ./data/YourSpecies/vis/geneCluster/ folder
    os.system('mv *.cpk ..; \
               mv coreGenomeTree.json ../vis/;\
               mv strainMetainfo.json ../vis/;\
               mv geneGainLossEvent.json ../vis/;\
               mv GC*.aln ../vis/geneCluster/;\
               mv GC*_tree.json ../vis/geneCluster/;')