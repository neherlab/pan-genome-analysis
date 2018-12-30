import os, sys, csv, json, glob, shutil
from collections import defaultdict, Counter, OrderedDict
from ete2 import Tree
from sf_miscellaneous import write_json
import math
import pandas as pd
import numpy as np


class Metadata(object):
    """docstring for Metadata"""
    def __init__(self, infile, data_description):
        self.data_description = pd.read_csv(data_description, sep='\t')
        self.data = pd.read_csv(infile, sep='\t')
        self.headers = self.data.columns
        for col in self.headers:
            try:
                self.data[col] = self.data[col].astype('float')
            except:
                pass

        self.strains = list(self.data.loc[:,'accession'])
        self.convert_threshold_data()


    def convert_threshold_data(self):
        for cat, dt in self.data_description.iterrows():
            if dt.loc["data_type"]=='mixed_continuous':
                self.convert_threshold_column(dt.loc['meta_category'])


    def convert_threshold_column(self, col, replace_patterns=('>=', '<=','=>', '=<', '>', '<', '=')):
        def is_numeric(s):
            try:
                _ = float(s)
                return True
            except:
                return False

        if True: #try:
            print(col)
            tmp = self.data[col]
            tmp_converted = []
            for x in tmp:
                if is_numeric(x):
                    tmp_converted.append(float(x))
                elif any([is_numeric(x.replace(t,"")) for t in replace_patterns]):
                    for t in replace_patterns:
                        if is_numeric(x.replace(t,"")):
                            tmp_converted.append(float(x.replace(t,"")))
                            continue
                else:
                    tmp_converted.append(np.nan)
            self.data[col] = tmp_converted
        # except:
        #     print("threshold conversion didn't work")

    def to_dict(self):
        return {v['accession']:v for v in self.data.to_dict('records')}


def metadata_load(path, infile):
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
            if len(icsv_line)==0:
                continue
            #assign unknown to empty meta-item
            icsv_line= ["unknown" if len(i)==0 else i for i in icsv_line ]
            accession=str(icsv_line[0])
            strain_meta_dict[ accession ] = icsv_line

            #store metadata belonging to each metadata_type (metadata_type as key)
            #headers[1:]: not consider accession as metadata
            metatable_strain_dt={}
            for header, meta_element in zip(headers, icsv_line):
                if len(meta_element)==0:
                    meta_element="unknown"
                metatable_strain_dt[header]=meta_element
                if header!='accession' and header!='strain':
                    metajson_dict[header].append(meta_element)
            #append metatable_strain_dt for each strain
            metatable_strains_json_dict["data"].append(metatable_strain_dt)

    return (headers, metajson_dict, metatable_strains_json_dict,
            strain_meta_dict)

def metadata_process(path, infile):

    (headers, metajson_dict, metatable_strains_json_dict,
            strain_meta_dict) = metadata_load(path, infile)
    #filter the redundant metadata for each metadata_type
    for k,v in metajson_dict.items():
        meta_list=list(set(v));meta_list.sort()
        metajson_dict[k]=meta_list
    output_path= ''.join([path,'geneCluster/'])
    with open(output_path+'strainMetainfo.json', 'wb') as strain_metadata_json:
        strain_metadata_json.write(json.dumps(metatable_strains_json_dict))
    return strain_meta_dict, headers, metajson_dict

def core_tree_to_json( node, path, metadata_process_result, strain_list):
    ## tree json with meta-info labels

    strain_meta_dict, headers= metadata_process_result
    node.name = node.name.replace("'", '')
    core_tree_dict = OrderedDict()
    core_tree_dict["name"]=node.name
    core_tree_dict["clade"]=strain_list.index(node.name)

    if node.children: ##branchset
        core_tree_dict["children"] = []
        for child in node.children: ## recursively
            core_tree_dict["children"].append(core_tree_to_json(child, path, metadata_process_result, strain_list))
    core_tree_dict["branch_length"]=float(node.dist)
    if not core_tree_dict["name"].startswith('NODE_'):
        ## accession, strain, antibiotics, collection_date, country, host
        accession=core_tree_dict["name"]
        core_tree_dict['attr']={}
        for index_head, head in enumerate(headers):
            if index_head!=0: # skip the accession head
                core_tree_dict['attr'][head]=strain_meta_dict[accession][index_head]
    return core_tree_dict

def process_mixed_continuous(meta_detail):
    #process mixed_continuous data: remove all the signs
    processed_elems=[]
    for raw_elem in meta_detail:
        new_elem=raw_elem.replace(' ','')
        for replace_str in ('>=','<=','>','<','='):
            if replace_str in new_elem:
                new_elem=new_elem.replace(replace_str,'')
        if new_elem.startswith('.'):
            new_elem='0'+new_elem
        if new_elem.endswith('/'):#<0.06/
            new_elem=new_elem.replace('/','')
        if '/' in new_elem:
            numer,denom=new_elem.split('/')
            new_elem=str(round(float(numer)/float(denom),3))
        if new_elem.replace('.','').isdigit():
            #new_elem=round(math.log(float(new_elem),2),2)
           new_elem = float(new_elem)
           if new_elem <= 0:
                new_elem=round(new_elem,2)
           else:
                new_elem=round(math.log(new_elem,2),2)

        processed_elems.append([raw_elem,new_elem])
    return processed_elems

def process_metajson(path, meta_data_config, metajson_dict):
    """ """
    metajson_exp={"color_options":{ }}
    meta_display_choice_dt={}
    meta_display_order=[]
    association_columns=[]
    if meta_data_config!='':
        with open(meta_data_config,'r') as meta_tidy_file:
            next(meta_tidy_file)
            ## meta_category data_type display and association
            for iline in meta_tidy_file:
                #meta_category  data_type  display  associate   log_scale
                meta_category, data_type, display, associate= iline.rstrip().split('\t')[:4]
                meta_display_order.append(meta_category)
                meta_display_choice_dt[meta_category]=(data_type, display)
                if associate=='yes':
                    association_columns.append(meta_category)
        metajson_exp['meta_display_order']=meta_display_order
    else:
        metajson_exp['meta_display_order']=metajson_dict.keys()

    ## write meta-dict.js in metajson
    for metatype in metajson_dict:
        metajson_exp['color_options'][metatype]={"menuItem":metatype}
        if metatype=='country':
            metajson_exp['color_options'][metatype]["type"]="discrete"
        elif metatype=='collection_date':#TODO create num_date #"menuItem":"date",
            metajson_exp['color_options'][metatype]["type"]="continuous"
        elif metatype=='host':
            metajson_exp['color_options'][metatype]["type"]="discrete"
        else:
            print(meta_display_choice_dt)
            if len(meta_display_choice_dt)!=0:
                if metatype in meta_display_choice_dt and meta_display_choice_dt[metatype][1]=='yes': #display
                    metajson_exp['color_options'][metatype]["type"]= meta_display_choice_dt[metatype][0]
                    if metajson_exp['color_options'][metatype]["type"]=='mixed_continuous':
                        metajson_dict[metatype]=process_mixed_continuous(metajson_dict[metatype])
                else:
                    metajson_exp['color_options'][metatype]["display"]='no'
            else:# infer the type
                valid_elem_set=set(metajson_dict[metatype])-set(['unknown'])
                num_isdigit=0
                sign_used=0
                coloring_type=''

                for elem in valid_elem_set:
                    for replace_str in ('>=','<=','>','<','='):
                        #"10.2": remove the 1st decimal point before isdigit()  (10..2 not)
                        if replace_str in elem:
                            sign_used+=1
                            elem=elem.replace(replace_str,'')
                        if elem.replace(' ','').replace('.', '', 1).isdigit():
                            num_isdigit+=1
                            break
                if num_isdigit==len(valid_elem_set):
                    if sign_used==0:
                        coloring_type="continuous"
                    else:
                        coloring_type="mixed_continuous"
                        metajson_dict[metatype]=process_mixed_continuous(metajson_dict[metatype])
                else:
                    coloring_type="discrete"
                metajson_exp['color_options'][metatype]["type"]=coloring_type
                print(metatype, ': undefined coloring type is now set to %s.'%coloring_type)

    with open(''.join([path,'metaConfiguration.js']),'wb') as meta_js_out:
        if len(metajson_dict['organism'])<=1:
            del metajson_dict['organism']
            if 'organism' in metajson_exp['meta_display_order']:
                metajson_exp['meta_display_order'].remove('organism')
                del metajson_exp['color_options']['organism']
        meta_js_content='var meta_details=%s,\nmeta_display=%s;\n'%(json.dumps(metajson_dict),json.dumps(metajson_exp))
        meta_js_content='%svar association_columns=%s;'%(meta_js_content,json.dumps(association_columns))
        meta_js_out.write(meta_js_content)

def json_parser( path, folders_dict, fpaths_dict, meta_info_file_path,
                 meta_data_config, clean_temporary_files ):
    """ create json files for web-visualiaztion
        input: strain_tree.nwk, metainfo.tsv
        output: json files for core gene SNP tree and strain metadata table
    """
    metaFile= '%s%s'%(path,'metainfo.tsv')
    ## create a link of user-provided meta_info_file
    if meta_info_file_path !='none':
        shutil.copy(meta_info_file_path, metaFile)

    main_data_path= path
    clustering_path= folders_dict['clustering_path']
    output_path= folders_dict['cluster_seq_path']
    vis_json_path= folders_dict['vis_json_path']
    vis_cluster_path=  folders_dict['vis_cluster_path']
    tree = Tree(output_path+'strain_tree.nwk',format=1)
    strain_list=[node.name for node in tree.traverse("preorder")]

    ## create strain tree json file
    strain_meta_dict, headers, metajson_dict = metadata_process(path, metaFile)
    metadata_process_result=(strain_meta_dict, headers)
    coreTree_dict=core_tree_to_json(tree, path, metadata_process_result, strain_list)
    coreTree_jsonString=json.dumps(coreTree_dict)
    with open(output_path+'coreGenomeTree.json', 'wb') as core_tree_json:
        core_tree_json.write(coreTree_jsonString)

    ## process meta json
    process_metajson(path, meta_data_config, metajson_dict)
    shutil.move(path+"metaConfiguration.js",vis_json_path+"metaConfiguration.js")

    ## Data organization
    ## Move: visualization-related files into ./vis/geneCluster/ folder
    #os.system('ln -sf %s/*.cpk %s/../'%(output_path,output_path))

    os.chdir(output_path)
    for f in ['coreGenomeTree.json', 'strainMetainfo.json']:
        try:
            shutil.move(f,vis_json_path+f)
        except:
            print("can't move ",f)

    shutil.copy('strain_tree.nwk', vis_json_path+'/strain_tree.nwk')
    for f in glob.glob('*_tree.json') + glob.glob('*.nwk') + glob.glob('*_aln*.fa') + glob.glob('*patterns.json'):
        shutil.move(f,vis_cluster_path+f)
    try:
        shutil.move(clustering_path+'allclusters_final.tsv', main_data_path+'allclusters_final.tsv')
    except:
        print("can't move allclusters_final.tsv")


    ## gzip aln files
    os.system('find '+vis_cluster_path+' -name \*.fa | xargs gzip -f ')
    if clean_temporary_files:
        # clean up record folders
        os.chdir(main_data_path);

        print('clean up temporary files (temporary core gene and post-processed cluster records, etc.)\n')
        tmp_core= folders_dict['tmp_core_seq_path']
        cluster_seq_path= folders_dict['cluster_seq_path']
        deleted= cluster_seq_path+'deleted_clusters/'
        split_long= cluster_seq_path+'update_long_branch_splits/'
        resolve_peak= cluster_seq_path+'update_uncluster_splits/'
        for dirname in [tmp_core, deleted, split_long, resolve_peak]:
            try:
                shutil.rmtree(dirname)
            except:
                print("{} can't be deleted".format(dirname))

        # clean up files
        for key, fpath in fpaths_dict.items():
            if key in ['cluster_fpath','cluster_final_fpath','cluster_cpk_final_fpath']:
                continue
            try:
                os.remove(fpath)
            except:
                print("{} can't be deleted".format(fpath))

    print('Pan-genome analysis is successfully accomplished, the results can be transferred to the local server for panX data visualization and exploration via link-to-server.py in the main folder.')
