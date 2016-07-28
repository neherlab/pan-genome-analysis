import os,glob,sys,time,shutil; import numpy as np
from itertools import izip; from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from ete2 import Tree
from SF00miscellaneous import times, read_fasta, load_pickle, write_pickle, write_in_fa, write_json
sys.path.append('./scripts/')

def make_dir(dname):
    import os
    if not os.path.isdir(dname):
        try:
            os.makedirs(dname)
        except OSError as e:
            print "Cannot create run_dir",e

def remove_dir(dname):
    import os, shutil
    if os.path.isdir(dname):
        import shutil
        shutil.rmtree(dname)

def tree_to_json(node, extra_attr = []):
    tree_json = {}
    str_attr = ['country','region','clade','strain', 'date', 'muts']
    num_attr = ['xvalue', 'yvalue', 'tvalue', 'num_date']
    #if hasattr(node, 'name'):
    #    tree_json['strain'] = node.name

    for prop in str_attr:
        if hasattr(node, prop):
            tree_json[prop] = node.__getattribute__(prop)
    for prop in num_attr:
        if hasattr(node, prop):
            try:
                tree_json[prop] = round(node.__getattribute__(prop),5)
            except:
                print "cannot round:", node.__getattribute__(prop), "assigned as is"
                tree_json[prop] = node.__getattribute__(prop)

    for prop in extra_attr:
        if len(prop)==2 and callable(prop[1]):
            if hasattr(node, prop[0]):
                tree_json[prop] = prop[1](node.__getattribute__(prop[0]))
        else:
            if hasattr(node, prop):
                tree_json[prop] = node.__getattribute__(prop)

    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            tree_json["children"].append(tree_to_json(ch, extra_attr))
    return tree_json

def make_strains_unique(aln):
    name_translation = {}
    for si,seq in enumerate(aln):
        name_translation[seq.id] = 'tax_'+str(si+1)
        seq.id = name_translation[seq.id]
        seq.name = seq.id
    return name_translation

def restore_strain_name(name_translation, seqs):
    reverse_translation = {v:k for k,v in name_translation.iteritems()}
    for seq in seqs:
        if seq.name not in reverse_translation:
            print("item not found",seq.name)
            continue
        if hasattr(seq, 'id'):
            seq.id = reverse_translation[seq.name]
            #seq.id = seq.id.replace('|','-',1)
        seq.name = reverse_translation[seq.name]

def calc_af(aln, alpha):
    aln_array = np.array(aln)
    af = np.zeros((len(alpha), aln_array.shape[1]))
    for ai, state in enumerate(alpha):
        af[ai] += (aln_array==state).mean(axis=0)
    af[-1] = 1.0 - af[:-1].sum(axis=0)
    return af

nuc_alpha = 'ACGT-N'; aa_alpha = 'ACDEFGHIKLMNPQRSTVWY*-X'
def resolve_polytomies(tree):
    for node in tree.get_nonterminals('preorder'):
        node.confidence = None
        if len(node.clades)>2:
            n = len(node.clades)
            children = list(node.clades)
            node.clades = []
            node.split(branch_length=1e-5)
            if n>3:
                node.clades[0].clades = children[:len(children)//2]
                node.clades[1].clades = children[len(children)//2:]
                for c in node.clades:
                    c.name=''
                    c.confidence = None
            else:
                node.clades[0] = children[0]
                node.clades[1].clades = children[1:]
                node.clades[1].confidence = None
                node.clades[1].name = None

def pad_nucleotide_sequences(aln_aa, seq_nuc):
    '''
    introduce gaps of 3 (---) into nucleotide sequences corresponding to aligned DNA sequences.

    Parameters:
    - aln_aa: amino acid alignment
    - seq_nuc: unaligned nucleotide sequences.

    Returns:
    - aligned nucleotide sequences with all gaps length 3
    '''
    aln_nuc = MultipleSeqAlignment([])
    for aa_seq  in aln_aa:
        try:
            tmp_nuc_seq = str(seq_nuc[aa_seq.id].seq)
        except KeyError as e:
            print aa_seq.id
            print 'Key not found, continue with next sequence'
            continue

        tmpseq = ''
        nuc_pos = 0
        for aa in aa_seq:
            if aa=='-':
                tmpseq+='---'
            else:
                tmpseq+=tmp_nuc_seq[nuc_pos:(nuc_pos+3)]
                nuc_pos+=3

        aln_nuc.append(SeqRecord(seq=Seq(tmpseq),id=aa_seq.id))

    return aln_nuc

def polytomies_midpointRooting(infileName, outfileName):
    # use ete2 to solve polytomies and midpoint rooting
    newickString=open(infileName, 'rb').readline().rstrip()
    tree = Tree(newickString,format=1);
    tree.resolve_polytomy(recursive=True)
    try:
        tree.set_outgroup( tree.get_midpoint_outgroup() );
    except:
        print 'can not conduct midpoint rooting'
    tree.ladderize()
    print 'ladderized'

    ## adding the missing node.name
    #for ind, node in enumerate(tree.traverse("postorder")):
    for ind, node in enumerate(tree.iter_descendants("postorder")):
        if node.name=='': node.name='%s%s'%('NODE_0',ind);

    with open('./'+outfileName, 'wb') as outfile:
        outfile.write(tree.write(format=1))

class mpm_tree(object):
    '''
    class that aligns a set of sequences and infers a tree
    '''

    def __init__(self, alignment_fname, **kwarks):
        self.fname_prefix=alignment_fname.split('geneCluster/')[1].split('.nu')[0]
        self.seqs = {x.id:x for x in SeqIO.parse(alignment_fname, 'fasta')}
        if 'run_dir' not in kwarks:
            import random
            self.run_dir = '_'.join(['temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
        else:
            self.run_dir = kwarks['run_dir']
        self.nuc=True

    def codon_align(self, alignment_tool="mafft", prune=True):
        '''
        takes a nucleotide alignment, translates it, aligns the amino acids, pads the gaps
        note that this suppresses any compensated frameshift mutations

        Parameters:
        - alignment_tool: ['mafft', 'muscle'] the commandline tool to use
        '''
        make_dir(self.run_dir)
        os.chdir(self.run_dir)

        # translate
        aa_seqs = {}
        for seq in self.seqs.values():
            tempseq = seq.seq.translate()
            # use only sequences that translate with out trouble
            if '*' not in str(tempseq)[:-1] or prune==False:
                aa_seqs[seq.id]=SeqRecord(tempseq,id=seq.id)
            else:
                print(seq.id,"has premature stops, discarding")

        tmpfname = 'temp_in.fasta'
        SeqIO.write(aa_seqs.values(), tmpfname,'fasta')

        if alignment_tool=='muscle':
            from Bio.Align.Applications import MuscleCommandline
            cline = MuscleCommandline(input=tmpfname, out=tmpfname[:-5]+'aligned.fasta')
            cline()
            aln_aa = AlignIO.read(tmpfname[:-5]+'aligned.fasta', "fasta")
        elif alignment_tool=='mafft':
            from Bio.Align.Applications import MafftCommandline
            from StringIO import StringIO
            mafft_cline = MafftCommandline(input=tmpfname)
            stdout, stderr = mafft_cline()
            aln_aa = AlignIO.read(StringIO(stdout), "fasta")
        else:
            print 'Alignment tool not supported:'+alignment_tool
            return

        #generate nucleotide alignment
        self.aln = pad_nucleotide_sequences(aln_aa, self.seqs)
        self.sequence_lookup = {seq.id:seq for seq in self.aln}
        os.chdir('..')
        remove_dir(self.run_dir)


    def align(self):
        '''
        align sequencences in self.seqs using mafft
        '''
        make_dir(self.run_dir)
        os.chdir(self.run_dir)

        SeqIO.write(self.seqs.values(), "temp_in.fasta", "fasta")
        os.system("mafft --anysymbol temp_in.fasta > temp_out.fasta")

        self.aln = AlignIO.read('temp_out.fasta', 'fasta')
        self.sequence_lookup = {seq.id:seq for seq in self.aln}
        os.chdir('..')
        remove_dir(self.run_dir)


    def build(self, root='midpoint', raxml=True, raxml_time_limit=0.5):
        '''
        build a phylogenetic tree using fasttree and raxML (optional)
        based on nextflu tree building pipeline
        '''
        from treetime import io
        from treetime import utils
        import subprocess
        make_dir(self.run_dir)
        os.chdir(self.run_dir)
        AlignIO.write(self.aln, 'origin.fasta', 'fasta')
        name_translation = make_strains_unique(self.aln)
        AlignIO.write(self.aln, 'temp.fasta', 'fasta')

        tree_cmd = ["fasttree"]
        if self.nuc: tree_cmd.append("-nt")
        tree_cmd.append("temp.fasta")
        tree_cmd.append(">")
        tree_cmd.append("initial_tree.newick")
        os.system(" ".join(tree_cmd))

        out_fname = "tree_infer.newick"

        ## only for tree with >3 strains
        if raxml and len(set([x.id for x in SeqIO.parse('temp.fasta', 'fasta')]))>3:
            if raxml_time_limit>0:
                tmp_tree = Phylo.read('initial_tree.newick','newick')
                resolve_iter = 0
                resolve_polytomies(tmp_tree)
                while (not tmp_tree.is_bifurcating()) and (resolve_iter<10):
                    resolve_iter+=1
                    resolve_polytomies(tmp_tree)
                Phylo.write(tmp_tree,'initial_tree.newick', 'newick')
                AlignIO.write(self.aln,"temp.phyx", "phylip-relaxed")
                print( "RAxML tree optimization with time limit", raxml_time_limit,  "hours")
                # using exec to be able to kill process
                end_time = time.time() + int(raxml_time_limit*3600)
                process = subprocess.Popen("exec raxml -f d -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
                while (time.time() < end_time):
                    if os.path.isfile('RAxML_result.topology'):
                        break
                    time.sleep(10)
                process.terminate()

                checkpoint_files = glob.glob("RAxML_checkpoint*")
                if os.path.isfile('RAxML_result.topology'):
                    checkpoint_files.append('RAxML_result.topology')
                if len(checkpoint_files) > 0:
                    last_tree_file = checkpoint_files[-1]
                    shutil.copy(last_tree_file, 'raxml_tree.newick')
                else:
                    shutil.copy("initial_tree.newick", 'raxml_tree.newick')
            else:
                shutil.copy("initial_tree.newick", 'raxml_tree.newick')

            try:
                print("RAxML branch length optimization")
                os.system("raxml -f e -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick")
                shutil.copy('RAxML_result.branches', out_fname)
            except:
                print("RAxML branch length optimization failed")
                shutil.copy('raxml_tree.newick', out_fname)
        else:
            #shutil.copy('initial_tree.newick', out_fname)
            polytomies_midpointRooting('initial_tree.newick',out_fname)

        # load the resulting tree as a treetime instance
        from treetime.gtr import GTR
        self.gtr = GTR.standard()
        self.tt = io.treetime_from_newick(self.gtr, out_fname)
        io.set_seqs_to_leaves(self.tt, self.aln)

        # provide short cut to tree and revert names that conflicted with newick format
        self.tree = self.tt.tree
        restore_strain_name(name_translation, self.aln)
        restore_strain_name(name_translation, self.tree.get_terminals())

        os.chdir('..')
        remove_dir(self.run_dir)
        self.is_timetree=False


    def ancestral(self, translate_tree = False):
        '''
        infer ancestral nucleotide sequences using maximum likelihood
        and translate the resulting sequences (+ terminals) to amino acids
        '''
        self.tt.set_additional_tree_params()
        try:
            self.tt.reconstruct_anc(method='ml')
        except:
            print "trouble at self.tt.reconstruct_anc(method='ml')"

        if translate_tree:
            for node in self.tree.find_clades():
                node.aa_sequence = np.fromstring(str(self.translate_seq("".join(node.sequence))), dtype='S1')

    def refine(self):
        '''
        determine mutations on each branch and attach as string to the branches
        '''
        for node in self.tree.find_clades():
            if node.name.startswith('NODE_')==False:
                node.ann=node.name

            if node.up is not None:
                node.muts = ",".join(["".join(map(str, x)) for x in node.mutations if '-' not in x])
                node.aa_muts = ",".join([anc+str(pos+1)+der for pos, (anc, der)
                                    in enumerate(zip(node.up.aa_sequence, node.aa_sequence))
                                    if anc!=der and '-' not in anc and '-' not in der])


    def translate_seq(self, seq):
        '''
        custom translation sequence that handles gaps
        '''
        if type(seq)!=str:
            str_seq = str(seq.seq)
        else:
            str_seq = seq
        try:
            # soon not needed as future biopython version will translate --- into -
            tmp_seq = Seq(str(Seq(str_seq.replace('---', 'NNN')).translate()).replace('X','-'))
        except:
            tmp_seq = Seq(str(Seq(str_seq.replace('-', 'N')).translate()).replace('X','-'))
            print("Trouble translating",seq)
        return tmp_seq


    def translate(self):
        '''
        translate the nucleotide alignment to an amino acid alignment
        '''
        aa_seqs = []
        for seq in self.aln:
            aa_seqs.append(SeqRecord(seq=self.translate_seq(seq), id=seq.id,
                                     name=seq.name, description=seq.description))
        self.aa_aln = MultipleSeqAlignment(aa_seqs)


    def diversity_statistics(self):
        ''' calculate alignment entropy of nucleotide alignments '''
        TINY = 1e-10
        if not hasattr(self, "aln"):
            print("calculate alignment first")
            return
        self.af = calc_af(self.aln, nuc_alpha)
        tmp_af = self.af[:-2]/self.af[:-2].sum(axis=0)
        self.entropy = np.mean(-(tmp_af*np.log(tmp_af+TINY)).sum(axis=0))
        self.diversity = np.mean(1.0-(tmp_af**2).sum(axis=0))


    def mutations_to_branch(self):
        self.mut_to_branch = defaultdict(list)
        for node in self.tree.find_clades():
            if node.up is not None:
                for mut in node.mutations:
                    self.mut_to_branch[mut].append(node)


    def export(self, path = '', extra_attr = ['aa_muts','ann','branch_length','name','longName']):
        ## write tree
        Phylo.write(self.tree, path+self.fname_prefix+'.nwk', 'newick')


        ## processing node name
        for node in self.tree.get_terminals():
            node.name = node.ann.split('|')[0]
            node.longName = node.ann.split(':')[0]

        ## write tree json
        timetree_fname = path+self.fname_prefix+'.tree.json'
        tree_json = tree_to_json(self.tree.root, extra_attr=extra_attr)
        write_json(tree_json, timetree_fname, indent=None)

        ## msa compatible
        for i_aln in self.aln:
            i_aln.id=i_aln.id.replace('|','-',1)
        for i_aa_aln in self.aa_aln:
            i_aa_aln.id=i_aa_aln.id.replace('|','-',1)

        AlignIO.write(self.aa_aln, path+self.fname_prefix+'.aa.aln', 'fasta')
        AlignIO.write(self.aln, path+self.fname_prefix+'.nu.aln', 'fasta')


        ## write seq json
        elems = {}
        for node in self.tree.find_clades():
            if hasattr(node, "sequence"):
                if hasattr(node, "longName")==False:
                    node.longName=node.name
                elems[node.longName] = {}
                nuc_dt= {pos:state for pos, (state, ancstate) in
                                enumerate(izip(node.sequence.tostring(), self.tree.root.sequence.tostring())) if state!=ancstate}
                nodeseq=node.sequence.tostring();nodeseq_len=len(nodeseq)
                elems[node.longName]['nuc']=nuc_dt

        elems['root'] = {}
        elems['root']['nuc'] = self.tree.root.sequence.tostring()

        self.sequences_fname=path+self.fname_prefix+'.seq.json'
        write_json(elems, self.sequences_fname, indent=None)


################################################################################
### functions to run the tree building and alignment routines
################################################################################

def multips(function_in_use, file_path, parallel, fa_files):
    """ running multiple threads """
    from multiprocessing import Process
    fa_files_len_part=len(fa_files)/parallel
    start = time.time(); procs = []

    for i in range(0, parallel):
        if i!=parallel-1:
            p = Process( target=function_in_use, args=(i,file_path, fa_files[ i*fa_files_len_part : (i+1)*fa_files_len_part ], ) )
        else: # the left-overs
            p = Process( target=function_in_use, args=(i,file_path, fa_files[ i*fa_files_len_part : len(fa_files) ], ) )
        p.Daemon = True;p.start();procs.append(p)
    for p in procs: p.join()

    print 'time-SF06*py:', times(start)

def align_and_makeTree(thread, alignFile_path, fa_files_list):
    for gene_cluster_nu_filename in fa_files_list:
        # extract GC_00002 from path/GC_00002.aln
        clusterID = gene_cluster_nu_filename.split('/')[-1].split('.')[0]
        start = time.time();
        geneDiversity_file = open(alignFile_path+'gene_diversity.txt', 'a')
        if len( read_fasta(gene_cluster_nu_filename) )==1: # nothing to do for singletons
            ## nu.aln
            shutil.copy(gene_cluster_nu_filename,
                        gene_cluster_nu_filename.replace('.nu.fa','.nu.aln'))
            ## aa.aln
            shutil.copy(gene_cluster_nu_filename.replace('.nu.fa','.aa.fa'),
                        gene_cluster_nu_filename.replace('.nu.fa','.aa.aln'))

            geneDiversity_file.write('%s\t%s\n'%(clusterID,'0.0'))

        else: # align and build tree
            print gene_cluster_nu_filename
            myTree = mpm_tree(gene_cluster_nu_filename)
            myTree.codon_align()
            myTree.translate()
            myTree.build(raxml=False)
            myTree.ancestral(translate_tree=True)
            myTree.refine()
            myTree.export(path=alignFile_path)
            myTree.diversity_statistics()
            diversity=myTree.diversity
            gene_diversity_values='{0:.3f}'.format(diversity)
            geneDiversity_file.write('%s\t%s\n'%(clusterID,gene_diversity_values))


def cluster_align_makeTree( path, species, parallel ):
    """
    create gene clusters as nucleotide/ amino_acid fasta files
    and build individual gene trees based on nu.fa files
    """
    def create_geneCluster_fa():
        """ dict storing amino_acid Id/Seq from '.fna' files
            input: '.fna', '.fastaCpk', '-orthamcl-allclusters.cpk'
            output:
        """
        fna_path=path+'protein_fna/'
        geneAAcidDt=defaultdict(list)
        for ifna in glob.glob(fna_path+"*.fna"):
            geneAAcidDt.update(read_fasta(ifna))

        ## dict storing nucleotide Id/Seq from '.fastaCpk' files
        istrain_cpk={}; strain_list= load_pickle(path+'strain_list.cpk');
        for istrain in strain_list:
            istrain_cpk[istrain]=load_pickle(path+istrain+'.fastaCpk')

        ## load gene cluster cpk file
        geneCluster_path=fna_path+'diamond_matches/'
        diamond_geneCluster_dt=load_pickle(geneCluster_path+species+'-orthamcl-allclusters.cpk')

        ## create cluster-genes fasta files
        fasta_path=path+'geneCluster/'; os.system('mkdir '+fasta_path);
        ## diamond_geneCluster_dt: {clusterID:[ count_strains,[memb1,...],count_genes }
        for clusterID, gene in diamond_geneCluster_dt.iteritems():
            ## geneCluster file name
            gene_cluster_nu_filename="%s%s"%(clusterID,'.nu.fa')
            gene_cluster_aa_filename="%s%s"%(clusterID,'.aa.fa')
            gene_cluster_nu_write=open( fasta_path+gene_cluster_nu_filename, 'wb')
            gene_cluster_aa_write=open( fasta_path+gene_cluster_aa_filename, 'wb')
            ## write nucleotide/amino_acid sequences into geneCluster files
            for gene_memb in gene[1]:
                ## gene_name format: strain_1|1-20:25
                ##! start_pos in ".fna" has been +1, so -1 for Bio_SeqIO gbk process
                strain_name,full_ann= gene_memb.split('|');
                contig_index=full_ann.split('-')[0];
                st, en= full_ann.split('-')[1].split(':'); st= str(int(st)-1)
                gene_memb_seq=str(istrain_cpk[strain_name][strain_name+'-'+contig_index+'-'+st+'-'+en]);
                write_in_fa(gene_cluster_nu_write,gene_memb, gene_memb_seq )
                #print gene_memb, geneAAcidDt[gene_memb]
                write_in_fa(gene_cluster_aa_write,gene_memb, geneAAcidDt[gene_memb])
            gene_cluster_nu_write.close(); gene_cluster_aa_write.close();

    create_geneCluster_fa()

    ## align, build_tree, make_geneTree_json
    fasta_path = path+'geneCluster/'
    if os.path.exists(fasta_path+'gene_diversity.txt'):
        os.system('rm '+fasta_path+'gene_diversity.txt')
    fa_files=glob.glob(fasta_path+"*.nu.fa")
    multips(align_and_makeTree, fasta_path, parallel, fa_files)

    ## write gene_diversity_Dt cpk file
    def write_gene_diversity_cpk():
        with open(fasta_path+'gene_diversity.txt', 'rb') as infile:
            write_pickle(fasta_path+'gene_diversity.cpk',
                         { i.rstrip().split('\t')[0]:i.rstrip().split('\t')[1] for i in infile})



################################################################################
#### cluster post processing and paralog splitting
################################################################################

def split_cluster(tree, max_branch_length=0.5, max_paralogs=30):
    '''
    linear regression to determine which clusters to split
    return tree/false depending on whether the cluster should be split or not
    '''
    # determine the optimal split
    best_split = find_best_split(tree)
    # evaluate linear discriminator
    return best_split.branch_length/max_branch_length + float(len(best_split.para_nodes))/max_paralogs > 1.0


def find_best_split(tree):
    '''
    iterate over all branches in the tree and find the branch with the largest
    intersection of node names up- and down-stream.
    '''

    #make sets of unique child node names for each node
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            node.leafs = set([node.name.split('|')[0]])
        else:
            node.leafs = set.union(*[c.leafs for c in node])

    # make sets of unique names tips *NOT* in the clade defined by child
    best_split = [None, 0, 0]
    tree.root.not_leafs=set()
    for node in tree.get_nonterminals(order='preorder'):
        for child in node:
            child.not_leafs = set(node.not_leafs).union(*[c.leafs for c in node if c!=child])
            child.para_nodes = set.intersection(child.not_leafs, child.leafs)

            if len(child.para_nodes)>best_split[1]:
                best_split = [child, len(child.para_nodes), child.branch_length]
            elif len(child.para_nodes)==best_split[1] and child.branch_length>best_split[2]:
                best_split = [child, len(child.para_nodes), child.branch_length]
    return best_split[0]


def explore_paralogs(path, nstrains, branch_length_cutoff=500, paralog_cutoff=0.30, plot=False):
    '''
    gather paralog statistics for all trees and plot if desired
    parameters:
    branch_length_cutoff -- cutoff used to determined whether or not to split cluster
                            measured in units of median branch length in the tree (at least 0.01).
                            (defaults to large value 500, i.e. 500 times the median branch length)
    paralog_cutoff       -- cutoff for paralog splitting as fraction of total strains.
                            (default 0.3 -- that is 30%)
    '''
    fname_list = glob.glob(data_path+'geneCluster/*nwk')
    paralog_stat = []; paralog_split_list = [];
    for fi,fname in enumerate(fname_list):
        #if fi%10==0: print('paralog check:', fname.split('/')[-1], fi, 'out of', len(fname_list))
        tree = Phylo.read(fname, 'newick')
        best_split = find_best_split(tree)
        median_branch_length = np.max(0.01, np.median([n.branch_length for n in tree.find_clades()]))
        if best_split is None:
            pass;#paralog_stat.append([fname, 0, 0])
        else:
            paralog_stat.append([fname, best_split.branch_length/median_branch_length, len(best_split.para_nodes)])
            do_split = split_cluster(tree,
                                     max_branch_length = branch_length_cutoff*median_branch_length,
                                     max_paralogs = paralog_cutoff*nstrains)
            if do_split:
                print('will split:', fname, tree.count_terminals(),
                      len(best_split.para_nodes), best_split.branch_length)
                paralog_split_list.append( fname, best_split )

    def plot_paralogs():
        '''
        plot branch length against # of paralogs across trees
        '''
        import matplotlib.pyplot as plt
        plt.ion()
        plt.figure()
        plt.scatter([x[1] for x in paralog_stat], [x[2] for x in paralog_stat])
        plt.ylabel('# paralogs')
        plt.xlabel('branch length')
        plt.savefig(data_path+'explore_paralogs.pdf')

    if plot: plot_paralogs()

    return paralog_split_list


def create_split_cluster_files(parallel, file_path, fname, new_fa_files_list,
                               gene_list1, gene_list2, diamond_geneCluster_dt):
    """
    TODO: add comments, explain function and input
    """
    orgin_nwk_name = fname.split('/')[-1]
    clusterID = orgin_nwk_name.replace('.nwk','')
    origin_cluster_nu_fa = orgin_nwk_name.replace('nwk','nu.fa')
    origin_cluster_aa_fa = orgin_nwk_name.replace('nwk','aa.fa')

    del diamond_geneCluster_dt[clusterID]

    ## write new cluster fa files
    origin_nu_fa_dt = read_fasta(file_path+origin_cluster_nu_fa)
    origin_aa_fa_dt = read_fasta(file_path+origin_cluster_aa_fa)
    sgs_index=0

    for split_gene_list in list(gene_list1), list(gene_list2):
        sgs_index+=1
        newClusterId="%s_%s"%(clusterID,sgs_index)
        gene_cluster_nu_filename="%s%s"%(newClusterId,'.nu.fa')
        gene_cluster_aa_filename="%s%s"%(newClusterId,'.aa.fa')
        gene_cluster_nu_write=open( file_path+gene_cluster_nu_filename, 'wb')
        gene_cluster_aa_write=open( file_path+gene_cluster_aa_filename, 'wb')

        new_fa_files_list.append(file_path+gene_cluster_nu_filename)

        ## write new split cluster files
        for gene_memb in split_gene_list:
            write_in_fa(gene_cluster_nu_write, gene_memb, origin_nu_fa_dt[gene_memb])
            write_in_fa(gene_cluster_aa_write, gene_memb, origin_aa_fa_dt[gene_memb])
        gene_cluster_nu_write.close(); gene_cluster_aa_write.close();

        diamond_geneCluster_dt[ newClusterId ] = [0,[],0]
        ## num_stains
        diamond_geneCluster_dt[ newClusterId ][0]=len(dict(Counter([ ig.split('|')[0] for ig in split_gene_list])).keys())
        ## num_genes
        diamond_geneCluster_dt[ newClusterId ][2]=len(dict(Counter([ ig for ig in split_gene_list])).keys())
        ## gene members
        diamond_geneCluster_dt[ newClusterId ][1]=[ ig for ig in split_gene_list ]


def update_gene_cluster(path, species, diamond_geneCluster_dt ):
    ## update gene cluster pickled file
    cluster_path = path+'protein_fna/diamond_matches/'
    write_pickle(cluster_path+species+'-orthamcl-allclusters_final.cpk',diamond_geneCluster_dt)


def postprocess_paralogs(parallel, path, species, nstrains, branch_length_cutoff=500, paralog_cutoff=0.5, plot=False):
    """
    splitting paralogs, discarding old gene clusters and creating new clusters of split paralogs
    input:
        TODO
    """

    ## exploring paralogs, default: False (not explore and plot), otherwise figure with statistics will saved
    if plot==True:
        explore_paralogs(path, nstrains, branch_length_cutoff=500, paralog_cutoff=0.5, plot=plot)

    file_path = path+'geneCluster/'
    cluster_path= path+'protein_fna/diamond_matches/'
    fname_list = glob.glob(file_path+'*nwk')
    diamond_geneCluster_dt=load_pickle(cluster_path+species+'-orthamcl-allclusters.cpk')
    new_fa_files_list=[]

    for fname in fname_list:
        tree = Phylo.read(fname, 'newick')
        best_split = find_best_split(tree)

        # determine median branch length used as scale for the paralog splitting criterion.
        median_branch_length = np.max(0.01, np.median([n.branch_length for n in tree.find_clades()]))
        if best_split is not None:
            do_split = split_cluster(tree,
                                     max_branch_length = branch_length_cutoff*median_branch_length,
                                     max_paralogs = paralog_cutoff*nstrains)
            if do_split:
                print('will split:', fname, tree.count_terminals(), len(best_split.para_nodes), best_split.branch_length)

                all_genes = set([n.name for n in tree.get_terminals()])
                gene_list1 = set([n.name for n in best_split.get_terminals()])
                gene_list2 = all_genes.difference(gene_list1)
                #print all_genes, gene_list1, gene_list2

                create_split_cluster_files( parallel, file_path, fname, new_fa_files_list, gene_list1, gene_list2, diamond_geneCluster_dt)


    print '---------------', new_fa_files_list
    ## make new aln and tree
    multips(align_and_makeTree, file_path, parallel, new_fa_files_list)

    ## write gene_diversity_Dt cpk file
    with open(file_path+'gene_diversity.txt', 'rb') as infile:
        write_pickle(file_path+'gene_diversity.cpk',{ i.rstrip().split('\t')[0]:i.rstrip().split('\t')[1] for i in infile})

    ## remove old gene cluster and create new split cluster
    update_gene_cluster(path, species, diamond_geneCluster_dt )
