import os,glob,sys,time,shutil; import numpy as np
from itertools import izip; from collections import defaultdict, Counter
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq; from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from ete2 import Tree
from SF00_miscellaneous import times, read_fasta, load_pickle, write_pickle, write_in_fa, write_json
sys.path.append('./scripts/')
sys.setrecursionlimit(2000)

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
                tree_json[prop] = round(node.__getattribute__(prop),7)
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
    from ete2 import Tree
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
        self.fname_prefix=alignment_fname.split('geneCluster/')[1].split('.fna')[0]
        temp_prefix=alignment_fname.split('/geneCluster/')[0].split('/')[-1]
        self.seqs = {x.id:x for x in SeqIO.parse(alignment_fname, 'fasta')}
        if 'run_dir' not in kwarks:
            import random
            self.run_dir = '_'.join([temp_prefix, 'temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,100000000))])
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
            os.system("mafft --amino temp_in.fasta 1> temp_out.fasta 2> mafft.log")
            aln_aa = AlignIO.read('temp_out.fasta', "fasta")
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
        os.system("mafft --anysymbol temp_in.fasta 1> temp_out.fasta 2> mafft.log")

        self.aln = AlignIO.read('temp_out.fasta', 'fasta')
        self.sequence_lookup = {seq.id:seq for seq in self.aln}
        os.chdir('..')
        remove_dir(self.run_dir)


    def build(self, root='midpoint', raxml=True, raxml_time_limit=0.5, treetime_used=True):
        '''
        build a phylogenetic tree using fasttree and raxML (optional)
        based on nextflu tree building pipeline
        '''
        import subprocess
        make_dir(self.run_dir)
        os.chdir(self.run_dir)
        AlignIO.write(self.aln, 'origin.fasta', 'fasta')
        name_translation = make_strains_unique(self.aln)
        AlignIO.write(self.aln, 'temp.fasta', 'fasta')

        tree_cmd = ["fasttree"]
        if self.nuc: tree_cmd.append("-nt")
        tree_cmd.append("temp.fasta")
        tree_cmd.append("1>")
        tree_cmd.append("initial_tree.newick")
        tree_cmd.append("2>")
        tree_cmd.append("fasttree.log")
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

        if treetime_used:
            # load the resulting tree as a treetime instance
            from treetime import TreeAnc
            self.tt = TreeAnc(tree=out_fname, aln=self.aln, gtr='Jukes-Cantor', verbose=0)

            # provide short cut to tree and revert names that conflicted with newick format
            self.tree = self.tt.tree
            self.tree.root.branch_length=0.0001
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
        is_valid = self.af[:-2].sum(axis=0)>0.5
        tmp_af = self.af[:-2,is_valid]/self.af[:-2,is_valid].sum(axis=0)
        self.entropy = np.mean(-(tmp_af*np.log(tmp_af+TINY)).sum(axis=0))
        self.diversity = np.mean(1.0-(tmp_af**2).sum(axis=0))


    def mutations_to_branch(self):
        self.mut_to_branch = defaultdict(list)
        for node in self.tree.find_clades():
            if node.up is not None:
                for mut in node.mutations:
                    self.mut_to_branch[mut].append(node)


    def export(self, path = '', extra_attr = ['aa_muts','ann','branch_length','name','longName'], RNA_specific=False):
        ## write tree
        Phylo.write(self.tree, path+self.fname_prefix+'.nwk', 'newick')

        ## processing node name
        for node in self.tree.get_terminals():
            node.name = node.ann.split('|')[0]
            node.longName = node.ann.split('-')[0]

        ## write tree json
        for n in self.tree.root.find_clades():
            if n.branch_length<1e-6:
                n.branch_length = 1e-6
        timetree_fname = path+self.fname_prefix+'_tree.json'
        tree_json = tree_to_json(self.tree.root, extra_attr=extra_attr)
        write_json(tree_json, timetree_fname, indent=None)

        ## msa compatible
        for i_aln in self.aln:
            i_aln.id=i_aln.id.replace('|','-',1)
        
        AlignIO.write(self.aln, path+self.fname_prefix+'_na.aln', 'fasta')

        if RNA_specific==False:
            for i_aa_aln in self.aa_aln:
                i_aa_aln.id=i_aa_aln.id.replace('|','-',1)

            AlignIO.write(self.aa_aln, path+self.fname_prefix+'_aa.aln', 'fasta')

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

        self.sequences_fname=path+self.fname_prefix+'_seq.json'
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
        try:
            # extract GC_00002 from path/GC_00002.aln
            clusterID = gene_cluster_nu_filename.split('/')[-1].split('.')[0]
            start = time.time();
            geneDiversity_file = open(alignFile_path+'gene_diversity.txt', 'a')
            if len( read_fasta(gene_cluster_nu_filename) )==1: # nothing to do for singletons
                ## na.aln
                gene_cluster_nu_aln_filename= gene_cluster_nu_filename.replace('.fna','_na.aln')
                ## geneSeqID separator '|' is replaced by '-' for msa viewer compatibility
                with open(gene_cluster_nu_aln_filename,'wb') as write_file:
                    for SeqID, Sequence in read_fasta(gene_cluster_nu_filename).iteritems():
                        write_in_fa(write_file, SeqID.replace('|','-'), Sequence)

                ## aa.aln
                gene_cluster_aa_filename= gene_cluster_nu_filename.replace('.fna','.faa')
                gene_cluster_aa_aln_filename= gene_cluster_nu_filename.replace('.fna','_aa.aln')
                ## geneSeqID separator '|' is replaced by '-' for msa viewer compatibility
                with open(gene_cluster_aa_aln_filename,'wb') as write_file:
                    for SeqID, Sequence in read_fasta(gene_cluster_aa_filename).iteritems():
                        write_in_fa(write_file, SeqID.replace('|','-'), Sequence)

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
        except:
            print("Aligning and tree building of %s failed"%gene_cluster_nu_filename)
            print(myTree.tree)


def cluster_align_makeTree( path, parallel ):
    """
    create gene clusters as nucleotide/ amino_acid fasta files
    and build individual gene trees based on fna files
    """
    def create_geneCluster_fa():
        """ dict storing amino_acid Id/Seq from '.faa' files
            input: '.faa', '_gene_nuc_dict.cpk', '-orthamcl-allclusters.cpk'
            output:
        """
        ## make sure the geneCluster folder is empty
        if os.path.isdir(path+'geneCluster/')==True:
            print 'remove previous folder: ',path+'geneCluster/'
            os.system('rm -rf %s'%(path+'geneCluster/'))

        faa_path=path+'protein_faa/'
        ## dict storing all genes' translation
        gene_aa_dict=defaultdict(list)
        for ifaa in glob.glob(faa_path+"*.faa"):
            gene_aa_dict.update(read_fasta(ifaa))

        ## dict storing nucleotide Id/Seq from '_gene_nuc_dict.cpk' files
        istrain_cpk={}; strain_list= load_pickle(path+'strain_list.cpk');
        nucleotide_dict_path= '%s%s'%(path,'nucleotide_fna/')
        for istrain in strain_list:
            istrain_cpk[istrain]=load_pickle(nucleotide_dict_path+istrain+'_gene_nuc_dict.cpk')

        ## load gene cluster cpk file
        geneCluster_path=faa_path+'diamond_matches/'
        diamond_geneCluster_dt=load_pickle(geneCluster_path+'orthamcl-allclusters.cpk')

        ## load geneID_to_geneSeqID geneSeqID cpk file
        geneID_to_geneSeqID_dict=load_pickle(path+'geneID_to_geneSeqID.cpk')

        ## create cluster-genes fasta files
        fasta_path=path+'geneCluster/'; os.system('mkdir '+fasta_path)
        ## diamond_geneCluster_dt: {clusterID:[ count_strains,[memb1,...],count_genes }
        for clusterID, gene in diamond_geneCluster_dt.iteritems():
            ## geneCluster file name
            gene_cluster_nu_filename="%s%s"%(clusterID,'.fna')
            gene_cluster_aa_filename="%s%s"%(clusterID,'.faa')
            gene_cluster_nu_write=open( fasta_path+gene_cluster_nu_filename, 'wb')
            gene_cluster_aa_write=open( fasta_path+gene_cluster_aa_filename, 'wb')
            ## write nucleotide/amino_acid sequences into geneCluster files
            for gene_memb in gene[1]:
                ## gene_name format: strain_1|locusTag
                strain_name= gene_memb.split('|')[0]
                gene_memb_seq=str(istrain_cpk[strain_name][gene_memb])
                geneSeqID=geneID_to_geneSeqID_dict[gene_memb]
                write_in_fa(gene_cluster_nu_write, geneSeqID, gene_memb_seq )
                write_in_fa(gene_cluster_aa_write,geneSeqID, gene_aa_dict[gene_memb])
            gene_cluster_nu_write.close(); gene_cluster_aa_write.close();

    create_geneCluster_fa()

    ## align, build_tree, make_geneTree_json
    fasta_path = path+'geneCluster/'
    if os.path.exists(fasta_path+'gene_diversity.txt'):
        os.system('rm '+fasta_path+'gene_diversity.txt')
    fa_files=glob.glob(fasta_path+"*.fna")
    multips(align_and_makeTree, fasta_path, parallel, fa_files)


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


def create_split_cluster_files(file_path, fname,
                               gene_list1, gene_list2, diamond_geneCluster_dt):
    """
    delete the old cluster and create two new clusters
    params:
        new_fa_files: list to which new file names are appeneded
        gene_list1/2: lists containing the genes in the new split clusters
        diamond_geneCluster_dt: cluster dictionary to be updated
    """
    orgin_nwk_name = fname.split('/')[-1]
    clusterID = orgin_nwk_name.replace('.nwk','')
    origin_cluster_nu_fa = orgin_nwk_name.replace('nwk','fna')
    origin_cluster_aa_fa = orgin_nwk_name.replace('nwk','faa')

    split_fa_files_set=set()
    #print 'xxxx', clusterID
    try:
        print('deleting:',orgin_nwk_name,gene_list1,gene_list2, clusterID)
        del diamond_geneCluster_dt[clusterID]
    except:
        print("can't delete",orgin_nwk_name,gene_list1,gene_list2, clusterID)

    ## write new cluster fa files
    origin_nu_fa_dt = read_fasta(file_path+origin_cluster_nu_fa)
    origin_aa_fa_dt = read_fasta(file_path+origin_cluster_aa_fa)
    sgs_index=0

    ## split_gene_list has geneSeqID instead of geneID
    for split_gene_list in (list(gene_list1), list(gene_list2)):
        sgs_index+=1
        newClusterId="%s_%s"%(clusterID,sgs_index)
        gene_cluster_nu_filename="%s%s"%(newClusterId,'.fna')
        gene_cluster_aa_filename="%s%s"%(newClusterId,'.faa')
        gene_cluster_nu_write=open( file_path+gene_cluster_nu_filename, 'wb')
        gene_cluster_aa_write=open( file_path+gene_cluster_aa_filename, 'wb')

        split_fa_files_set |=  set([file_path+gene_cluster_nu_filename])

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
        diamond_geneCluster_dt[ newClusterId ][1]=[ ig.split('-')[0] for ig in split_gene_list ]
    return split_fa_files_set

def update_gene_cluster(path, diamond_geneCluster_dt ):
    ## update gene cluster pickled file
    cluster_path = path+'protein_faa/diamond_matches/'
    write_pickle(cluster_path+'orthamcl-allclusters_final.cpk',diamond_geneCluster_dt)

def update_diversity_cpk_file(path):
    ## write gene_diversity_Dt cpk file
    output_path = path+'geneCluster/'
    with open(output_path+'gene_diversity.txt', 'rb') as infile:
        write_pickle(output_path+'gene_diversity.cpk',{ i.rstrip().split('\t')[0]:i.rstrip().split('\t')[1] for i in infile})

def postprocess_paralogs_iterative(parallel, path, nstrains,
                         branch_length_cutoff=500, paralog_cutoff=0.5, plot=False):

    cluster_path= path+'protein_faa/diamond_matches/'
    diamond_geneCluster_dt=load_pickle(cluster_path+'orthamcl-allclusters.cpk')

    split_result= postprocess_paralogs( parallel, path, nstrains,
                                            diamond_geneCluster_dt,
                                            set(),
                                            branch_length_cutoff=branch_length_cutoff,
                                            paralog_cutoff=paralog_cutoff,
                                            plot=plot)
    n_split_clusters, new_fa_files_set = split_result
    iteration=0
    while(n_split_clusters):
        print('---- split a total of ',n_split_clusters, 'in iteration', iteration)
        split_result= postprocess_paralogs( parallel, path, nstrains,
                                                diamond_geneCluster_dt,
                                                new_fa_files_set,
                                                branch_length_cutoff=branch_length_cutoff,
                                                paralog_cutoff=paralog_cutoff,
                                                plot=plot)
        n_split_clusters, new_fa_files_set = split_result
        iteration+=1

    
    # output_path = path+'geneCluster/'
    # with open(output_path+'gene_diversity.txt', 'rb') as infile:
    #     write_pickle(output_path+'gene_diversity.cpk',{ i.rstrip().split('\t')[0]:i.rstrip().split('\t')[1] for i in infile})
    
    ## write gene_diversity_Dt cpk file
    update_diversity_cpk_file(path)

    ## remove old gene cluster and create new split cluster
    update_gene_cluster(path, diamond_geneCluster_dt )


def postprocess_paralogs(parallel, path, nstrains, diamond_geneCluster_dt,
                        new_fa_files_set, branch_length_cutoff=500, paralog_cutoff=0.5,
                        plot=False):
    """
    splitting paralogs, discarding old gene clusters and creating new clusters of split paralogs
    params:
        parallel: number of threads to use
        path:     path to data
        nstrains: total number of strains
        branch_length_cutoff: multiple of median branch length to split
                              (contribution in linear classifier, parameters to split_cluster)
        paralog_cutoff: fraction of nstrains required for splitting
        plot:      save figure with paralog statistics
    """

    ## exploring paralogs, default: False (not explore and plot), otherwise figure with statistics will saved
    if plot==True:
        explore_paralogs(path, nstrains, branch_length_cutoff=branch_length_cutoff,
                         paralog_cutoff=paralog_cutoff, plot=plot)

    file_path = path+'geneCluster/'

    if len(new_fa_files_set)==0:
        fname_list = glob.glob(file_path+'*nwk')
    else:
        fname_list = [ new_fa.replace('.fna','.nwk') for new_fa in new_fa_files_set ]
        print fname_list

    new_fa_files_set= set()
    n_split_clusters = 0

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

                new_fa_files = create_split_cluster_files(file_path, fname, gene_list1, gene_list2, diamond_geneCluster_dt)
                new_fa_files_set |= new_fa_files
                n_split_clusters+=1

    print 'new_split_fasta_files', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), new_fa_files_set
    ## make new aln and tree
    multips(align_and_makeTree, file_path, parallel, list(new_fa_files_set))

    return n_split_clusters, new_fa_files_set


def load_sorted_clusters(path):
    '''
    load gene clusters and sort 1st by abundance and then by clusterID
    '''
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    diamond_geneCluster_dt=load_pickle(geneClusterPath+'orthamcl-allclusters_final.cpk')
    from operator import itemgetter
    # sort by decreasing abundance (-v[0], minus to achieve decreasing)
    # followed by increasing clusterID GC_00001
    return sorted(diamond_geneCluster_dt.iteritems(),
                   key=lambda (k,v): (-itemgetter(0)(v),k), reverse=False)

#=============================================#
# postprocessing unclustered genes (peaks)     
#=============================================#

# #from SF06_2_unclustered_genes import find_and_merge_unclustered_genes
# from SF06_2_unclustered_genes import find_and_merge_unclustered_genes, cut_tree_from_merged_clusters

# def postprocess_unclustered_genes(n_threads, path, nstrains, window_size=5, strain_proportion=0.3 , sigma_scale=3):
#     diamond_geneCluster_dt=load_pickle(geneClusterPath+'orthamcl-allclusters_final.cpk')
#     find_and_merge_unclustered_genes(n_threads, path, nstrains, window_size, strain_proportion , sigma_scale)
#     cut_tree_from_merged_clusters(path)