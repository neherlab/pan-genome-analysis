import numpy as np
import resource, multiprocessing
import os, sys, glob, time, shutil
from itertools import izip
from collections import defaultdict, Counter
from ete2 import Tree
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from sf_miscellaneous import times, read_fasta, \
    load_pickle, write_pickle, write_in_fa, write_json, multips

sys.setrecursionlimit(50000)
nuc_alpha = 'ACGT-N'
aa_alpha = 'ACDEFGHIKLMNPQRSTVWY*-X'

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
                tmpseq='%s---'%tmpseq
            else:
                tmpseq='%s%s'%(tmpseq,tmp_nuc_seq[nuc_pos:(nuc_pos+3)])
                nuc_pos+=3

        aln_nuc.append(SeqRecord(seq=Seq(tmpseq),id=aa_seq.id))

    return aln_nuc

def polytomies_midpointRooting(infileName, outfileName, clusterID):
    # use ete2 to solve polytomies and midpoint rooting
    from ete2 import Tree
    newickString=open(infileName, 'rb').readline().rstrip()
    tree = Tree(newickString,format=1);
    tree.resolve_polytomy(recursive=True)
    try:
        tree.set_outgroup( tree.get_midpoint_outgroup() )
    except:
        print clusterID, ' can not conduct midpoint rooting'
    tree.ladderize()#print 'ladderized'

    ## adding the missing node.name
    #for ind, node in enumerate(tree.traverse("postorder")):
    for ind, node in enumerate(tree.iter_descendants("postorder")):
        if node.name=='': node.name='%s%s'%('NODE_0',ind);

    with open('./%s'%outfileName, 'wb') as outfile:
        outfile.write(tree.write(format=1))

class mpm_tree(object):
    '''
    class that aligns a set of sequences and infers a tree
    '''

    def __init__(self, cluster_seq_filepath, **kwarks):
        self.clusterID= cluster_seq_filepath.split('/')[-1].split('.fna')[0]
        if 'speciesID' in kwarks:
            folderID=kwarks['speciesID']
        else:
            folderID= cluster_seq_filepath.split('/')[-3]
        self.seqs = {x.id:x for x in SeqIO.parse(cluster_seq_filepath, 'fasta')}
        if 'run_dir' not in kwarks:
            import random
            #self.run_dir = '_'.join(['tmp', self.clusterID])
            self.run_dir = '_'.join([folderID, 'temp', time.strftime('%H%M%S',time.gmtime()), str(random.randint(0,100000000))])
        else:
            self.run_dir = kwarks['run_dir']
        self.nuc=True

    def codon_align(self, alignment_tool="mafft", prune=True, discard_premature_stops=False):
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
            tempseq = seq.seq.translate(table="Bacterial")
            # use only sequences that translate without trouble
            if not discard_premature_stops or '*' not in str(tempseq)[:-1] or prune==False:
                aa_seqs[seq.id]=SeqRecord(tempseq,id=seq.id)
            else:
                print(seq.id,"has premature stops, discarding")

        tmpfname = 'temp_in.fasta'
        SeqIO.write(aa_seqs.values(), tmpfname,'fasta')

        if alignment_tool=='mafft':
            os.system('mafft --reorder --amino --anysymbol temp_in.fasta 1> temp_out.fasta 2> mafft.log')
            aln_aa = AlignIO.read('temp_out.fasta', "fasta")
        elif alignment_tool=='muscle':
            from Bio.Align.Applications import MuscleCommandline
            cline = MuscleCommandline(input=tmpfname, out=tmpfname[:-5]+'aligned.fasta')
            cline()
            aln_aa = AlignIO.read(tmpfname[:-5]+'aligned.fasta', "fasta")
        else:
            print 'Alignment tool not supported:'+alignment_tool
            #return

        #generate nucleotide alignment
        self.aln = pad_nucleotide_sequences(aln_aa, self.seqs)
        os.chdir('..')
        remove_dir(self.run_dir)

    def align(self):
        '''
        align sequencences in self.seqs using mafft
        '''
        make_dir(self.run_dir)
        os.chdir(self.run_dir)

        SeqIO.write(self.seqs.values(), "temp_in.fasta", "fasta")
        os.system('mafft --reorder --anysymbol temp_in.fasta 1> temp_out.fasta 2> mafft.log')

        self.aln = AlignIO.read('temp_out.fasta', 'fasta')
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
        tree_cmd.append("temp.fasta 1> initial_tree.newick 2> fasttree.log")
        os.system(" ".join(tree_cmd))

        out_fname = "tree_infer.newick"

        if raxml==False:
            #shutil.copy('initial_tree.newick', out_fname)
            polytomies_midpointRooting('initial_tree.newick',out_fname, self.clusterID)
        elif len(set([x.id for x in SeqIO.parse('temp.fasta', 'fasta')]))>3:
        ## only for tree with >3 strains
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

        if treetime_used:
            # load the resulting tree as a treetime instance
            from treetime.treetime import TreeAnc
            self.tt = TreeAnc(tree=out_fname, aln=self.aln, gtr='Jukes-Cantor', verbose=0)
            # provide short cut to tree and revert names that conflicted with newick format
            self.tree = self.tt.tree
        else:
            self.tree = Phylo.read(out_fname,'newick')
        self.tree.root.branch_length=0.0001
        restore_strain_name(name_translation, self.aln)
        restore_strain_name(name_translation, self.tree.get_terminals())

        for node in self.tree.find_clades():
            if node.name is not None:
                if node.name.startswith('NODE_')==False:
                    node.ann=node.name
            else:
                node.name='NODE_0'

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
            tmp_seq = Seq(str(Seq(str_seq.replace('---', 'NNN')).translate(table="Bacterial")).replace('X','-'))
        except:
            tmp_seq = Seq(str(Seq(str_seq.replace('-', 'N')).translate(table="Bacterial")).replace('X','-'))
            print("Trouble translating", self.clusterID)#print("Trouble translating",seq)
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

    def mean_std_seqLen(self):
        """ returen mean and standard deviation of sequence lengths """
        seqLen_arr = np.array([ len(seq) for seq in self.seqs.values()])
        return np.mean(seqLen_arr, axis=0), np.std(seqLen_arr, axis=0)

    def paralogy_statistics(self):
        best_split = find_best_split(self.tree)
        return len(best_split.para_nodes), best_split.branch_length

    def diversity_statistics_nuc(self):
        ''' calculate alignment entropy of nucleotide alignments '''
        TINY = 1e-10
        if not hasattr(self, "aln"):
            print("calculate alignment first")
            return
        self.af_nuc = calc_af(self.aln, nuc_alpha)
        is_valid = self.af_nuc[:-2].sum(axis=0)>0.5
        tmp_af = self.af_nuc[:-2,is_valid]/self.af_nuc[:-2,is_valid].sum(axis=0)
        #self.entropy_nuc = np.mean(-(tmp_af*np.log(tmp_af+TINY)).sum(axis=0))
        self.diversity_nuc = np.mean(1.0-(tmp_af**2).sum(axis=0))

    def diversity_statistics_aa(self):
        ''' calculate alignment entropy of nucleotide alignments '''
        TINY = 1e-10
        if not hasattr(self, "aln"):
            print("calculate alignment first")
            return
        self.af_aa = calc_af(self.aa_aln, aa_alpha)
        is_valid = self.af_aa[:-2].sum(axis=0)>0.5
        tmp_af = self.af_aa[:-2,is_valid]/self.af_aa[:-2,is_valid].sum(axis=0)
        #self.entropy_aa = np.mean(-(tmp_af*np.log(tmp_af+TINY)).sum(axis=0))
        self.diversity_aa = np.mean(1.0-(tmp_af**2).sum(axis=0))

    def mutations_to_branch(self):
        self.mut_to_branch = defaultdict(list)
        for node in self.tree.find_clades():
            if node.up is not None:
                for mut in node.mutations:
                    self.mut_to_branch[mut].append(node)

    def reduce_alignments(self):
        for attr, aln, alpha, freq in [["aln_reduced", self.aln, nuc_alpha, self.af_nuc],
                                 ["aa_aln_reduced", self.aa_aln, aa_alpha, calc_af(self.aa_aln, aa_alpha)]]:
            try:
                consensus = np.array(list(alpha))[freq.argmax(axis=0)]
                aln_array = np.array(aln)
                aln_array[aln_array==consensus]='.'
                new_seqs = [SeqRecord(seq=Seq("".join(consensus)), name="consensus", id="consensus")]
                for si, seq in enumerate(aln):
                    new_seqs.append(SeqRecord(seq=Seq("".join(aln_array[si])), name=seq.name,
                                       id=seq.id, description=seq.description))
                self.__setattr__(attr, MultipleSeqAlignment(new_seqs))

            except:
                print("sf_geneCluster_align_MakeTree: aligment reduction failed")


    #def export(self, path = '', extra_attr = ['aa_muts','ann','branch_length','name','longName'], RNA_specific=False):
    def export(self, path = '', extra_attr = ['aa_muts','annotation','branch_length','name','accession'], RNA_specific=False):
        ## write tree
        Phylo.write(self.tree, path+self.clusterID+'.nwk', 'newick')

        ## processing node name
        for node in self.tree.get_terminals():
            #node.name = node.ann.split('|')[0]
            node.accession = node.ann.split('|')[0]
            #node.longName = node.ann.split('-')[0]
            node.name = node.ann.split('-')[0]
            #NZ_CP008870|HV97_RS21955-1-fabG_3-ketoacyl-ACP_reductase
            annotation=node.ann.split('-',2)
            if len(annotation)==3:
                node.annotation= annotation[2]
            else:
                node.annotation= annotation[0]

        ## write tree json
        for n in self.tree.root.find_clades():
            if n.branch_length<1e-6:
                n.branch_length = 1e-6
        timetree_fname = path+self.clusterID+'_tree.json'
        tree_json = tree_to_json(self.tree.root, extra_attr=extra_attr)
        write_json(tree_json, timetree_fname, indent=None)

        ## msa compatible
        self.reduce_alignments()
        for i_aln in self.aln:
            i_aln.id=i_aln.id.replace('|','-',1)

        AlignIO.write(self.aln, path+self.clusterID+'_na_aln.fa', 'fasta')
        AlignIO.write(self.aln_reduced, path+self.clusterID+'_na_aln_reduced.fa', 'fasta')

        if RNA_specific==False:
            for i_aa_aln in self.aa_aln:
                i_aa_aln.id=i_aa_aln.id.replace('|','-',1)

            AlignIO.write(self.aa_aln, path+self.clusterID+'_aa_aln.fa', 'fasta')
            AlignIO.write(self.aa_aln_reduced, path+self.clusterID+'_aa_aln_reduced.fa', 'fasta')

        ## write seq json
        write_seq_json=0
        if write_seq_json:
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

            self.sequences_fname=path+self.clusterID+'_seq.json'
            write_json(elems, self.sequences_fname, indent=None)


################################################################################
### functions to run the tree building and alignment routines
################################################################################

def align_and_makeTree( fna_file_list, alignFile_path, simple_tree):
    for gene_cluster_nu_filename in fna_file_list:
        try:
            # extract GC_00002 from path/GC_00002.aln
            clusterID = gene_cluster_nu_filename.split('/')[-1].split('.')[0]
            start = time.time();
            geneDiversity_file = open(alignFile_path+'gene_diversity.txt', 'a')
            if len( read_fasta(gene_cluster_nu_filename) )==1: # nothing to do for singletons
                ## na_aln.fa
                gene_cluster_nu_aln_filename= gene_cluster_nu_filename.replace('.fna','_na_aln.fa')
                ## geneSeqID separator '|' is replaced by '-' for msa viewer compatibility
                with open(gene_cluster_nu_aln_filename,'wb') as write_file:
                    for SeqID, Sequence in read_fasta(gene_cluster_nu_filename).iteritems():
                        write_in_fa(write_file, SeqID.replace('|','-'), Sequence)

                ## aa_aln.fa
                gene_cluster_aa_filename= gene_cluster_nu_filename.replace('.fna','.faa')
                gene_cluster_aa_aln_filename= gene_cluster_nu_filename.replace('.fna','_aa_aln.fa')
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
                if simple_tree==0:
                    myTree.build(raxml=False,treetime_used=True)
                    myTree.ancestral(translate_tree=True)
                    myTree.refine()
                else:
                    myTree.build(raxml=False,treetime_used=False)
                myTree.diversity_statistics_nuc()
                myTree.export(path=alignFile_path)
                #myTree.diversity_statistics_aa()
                #random_alnID=myTree.seqs.keys()[0].split('-')[0]
                diversity_nuc= round(myTree.diversity_nuc,3)#diversity_aa=round(myTree.diversity_aa,3)
                #bestSplit_paraNodes,bestSplit_branchLen = myTree.paralogy_statistics()
                #mean_seqLen, std_seqLen=  myTree.mean_std_seqLen()
                #mean_seqLen, std_seqLen= [ round(i,3) for i in mean_seqLen, std_seqLen ]
                geneDiversity_file.write('%s\t%s\n'%(clusterID,diversity_nuc))
                if 0:
                    cluster_correl_stats_file = open(alignFile_path+'cluster_correl_stats.txt', 'a')
                    cluster_correl_stats_file.write('%s\n'%'\t'.join([
                     str(i) for i in [clusterID, random_alnID, diversity_nuc, \
                        mean_seqLen, std_seqLen, bestSplit_paraNodes, bestSplit_branchLen ] ]))
        except:
            print("Aligning and tree building of %s failed"%gene_cluster_nu_filename)


def mem_check(flag):
    print(flag, ' memory usage: %.2f GB' % round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024.0**2),1))

def create_geneCluster_fa(path,folders_dict):
    """ dict storing amino_acid Id/Seq from '.faa' files
        input: '.faa', '_gene_nuc_dict.cpk', 'allclusters.cpk'
        output:
    """
    ## make sure the geneCluster folder is empty
    os.system('rm -rf %s'%(path+'geneCluster/'))

    clustering_path= folders_dict['clustering_path']
    geneCluster_dt= load_pickle(clustering_path+'allclusters.cpk')
    protein_path= folders_dict['protein_path']
    nucleotide_path= folders_dict['nucleotide_path']

    geneID_to_geneSeqID_dict=load_pickle(path+'geneID_to_geneSeqID.cpk')
    gene_aa_dict= load_pickle('%s%s'%(protein_path,'all_protein_seq.cpk'))
    gene_na_dict= load_pickle('%s%s'%(nucleotide_path,'all_nucleotide_seq.cpk'))

    ## create cluster-genes fasta files
    cluster_seqs_path=path+'geneCluster/'
    os.system('mkdir '+cluster_seqs_path)

    ## write nuc/aa sequences for each cluster
    for clusterID, gene in geneCluster_dt.iteritems():
        ## geneCluster file name
        gene_cluster_nu_filename="%s%s"%(clusterID,'.fna')
        gene_cluster_aa_filename="%s%s"%(clusterID,'.faa')
        with open( cluster_seqs_path+gene_cluster_nu_filename, 'wb') as gene_cluster_nu_write, \
            open( cluster_seqs_path+gene_cluster_aa_filename, 'wb') as gene_cluster_aa_write:
            ## write nucleotide/amino_acid sequences into geneCluster files
            for gene_memb in gene[1]:
                ## gene_name format: strain_1|locusTag
                strain_name= gene_memb.split('|')[0]
                geneSeqID=geneID_to_geneSeqID_dict[gene_memb]
                write_in_fa(gene_cluster_nu_write, geneSeqID, gene_na_dict[strain_name][gene_memb] )
                write_in_fa(gene_cluster_aa_write, geneSeqID, gene_aa_dict[strain_name][gene_memb])

def cluster_align_makeTree( path, folders_dict, parallel, disable_cluster_postprocessing, simple_tree):
    """
    create gene clusters as nucleotide/ amino_acid fasta files
    and build individual gene trees based on fna files
    """
    proc= multiprocessing.Process(target=create_geneCluster_fa, args=(path, folders_dict))
    proc.start(); proc.join()

    ## align, build_tree, make_geneTree_json
    cluster_seqs_path = path+'geneCluster/'
    if os.path.exists(cluster_seqs_path+'gene_diversity.txt'):
        os.system('rm '+cluster_seqs_path+'gene_diversity.txt')

    if 0:
        with open(cluster_seqs_path+'cluster_correl_stats.txt', 'wb') as cluster_correl_stats_file:
            cluster_correl_stats_file.write('%s\n'%'\t'.join(
                        ['clusterID', 'random_alnID', 'diversity_nuc', \
                        'mean_seqLen', 'std_seqLen', 'bestSplit_paraNodes', 'bestSplit_branchLen'
                        ]))

    fna_file_list=glob.glob(cluster_seqs_path+"*.fna")
    multips(align_and_makeTree, parallel, fna_file_list,
        cluster_seqs_path, simple_tree)

    ## if cluster_postprocessing skipped, rename allclusters.cpk as the final cluster file
    if disable_cluster_postprocessing==1:
        update_diversity_cpk(path)
        clustering_path= '%s%s'%(path,'protein_faa/diamond_matches/')
        os.system('cp %sallclusters.tsv %sallclusters_final.tsv'%(clustering_path,clustering_path))

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
            ## calcuate split branch length
            ## if the parent of the best split is the root, the branch_length is the sum of the two children of the root.
            child.split_bl = child.branch_length+ [c.branch_length for c in node if c!=child][0] if node==tree.root else 0
            more_para_nodes = len(child.para_nodes)>best_split[1]
            longer_branch  = len(child.para_nodes)==best_split[1] and child.split_bl>best_split[2]
            if more_para_nodes or longer_branch:
                best_split = [child, len(child.para_nodes), child.split_bl]

    return best_split[0]

def update_diversity_cpk(path):
    ## write gene_diversity_Dt cpk file
    output_path = path+'geneCluster/'
    with open(output_path+'gene_diversity.txt', 'rb') as infile:
        write_pickle(output_path+'gene_diversity.cpk',{ i.rstrip().split('\t')[0]:i.rstrip().split('\t')[1] for i in infile})

def update_geneCluster_cpk(path, geneCluster_dt):
    ## update gene cluster pickled file
    cluster_path = path+'protein_faa/diamond_matches/'
    write_pickle(cluster_path+'allclusters_postprocessed.cpk',geneCluster_dt)
    #write_final_cluster(path, geneCluster_dt)

def write_final_cluster(path):
    geneCluster_dt=load_sorted_clusters(path)
    outfileName='allclusters_final.tsv'
    with open(path+'protein_faa/diamond_matches/'+outfileName, 'wb') as outfile:
        for clusterID, cluster_stat in geneCluster_dt:
            outfile.write('\t'.join([gene for gene in cluster_stat[1]]))
            outfile.write('\n')

def load_sorted_clusters(path):
    '''
    load gene clusters and sort 1st by abundance and then by clusterID
    '''
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(geneClusterPath+'allclusters_postprocessed.cpk')
    from operator import itemgetter
    # sort by decreasing abundance (-v[0], minus to achieve decreasing)
    # followed by increasing strain count
    return sorted(geneCluster_dt.iteritems(),
                key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)
    #return sorted(geneCluster_dt.iteritems(),
    #            key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)
