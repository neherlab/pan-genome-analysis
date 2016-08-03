from __future__ import print_function, division
import os,sys,copy;import numpy as np
from treetime import treeanc as ta
from treetime.gtr import GTR
from treetime import io
from Bio import Phylo, AlignIO
from SF00miscellaneous import write_json
import sys


def infer_gene_gain_loss(path, rates = [1.0, 1.0]):
    # initialize GTR model with given parameters
    mu = np.sum(rates)
    gene_pi = np.array(rates)/mu
    gain_loss_model = GTR.custom(pi = gene_pi, mu=mu,
                           W=np.ones((2,2)),
                           alphabet = np.array(['0','1']))
    # add "unknown" state to profile
    gain_loss_model.profile_map['-'] = np.ones(2)
    root_dir = os.path.dirname(os.path.realpath(__file__))

    # define file names for pseudo alignment of presence/absence patterns as in 001001010110
    sep='/'
    fasta = sep.join([path.rstrip(sep), 'geneCluster', 'genePresence.aln'])
    # strain tree based on core gene SNPs
    nwk =  sep.join([path.rstrip(sep), 'geneCluster', 'tree_result.newick'])

    # instantiate treetime with custom GTR
    t = io.treetime_from_newick(gain_loss_model, nwk)
    # fix leaves names since Bio.Phylo interprets numeric leaf names as confidence
    for leaf in t.tree.get_terminals():
        if leaf.name is None:
            leaf.name = str(leaf.confidence)
    # load alignment and associate with tree leafs
    failedleaves = io.set_seqs_to_leaves(t, AlignIO.read(fasta, 'fasta'))
    if failedleaves > 0:
        print(' '.join(["Warning",str(failedleaves),"leaves have failed"]))
    
    print(t.tree.get_terminals()[10])
    print(t.tree.get_terminals()[10].sequence[9375:9390])

    t.reconstruct_anc(method='ml')
    
    print(t.tree.get_terminals()[10])
    print(t.tree.get_terminals()[10].sequence[9375:9390])

    for n in t.tree.find_clades():
        n.genepresence = n.sequence

    return t


def export_gain_loss(tree, path):
    '''
    '''
    # write final tree with internal node names as assigned by treetime
    sep='/'
    tree_fname = sep.join([path.rstrip(sep), 'geneCluster', 'tree_result.newick'])
    Phylo.write(tree.tree, tree_fname, 'newick')

    from collections import defaultdict
    gene_gain_loss_dict=defaultdict(str)
    for node in tree.tree.find_clades(order='preorder'):# order does not matter much here
            if node.up is None: continue
            #print(node.name ,len(node.geneevents),node.geneevents)
            gain_loss = [ str(int(ancestral)*2+int(derived))
                        for ancestral,derived in zip(node.up.genepresence, node.genepresence)]
            gene_gain_loss_dict[node.name]="".join(gain_loss)

    # export gene loss dict to json for visualization
    gene_loss_fname = sep.join([path.rstrip(sep), 'geneCluster', 'genePresence.json'])
    write_json(gene_gain_loss_dict, gene_loss_fname, indent=1) 


def process_gain_loss(path):
    tree = infer_gene_gain_loss(path)
    export_gain_loss(tree, path)
    
    
def create_visible_pattern_dictionary(tree):
    """
    create a sequence in all leaves such that each presence absence pattern occurs only once
    """
    #create a pattern dictionary  
    #patterndict = {pattern_tuple: [first position in pseudoalignment with pattern, number of genes with this pattern]}
    #clusterdict = {first position with pattern: number of genes with pattern}
    #initialize dictionaries
    tree.tree.patterndict = {}
    numstrains = len(tree.tree.get_terminals())
    corepattern = ('1',)*numstrains
    nullpattern = ('0',)*numstrains
    tree.tree.clusterdict = {}
    #create dictionaries
    numgenes = tree.tree.get_terminals()[0].genepresence.shape[0]
    for genenumber in range(numgenes):
        pattern=()
        for leaf in tree.tree.get_terminals():
            pattern = pattern + (leaf.genepresence[genenumber],)
        if pattern == nullpattern:
            print("Warning: There seems to be a nullpattern in the data! Check your presence absence pseudoalignment at pos", genenumber+1)
        if pattern in tree.tree.patterndict:
            tree.tree.patterndict[pattern][1] = tree.tree.patterndict[pattern][1]+1
            tree.tree.clusterdict[tree.tree.patterndict[pattern][0]] = tree.tree.patterndict[pattern][1]
        else:
            tree.tree.patterndict[pattern] = [genenumber,1]
            tree.tree.clusterdict[tree.tree.patterndict[pattern][0]] = tree.tree.patterndict[pattern][1]
               
    #thin sequence to unique pattern
    for node in tree.tree.find_clades():
        if hasattr(node, 'sequence'):
            if len(node.sequence) != numgenes:
                print ("Warning: Nonmatching number of genes in sequence")
            node.patternseq = node.sequence[sorted(tree.tree.clusterdict.keys())]
            # add the all zero pattern at the end of all pattern
            node.patternseq = np.append(node.patternseq,['0',])

    # add an artificial pattern of all zero
    tree.tree.patterndict[nullpattern] = [numgenes,0]
    tree.tree.clusterdict[tree.tree.patterndict[nullpattern][0]] = tree.tree.patterndict[nullpattern][1]
    tree.tree.pattern_abundance = [tree.tree.clusterdict[key] for key in sorted(tree.tree.clusterdict.keys())]
    #save the index of the first core pattern (in most cases this should be zero)
    tree.tree.corepattern_index = sorted(tree.tree.clusterdict.keys()).index(tree.tree.patterndict[corepattern][0])

def index2pattern(index,numstrains):
    pattern = [0] * numstrains
    for ind in index:
        pattern[ind] = 1
    return tuple(pattern)    
    
def create_ignoring_pattern_dictionary(tree,p = 0):
    """
    create a dictionary of pattern that correspond to extended core genes and extended unique genes
    these pattern will be ignored in the inference of gene gain/loss rates
    """
    #create a pattern dictionary  
    #patterndict = {pattern_tuple: [first position in pseudoalignment with pattern, number of genes with this pattern]}
    #initialize dictionaries
    import itertools
    tree.tree.unpatterndict = {}
    numstrains = len(tree.tree.get_terminals())
    if p == 0:
        p = numstrains/10
    corepattern = ('1',)*numstrains
    nullpattern = ('0',)*numstrains
    
    #all sets of indices for p or less of numstrains individuals
    myindices = iter(())
    for i in range(p):
        myindices = itertools.chain(myindices, itertools.combinations(range(numstrains),i+1))
        
    for indices in myindices:
        tree.tree.unpatterndict[index2pattern(indices,numstrains)] = [0]


#def extract_ignoring_pattern_dictionary(tree):
    #i will need sum([int(i) for i in pattern]) here
    

if __name__=='__main__':
    species= 'Papn'
    #path = '/ebio/ag-neher/share/users/wding/mpam/data/'+species+'/'
    path = '/home/franz/tmp/'
    tree = infer_gene_gain_loss(path)

    #outpath = '.'
    outpath = '/home/franz/tmp/'
    export_gain_loss(tree, outpath)