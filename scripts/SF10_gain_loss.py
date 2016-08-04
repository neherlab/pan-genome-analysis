from __future__ import print_function, division
import os,sys,copy;import numpy as np
from treetime import treeanc as ta
from treetime.gtr import GTR
from treetime import io
from Bio import Phylo, AlignIO
from SF00miscellaneous import write_json
import sys
from treetime import seq_utils


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

    t.tree.one_mutation = 1.0
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
    #patterndict = {pattern_tuple: [first position in pseudoalignment with pattern, number of genes with this pattern,indicator to include this pattern in the estimation]}
    #clusterdict = {first position with pattern: [number of genes with pattern,indicator to include gene in gtr inference]}
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
            tree.tree.clusterdict[tree.tree.patterndict[pattern][0]] = [tree.tree.patterndict[pattern][1],1]
        else:
            tree.tree.patterndict[pattern] = [genenumber,1,1]
            tree.tree.clusterdict[tree.tree.patterndict[pattern][0]] = [tree.tree.patterndict[pattern][1],1]
               
    #thin sequence to unique pattern
    for node in tree.tree.find_clades():
        if hasattr(node, 'sequence'):
            if len(node.sequence) != numgenes:
                print ("Warning: Nonmatching number of genes in sequence")
            node.patternseq = node.sequence[sorted(tree.tree.clusterdict.keys())]
            # add the all zero pattern at the end of all pattern
            node.patternseq = np.append(node.patternseq,['0',])

    # add an artificial pattern of all zero
    tree.tree.patterndict[nullpattern] = [numgenes,0,0]
    tree.tree.clusterdict[tree.tree.patterndict[nullpattern][0]] = [tree.tree.patterndict[nullpattern][1],0]
    #create lists for abundance of pattern and inclusion_flag
    tree.tree.pattern_abundance = [tree.tree.clusterdict[key][0] for key in sorted(tree.tree.clusterdict.keys())]
    tree.tree.pattern_include = [tree.tree.clusterdict[key][1] for key in sorted(tree.tree.clusterdict.keys())]
    #save the index of the first core pattern (in most cases this should be zero)
    tree.tree.corepattern_index = sorted(tree.tree.clusterdict.keys()).index(tree.tree.patterndict[corepattern][0])

def index2pattern(index,numstrains):
    pattern = [0] * numstrains
    for ind in index:
        pattern[ind] = 1
    return tuple(pattern)    


def index2pattern_reverse(index,numstrains):
    pattern = [1] * numstrains
    for ind in index:
        pattern[ind] = 0
    return tuple(pattern)  
    
def create_ignoring_pattern_dictionary(tree,p = 0):
    """
    create a dictionary of pattern that correspond to extended core genes and extended unique genes
    these pattern will be ignored in the inference of gene gain/loss rates
    """
    #create a pattern dictionary  
    #unpatterndict = {pattern_tuple: [first position in pseudoalignment with pattern, number of genes with this pattern]}
    #initialize dictionaries
    import itertools
    tree.tree.unpatterndict = {}
    numstrains = len(tree.tree.get_terminals())
    if p == 0:
        p = int(numstrains/10)
    corepattern = ('1',)*numstrains
    nullpattern = ('0',)*numstrains
    
    #all sets of indices for p or less of numstrains individuals
    myindices = iter(())
    for i in range(p):
        myindices = itertools.chain(myindices, itertools.combinations(range(numstrains),i+1))
        
    for indices in myindices:
        tree.tree.unpatterndict[index2pattern(indices,numstrains)] = [-1,0,0]
        tree.tree.unpatterndict[index2pattern_reverse(indices,numstrains)] = [-1,0,0]


def set_visible_pattern_to_ignore(tree,p = 0):
    
    numstrains = len(tree.tree.get_terminals())
    if p == 0:
        p = int(numstrains/10)
    for pattern in tree.tree.patterndict.keys():
        freq = sum([int(i) for i in pattern])
        if freq <= p or freq >= numstrains - p:
            tree.tree.patterndict[pattern][2] = 0
            tree.tree.clusterdict[tree.tree.patterndict[pattern][0]][1] = 0
    
    tree.tree.pattern_include = [tree.tree.clusterdict[key][1] for key in sorted(tree.tree.clusterdict.keys())]


def compute_lh(tree,verbose=0):
    """
    compute the likelihood for each gene presence pattern in the sequence
    """

    L = tree.tree.get_terminals()[0].sequence.shape[0]
    n_states = tree.gtr.alphabet.shape[0]
    if verbose > 2:
        print ("Walking up the tree, computing likelihoods...")
    for leaf in tree.tree.get_terminals():
        # in any case, set the profile
        leaf.profile = seq_utils.seq2prof(leaf.sequence, tree.gtr.profile_map)
        leaf.lh_prefactor = np.zeros(L)
    for node in tree.tree.get_nonterminals(order='postorder'): #leaves -> root
        # regardless of what was before, set the profile to ones
        node.lh_prefactor = np.zeros(L)
        node.profile = np.ones((L, n_states)) # we will multiply it
        for ch in node.clades:
            ch.seq_msg_to_parent = tree.gtr.propagate_profile(ch.profile,
                ch.branch_length,
                rotated=False, # use unrotated
                return_log=False) # raw prob to transfer prob up
            node.profile *= ch.seq_msg_to_parent
            node.lh_prefactor += ch.lh_prefactor
        pre = node.profile.sum(axis=1) #sum over nucleotide states

        node.profile = (node.profile.T/pre).T # normalize so that the sum is 1
        node.lh_prefactor += np.log(pre) # and store log-prefactor
    tree.tree.root.pattern_profile_lh = (np.log(tree.tree.root.profile).transpose() + tree.tree.root.lh_prefactor).transpose()

def change_gtr_parameters_forgainloss(tree,pi_present,mu):
    genepi = np.array([1.0-pi_present,pi_present])
    genepi /= genepi.sum()
    tree.gtr.Pi = np.diagflat(genepi)
    # change speed
    tree.gtr.mu = mu
    # flow matrix
    tree.gtr.W = np.ones((2,2))
    np.fill_diagonal(tree.gtr.W, - ((tree.gtr.W).sum(axis=0) - 1.))
    tree.gtr._check_fix_Q()
    # meanwhile tree.gtr._check_fix_Q() keeps mu        
    tree.gtr._eig()
    
            
def compute_totallh(tree,params,adjustcore = True,verbose = 0):
    """
    compute the total likelihood for all gene presence pattern in the sequence
    be careful: this function changes the gtr model 
    """
    # change the relation of genegain and geneloss rate and the speed
    pi_present = params[0]
    mymu = params[1]
    change_gtr_parameters_forgainloss(tree,pi_present,mymu)

    # this gives the log likelihoods for each pattern
    compute_lh(tree)
    tree.tree.root.pattern_lh =  np.log(np.sum(np.exp(tree.tree.root.pattern_profile_lh)*np.diag(tree.gtr.Pi),axis=1))
    # have to use nansum instead of sum since some pattern have nan likelihoods (maybe due to zero branch_length?)
    tree.tree.root.total_llh =  np.nansum(tree.tree.root.pattern_lh * np.array(tree.tree.pattern_abundance) * np.array(tree.tree.pattern_include))
    #adjust for pattern that should not be included
    ll_forsumofignored = np.nansum(np.exp( tree.tree.root.pattern_lh) * np.subtract(1,tree.tree.pattern_include))
    if verbose > 2:
        print("adjusting nullpattern")
    tree.tree.root.total_llh = tree.tree.root.total_llh - ( np.log(1.- ll_forsumofignored) * np.sum(np.array(tree.tree.pattern_abundance) * np.array(tree.tree.pattern_include)) )
    return tree.tree.root.total_llh * -1.

def set_seq_to_patternseq(tree):
    for node in tree.tree.find_clades():
        if hasattr(node, 'patternseq'):
            node.sequence = node.patternseq
        else:
            delattr(node,sequence)

def set_seq_to_genepresence(tree):
    for node in tree.tree.find_clades():
        if hasattr(node, 'genepresence'):
            node.sequence = node.genepresence
        else:
            delattr(node,sequence)
            
def plot_ll(filename,tree,mu =1.0):
    import matplotlib.pyplot as plt    
    xaxis = np.linspace(0.0001, 0.999, num=500)
    graph = [compute_totallh(tree,[x,mu]) for x in xaxis]
    plt.plot(xaxis, graph)
    plt.savefig(filename)
    
    
def plot_ll_mu(filename,tree,pi_present =0.5):
    import matplotlib.pyplot as plt    
    xaxis = np.linspace(0.01, 250, num=500)
    graph = [compute_totallh(tree,[pi_present,x]) for x in xaxis]
    plt.plot(xaxis, graph)
    plt.savefig(filename)

if __name__=='__main__':
    species= 'Papn'
    #path = '/ebio/ag-neher/share/users/wding/mpam/data/'+species+'/'
    path = '/home/franz/tmp/'
    tree = infer_gene_gain_loss(path)

    #outpath = '.'
    outpath = '/home/franz/tmp/'
    export_gain_loss(tree, outpath)