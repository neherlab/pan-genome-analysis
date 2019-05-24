from __future__ import print_function, division
import os,sys,copy;import numpy as np
from collections import defaultdict
sys.path.append('./')
from treetime import TreeAnc
from treetime import GTR
from treetime import seq_utils
from Bio import Phylo, AlignIO
from sf_miscellaneous import write_json, write_pickle
from sf_geneCluster_align_makeTree import load_sorted_clusters

def infer_gene_gain_loss(path, rates = [1.0, 1.0]):
    # initialize GTR model with default parameters
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
    nwk =  sep.join([path.rstrip(sep), 'geneCluster', 'strain_tree.nwk'])

    # instantiate treetime with custom GTR
    t = TreeAnc(nwk, gtr =gain_loss_model, verbose=2)
    # fix leaves names since Bio.Phylo interprets numeric leaf names as confidence
    for leaf in t.tree.get_terminals():
        if leaf.name is None:
            leaf.name = str(leaf.confidence)
    t.aln = fasta
    t.tree.root.branch_length=0.0001
    t.reconstruct_anc(method='ml')

    for n in t.tree.find_clades():
        n.genepresence = n.sequence

    return t


def export_gain_loss(tree, path, merged_gain_loss_output):
    '''
    '''
    # write final tree with internal node names as assigned by treetime
    sep='/'
    output_path= sep.join([path.rstrip(sep), 'geneCluster/'])
    events_dict_path= sep.join([ output_path, 'dt_geneEvents.cpk'])
    gene_pattern_dict_path= sep.join([ output_path, 'dt_genePattern.cpk'])

    tree_fname = sep.join([output_path, 'strain_tree.nwk'])
    Phylo.write(tree.tree, tree_fname, 'newick')


    gene_gain_loss_dict=defaultdict(str)
    preorder_strain_list= [] #store the preorder nodes as strain list
    for node in tree.tree.find_clades(order='preorder'):# order does not matter much here
            if node.up is None: continue
            #print(node.name ,len(node.geneevents),node.geneevents)
            gain_loss = [ str(int(ancestral)*2+int(derived))
                        for ancestral,derived in zip(node.up.genepresence, node.genepresence)]
            gene_gain_loss_dict[node.name]="".join(gain_loss)
            preorder_strain_list.append(node.name)

    gain_loss_array = np.array([[i for i in gain_loss_str]
                                for gain_loss_str in gene_gain_loss_dict.values()], dtype=int)
    # 1 and 2 are codes for gain/loss events
    events_array = ((gain_loss_array == 1) | (gain_loss_array == 2)).sum(axis=0)
    events_dict =  { index:event for index, event in enumerate(events_array) }

    write_pickle(events_dict_path, events_dict)

    if merged_gain_loss_output:
        ## export gene loss dict to json for visualization
        #gene_loss_fname = sep.join([ output_path, 'geneGainLossEvent.json'])
        #write_json(gene_gain_loss_dict, gene_loss_fname, indent=1)
        write_pickle(gene_pattern_dict_path, gene_gain_loss_dict)
    else:
        ## strainID as key, presence pattern as value (converted into np.array)
        sorted_genelist = load_sorted_clusters(path)
        strainID_keymap= {ind:k for ind, k in enumerate(preorder_strain_list)}
        #presence_arr= np.array([ np.fromstring(gene_gain_loss_dict[k], np.int8)-48 for k in preorder_strain_list])
        presence_arr= np.array([ np.array(gene_gain_loss_dict[k],'c') for k in preorder_strain_list])
        ## if true, write pattern dict instead of pattern string in a json file
        pattern_json_flag=False
        for ind, (clusterID, gene) in enumerate(sorted_genelist):
            pattern_fname='%s%s_patterns.json'%(output_path,clusterID)
            if pattern_json_flag:
                pattern_dt= { strainID_keymap[strain_ind]:str(patt) for strain_ind, patt in enumerate(presence_arr[:, ind])}
                write_json(pattern_dt, pattern_fname, indent=1)
            #print(preorder_strain_list,clusterID)
            #print(''.join([ str(patt) for patt in presence_arr[:, ind]]))
            with open(pattern_fname,'w') as write_pattern:
                write_pattern.write('{"patterns":"'+''.join([ str(patt) for patt in presence_arr[:, ind]])+'"}')


def process_gain_loss(path, merged_gain_loss_output):
    ##  infer gain/loss event
    tree = infer_gene_gain_loss(path)
    create_visible_pattern_dictionary(tree)
    set_seq_to_patternseq(tree)
    set_visible_pattern_to_ignore(tree,p=-1,mergeequalstrains=True)

    def myminimizer(c):
        return compute_totallh(tree,c)

    from scipy.optimize import minimize
    with np.errstate(divide='ignore'):
        try:
            res1 = minimize(myminimizer,[0.5,1.],method='L-BFGS-B',bounds = [(0.01,0.99),(0.1,100.)])
            success1 = res1.success
        except:
            res1 = type('', (), {})()
            res1.fun = np.inf
            success1 = False
        try:
            res2 = minimize(myminimizer,[0.2,1.],method='L-BFGS-B',bounds = [(0.01,0.99),(0.1,100.)])
            success2 = res2.success
        except:
            res2 = type('', (), {})()
            res2.fun = np.inf
            success2 = False
        try:
            res3 = minimize(myminimizer,[0.8,1.],method='L-BFGS-B',bounds = [(0.01,0.99),(0.1,100.)])
            success3 = res3.success
        except:
            res3 = type('', (), {})()
            res3.fun = np.inf
            success3 = False
    if (success1 or success2 or success3) == True:
        print('successfully estimated the gtr parameters. Reconstructing ancestral states...')
        #get the best of the three numerical estimates
        minimalindex = (res1.fun,res2.fun,res3.fun).index(min(res1.fun,res2.fun,res3.fun))
        res = (res1,res2,res3)[minimalindex]
        change_gtr_parameters_forgainloss(tree,res.x[0],res.x[1])
        print('estimated gain/loss rate: ',res.x[0],res.x[1])
        tree.reconstruct_anc(method='ml')
        export_gain_loss(tree,path,merged_gain_loss_output)
    else:
        print('Warning: failed to estimated the gtr parameters by ML.')
        #import ipdb;ipdb.set_trace()
        change_gtr_parameters_forgainloss(tree,0.5,1.0)
        tree.reconstruct_anc(method='ml')
        export_gain_loss(tree,path,merged_gain_loss_output)


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

    #thin sequence to unique pattern and save result to node.patternseq
    for node in tree.tree.find_clades():
        if hasattr(node, 'sequence'):
            if len(node.sequence) != numgenes:
                print ("Warning: Nonmatching number of genes in sequence")
            node.patternseq = node.sequence[sorted(tree.tree.clusterdict.keys())]
            # add the all zero pattern at the end of all pattern
            node.patternseq = np.append(node.patternseq,['0',])

    # add an artificial pattern of all zero (nullpattern)
    tree.tree.patterndict[nullpattern] = [numgenes,0,0]
    tree.tree.clusterdict[tree.tree.patterndict[nullpattern][0]] = [tree.tree.patterndict[nullpattern][1],0]
    #create lists for abundance of pattern and inclusion_flag, resp..
    tree.tree.pattern_abundance = [tree.tree.clusterdict[key][0] for key in sorted(tree.tree.clusterdict.keys())]
    tree.tree.pattern_include = [tree.tree.clusterdict[key][1] for key in sorted(tree.tree.clusterdict.keys())]
    #save the index of the first core pattern
    # check whether there is a corepattern (there should always be a corepattern, unless you are using single cell sequencing data.)
    if corepattern in tree.tree.patterndict:
        tree.tree.corepattern_index = sorted(tree.tree.clusterdict.keys()).index(tree.tree.patterndict[corepattern][0])

def index2pattern(index,numstrains):
    """
    transforms a set of indices to the pattern where only these are 1, all others are 0
    """
    pattern = [0] * numstrains
    for ind in index:
        pattern[ind] = 1
    return tuple(pattern)


def index2pattern_reverse(index,numstrains):
    """
    transforms a set of indices to the pattern where only these are 0, all others are 1
    """
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


def create_distance_matrix(tree):
    numstrains = len(tree.tree.get_terminals())
    tree.tree.distance_matrix = np.zeros([numstrains,numstrains])
    i = 0
    for leaf1 in tree.tree.get_terminals():
        j = 0
        for leaf2 in tree.tree.get_terminals():
            tree.tree.distance_matrix[i,j] = tree.tree.distance(leaf1,leaf2)
            j += 1
        i += 1

def merge_strains(distances,indices,mindist = 0.0):
    remain = set(indices)
    final = set()
    tempdelset = set()
    while len(remain) >0:
        i = remain.pop()
        final.add(i)
        for j in remain:
            if distances[i,j] <= mindist:
                tempdelset.add(j)
        remain.difference_update(tempdelset)
        tempdelset.clear
    return len(final)

def set_visible_pattern_to_ignore(tree,p = -1,mergeequalstrains = False,lowfreq = True, highfreq = True):
    """
    sets all pattern with at most p strains or at least numstrains-p strains to ignore
    """
    numstrains = len(tree.tree.get_terminals())
    if mergeequalstrains:
        if not hasattr(tree.tree,'distance_matrix'):
            create_distance_matrix(tree)
        numstrains = merge_strains(tree.tree.distance_matrix,np.array(range(numstrains)) )
    if p == -1:
        p = int(numstrains/10)
    for pattern in tree.tree.patterndict.keys():
        #freq = sum([int(i) for i in pattern])
        freq = pattern.count('1')
        if lowfreq:
            if mergeequalstrains:
                # indices of individuals in the pattern
                strainindices = np.where(np.array(pattern) == '1')[0]
                #no_observedstrains = len(strainindices)
                lowfreq = merge_strains(tree.tree.distance_matrix,strainindices)
            else:
                lowfreq = freq
            if lowfreq <= p:
                tree.tree.patterndict[pattern][2] = 0
                tree.tree.clusterdict[tree.tree.patterndict[pattern][0]][1] = 0
        if highfreq:
            if mergeequalstrains:
                # indices of individuals in the pattern
                strainindices = np.where(np.array(pattern) == '0')[0]
                #no_observedstrains = len(strainindices)
                highfreq = merge_strains(tree.tree.distance_matrix,strainindices)
            else:
                highfreq = numstrains -freq
            if highfreq <= p:
                tree.tree.patterndict[pattern][2] = 0
                tree.tree.clusterdict[tree.tree.patterndict[pattern][0]][1] = 0

    tree.tree.pattern_include = [tree.tree.clusterdict[key][1] for key in sorted(tree.tree.clusterdict.keys())]
    if sum(tree.tree.pattern_include) == 0:
        print('WARNING all pattern have been excluded, estimation of parameters is thus impossible')


def _check_seq_and_patternseq(tree):
    for leaf in tree.tree.get_terminals():
        if all(leaf.sequence != leaf.patternseq):
            print('WARNING: wrong pattern in ',leaf.name)
        else:
            print(leaf.name, ' is ok')

def compute_lh(tree,verbose=0):
    """
    compute the likelihood for each gene presence pattern in the sequence given the gtr parameters
    """

    min_branch_length = 1e-10
    L = tree.tree.get_terminals()[0].sequence.shape[0]
    n_states = tree.gtr.alphabet.shape[0]
    if verbose > 2:
        print ("Walking up the tree, computing likelihoods for the pattern in the leaves...")
    for leaf in tree.tree.get_terminals():
        # in any case, set the profile
        leaf.profile = seq_utils.seq2prof(leaf.sequence, tree.gtr.profile_map)
        leaf.lh_prefactor = np.zeros(L)
    for node in tree.tree.get_nonterminals(order='postorder'): #leaves -> root
        # regardless of what was before, set the profile to ones
        node.lh_prefactor = np.zeros(L)
        node.profile = np.ones((L, n_states)) # this has to be ones in each entry -> we will multiply it
        for ch in node.clades:
            ch.seq_msg_to_parent = tree.gtr.propagate_profile(ch.profile,
                max(ch.branch_length, min_branch_length),
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
    tree.gtr.Pi = genepi
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
    compute the total likelihood for all genes with presence pattern set to include
    conditioned on not observing pattern set to not include (e.g. the nullpattern)
    be careful: this function changes the gtr model
    """
    # change the relation of genegain and geneloss rate and the speed
    pi_present = params[0]
    mymu = params[1]
    change_gtr_parameters_forgainloss(tree,pi_present,mymu)

    # this gives the log likelihoods for each pattern
    compute_lh(tree)
    tree.tree.root.pattern_lh =  np.log(np.sum(np.exp(tree.tree.root.pattern_profile_lh)*tree.gtr.Pi,axis=1))
    #compute the likelihood of all genes with included pattern
    tree.tree.root.total_llh =  np.sum(tree.tree.root.pattern_lh * np.array(tree.tree.pattern_abundance) * np.array(tree.tree.pattern_include))
    #adjust for pattern that should not be included
    ll_forsumofignored = np.sum(np.exp( tree.tree.root.pattern_lh) * np.subtract(1,tree.tree.pattern_include))
    if verbose > 2:
        print("adjusting for all pattern that have been set to pattern_include == 0")
    tree.tree.root.total_llh = tree.tree.root.total_llh - ( np.log(1.- ll_forsumofignored) * np.sum(np.array(tree.tree.pattern_abundance) * np.array(tree.tree.pattern_include)) )
    if verbose > 3:
        print("totalLH:", pi_present, mymu, tree.tree.root.total_llh)
    if np.isnan(tree.tree.root.total_llh) or np.isinf(tree.tree.root.total_llh):
        return 1e50
    else:
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
    plt.close()

def plot_ll_mu(filename,tree,pi_present =0.5,mu_max = 10):
    import matplotlib.pyplot as plt
    xaxis = np.linspace(0.01, mu_max, num=500)
    graph = [compute_totallh(tree,[pi_present,x]) for x in xaxis]
    plt.plot(xaxis, graph)
    plt.savefig(filename)
    plt.close()
