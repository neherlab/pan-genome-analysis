from __future__ import print_function, division
import os,sys,copy;import numpy as np
from treetime import treeanc as ta
from treetime.gtr import GTR
from treetime import io
from Bio import Phylo, AlignIO
from SF00miscellaneous import write_json, write_pickle
import sys


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
    nwk =  sep.join([path.rstrip(sep), 'geneCluster', 'tree_result.newick'])

    # instantiate treetime with custom GTR
    t = io.treetime_from_newick(gain_loss_model, nwk)
    # fix leaves names since Bio.Phylo interprets numeric leaf names as confidence
    for leaf in t.tree.get_terminals():
        if leaf.name is None:
            leaf.name = str(leaf.confidence)
    # load alignment and associate with tree leafs
    io.set_seqs_to_leaves(t, AlignIO.read(fasta, 'fasta'))

    t.reconstruct_anc(method='ml')

    for n in t.tree.find_clades():
        n.genepresence = n.sequence

    return t


def export_gain_loss(tree, path, species):
    '''
    '''
    # write final tree with internal node names as assigned by treetime
    sep='/'
    output_path= sep.join([path.rstrip(sep), 'geneCluster'])
    tree_fname = sep.join([output_path, 'tree_result.newick'])
    Phylo.write(tree.tree, tree_fname, 'newick')

    from collections import defaultdict
    gene_gain_loss_dict=defaultdict(str)
    for node in tree.tree.find_clades(order='preorder'):# order does not matter much here
            if node.up is None: continue
            #print(node.name ,len(node.geneevents),node.geneevents)
            gain_loss = [ str(int(ancestral)*2+int(derived))
                        for ancestral,derived in zip(node.up.genepresence, node.genepresence)]
            gene_gain_loss_dict[node.name]="".join(gain_loss)

    gain_loss_array = np.array([[i for i in gain_loss_str]
                                for gain_loss_str in gene_gain_loss_dict.values()], dtype=int)
    # 1 and 2 are codes for gain/loss events
    events_array = ((gain_loss_array == 1) | (gain_loss_array == 2)).sum(axis=0)
    events_dict =  { index:event for index, event in enumerate(events_array) }
    events_dict_path= sep.join([output_path, 'dt_geneEvents.cpk'])
    write_pickle(events_dict_path, events_dict)

    # export gene loss dict to json for visualization
    gene_loss_fname = sep.join([ output_path, species+'-genePresence.json'])
    write_json(gene_gain_loss_dict, gene_loss_fname, indent=1)

def process_gain_loss(path, species):
    tree = infer_gene_gain_loss(path)
    export_gain_loss(tree, path, species)