import os, sys, time, random, glob, subprocess, shutil
from sf_miscellaneous import times
from ete2 import Tree
from sf_geneCluster_align_makeTree import mpm_tree

def resolve_polytomies(infileName, outfileName):
    newickString=open(infileName, 'rb').readline().rstrip().replace('[&R] ', '')
    tree = Tree(newickString);
    tree.resolve_polytomy(recursive=True)
    with open(outfileName, 'wb') as outfile:
        outfile.write(tree.write(format=1))

def midpointRooting(infileName, outfileName):
    """ using ete2 for mid-point rooting """
    newickString=open(infileName, 'rb').readline().rstrip().replace('[&R] ', '')
    tree = Tree(newickString);
    #tree.resolve_polytomy(recursive=True)
    if tree.get_midpoint_outgroup()!=None:
        tree.set_outgroup( tree.get_midpoint_outgroup() )
    tree.ladderize()
    with open(outfileName, 'wb') as outfile:
        outfile.write(tree.write(format=1))

def aln_to_Newick(path, raxml_timelimit, raxml_path, threads, scratch=''):
    """ function: build core gene SNP tree using SNP alignment
        input: SNP_whole_matrix.aln
        output: tree_result.newick
    """
    output_path = (scratch.rstrip('/')+'/' if len(scratch) else '') + '_'.join(['temp_coretree', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
    print(path, '../%s/SNP_whole_matrix.aln'%output_path, glob.glob('../%s/SNP_whole_matrix.aln'%output_path))
    cwd = os.getcwd()
    os.system('mkdir %s'%output_path)
    os.system('ln -sf %s/SNP_whole_matrix.aln %s'%(path+'geneCluster',output_path))
    os.chdir(output_path)

    ## run fasttree
    start = time.time();

    os.system('fasttree -gtr -nt -gamma -nosupport -mlacc 2 -slownni SNP_whole_matrix.aln > initial_tree.newick0') ;
    print ' fasttree time-cost:', times(start)

    resolve_polytomies('initial_tree.newick0','initial_tree.newick')

    ## run raxml
    start = time.time();
    out_fname = "tree_infer.newick"
    if raxml_timelimit>0:
        print '%s%d%s'%('RAxML tree optimization within the timelimit of ',raxml_timelimit, ' minutes')
        # exec for killing process
        end_time = time.time() + int(raxml_timelimit*60) #
        raxml_program= 'raxml' if raxml_path=='' else raxml_path
        process = subprocess.Popen('exec '+raxml_program+' -f d -T '+str(threads)+' -j -s SNP_whole_matrix.aln -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick', shell=True)
        while (time.time() < end_time):
            if os.path.isfile('RAxML_result.topology'):
                break
            time.sleep(10)
        process.terminate()

        checkpoint_files = glob.glob('RAxML_checkpoint*')
        if os.path.isfile('RAxML_result.topology'):
            checkpoint_files.append('RAxML_result.topology')
        if len(checkpoint_files) > 0:
            last_tree_file = checkpoint_files[-1]
            shutil.copy(last_tree_file, 'raxml_tree.newick')
        else:
            shutil.copy('initial_tree.newick', 'raxml_tree.newick')
    else:
        shutil.copy('initial_tree.newick', 'raxml_tree.newick')

    print 'RAxML branch length optimization and rooting'
    os.system(raxml_program+' -f e -T 6 -s SNP_whole_matrix.aln -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick')
    shutil.copy('RAxML_result.branches', out_fname)

    print ' raxml time-cost:', times(start)
    midpointRooting(out_fname,'tree_result.newick')
    shutil.copy('tree_result.newick', '%s/tree_result.newick'%(path+'geneCluster') )
    os.chdir(cwd)
    if not 'scratch' in output_path:
        os.system('rm -r %s'%output_path)

    # if 0: phyloTree based visualization
    #     simple_tree=0
    #     ##new tree layout
    #     cluster_seqs_path = '%s%s'%(path,'geneCluster/')
    #     myTree = mpm_tree('%s%s'%(cluster_seqs_path,'SNP_whole_matrix.aln'))
    #     myTree.align()
    #     myTree.aa_aln=None
    #     if simple_tree==0:
    #         myTree.build(raxml=True,treetime_used=True)
    #         myTree.ancestral(translate_tree=True)
    #         myTree.refine()
    #     else:
    #         myTree.build(raxml=True,treetime_used=False)
    #     myTree.layout()
    #     myTree.export(path=cluster_seqs_path)
