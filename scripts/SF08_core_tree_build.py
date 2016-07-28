import os;from SF00miscellaneous import times
from ete2 import Tree

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
    
def aln_to_Newick(path, raxml_timelimit, threads):
    """ function: build core gene SNP tree using SNP alignment
        input: SNP_whole_matrix.aln 
        output: tree_result.newick
    """
    import sys, time, random, glob, subprocess, shutil
    output_path = '_'.join(['temp_coretree', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
    os.system('mkdir %s'%output_path)
    os.system('ln -sf ../%s/SNP_whole_matrix.aln %s'%(path+'geneCluster',output_path))
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
        process = subprocess.Popen('exec raxml -f d -T '+str(threads)+' -j -s SNP_whole_matrix.aln -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick', shell=True)
        while (time.time() < end_time):
            if os.path.isfile('RAxML_result.topology'):
                break
            time.sleep(10)
        process.terminate()

        checkpoint_files = [file for file in glob.glob('RAxML_checkpoint*')]
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
    os.system('raxml -f e -T 6 -s SNP_whole_matrix.aln -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick')
    shutil.copy('RAxML_result.branches', out_fname)

    print ' raxml time-cost:', times(start)
    midpointRooting(out_fname,'tree_result.newick')
    shutil.copy('tree_result.newick', '../%s/tree_result.newick'%(path+'geneCluster') )
    os.chdir('../')
    os.system('rm -r %s'%output_path)