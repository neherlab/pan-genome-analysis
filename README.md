# Microbial pan-genome analysis toolkit (MPA)
Author: Wei Ding and Richard Neher

Overview:
MPAM is based on an automated pan-genome identification pipeline that determines clusters of orthologous genes. The pipeline starts with a set of annotated sequences (e.g. NCBI RefSeq) of a bacterial species.
The genomes are split into individual genes and all genes from all strains are compared to each other via the fast protein alignment tool [DIAMOND](http://www.nature.com/nmeth/journal/v12/n1/full/nmeth.3176.html) and then clustered into orthologous groups using [orthAgogue](https://code.google.com/archive/p/orthagogue/) and MCL. After the construction of gene clusters, genes within clusters are aligned and the corresponding phylogenetic tree is computed, with mutations mapped into each tree and various summary statistics calculated.

1. All dependencies:
  1.1 Required softwares:
	(1) DIAMOND (fast protein alignment tool)
	Install: (source: https://github.com/bbuchfink/diamond)
	wget http://github.com/bbuchfink/diamond/releases/download/v0.7.12/diamond-linux64.tar.gz
	tar xzf diamond-linux64.tar.gz

    (2) orthAgogue
       ("a tool for high speed estimation of homology relations within and between species in massive data sets.")
       Packages and software requirements:
       Please install including all the required dependencies as specified here:
       https://code.google.com/archive/p/orthagogue/ 

    (3) MCL (Markov Cluster Algorithm)
       Install: (Homepage: http://micans.org/mcl/)
       sudo apt-get install mcl

    (4) mafft (multiple alignment program)
        Download and install:
        http://mafft.cbrc.jp/alignment/software/linux.html

    (5) fasttree
      sudo apt-get update
      sudo apt-get install fasttree

	(6) raxml
      Download and install:
      https://github.com/stamatak/standard-RAxML

  1.2 Required python packages:
    pip install numpy scipy biopython ete2
    !!treetime/nextstrain

2. How to run:
    sh run.sh
    Description:
    This calls run-pipeline.py to run each step using scripts located in folder ./scripts/
    run-pipeline.py [-h] -fn folder_name -sl strain_list
                       [-st steps [steps ...]] [-rt raxml_max_time]
                       [-t threads] [-bp blast_file_path]

    mandatory parameters: -fn folder_name / -sl strain_list / [-st steps [steps ...]]
    NOTICE: strain_list format should be species_name+'-RefSeq', e.g.: Saureus-RefSeq.txt
    Example: python ./scripts/run-pipeline.py  -fn /ebio/ag-neher/share/users/wding/mpam/data/Pat3 -sl Pat3-RefSeq.txt -st 11 -t 64 > Pat3-11.log 2>&1

