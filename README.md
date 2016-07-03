# Microbial pan-genome analysis toolkit (MPA)
Author: Wei Ding and Richard Neher

Overview:
MPAM is based on an automated pan-genome identification pipeline that determines clusters of orthologous genes. The pipeline starts with a set of annotated sequences (e.g. NCBI RefSeq) of a bacterial species.
The genomes are split into individual genes and all genes from all strains are compared to each other via the fast protein alignment tool [DIAMOND](http://www.nature.com/nmeth/journal/v12/n1/full/nmeth.3176.html) and then clustered into orthologous groups using [orthAgogue](https://code.google.com/archive/p/orthagogue/) and MCL. After the construction of gene clusters, genes within clusters are aligned and the corresponding phylogenetic tree is computed, with mutations mapped into each tree and various summary statistics calculated.

1. Dependencies:

  1.1 Required software:

    * DIAMOND (fast protein alignment tool)
      - Install: (source: https://github.com/bbuchfink/diamond)
      - wget http://github.com/bbuchfink/diamond/releases/download/v0.7.12/diamond-linux64.tar.gz
      - tar xzf diamond-linux64.tar.gz

    * orthAgogue: Please install including all the required dependencies as specified [here] (https://code.google.com/archive/p/orthagogue/)

    * [MCL Markov Cluster Algorithm](http://micans.org/mcl/)
      - sudo apt-get install mcl

    * mafft (multiple alignment program)
      - Download and install from http://mafft.cbrcj.p/alignment/software/linux.html
      - OR sudo apt-get install mafft

    * [fasttree](http://www.microbesonline.org/fasttree/)
      - sudo apt-get install fasttree

    * [raxml](https://github.com/stamatak/standard-RAxML)
      - sudo apt-get install raxml

  1.2 Required python packages:
     - pip install numpy scipy biopython ete2
     - [treetime](http://github.com/neherlab/treetime)

2. How to run:
  - sh run.sh
```    
    Description:
    This calls run-pipeline.py to run each step using scripts located in folder ./scripts/
    run-pipeline.py [-h] -fn folder_name -sl strain_list
                       [-st steps [steps ...]] [-rt raxml_max_time]
                       [-t threads] [-bp blast_file_path]

    mandatory parameters: -fn folder_name / -sl strain_list / [-st steps [steps ...]]
    NOTICE: strain_list format should be species_name+'-RefSeq', e.g.: Saureus-RefSeq.txt
    Example: python ./scripts/run-pipeline.py  -fn /ebio/ag-neher/share/users/wding/mpam/data/Pat3 -sl Pat3-RefSeq.txt -st 11 -t 64 > Pat3-11.log 2>&1
```
