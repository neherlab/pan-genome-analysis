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
    Example: python ./scripts/run-pipeline.py  -fn /ebio/ag-neher/share/users/wding/mpam/data/TestSet -sl TestSet-RefSeq.txt -st 1 2 3 4 5 6 7 8 9 10 11 -t 64 > TestSet.log 2>&1
```

##**Step-by-Step tutorial:**<br />
In `data/TestSet`, you will find a small set of four *Pseudomonas aeruginosa* genomes that is used in this tutorial. Your own data should also reside in such a folder within `data/` -- we will refer to this folder as *run directory* below. The name of the run directory is used as a species name in down-stream analysis.
To run `pan-genome-analysis` pipeline, you need to execute a series of steps that can be started using the `run.sh` script
To run ...
```
python ./scripts/run-pipeline -fn data/TestSet -sl TestSet-RefSeq.txt 1 2
```
<br />

**Step01: specify the set of strains**<br />
The pipeline can either down load sequences from GenBank or run on genomes you provide. You need to provide a file within the run directory that contains a list of NCBI RefSeq accession numbers for fetching GenBank files or the names of the files (without file ending) provided.<br />
If using own GenBank files, step02 can be skipped and corresponding GenBank files should be placed in the same folder where strain list is located.<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
TestSet-RefSeq.txt (accession list for download RefSeq strains)<br />
- Output:<br />
In folder `./data/TestSet/:`<br />
strain_list.cpk (cPickled file for the strain list )<br />

**Step02: download RefSeq GenBank (*.gbk) file**<br />
If using own GenBank files, step02 can be skipped. Otherwise, run 
to fetch NCBI RefSeq GenBank (\*.gbk) file from strain list<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
strain_list.cpk<br />
- Output:<br />
In folder `./data/TestSet/:`<br />
\*.gbk files<br />

**Step03: extract gene sequences from GenBank (*.gbk) file**<br />
Extract gene sequences in GenBank (\*.gbk) file for preparing nucleotide sequence (\*.fastaCpk) for gene cluster and amino acid sequences for Diamond input (\*.fna)<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
\*.gbk file<br />
- Output:<br />
In folder `./data/TestSet/:`<br />
\*.fastaCpk file (nucleotide sequences)<br />
In folder `./data/TestSet/protein_fna:`<br />
\*.fna file (amino acid sequences for DIAMOND input)<br />

**Step04: extract metadata from GenBank (\*.gbk) file (Alternative: use manually curated metadata table)**<br />
Extracting meta-information ( E.g.: country,collection_date, host, strain) or provide a simple tab-separated values (TSV) table.<br />

| strain        | location      | ...   |
| ------------- |:-------------:| -----:|
| NC_01         | Germany       | ...   |
| NC_02         | Switzerland   | ...   |

- Input:<br />
In folder `./data/TestSet/:`<br />
\*.gbk file<br />
- Output:<br />
In folder `./data/TestSet/:`<br />
metainfo_curated.txt and meta-dict-TestSet.js (metadata for visualization)<br />

**Step05: compute gene clusters**<br />
Conduct all-against-all protein sequences comparison by Diamond and cluster genes using Orthagogue and MCL<br />
- Input:<br />
In folder `./data/TestSet/protein_fna/:`<br />
\*.fna file<br />
- Output:<br />
In folder `./data/TestSet/protein_fna/diamond_matches/:`<br />
TestSet-orthamcl-allclusters.cpk (dictionary for gene clusters)<br />
diamond_geneCluster_dt: {clusterID:[ count_strains,[memb1,...],count_genes }<br />

**Step06: build alignments, gene trees from gene clusters and split paralogs**<br />
Load nucleotide sequences in gene clusters, construct nucleotide and amino acid alignment, build a gene tree based on nucleotide alignment, split paralogs and export the gene tree in json file for visualization<br />
- Input:<br />
In folder `./data/TestSet/protein_fna/diamond_matches/:`<br />
TestSet-orthamcl-allclusters.cpk file<br />
- Output:<br />
In folder `./data/TestSet/protein_fna/diamond_matches/:`<br />
TestSet-orthamcl-allclusters_final.cpk ( final gene clusters)<br />
In folder `./data/TestSet/geneCluster/:`<br />
GC\*.nu.fa (nucleotide fasta)<br />
GC\*.nu.aln (nucleotide alignment)<br />
GC\*.aa.fa (amino acid fasta)<br />
GC\*.aa.aln (amino acid alignment)<br />
GC\*.tree.json (gene tree in json file)<br />

**Step07: construct core gene SNP matrix**<br />
Call SNPs in strictly core genes (without no gene duplication) and build SNP matrix for strain tree<br />
- Input:<br />
In folder `./data/TestSet/protein_fna/diamond_matches/:`<br />
TestSet-orthamcl-allclusters_final.cpk file<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/:`<br />
SNP_whole_matrix.aln (SNP matrix as pseudo alignment)<br />
snp_pos.cpk (snp positions)<br />

**Step08:  build the strain tree using core gene SNPs**<br />
Use fasttree to build core genome phylogeny and further refine it by RAxML<br />
- Input:<br />
In folder `./data/TestSet/geneCluster/:`<br />
SNP_whole_matrix.aln<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/:`<br />
tree_result.newick<br />

**Step09: infer gene gain and loss event**<br />
Use ancestral reconstruction algorithm (treetime) to conduct gain and loss events inference<br />
- Input:<br />
In folder `./data/TestSet/geneCluster/:`<br />
TestSet-orthamcl-allclusters_final.cpk file<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/:`<br />
TestSet-genePresence.json (gene gain/loss event)<br />
dt_geneEvents.cpk (number of gene gain/loss events)<br />
tree_result.newick (final strain tree with inner nodes)<br />

**Step10: export gene cluster json file**<br />
Export json file for gene cluster datatable visualization<br />
- Input:<br />
In folder `./data/TestSet/geneCluster/:`<br />
TestSet-orthamcl-allclusters_final.cpk file (gene cluster dictionary)<br />
gene_diversity.cpk (diversity for each gene cluster)<br />
locusTag_to_geneId.cpk (locus_tag for each gene)<br />
dt_geneEvents.cpk (gain/loss event count)<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/`<br />
TestSet-geneCluster.json (gene cluster json for datatable visualization)<br />

**Step11: export tree and metadata json file**<br />
Export json files for tree and metadata visualization<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
metainfo_curated.txt (metadata table)<br />
In folder `./data/TestSet/geneCluster/:`<br />
tree_result.newick (strain tree)<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/`<br />
TestSet-tree-tnt-version.json (strain tree visualization)<br />
TestSet-tnt-nodeAttri-dataTable.json (strain metadata table visualization)
