# panX: Microbial pan-genome analysis and exploration toolkit
Author: Wei Ding, Franz Baumdicker and Richard Neher

Overview:
panX is based on an automated pan-genome identification pipeline that determines clusters of orthologous genes. The pipeline starts with a set of annotated sequences (e.g. NCBI RefSeq) of a bacterial species.
The genomes are split into individual genes and all genes from all strains are compared to each other via the fast protein alignment tool [DIAMOND](http://www.nature.com/nmeth/journal/v12/n1/full/nmeth.3176.html) and then clustered into orthologous groups using MCL and panX post-processing procedures. After the construction of gene clusters, genes within clusters are aligned and the corresponding phylogenetic tree is computed, with mutations mapped into each tree and various summary statistics calculated.

Quick start:
```git clone https://github.com/neherlab/pan-genome-analysis.git```

Enter the folder ```pan-genome-analysis```:
```git submodule update --init```

Install dependencies and then run the test:
```sh run-TestSet.sh```

1. Dependencies:
  1. Required software:
    * [DIAMOND](https://github.com/bbuchfink/diamond) (fast protein alignment tool)
      - located in ./tools/diamond
    * [MCL Markov Cluster Algorithm](http://micans.org/mcl/)
      - sudo apt-get install mcl
    * mafft (multiple alignment program)
      - sudo apt-get install mafft
      - Or: download and install from http://mafft.cbrcj.p/alignment/software/linux.html
    * [fasttree](http://www.microbesonline.org/fasttree/)
      - sudo apt-get install fasttree
    * [raxml](https://github.com/stamatak/standard-RAxML)
      - sudo apt-get install raxml

  2. Required python packages:
     - pip install numpy scipy biopython ete2
     - [treetime](http://github.com/neherlab/treetime) (please fetch treetime using the following two commands)
     ```
       git submodule init
       git submodule update
     ```

2. How to run:
  - sh run.sh
  ```
    Description:
    This calls panX.py to run each step using scripts located in folder ./scripts/
    panX.py [-h] -fn folder_name -sl strain_list
                       [-st steps [steps ...]] [-rt raxml_max_time]
                       [-t threads] [-bp blast_file_path]

    mandatory parameters: -fn folder_name / -sl strain_list / [-st steps [steps ...]]
    NOTICE: strain_list format should be species_name+'-RefSeq', e.g.: Saureus-RefSeq.txt
    Example: ./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -t 64 > TestSet.log 2>&1
  ```
  The result will be a number of files that contain the all the information necessary for visualizing the pan-genome in the browser using [pan-genome-visualization](https://github.com/neherlab/pan-genome-visualization).
  ```
  ./data
      YourSpecies               # folder specific to the your pan genome
        - YourSpecies-RefSeq.txt    # INPUT: GenBank accession numbers
        - input_GenBank              # INPUT: genomes in GenBank format
          - strain1.gbk
          - strain2.gbk
          ...
        - vis
          - geneCluster.json
          - strainMetainfo.json
          - coreGenomeTree.json
          - geneCluster/       # Folder contain orthologous clusters
            - GC000001*_na_aln.fa
            - GC000001*_aa_aln.fa
            - GC000001*_tree.json
            - GC000001*_patterns.json
            ...
  ```
  In which step of the analysis different files and directories are produced is described in more detail below.


##**Step-by-Step tutorial:**<br />
In `data/TestSet`, you will find a small set of five *Mycoplasma genitalium* genomes that is used in this tutorial. Your own data should also reside in such a folder within `data/` -- we will refer to this folder as *run directory* below. The name of the run directory is used as a species name in down-stream analysis.
To run `pan-genome-analysis` pipeline, you need to execute a series of steps that can be started using the `run-TestSet.sh` script
To run ...
```
python ./scripts/run-pipeline -fn data/TestSet -sl TestSet-RefSeq.txt -t 32
```
All steps can be run in order by omitting the `-st` option, whereas using `-st 5 6` will specify the analysis steps you want to run.
<br />

**Step01: specify the set of strains**<br />
The pipeline can either download sequences from GenBank or run on genomes you provide. You need to provide a file within the run directory that contains a list of NCBI RefSeq accession numbers for fetching GenBank files or the names of the files (without file ending) provided.<br />
If using own GenBank files, step02 can be skipped and corresponding GenBank files should be placed in the same folder where strain list is located.<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
TestSet-RefSeq.txt<br />
- Output:<br />
In folder `./data/TestSet/:`<br />
strain_list.cpk (cPickled file for the strain list )<br />

**Step02: download RefSeq GenBank (*.gbk) file**<br />
If using own GenBank files, step02 can be skipped. Otherwise, run
to fetch NCBI RefSeq GenBank (\*.gbk) file from strain list<br />
- Input: In folder `./data/TestSet/:`<br />
strain_list.cpk<br />
- Output: In folder `./data/TestSet/:`<br />
\*.gbk files<br />

**Step03: extract gene sequences from GenBank (*.gbk) file**<br />
Extract gene sequences in GenBank (\*.gbk) file for preparing nucleotide sequence (\*_genes_dict.cpk) for gene cluster and amino acid sequences for Diamond input (\*.faa)<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
\*.gbk file<br />
- Output:<br />
In folder `./data/TestSet/nucleotide_fna:`<br />
\*.fna file (nucleotide sequences)<br />
In folder `./data/TestSet/protein_faa:`<br />
\*.faa file (amino acid sequences for DIAMOND input)<br />

**Step04: extract metadata from GenBank (\*.gbk) file (Alternative: use manually curated metadata table)**<br />
Extracting meta-information ( E.g.: country,collection_date, host, strain) or provide a simple tab-separated values (TSV) table.<br />

| strain | location    | ... |
| -------|:-----------:| ---:|
| NC_01  | Germany     | ... |
| NC_02  | Switzerland | ... |

- Input:<br />
In folder `./data/TestSet/:`<br />
\*.gbk file<br />
- Output:<br />
In folder `./data/TestSet/:`<br />
metainfo.tsv  (metadata for visualization)<br />

**Step05: compute gene clusters**<br />
Conduct all-against-all protein sequences comparison by Diamond and cluster genes using MCL<br />
- Input:<br />
In folder `./data/TestSet/protein_faa/:`<br />
\*.faa file<br />
- Output:<br />
In folder `./data/TestSet/protein_faa/diamond_matches/:`<br />
allclusters.cpk (dictionary for gene clusters)<br />
diamond_geneCluster_dt: {clusterID:[ count_strains,[memb1,...],count_genes }<br />

**Step06: build alignments, gene trees from gene clusters and run phylogeny-based post-processing**<br />
Load nucleotide sequences in gene clusters, construct nucleotide and amino acid alignment, build a gene tree based on nucleotide alignment, split paralogs and export the gene tree in json file for visualization<br />
- Input:<br />
In folder `./data/TestSet/protein_faa/diamond_matches/:`<br />
allclusters.cpk file<br />
- Output:<br />
In folder `./data/TestSet/protein_faa/diamond_matches/:`<br />
allclusters_final.tsv ( final gene clusters)<br />
In folder `./data/TestSet/geneCluster/:`<br />
GC\*.fna (nucleotide fasta)<br />
GC\*_na_aln.fa (nucleotide alignment)<br />
GC\*.faa (amino acid fasta)<br />
GC\*_aa_aln.fa (amino acid alignment)<br />
GC\*_tree.json (gene tree in json file)<br />

**Step07: construct core gene SNP matrix**<br />
Call SNPs in strictly core genes (without no gene duplication) and build SNP matrix for strain tree<br />
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
- Output:<br />
In folder `./data/TestSet/geneCluster/:`<br />
genePresence.aln  (gene presence and absence pattern)<br />
dt_genePresence.cpk (gene presence and absence pattern as cPickled file)<br />
dt_geneEvents.cpk (number of gene gain/loss events)<br />
GC000\*_patterns.json (gene gain/loss pattern for each gene cluster)<br />

**Step10: export gene cluster json file**<br />
Export json file for gene cluster datatable visualization<br />
In folder `./data/TestSet/geneCluster/:`<br />
gene_diversity.cpk (diversity for each gene cluster)<br />
dt_geneEvents.cpk (gain/loss event count)<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/`<br />
geneCluster.json (gene cluster json for datatable visualization)<br />

**Step11: export tree and metadata json file**<br />
Export json files for tree and metadata visualization<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
metainfo.tsv (metadata table)<br />
In folder `./data/TestSet/geneCluster/:`<br />
tree_result.newick (strain tree)<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/`<br />
coreGenomeTree.json (strain tree visualization)<br />
strainMetainfo.json (strain metadata table visualization)
- Data collection for visualization (sending data to server)
In folder `./data/TestSet/vis/`<br />
geneCluster.json
coreGenomeTree.json
strainMetainfo.json
In folder `./data/TestSet/vis/geneCluster/`<br />
GC000\*_na_aln.fa
GC000\*_aa_aln.fa
GC000\*_tree.json
GC000\*_patterns.json
