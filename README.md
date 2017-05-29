# panX: microbial pan-genome analysis and exploration
Author: Wei Ding, Franz Baumdicker and Richard Neher

Overview:
[**panX**](http://pangenome.de) is a software package for pan-genome analysis, interactive visualization and exploration. The analysis pipeline is based on [DIAMOND](https://github.com/bbuchfink/diamond) ([Buchfink et al. 2015 Nature Methods](http://www.nature.com/nmeth/journal/v12/n1/full/nmeth.3176.html)), MCL and phylogeny-aware post-processing, which takes a set of annotated bacterial strains as input (e.g. NCBI RefSeq records or user's own data in **GenBank** format).

Alls genes from all strains are compared to each other via the fast protein alignment tool DIAMOND and then clustered into orthologous groups using MCL and adaptive phylogenetic post-processing, which split distantly related genes and paralogs if necessary. For each gene cluster, corresponding alignment and phylogeny are constructed. All core gene SNPs are then used to build strain/species phylogeny.

The results can be interactively explored using a [**powerful web-based visualization application**](https://github.com/neherlab/pan-genome-visualization) (either hosted by web server or used locally on desktop). The web application integrates various interconnected viewers (pan-genome statistical charts, gene cluster table, alignment, comparative phylogenies, metadata table) and allows rapid search and filter of gene clusters by gene name, annotation, duplication, diversity, gene gain/loss events, etc. Strain-specific **metadata** are integrated into strain phylogeny such that genes related to adaptation, antibiotic resistance, virulence etc can be readily identified.

### **Quick start:**

`git clone https://github.com/neherlab/pan-genome-analysis.git`

Enter the folder `pan-genome-analysis`:
`git submodule update --init`

Install dependencies and then run the test:
`sh run-TestSet.sh`

The results can be explored via our interactive [**pan-genome-visualization**](https://github.com/neherlab/pan-genome-visualization) application.

### **Pipeline overview:**

![panX](/panX-pipeline.png)

1. Dependencies:
  1. Required software:
      * [DIAMOND](https://github.com/bbuchfink/diamond) (already located in `./tools/diamond`)
      * [MCL](http://micans.org/mcl/)
        `sudo apt-get install mcl`
      * [mafft](http://mafft.cbrc.jp/alignment/software/)
        `sudo apt-get install mafft`
      * [fasttree](http://www.microbesonline.org/fasttree/)
        `sudo apt-get install fasttree`
      * [raxml](https://github.com/stamatak/standard-RAxML)
        `sudo apt-get install raxml`

  2. Required python packages:
      - `pip install numpy scipy biopython ete2`
      - [treetime](http://github.com/neherlab/treetime):
      `git submodule update --init`

2. How to run:

    To run the test set: ` sh run-TestSet.sh `

    In `data/TestSet`, you will find a small set of five *Mycoplasma genitalium* genomes that is used in this tutorial. Your own data should also reside in such a folder within `data/` -- we will refer to this folder as *run directory* below. The name of the run directory is used as a species name in down-stream analysis.

    All steps can be run in order by omitting the `-st` option, whereas using `-st 5 6` will specify the analysis steps. `-t ` sets the number of CPU cores.
    <br />
    ```
    ./panX.py -fn data/TestSet -sl TestSet-RefSeq.txt -t 32 > TestSet.log 2> TestSet.err
    ```

    This calls panX.py to run each step using scripts located in folder ./scripts/
    ```
    ./panX.py [-h] -fn folder_name -sl strain_list
                       [-st steps [steps ...]] [-rt raxml_max_time]
                       [-t threads] [-bp blast_file_path]

    Mandatory parameters: -fn folder_name / -sl strain_list
    NOTICE: strain_list format should be species_name+'-RefSeq', e.g.: Saureus-RefSeq.txt
    Example: ./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -t 32 > TestSet.log 2> TestSet.err
    ```

    **Useful options**:

      Soft core-gene:

        -cg    core-genome threshold [e.g.: 0.7] percentage of strains used to decide whether a gene is core
      Large dataset (use divide-and-conquer(DC) strategy which scales approximately linearly with the number of genomes):

        -dmdc  apply DC strategy to run DIAMOND on subsets and then combine the results
        -dcs   subset size used in DC strategy [default:50]

        E.g.: ./panX.py -dmdc -dcs 50 -fn ...


    The results contain files required for visualizing the pan-genome using [pan-genome-visualization](https://github.com/neherlab/pan-genome-visualization).
    ```
    ./data
        YourSpecies               # folder specific to the your pan genome
          - YourSpecies-RefSeq.txt    # INPUT: GenBank accession numbers
          - input_GenBank              # INPUT: genomes in GenBank format
            - strain1.gbk
            - strain2.gbk
            ...
          - vis
            - geneCluster.json       # for clusters table: gene clusters and their summary statistics
            - strainMetainfo.json    # for metadata table: strain-associated metadata
            - coreGenomeTree.json    # core genome SNP tree (json file)
            - strain_tree.nwk        # core genome SNP tree (newick file)
            - geneCluster/           # folder contain orthologous clusters
              - GC00000001_na_aln.fa
              - GC00000001_aa_aln.fa
              - GC00000001_na_aln_reduced.fa
              - GC00000001_aa_aln_reduced.fa
              - GC00000001_tree.json
              - GC00000001_patterns.json
    ```
    In which step different files and directories are produced is described in more details in step-tutorials.txt.

