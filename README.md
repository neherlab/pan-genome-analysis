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
            - geneCluster.json
            - strainMetainfo.json
            - coreGenomeTree.json
            - geneCluster/       # Folder contain orthologous clusters
              - GC00000001_na_aln.fa
              - GC00000001_aa_aln.fa
              - GC00000001_na_aln_reduced.fa
              - GC00000001_aa_aln_reduced.fa
              - GC00000001_tree.json
              - GC00000001_patterns.json
    ```
    In which step different files and directories are produced is described in more details below.


### **Step-by-Step tutorial:**<br />


**Step01: specify the set of strains**<br />
Load the strain file within the run directory which contains a list of NCBI RefSeq accession numbers or names of own GenBank files (without file ending).<br />

**Step03: extract gene sequences from GenBank (*.gbk) file**<br />
Extract genes from GenBank (\*.gbk) file as nucleotide and amino acid sequences<br />
- Input:<br />
In folder `./data/TestSet/:`<br />
\*.gbk file<br />
- Output:<br />
In folder `./data/TestSet/nucleotide_fna:`<br />
\*.fna file (nucleotide sequences)<br />
In folder `./data/TestSet/protein_faa:`<br />
\*.faa file (amino acid sequences)<br />

**Step04: extract metadata from GenBank (\*.gbk) file (Alternative: provide manually curated metadata table)**<br />
Extracting metadata ( E.g.: country, collection_date, host, strain) or provide a tab-separated values (TSV) file.<br />

| strain | location    | host age| serotype | benzylpenicillin MIC (ug/mL) |... |
| -------|:-----------:| -------:| --------:| ----------------------------:|---:|
| NC_01  | Germany     | 35      | 23A      | 0.016                        |... |
| NC_02  | Switzerland | 66      | 23B      | 4                            |... |

- Input:<br />
In folder `./data/TestSet/:`<br />
\*.gbk file<br />
- Output:<br />
In folder `./data/TestSet/:`<br />
metainfo.tsv  (metadata for visualization)<br />

**Step05: compute gene clusters**<br />
all-against-all protein sequences comparison by Diamond and clustering of genes using MCL<br />
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
Call SNPs in strictly core genes (without gene duplication) and build SNP matrix for strain tree<br />
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
Use ancestral reconstruction algorithm (treetime) to infer gain and loss events<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/:`<br />
genePresence.aln  (gene presence and absence pattern)<br />
GC000\*_patterns.json (gene gain/loss pattern for each gene cluster)<br />

**Step10: export gene cluster json file**<br />
Export json file for gene cluster datatable visualization<br />
In folder `./data/TestSet/geneCluster/:`<br />
- Output:<br />
In folder `./data/TestSet/geneCluster/`<br />
geneCluster.json (gene cluster json for datatable visualization)<br />

**Step11: export tree and metadata json file**<br />
Export json files for strain tree and metadata visualization<br />
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
