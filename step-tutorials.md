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
strain_tree.nwk<br />

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
strain_tree.nwk (strain tree)<br />
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