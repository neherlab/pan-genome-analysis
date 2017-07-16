## Table of contents
  * [Mandatory parameters](#mandatory-parameters)
  * [Specify steps and threads](#specify-steps-and-threads)
  * [Clustering](#clustering)
    + [Parameters used to run Diamond](#parameters-used-to-run-diamond)
    + [Running diamond with divide-and-conquer algorithm for large dataset](#running-diamond-with-divide-and-conquer-algorithm-for-large-dataset)
    + [Other clustering parameters](#other-clustering-parameters)
  * [Cluster post-processing](#cluster-post-processing)
    + [Split long branches](#split-long-branches)
    + [Split paralogy](#split-paralogy)
    + [Resolve unclustered genes](#resolve-unclustered-genes)
  * [Metadata processing](#metadata-processing)
  * [Core genome cutoff](#core-genome-cutoff)
  * [Association inference](#association-inference)
  * [Other options](#other-options)

```
Usage example: ./panX -fn data/TestSet -sl TestSet-RefSeq.txt 1>TestSet.log 2>TestSet.err
For help: ./panX -h
```
##### Mandatory parameters
  - -fn --folder_name

    the absolute path for project folder
  - -sl --strain_list

    the file name for strain list (e.g.: Pa-RefSeq.txt)

##### Specify steps and threads
panX runs through all steps by default, which contain pan-genome analysis and data preparation for visualization.
Alternatively, user can run specific steps.
  - -st --steps

    e.g.: -st 5 6 (It supposes that steps before 5 have been already finished.)
  - -t --threads

    number of threads (default:1)

#### Clustering
panX generates orthologous gene clusters by itself and can alternatively use clustering output from other pan-genome tools. To ingest the output from roary, orthofinder, or other tools, one can specify the path of output from these tools via:
  - -rp --roary_file_path

    the absolute path of roary result (e.g.: /path/roary.out) ,
  - -op --orthofinder_file_path

    the absolute path of orthofinder result (e.g.: /path/orthofinder.out)

For other tools, panX expects a tab-separated values (TSV) file containing each gene cluster per row and each gene separated by | between strain name and locus_tag (strain01|gene01)
  - -otp

    the absolute path of result from other orthology inference tool  (e.g.: /path/other_tool.out)

##### Parameters used to run diamond
[DIAMOND](https://github.com/bbuchfink/diamond) [(Buchfink et al. 2015 Nature Methods)] (http://www.nature.com/nmeth/journal/v12/n1/full/nmeth.3176.html)

  - -dmp --diamond_path

    alternative diamond path provided by user
  - -dme --diamond_evalue (default:0.001)

    e-value threshold
  - -dmt --diamond_max_target_seqs (default:600)

    maximum number of target sequences per query
    estimation: #strain * #max_duplication
    (50*10=500, max_duplication is the maximum of the counts of gene duplication in all strains)
  - -dmi --diamond_identity (default:0)

    sequence identity threshold to report an alignment
  - -dmqc --diamond_query_cover (default:0)

    query sequence coverage threshold to report an alignment
  - -dmsc --diamond_subject_cover (default:0)

    subject sequence coverage threshold to report an alignment

##### Running diamond with divide-and-conquer algorithm for large dataset
  - -dmdc --diamond_divide_conquer

    using divide-and-conquer(DC) algorithm
  - -dcs --subset_size (default:50)

    subset size (number of strains in a subset) used in divide-and-conquer(DC) algorithm
  - -dmsi --diamond_identity_subproblem (default:90)

    subproblem alignment: sequence identity threshold to report an alignment
  - -dmsqc --diamond_query_cover_subproblem (default:90)

    subproblem alignment: query sequence coverage threshold to report an alignment
  - -dmssc --diamond_subject_cover_subproblem (default:90)

    subproblem alignment: subject sequence coverage threshold to report an alignment

panX can also use a matrix of all-against-all blast hits to generate clusters with MCL.
  - -bp --blast_file_path

    the absolute path for blast result (e.g.: /path/blast.out)

##### Other clustering parameters
MCL parameter
  - -imcl --mcl_inflation (default:1.5)

  inflation parameter (this parameter affects granularity)

RNA clustering
  - -rna --enable_RNA_clustering cluster rRNAs

Running Blastn on RNAs
  - -bmt --blastn_RNA_max_target_seqs (default:100)

    the maximum number of target sequences per query; estimation: #strain * #max_duplication

Disable postprocessing
 - -np --disable_cluster_postprocessing

    Not apply postprocessing

#### Cluster post-processing
##### Split long branches
  For each gene tree, panX scans the branches longer than the diversity cutoff and splits the tree from those branches into sub-clusters.
  The diversity cutoff is defined by the formula: (0.1+fcd*core_diversity)/(1+fcd*core_diversity)
  - -fcd --factor_core_diversity (default:2.0)

    factor used to define the diversity cutoff for splitting long branches on gene tree of a cluster

  Alternatively, the diversity cutoff can be specified explicitly with the argument -slb, which is mutually exclusive from -fcd .
  If this value is not given (by default), panX uses the diversity cutoff calculated from core genome diversity as described in above-mentioned formula.
  - -slb --split_long_branch_cutoff (default: '')

    using the diversity cutoff specified by user

##### Split paralogy
Relatively recent gene duplications (paralogs) are not considered when only applying long branch splitting.
panX searches paralogy evidence via two tree traversals by taking both paralogy score and branch length into consideration. In order to identify these paralogs and split them into orthologs, a paralogy score is computed for each branch on each gene tree, which is the count of strains having genes on both sides of the branch.

The splitting procedure can be summarized in these conditions:
branch_length/paralog_branch_cutoff + paralog_score/(1.5 * nstrains)>1 and paralog_score > 1
panX uses the defined diversity cutoff as paralog branch_cutoff by default.
The paralogy splitting is iteratively conducted until no paralogs can be found.

Alternatively, paralog branch_cutoff can be specified with the argument -pbc instead of using the defined diversity cutoff.
  - -pbc --paralog_branch_cutoff' (default:'')

    branch_length cutoff for paralogy splitting provided by user

##### Resolve unclustered genes
  (peaks in gene length distribution)

  From our experience, there can be sometimes a small number of unclustered genes, despite of their sequence similarity.
  panX spots those sequences by inspecting the peaks in the gene length distribution (average length of genes from each cluster), statistically speaking, by comparing the length distribution to smoothed length distribution.

  Smoothed length distribution is calculated based on the convolution with defined sliding window by formula:
      smoothed_length_distribution = convolve(cluster_length_distribution, **window**)
    The signal is detected when:
    (cluster_length_distribution - smoothed_length_distribution) > maximum(**strain_proportion**\*nstrains, **sigma_scale**\*squareRoot(smoothed_length_distribution))
  The related parameters are:
  - -ws --window_size_smoothed (default:5)

    postprocess_unclustered_genes: window size for smoothed cluster length distribution
  - -spr --strain_proportion (default:0.3)

    postprocess_unclustered_genes: strain proportion
  - -ss --sigma_scale  (default:3)

    postprocess_unclustered_genes: sigma scale
  For each peak, all genes in these clusters are aligned, with corresponding phylogeny being inferred.
  Sub-clusters are then split using long branch splitting principle as described above.

#### Metadata processing
  panX parses metadata from GenBank file by default.
  To link customized metadata to the strains, a TSV table needs to be provided. [(example)](https://github.com/neherlab/pan-genome-analysis/blob/master/metadata/metainfo.tsv)
  This will only use the meta-information given by user.
  - -mi --metainfo_fpath

    the absolute path for meta_information file (e.g.: /path/meta.out)

For metadata such as drug concentration measurements, which user might want to have the branch and presence/absence association to be visualized, an extra TSV table is required to specify what metadata types are used for above-mentioned association inference. [(example)](https://github.com/neherlab/pan-genome-analysis/blob/master/metadata/meta_config.tsv)
  - -mtf --meta_tidy_fpath

    the absolute path for pre-defined metadata structure (discrete/continuous data type, etc.)

#### Core genome cutoff
  Core genome can be defined either as **strictly core** or **soft core**.
  Strictly core genome contain genes present in **all** strains;
  soft core genome includes genes shared by the **majority** of strains (e.g.: 80% of strains).
  Correspondingly, setting core genome threshold affects core genome diversity and strain tree based on core gene SNPs, specially when incomplete genomes or outliers are involved in the dataset.
  - -cg --core_genome_threshold (default:1.0)

    percentage of strains used to decide whether a gene is core (1.0 for strictly core gene; < 1.0 for soft core genes)

  Core gene strain constraint
  If incomplete genomes are included in the analysis, user can provide an extract list to ensure that soft core genes must be found in all strains in the list.
  - -csf --core_gene_strain_fpath

    the absolute path for user-provided subset of strains for restricting core genes

#### Association inference
  - -iba --infer_branch_association

    infer branch association

  - -mtf  --meta_tidy_fpath

  Example: [meta_config.tsv](https://github.com/neherlab/pan-genome-analysis/blob/master/metadata/meta_config.tsv)

  E.g.: ./panX.py -iba -mtf ./data/yourSpecies/meta_config.tsv -fn ...

#### Other options
  - -rxm --raxml_path

    absolute path of raxml
  - -ct --clean_temporary_files

    remove all temporary files
  - -sitr --simple_tree

    simple tree: does not apply treetime for ancestral inference
  - -slt --store_locus_tag

    store locus_tags in a separate file instead of saving locus_tags in gene cluster json for large dataset
  - -rlt --raw_locus_tag

    use raw locus_tag from GenBank instead of strain_ID + locus_tag
  - -otc --optional_table_column

    add customized column in gene cluster json file for visualization
  - -mglo --merged_gain_loss_output

    not split gene presence/absence and gain/loss pattern into separate files for each cluster
  - -ngbk --gbk_present

    use nucleotide/amino acid sequence files (fna/faa) when no genBank files given (this option does not consider annotations)
  - -rt --raxml_max_time

    RAxML tree optimization: maximal running time (minutes, default:30min)
  - -dgl --disable_gain_loss

    disable enable gene gain and loss inference (not recommended)
