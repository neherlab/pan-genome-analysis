# Defects of the input data

The inputs of `panx-data` are the pipeline output of the 2022 run, the 2018 archive of the older datasets, and objects of the old bucket `pangenome.org-data-master` (see `README.md`). This file lists every defect of these inputs that the build or its checks found, and how the build handles it. Paths are relative to the dataset folder.

## Pipeline output

- **`Pseudomonadales`, 179 gene trees**: 279 tips name `NZ_CP020737`, which is no strain of the dataset, in the 2018 archive and in the old bucket. The release build prunes these tips and joins the branch lengths around them. A tree with no strain tip left, such as the tree of `GC00000157_9`, whose tips are all `NZ_CP020737` while its cluster has one gene, is left out of the release. `release-check` lists each pruned tree as `removed-tips` and each left-out tree as `removed-tree`
- **`S_pneumoniae616`, `geneCluster/GC00000001_21_aa_aln_reduced.fa`**: the 2018 archive holds the file twice, as a plain member and as a `.gz` member; the plain member holds gzip bytes. The build reads the `.gz` member
- **`Achromobacter_xylosoxidans`, `metainfo.tsv`**: the 2022 pipeline output starts each row with an extra column, the versioned accession (`NZ_LN890477.1`), and names the accession column `accessions`. The old bucket holds a corrected file of 2022-09-09 without the extra column, like the `metainfo.tsv` of 211 other datasets that S3 replaced on that day. The build reads the files of 2022-09-09
- **Files of no cluster in `geneCluster/`**: files whose cluster id no row of `geneCluster.json` names: `Neisseria_gonorrhoeae_ECDC` 660, `S_pneumoniae616` 1,178, `P_marinus_meta` 84,512, `Enterobacteriales` 16,772, `Pseudomonadales` 121,903, `Vibrionales` 31,706. The release build ignores them
- **`P_marinus_meta`, `vis/geneCluster/npy_bk/` and `npy_sort/`**: 10,994 NumPy working files of the pipeline in the 2018 archive, which the old bucket does not serve. The build leaves them out
- **`Neisseria_gonorrhoeae_ECDC`, strains**: the rows of `strainMetainfo.json` and the tips of the strain tree have no `strain` field (all 1,054 strains). The release writes the strain name as `null`
- **Older datasets, `mixed_continuous` metadata columns**: the pairs `[text, log2 number]` of `meta_details` include pairs without a number, such as `["unknown", "unknown"]` in the MIC columns of `Neisseria_gonorrhoeae_ECDC`. The release leaves these pairs out of the column's `numbers`

## Old bucket

The 2022 upload into `pangenome.org-data-master` left files out or broke them. `verify` compares the rebuild of the old layout with the listings of the bucket and writes each difference to `out/deviations.tsv`. The build takes every file from the archives, so a release has the correct files.

- **Files that the bucket lacks**: 16,054 files: `Staphylococcus_pseudintermedius` 8,098, `Bordetella_pertussis` 7,951, `Bacillus_velezensis` 2, `Phaeobacter_inhibens` 2, and `Vibrionales` 1 (`geneCluster/GC00003645_5_aa_aln.fa`). The clusters of these files fail to load on the old site
- **Files that the bucket holds empty or shortened**: 17 objects of 13 clusters. The bucket holds the plain file, the `.gz` copy, or both:

  | Dataset                        | Cluster file in `geneCluster/`               | Object in the bucket              |
  | ------------------------------ | -------------------------------------------- | --------------------------------- |
  | `Bifidobacterium_bifidum`      | `933b5cfc1bc8580a_refined_aa_aln_reduced.fa` | empty                             |
  | `Bordetella_pertussis`         | `70bf7447c8eaf370_refined_na_aln_reduced.fa` | plain and `.gz` shortened         |
  | `Bordetella_pertussis`         | `7a393f704311812d_refined_na_aln_reduced.fa` | empty                             |
  | `Citrobacter_portucalensis`    | `a3eee768a9b3f979_refined_na_aln.fa`         | empty, `.gz` shortened            |
  | `Faecalibacterium_prausnitzii` | `87ee552fad6a0db4_refined_patterns.json`     | `.gz` shortened                   |
  | `Flavobacterium_psychrophilum` | `d85b431834601aed_refined_na_aln_reduced.fa` | `.gz` shortened                   |
  | `Glaesserella_parasuis`        | `dea6af21e5e5dc62_refined_na_aln.fa`         | `.gz` shortened                   |
  | `Lactobacillus_johnsonii`      | `e11c2b7561fc0764_refined.nwk`               | empty, `.gz` shortened            |
  | `Mycobacterium_intracellulare` | `6d682de0098b0ee9_refined_na_aln_reduced.fa` | `.gz` shortened                   |
  | `P_marinus_meta`               | `GC00003060_r1_1_tree.json`                  | `.gz` shortened (46 of 217 bytes) |
  | `S_pneumoniae616`              | `GC00000001_21_aa_aln_reduced.fa`            | shortened: holds gzip bytes       |
  | `Salmonella_bongori`           | `52638a9ac3895224_refined_na_aln.fa`         | empty, `.gz` shortened            |
  | `Xanthomonas_oryzae`           | `d8998e8bfca5c2d6_refined_aa_aln_reduced.fa` | `.gz` shortened                   |

- **`S_pneumoniae616`, `geneCluster/GC00000034_2_na_aln_reduced.fa`**: the bucket holds a version newer than the 2018 archive with 4 records of genomes of other species (`NZ_CP029093`, `NZ_CP027169`, `NZ_CP020560`, `NC_009656`) and none of the 17 records of the cluster's full alignment. The build takes the version of the 2018 archive, the only exception to the rule that a newer file of the bucket replaces the archive's
- **Stale alignment archives**: 424 old `*_gene_alignments.tar.gz` objects that no page links. The build does not read them

## Findings of `release-check`

`release-check` compares a release with the inputs. For the release `20261111T161409Z` of all 270 datasets, it found:

- **`Pseudomonadales`, gene trees with tips of no strain**: 179 clusters whose trees the build pruned (`removed-tips`, 1 to 6 tips each), and 122 clusters whose trees the build left out because no strain tip remained (`removed-tree`), such as `GC00013598`. These are the trees of the `NZ_CP020737` tips above
- **`organism` column of `metainfo.tsv`**: in 33 datasets, the file of 2022-09-09 in the old bucket has an `organism` column that `metadata.tsv` of the release does not have, because `meta_display_order` of the dataset's `metaConfiguration.js` does not name it, and the release keeps only the metadata columns that the viewer shows. The datasets: `Acinetobacter_indicus`, `Acinetobacter_seifertii`, `Aeromonas_caviae`, `Agrobacterium_tumefaciens`, `Bacillus_altitudinis`, `Bacillus_safensis`, `Bartonella_krasnovii`, `Bordetella_hinzii`, `Bordetella_parapertussis`, `Citrobacter_portucalensis`, `Clostridium_butyricum`, `Dolosigranulum_pigrum`, `Edwardsiella_piscicida`, `Enterobacter_kobei`, `Enterobacter_roggenkampii`, `Escherichia_fergusonii`, `Faecalibacterium_prausnitzii`, `Histophilus_somni`, `Klebsiella_grimontii`, `Lactococcus_cremoris`, `Leclercia_adecarboxylata`, `Mammaliicoccus_sciuri`, `Metamycoplasma_hominis`, `Pectobacterium_brasiliense`, `Pediococcus_acidilactici`, `Proteus_mirabilis`, `Proteus_terrae`, `Shewanella_algae`, `Staphylococcus_argenteus`, `Streptomyces_clavuligerus`, `Weissella_cibaria`, `Xylella_taiwanensis`, `Yersinia_ruckeri`
- **`S_pneumoniae616`, `geneCluster/GC00000001_21_aa_aln_reduced.fa`**: the plain and `.gz` members of the 2018 archive differ in content (see above); the build uses the `.gz` member
