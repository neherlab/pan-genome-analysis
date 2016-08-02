#Roadmap
  - harmonize file names:
    - pickle: .cpk
    - nucleotide: .fna
    - amino acid: .faa
    - aligned nucleotide: *_na.aln
    - aligned amino acid: *_aa.aln
    - tree: *.nwk
    - tree as json: *_tree.json
    - seq as json: *_seq.json
  - rework folder structure with contained folder for visualization and input
    - adjust js and python
  - clean-up step removing tmp files
  - catch logs and put them into a log folder
  - rename visualization output
    - dataset/YourSpecies/
                          geneCluster.json
                          genePresence.json
                          strainMetainfo.json
                          coreGenomeTree.json
                          geneCluster/
                              GC_0000001.json
                              ...
  - column for gene names
  - blastn all-against-all for rRNAs and other non-translated features such as tRNAs
  - mapping
  - adjust tree vis for >200 nodes
