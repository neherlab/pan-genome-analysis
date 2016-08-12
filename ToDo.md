#Roadmap
  x harmonize file names:
    x pickle: .cpk
    x nucleotide: .fna
    x amino acid: .faa
    x aligned nucleotide: *_na.aln
    x aligned amino acid: *_aa.aln
    x tree: *.nwk
    x tree as json: *_tree.json
    x - seq as json: *_seq.json
  x rework folder structure with contained folder for visualization and input
    x adjust js and python
  x rename visualization output
    x dataset/YourSpecies/
                          geneCluster.json
                          genePresence.json
                          strainMetainfo.json
                          coreGenomeTree.json
                          geneCluster/
                              GC_0000001.json
                              ...
  
  x column for gene names
  - clean-up step removing tmp files
  - catch logs and put them into a log folder
  - mapping
  
  - blastn all-against-all for rRNAs and other non-translated features such as tRNAs
  - adjust tree vis for >200 nodes
