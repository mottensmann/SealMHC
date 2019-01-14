# `SealMHC`

**Major histocompatibility complex characterisation in Antarctic fur seals** 
_(Master´s project)_

___

### Documentation
All analyses in processing steps are documented in _Rmarkdown_ files that can be found at the top level of this repository

### Dependencies

Analyses make use of both `R` and `Python` scripts.
In order to to repeat the clustering of Illumina MiSeq reads the pipeline [mhc_cluster](https://github.com/mottensmann/mhc_cluster/) needs to be set up as well.

### Raw data

*File sizes of raw sequence data are too large to integrate them within this repository*

#### Fur seal transcriptome assembly
* [Download](http://ramadda.nerc-bas.ac.uk/repository/entry/show/Polar+Data+Centre/NERC-BAS+Datasets/Genomics/Transcriptomes/Arctocephalus_gazella?entryid=synth%3A2d2268fe-907c-45b0-a493-0a6cab8642e6%3AL1RyYW5zY3JpcHRvbWVzL0FyY3RvY2VwaGFsdXNfZ2F6ZWxsYQ%3D%3D) and save as "blast/arc_gaz_transcriptome.fasta"

#### Fur seal genome assembly
* [Download](http://datadryad.org/resource/doi:10.5061/dryad.8kn8c) and save as "blast/arc_gaz_genome.fasta"

* Raw DQB/DRB forward and reverse reads are available upon request
