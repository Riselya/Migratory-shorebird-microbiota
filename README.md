# Migratory-shorebird-microbiota
Data and analysis to go with our article "Gut microbiota of a long-distance migrant demonstrates resistance against environmental microbe incursions", By Alice Risely, David Waite, Bea Ujvari, Marcel Klaassen and Bethany Hoye. 
Molecular Ecology 2017 http://onlinelibrary.wiley.com/doi/10.1111/mec.14326/full

We do not include sequences here, which can be downloaded from NCBI (SRA accession number SRP106581). This dataset includes only the OTU table generated from sequence analysis within Mothur. 

This article describes changes to gut microbiota dynamics in a population of red-necked stint occupying two sites in Victoria, Australia. The two sites are coded as WTP and Flinders. The sub-population at Flinders were captured over three time points (September, January and March). In addition, other important variables examined are bird age and bird ID (some individuals were recaptured).

Please note that I used the Bioconductor package 'Phyloseq' in this analysis, and this will need to be downloaded from Bioconductor prior to analysis (this is included in the code). In addition, I use SourceTracker package to estimate stint OTU sourcing from environmental samples.

This project contains the files:

susceptibility.R = The R code which includes all analyses presented in the article, in order of presentation (Figs 1 - 6). If you download all files in this repository and run this code, then everything should work. Please note the only analyses not included here are those that were created using LEFSE GALAXY, not R (Fig 3b and 4c).

Source files for OTU analysis:

susceptibility.metadata.csv = this contains the meta-data for each sample analysed

susceptibility.shared = OTU table

susceptibility.taxonomy = Taxonomic information for each OTU

susceptibility.tre = Phylogenetic tree file

miseqR.R = source code downloaded from one of the workflows I followed.

SOURCETRACKER ANLYSIS

Files for SourceTracker analysis (code for this analysis embedded in susceptibility.R script):

Sourcetracker.R = This is the source code for the Bayesian analysis. Please note these analyses take several hours to run.

count_table_filtered.txt = OTU table whereby all OTUs that were present in the negative control were filtered out so that they did not bias analyses (see Methods section in article)

metadata_wtp_source.csv = Metadata in SourceTracker format where the source is the WTP sediment samples

metadata_flind_source.csv = Metadata in SourceTracker format where the source is the Flinders sediment samples



