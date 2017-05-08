# Migratory-shorebird-microbiota
Data and analysis to go with our article "Gut microbiota of a long-distance migrant demonstrates resistance against environmental microbe incursions", By Alice Risely, David Waite, Bea Ujvari, Marcel Klaassen and Bethany Hoye.

We do not include sequences here, which can be downloaded from NCBI (SRA accession number SRP106581). This dataset includes only the OTU table generated from sequence analysis within Mothur.

This article describes changes to gut microbiota dynamics in a population of red-necked stint occupying two sites in Victoria, Australia. The two sites are coded as WTP and Flinders. The sub-population at Flinders were captured over three time points (September, January and March). In addition, other important variables examined are bird age and bird ID (some individuals were recaptured).

Please note that I used the Bioconductor package 'Phyloseq' in this analysis, and this will need to be downloaded from Bioconductor prior to analysis (this is included in the code, but sometimes can be problematic to download).

This project contains the files:

Risely_et_al_shorebird_microbiota.R = The R code which includes all analyses presented in the article. 

Source files for OTU analysis:

clean.metadata.csv = this contains the meta-data for each sample analysed
phylo.shared = OTU table
phylo.taxonomy = Taxonomic information for each OTU

miseqR.R = source code downloaded from one of the workflows I followed.


Files for SourceTracker analysis:

count_table_filtered.txt = OTU table whereby all OTUs that were present in the negative control were filtered out so that they did not bias analyses (see Methods section in article)

metadata_wtp_source.csv = Metadata in SourceTracker format where the source is the WTP sediment samples
metadata_flind_source.csv = Metadata in SourceTracker format where the source is the Flinders sediment samples
