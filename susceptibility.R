
# R script for analysis included in the publication "Gut microbiota of long-distance migrant demonstrates resistance against environemtnal microbe incursions"

# Amplicon sequences can be downloaded from NCBI. See publication for details.
# I have tried to annoate what I'm doing as much as possible without going into too much detail.

# Using R version 3.4.1

# email Alice Risely riselya@gmail.com for queries

# workflows used for analysis:
# http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# https://f1000research.com/articles/5-1492/v1


# To download phyloseq, run the following two lines of code:

# source('http://bioconductor.org/biocLite.R')
# biocLite('phyloseq')

library(phyloseq)

##other packages you may need, althoguh I think a few are included with the phyloseq package.

library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(ape)
library(gridExtra)
library(ade4)
library(plyr)
library(tidyr)

##some extra functions, by Michelle Berry, downloaded from second workflow link above

source("miseqR.R")

theme_set(theme_bw())

################################# DATA MANIPULATION

## data already trimmed to those over 10 sequences within Mothur to create this OTU table
## data includes bird, env and negative control samples


# Assign variables for imported data

sharedfile = "susceptibility.shared"
taxfile = "susceptibility.taxonomy"

##import mothur shared and tax files

mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

##import metadata
map <- read.csv("susceptibility.metadata.csv", header=T, row.names=1)

##make meta data into phyloseq format

map <- sample_data(map)
str(map)

head(map)
tail(map)

##order some of the metadata variables

map$month<-factor(map$month, level=c("September","December","January","March"),ordered = T)
map$age<-factor(map$age, level=c("Adult","Second year","Unknown","Juvenile"))
map$type<-factor(map$type, level=c("stint","env","neg"))
levels(map$age)
levels(map$month)

##merge the metadata into the phyloseq object

moth_merge <- merge_phyloseq(mothur_data, map)

##import phylogenetic tree

tree<-read.tree("susceptibility.tre")

##make tree into phyloseq format

tree<-phy_tree(tree)

##merge tree into current phyloseq object so all object now contained in moth_merge

moth_merge <-merge_phyloseq(moth_merge, tree)

moth_merge
##everything looks fine and no errors!

##rename columns to relevant taxonomic group

colnames(tax_table(moth_merge))

colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")

stint_all1 <- moth_merge %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )

stint_all1
sample_data(stint_all1)
##a few OTUs removed


###now to look at what is in the negative control and remove those OTUs from the rest of the samples

################################################################look at read library


sample_sum_df <- data.frame(sum = sample_sums(stint_all1))
plot(sample_sum_df$sum)


##9185 lowest read of ~7000


ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  ylab("Frequency")+
  theme(axis.title.y = element_blank())

##############################################################################################################################


##exclude 9185 as lowest read count
##exclude of 8688 as mostly unclassified taxa (can skip this step and exclude later)

stint_all1<-stint_all1 %>%
  subset_samples(Group != "8688")
stint_all1<-stint_all1 %>%
  subset_samples(Group != "9185")


stint_all1<-prune_taxa(taxa_sums(stint_all1)>0, stint_all1)

##exclude OTUs that present in negative control sample at over 5 reads

neg_control<-stint_all1 %>% subset_samples(type == "neg")
neg_control <- prune_taxa(taxa_sums(neg_control) > 5, neg_control)
neg_control

otu_table(neg_control)


##40 OTUs which have over 1 reads

##we should get rid of these

badtaxa<-taxa_names(neg_control)
alltaxa<-taxa_names(stint_all1)
alltaxa1 <- alltaxa[!(alltaxa %in% badtaxa)]

stint_all1 = prune_taxa(alltaxa1, stint_all1)
sample_data(stint_all1)

##get rid of negative control altogether and any with low reads

sample_sums(stint_all1)

stint_all1<-prune_samples(sample_sums(stint_all1)>=60, stint_all1)

sample_sums(stint_all1)

# delete replicates

stint_all1<-stint_all1 %>%
  subset_samples(remove_replicate == "No")


##scale abundance to read depth


smin <- min(sample_sums(stint_all1)) ##9809
smean <- mean(sample_sums(stint_all1)) ##37993
smax <- max(sample_sums(stint_all1)) ##85993


stint_all <- stint_all1 %>%
  scale_reads(n=9809) ##scale to smallest library size

# can also rarefy, makes no difference:
# stint_all <-rarefy_even_depth(stint_all1, sample.size = min(sample_sums(stint_all1)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


################################ ANALYSIS #########################################################################

## TABLE 1

##make object for just for stint, excluding environmental samples

stint<-stint_all %>%
  subset_samples(type == "stint")

##prune taxa

stint <- prune_taxa(taxa_sums(stint) > 0, stint)

############## work out prevalence and abundance of common OTUs across all stint

prev0 = apply(X = otu_table(stint),
              MARGIN = ifelse(taxa_are_rows(stint), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(stint),
                    tax_table(stint))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)] ###only keep those phyla which occur in at least 5 samples
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))

prevdf$Prevalence<-(prevdf$Prevalence/85)*100
prevdf$rel_abund<-(prevdf$TotalAbundance/(sum(prevdf$TotalAbundance)))*100

head(prevdf)

#order by prevalence and then abundance

prevdf<-prevdf[order(-prevdf$Prevalence,-prevdf$rel_abund),]
#reorder columns
prevdf<-prevdf[,c(1,2,9,3,4,5,6,7,8)]

#TABLE 1
head(prevdf,10)

###Fig. S2

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(stint)
prevalenceThreshold

stint1 = prune_taxa((prev0 > prevalenceThreshold), stint)


ggplot(prevdf, aes(TotalAbundance, Prevalence)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)


############################################ TABLE 2

###object for environmental samples
env<-stint_all %>%
  subset_samples(type == "env")
env<-prune_taxa(taxa_sums(env) > 0, env)

##Just flinders env

env_flind<-env%>% subset_samples(site == "Flinders")
env_flind <- prune_taxa(taxa_sums(env_flind) > 0, env_flind)

##just wto env samples
env_wtp<-env%>% subset_samples(site == "WTP")
env_wtp <- prune_taxa(taxa_sums(env_wtp) > 0, env_wtp)


################### WTP most common genera

prev0 = apply(X = otu_table(env_wtp),
              MARGIN = ifelse(taxa_are_rows(env_wtp), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf_env_wtp = data.frame(Prevalence = prev0,
                            TotalAbundance = taxa_sums(env_wtp),
                            tax_table(env_wtp))

prevdf_env_wtp$rel_abund<-(prevdf_env_wtp$TotalAbundance/(sum(prevdf_env_wtp$TotalAbundance)))*100
prevdf_env_wtp<-prevdf_env_wtp[,c(1,2,9,3,4,5,6,7,8)]

prevdf_env_wtp<-prevdf_env_wtp[order(-prevdf_env_wtp$rel_abund),]
head(prevdf_env_wtp)

################just flinders

prev0 = apply(X = otu_table(env_flind),
              MARGIN = ifelse(taxa_are_rows(env_flind), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf_env_flind = data.frame(Prevalence = prev0,
                              TotalAbundance = taxa_sums(env_flind),
                              tax_table(env_flind))


prevdf_env_flind$rel_abund<-(prevdf_env_flind$TotalAbundance/(sum(prevdf_env_flind$TotalAbundance)))*100
prevdf_env_flind<-prevdf_env_flind[,c(1,2,9,3,4,5,6,7,8)]

prevdf_env_flind<-prevdf_env_flind[order(-prevdf_env_flind$rel_abund),]
head(prevdf_env_flind)



##########################################################################

#FIG 1a)

stint_phylum <- stint_all %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


stint_phylum$Phylum<-factor(stint_phylum$Phylum)


unique(stint_phylum$Phylum)


stint_phylum$Phylum<-factor(stint_phylum$Phylum, level = c("Actinobacteria", "Proteobacteria","Firmicutes",           
                                                           "Bacteria_unclassified", "Tenericutes","Bacteroidetes","Cyanobacteria",        
                                                           "Spirochaetae" ,"Deferribacteres","Fusobacteria","Acidobacteria",
                                                           "Chloroflexi","Parcubacteria","Planctomycetes","TM6","Verrucimicrobia"))
####


colors <- c("navy","royalblue1","darkgreen","cadetblue2","blue2","mediumseagreen","deeppink","gold","darkorange","firebrick","black","darkcyan","grey","lightpink","springgreen","maroon")


names(colors) <- levels(stint_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

##Fig 1a)

ggplot(stint_phylum, aes(x = type, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+colScale+theme_bw()

######################################################################

##FIG 1b) ##diversity

min_lib <- min(sample_sums(stint_all))

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(stint_all)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(stint_all)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(stint_all)


###

set.seed(3)

for (i in 1:100) {
  # Subsample
  r <-  rarefy_even_depth(stint_all, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed"))) 
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

s <- data.frame(sample_data(stint_all))

alphadiv_rich <- merge(rich_stats, s, by = "row.names") 
alphadiv_even <- merge(even_stats, s, by = "row.names") 

ggplot(alphadiv_rich, aes(x = site, y = mean))+
  geom_jitter(aes(shape = type, col=site), width=0.2)+
  scale_color_manual(values=c("black","darkgray"))+
  ylab("Number of OTUs")

############################################### Fig 1c)

cbPalette <- c("black", "darkgrey")

set.seed(1)
stint_nmds <- ordinate(
  physeq = stint_all, 
  method = "NMDS", 
  distance = "unifrac"
)
str(stint_nmds)

plot_ordination( physeq = stint_all,
                 ordination = stint_nmds,
                 color = "site")+
  
  geom_point(aes(color = site), alpha = 0.7, size = 2)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_manual(values = cbPalette)+theme(legend.position="none")+ylab("Axis 2")+xlab("Axis 1")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))


#############################################################################################################

################ SOURCETRACKER ANALYSIS

### see http://www.nature.com/nmeth/journal/v8/n9/abs/nmeth.1650.html for methods

# This analysis (Fig2) estimates how much birds from each site source their OTUs from sediment, within a Bayesian framework.
# Note this is a completely seperate analysis, not within Phyloseq. 
# However, I have placed it within this script to be chronologically consistent with paper. 
# This analysis takes a long time (overnight). Therefore can skip.

##Start with using WTP env samples as a source and then comparing birds from each site to WTP:

metadata<- read.csv("metadata_wtp_source.csv",  header=T, row.names=1)
metadata<-metadata[order(metadata$Group),]
str(metadata)

# load filtered OTU table. This OTU table is different to the one in previous analyses, 
# as any OTU which was found in the negative control has been excluded completely.
# It has also been reformatted to QIIME format for Sourcetracker

otus<- read.table('count_table_filtered.txt', sep = '\t', header = T, row.names = 1)
otus<-otus[,1:3102]

#convert to matrix
otus<-as.matrix(otus)

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
##example data
#if(is.element('Description',colnames(metadata))) desc <- metadata$Description
if(is.element('Group',colnames(metadata))) 
  desc <- metadata$Group

# load SourceTracker package
# available frmo https://github.com/danknights/sourcetracker ("example.r")
source('SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix], rarefaction_depth = 1000)

# Estimate source proportions in data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

sourcetracker1<-data.frame(results$proportions)

summary(sourcetracker1$soil)
summary(sourcetracker1$Flinders)
summary(sourcetracker1$bird)
sourcetracker1<-sourcetracker1[order(-sourcetracker1$soil),]

#results
#write.csv(sourcetracker1,"results_wtp_source_.csv")

##REPEAT USING FLINDERS ENV SAMPLES AS A SOURCE:
#metadata<- read.csv("metadata_FLIND_source.csv",  header=T, row.names=1)

## AND GO FROM LINE 443

## For comparing birds between sites (top arrow in Fig 2a), can use either of the following metadata files (results same):

#metadata<- read.csv("metadata_bird_bird_wtp.source.csv",  header=T, row.names=1)
#metadata<- read.csv("metadata_bird_bird_flind.source.csv",  header=T, row.names=1)

##compare two env samples to each other (bottom arrow of Fig. 2a)
#metadata<- read.csv("metadata_env_env.csv",  header=T, row.names=1)

########################################################################################################################


###############FIG 3: comparison between sites

#Fig 3a) 
##########################  NMDS  ###########################################

set.seed(1)
stint_nmds <- ordinate(
  physeq = stint, 
  method = "NMDS", 
  distance = "bray"
)
stint_nmds

plot_ordination(
  physeq = stint,
  ordination = stint_nmds)+
  
  
  geom_point(aes(fill = site), alpha = 1, size = 5, pch=21)+
  scale_shape_manual(values=c(17,16))+theme_bw()+
  scale_fill_manual(values=c("orangered", "cadetblue3"))+ theme(legend.position="none")+
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title=element_text(size=16))

stint_bray <- phyloseq::distance(stint, method = "bray")
stint_unifrac <- phyloseq::distance(stint, method = "unifrac")
str(stint_bray)

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint))

# Adonis test
adonis(stint_bray ~ site, data = sampledf)
adonis(stint_unifrac ~ site, data = sampledf)

################################################################


## Fig 3b) created in Lefse software, not R. 

# Fig 3c) diversity between stint occupying different sites

min_lib <- min(sample_sums(stint))

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(stint)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(stint)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(stint)


###

set.seed(3)

for (i in 1:100) {
  # Subsample
  r <-  stint
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed"))) 
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

s <- data.frame(sample_data(stint))

alphadiv_rich <- merge(rich_stats, s, by = "row.names") 
alphadiv_even <- merge(even_stats, s, by = "row.names") 

###Plot

ggplot(alphadiv_rich, aes(x = site, y = mean))+geom_boxplot()+geom_jitter(width=0.2)+
  ylab("Number of OTUs")

#################### SHANNON INDEX 

set.seed(3)

for (i in 1:100) {
  # Subsample
  r <-  stint
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon"))) 
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

s <- data.frame(sample_data(stint))

alphadiv_rich <- merge(rich_stats, s, by = "row.names") 
alphadiv_even <- merge(even_stats, s, by = "row.names") 

#Plot

ggplot(alphadiv_rich, aes(x = site, y = mean))+geom_boxplot()+geom_jitter(width=0.2)+
  ylab("Shannon Index")

#############################################################################################

#Fig 3d)

ggplot(stint_phylum, aes(x = site, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+
  scale_fill_manual(values = colors)

############################################################################################


#######FIG 4: COMPARISON OF MIGRANT (ADULT) AND RESIDENT (SECOND YEAR) STINT IN SEPTEMBER, AT FLINDERS.

stint_sep<-stint %>%
  subset_samples(month == "September")

stint_sep <- prune_taxa(taxa_sums(stint_sep) > 0, stint_sep)


#### Fig 4a)

set.seed(1)

stint_sep_nmds <- ordinate(
  physeq = stint_sep, 
  method = "NMDS", 
  distance = "bray"
)
stint_sep_nmds

plot_ordination(
  physeq = stint_sep,
  ordination = stint_sep_nmds) + 
  geom_point(aes(fill = age), alpha = 1, size = 5, pch=21)+
  scale_shape_manual(values=c(17,16))+theme_bw()+
  scale_color_manual(values=c("red", "cadetblue3"))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title=element_text(size=16))+
  ggtitle("Ordination of migrant and resident stint, using Bray-Curtis distances")

####using unifrac

set.seed(1)

stint_sep_nmds <- ordinate(
  physeq = stint_sep, 
  method = "NMDS", 
  distance = "unifrac"
)
stint_sep_nmds

plot_ordination(
  physeq = stint_sep,
  ordination = stint_sep_nmds) + 
  geom_point(aes(fill = age), alpha = 1, size = 5, pch=21)+
  scale_shape_manual(values=c(17,16))+theme_bw()+
  scale_color_manual(values=c("red", "cadetblue3"))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title=element_text(size=16))+
  ggtitle("Ordination of migrant and resident stint, using Unifrac distances")

################## stats 

# Calculate bray curtis distance matrix
stint_sep_bray <- phyloseq::distance(stint_sep, method = "bray")
stint_sep_unifrac <- phyloseq::distance(stint_sep, method = "unifrac")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint_sep))

# Adonis test
adonis(stint_sep_bray ~ age, data = sampledf)
adonis(stint_sep_unifrac ~ age, data = sampledf)


#######################################################################

###Fig 4b)

stint_sep_phylum <- stint_sep %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


unique(stint_sep_phylum$Phylum)

stint_sep_phylum$Phylum<-factor(stint_sep_phylum$Phylum)

stint_sep_phylum$Phylum<-factor(stint_sep_phylum$Phylum, level = c("Actinobacteria", "Proteobacteria","Firmicutes",           
                                                           "Bacteria_unclassified", "Tenericutes","Bacteroidetes","Cyanobacteria",        
                                                           "Spirochaetae" ,"Deferribacteres","Fusobacteria"))
##############################################


colors <- c("navy","royalblue1","darkgreen","cadetblue2","blue2","mediumseagreen","deeppink","gold","darkorange","firebrick")


names(colors) <- levels(stint_sep_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

ggplot(stint_sep_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+ colScale+ facet_wrap(~age, scales="free", nrow=1, ncol=2)

####################################################

##Fig 4d)

min_lib <- min(sample_sums(stint_sep))

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(stint_sep)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(stint_sep)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(stint_sep)


###

set.seed(3)

for (i in 1:100) {
  # Subsample
  r <-  stint_sep
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed"))) 
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

s <- data.frame(sample_data(stint_sep))

alphadiv_rich <- merge(rich_stats, s, by = "row.names") 
alphadiv_even <- merge(even_stats, s, by = "row.names") 

###Plot

ggplot(alphadiv_rich, aes(x = age, y = mean))+geom_boxplot()+geom_jitter(width=0.2)+
  ylab("Number of OTUs")+ggtitle("OTU richness for migrants (adults) and residents (second years) in September")


####################################################

## Figure 5

## object for just flinders birds

stint_flind<-stint %>%
  subset_samples(site == "Flinders")

stint_flind <- prune_taxa(taxa_sums(stint_flind) > 0, stint_flind)



set.seed(1)

stint_flind_nmds <- ordinate(
  physeq = stint_flind, 
  method = "NMDS", 
  distance = "unifrac"
)

axes<-data.frame(stint_flind_nmds$points)


Palette <- c("red","cadetblue3","darkgrey", "darkorange")

plot_ordination(
  physeq = stint_flind,
  ordination = stint_flind_nmds)+
  
  geom_point(aes(fill = age, shape = month), alpha = 1, size = 5) +
  theme_bw()+scale_fill_manual(values=Palette)+scale_shape_manual(values=c(21,22,23))+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_path(aes(group = bird_id),col = "darkgrey", arrow=arrow(length=unit(0.5,"cm")))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))+ theme(legend.position="none")

###################################################################################

##Fig 6

Palette <- c("red","cyan3","grey30", "darkorange")

alphadiv_rich<-subset(alphadiv_rich, site == "Flinders")

ggplot(alphadiv_rich, aes(x = month, y = mean)) +
  geom_boxplot(width = 0.4, fill="lightgrey") +
  ylab("Observed richness")+xlab("")+
  geom_jitter(aes(fill=age),width = 0.1, size = 3, alpha = 0.7, pch=21)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(aes(group = bird_id, col=age), size = 0.5, linetype = "longdash")+theme(axis.text=element_text(size=18), axis.title=element_text(size=18))+
  scale_fill_manual(values = Palette)+theme(legend.position="none")+scale_color_manual(values = Palette)

#########mixed model


library(nlme)

model<-lme(mean ~ age + month, random = ~ 1|bird_id , data = alphadiv_rich)
summary(model)

#################################### END ###############################

