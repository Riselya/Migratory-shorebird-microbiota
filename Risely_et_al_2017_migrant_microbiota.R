


##R script for analysis included in the publication "Gut microbiota of long-distance migrant demonstrates resistance against environemtnal microbe incursions"

## Amplicon sequences can be downloaded from NCBI. See publication for details.

##workflows used for analysis:
#https://f1000research.com/articles/5-1492/v1
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

library(phyloseq)

##other packages you may need, althoguh I think a few are included with the phyloseq package so you may not need to load them

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

source("miseqR.R")

theme_set(theme_bw())


################# USING ONLY OTUS WITH OVER 10 SEQUENCES (AND WITH PHYLO TREE) ##################################################################

##data already trimmed to those over 10 sequences within Mothur to create this OTU table


# Assign variables for imported data

sharedfile = "phylo.shared"
taxfile = "phylo.taxonomy"

##import mothur shared and tax files

mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

##import metadata
map <- read.csv("clean.metadata.csv", header=T, row.names=1)

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

tree<-read.tree("Riseley.otus.pick.rerooted.tre")

##make tree into phyloseq format

tree<-phy_tree(tree)

##merge tree into current phyloseq object so all object now contained in moth_merge

moth_merge <-merge_phyloseq(moth_merge, tree)

moth_merge
sample_data(moth_merge)
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


####get rid of all sigletons

#stint_all1<-prune_taxa(taxa_sums(stint_all1)>1, stint_all1)

#sample_data(stint_all1)


##get rid of 9185 as lowest read count
##get rid of 8688 as mostly unclassified taxa... not sure why

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


##scale abundance to read depth



smin <- min(sample_sums(stint_all1)) ##9809
smean <- mean(sample_sums(stint_all1)) ##37993
smax <- max(sample_sums(stint_all1)) ##85993


stint_all <- stint_all1 %>%
  scale_reads(n=9809) ##scale to smallest library size


#########################################################################################################


##make object for just stint, excluding environmental samples

stint<-stint_all %>%
  subset_samples(type == "stint")

stint<-stint %>%
  subset_samples(remove_replicate == "No")

##prune taxa

stint <- prune_taxa(taxa_sums(stint) > 0, stint)

################################

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

##################################


##just flinders birds

stint_flind<-stint %>%
  subset_samples(site == "Flinders")

stint_flind <- prune_taxa(taxa_sums(stint_flind) > 0, stint_flind)

###just WTP birds

stint_wtp<-stint %>%
  subset_samples(site == "WTP")

stint_wtp <- prune_taxa(taxa_sums(stint_wtp) > 0, stint_wtp)

##just september birds (migrant/non-migrant analysis)

stint_sep<-stint %>%
  subset_samples(month == "September")

stint_sep <- prune_taxa(taxa_sums(stint_sep) > 0, stint_sep)

##make object to look at just recaptures


recaptures<-stint%>%
  subset_samples( Recap == "Yes")
recaptures <- prune_taxa(taxa_sums(recaptures) > 0, recaptures)

#####stint without repeat captures

stint_norecaps<-stint_flind %>%
  subset_samples(recap_delete == "No")

stint_norecaps<- prune_taxa(taxa_sums(stint_norecaps) > 0, stint_norecaps)


##########################################  Prevalence  ##############################################


prev0 = apply(X = otu_table(stint),
              MARGIN = ifelse(taxa_are_rows(stint), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(stint),
                    tax_table(stint))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))

prevdf$Prevalence<-(prevdf$Prevalence/85)*100


head(prevdf)
str(prevdf)
summary(prevdf$Prevalence)
###only keep those phyla which occur in at least 5 samples

prevdf1


# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(stint)
prevalenceThreshold

stint1 = prune_taxa((prev0 > prevalenceThreshold), stint)
stint1

ggplot(prevdf, aes(TotalAbundance, Prevalence)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)

#write.csv(prevdf, "rerun.prevalence.csv")



#####################################make prevalence tables for environmental samples

prev0 = apply(X = otu_table(env),
              MARGIN = ifelse(taxa_are_rows(env), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf_env = data.frame(Prevalence = prev0,
                        TotalAbundance = taxa_sums(env),
                        tax_table(env))


prevdf_env<-prevdf_env[order(-prevdf_env$TotalAbundance),]
head(prevdf_env)

##########just WTP

prevdf_env_wtp$sampletype<-"Env WTP"

prev0 = apply(X = otu_table(env_wtp),
              MARGIN = ifelse(taxa_are_rows(env_wtp), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf_env_wtp = data.frame(Prevalence = prev0,
                            TotalAbundance = taxa_sums(env_wtp),
                            tax_table(env_wtp))

prevdf_env_wtp<-prevdf_env_wtp[order(-prevdf_env_wtp$TotalAbundance),]
head(prevdf_env_wtp)

################just flinders



prev0 = apply(X = otu_table(env_flind),
              MARGIN = ifelse(taxa_are_rows(env_flind), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf_env_flind = data.frame(Prevalence = prev0,
                              TotalAbundance = taxa_sums(env_flind),
                              tax_table(env_flind))

prevdf_env_flind<-prevdf_env_flind[order(-prevdf_env_flind$TotalAbundance),]
head(prevdf_env_flind)


#################################Core microbiome of adults and second years in september


stint_sep_sy<-stint_sep %>%
  subset_samples(age == "Second year")

stint_sep_sy <- prune_taxa(taxa_sums(stint_sep_sy) > 0, stint_sep_sy)

prev0 = apply(X = otu_table(stint_sep_sy),
              MARGIN = ifelse(taxa_are_rows(stint_sep_sy), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prev_sy = data.frame(Prevalence = prev0,
                     TotalAbundance = taxa_sums(stint_sep_sy),
                     tax_table(stint_sep_sy))
prev_sy$sampletype<-"Second year"

##order by OTU prevalence and abundance

prev_sy<-prev_sy[order(-prev_sy$Prevalence,-prev_sy$TotalAbundance),]
head(prev_sy,20)

##Helicobacter  OTU 1
##Cetobacterium OTU 3
##Cetobacterium OTU 12
##Gammaproteobacteria OTU 15
##Fusobacterium OTU 7


####adults

stint_sep_adults<-stint_sep %>%
  subset_samples(age == "Adult")

stint_sep_adults <- prune_taxa(taxa_sums(stint_sep_adults) > 0, stint_sep_adults)

prev0 = apply(X = otu_table(stint_sep_adults),
              MARGIN = ifelse(taxa_are_rows(stint_sep_adults), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prev_adults = data.frame(Prevalence = prev0,
                         TotalAbundance = taxa_sums(stint_sep_adults),
                         tax_table(stint_sep_adults))
prev_adults$sampletype<-"Adults"


prev_adults<-prev_adults[order(-prev_adults$Prevalence,-prev_adults$TotalAbundance),]
head(prev_adults,10)

write.csv(prev_adults, "prev_adults.csv")


################################################################



####################################TAXANOMIC GRAPHS #########################################
###create function to make mean

merge_samples_mean <- function(physeq, group){
  
  
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  
  
  
  merged <- merge_samples(physeq, group)
  
  
  
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  
  
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}

##is there a difference in taxonomic distribution between sites?

taxGlomRank = "Genus"
length(get_taxa_unique(stint, taxonomic.rank = taxGlomRank))


h1 = 0.5
stint_tree = tip_glom(stint, h = h1)

stint_tree1<-merge_samples_mean(stint_tree, group = "site")

sample_data(stint_tree1)$site <- factor(sample_data(stint_tree1)$site)


plot_tree(stint_tree1, nodelabf=nodeplotblank,size = "abundance", color = "site")+
  guides(colour=guide_legend(override.aes=list(size=10)))+scale_size_area(max_size = 15)

plot_tree(stint_tree1, nodelabf=nodeplotblank,size = "abundance", color = "Phylum")+
  guides(colour=guide_legend(override.aes=list(size=10)))+scale_size_area(max_size = 15)



plot_tree(stint_tree, nodelabf=nodeplotblank, size = "abundance", color = "Phylum")+
  guides(colour=guide_legend(override.aes=list(size=10)))+
  scale_size_area()
tax_table(stint_tree)


############################################################### both env and stint (and neg control)



h1 = 0.5
all_tree = tip_glom(stint_all, h = h1)
tree_all1<-merge_samples_mean(all_tree, group = "type")

##tree which differentiates between stint and enviro samples. Black is stint and blue is enviro (average abundances of bacteria in that lineage)
plot_tree(tree_all1, nodelabf=nodeplotblank, color = "type",size = "abundance")+
  guides(colour=guide_legend(override.aes=list(size=10)))+scale_size_area(max_size = 15)

plot_tree(tree_all1, nodelabf=nodeplotblank, color = "Class",size = "abundance")+
  guides(colour=guide_legend(override.aes=list(size=10)))+scale_size_area(max_size = 15)
sample_data(tree_all1)





#############################################stacked barplots #################################################################

stint_phylum <- stint %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

str(stint_phylum)
head(stint_phylum)

unique(stint_phylum$Phylum)

stint_phylum$Phylum<-factor(stint_phylum$Phylum)




summary<- stint_phylum%>%
  group_by(Phylum)%>%
  summarise(Mean=mean(Abundance), Max=max(Abundance), Min=min(Abundance), Median=median(Abundance), Std=sd(Abundance))

summary<-summary[order(-summary$Mean),]

summary$Phylum

sum(summary$Mean)

summary$percent<-(summary$Mean/10200)*100


stint_phylum$Phylum<-factor(stint_phylum$Phylum, level = c("Actinobacteria", "Proteobacteria","Firmicutes",           
                                                           "Bacteria_unclassified", "Tenericutes","Bacteroidetes","Cyanobacteria",        
                                                           "Spirochaetae" ,"Deferribacteres","Fusobacteria"))
##############################################


stint_phylum$Sample<-factor(stint_phylum$Sample, levels = stint_phylum$Sample[order(stint_phylum$month, stint_phylum$age, stint_phylum$Phylum, stint_phylum$Abundance )])

#stint_phylum$Sample<-factor(stint_phylum$Sample, levels = stint_phylum$Sample[order(stint_phylum$weight)])

colors <- c("navy","royalblue1","darkgreen","cadetblue2","blue2","mediumseagreen","deeppink","gold","darkorange","firebrick")


names(colors) <- levels(stint_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

ggplot(stint_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+
  colScale



ggplot(stint_phylum, aes(x = month, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+
  colScale

sample_df<-sample_data(stint)
unique(sample_df$bird_id)


##by site - super similar ##################################################################

stint_phylum <- stint %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum


stint_phylum$Phylum<-factor(stint_phylum$Phylum, level = c("Actinobacteria", "Proteobacteria","Firmicutes",           
                                                           "Bacteria_unclassified", "Tenericutes","Bacteroidetes","Cyanobacteria",        
                                                           "Spirochaetae" ,"Deferribacteres","Fusobacteria"))




colors <- c("navy","royalblue1","darkgreen","cadetblue2","blue2","mediumseagreen","deeppink","gold","darkorange","firebrick")


names(colors) <- levels(stint_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)


ggplot(stint_phylum, aes(x = site, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+
  scale_fill_manual(values = colors)

################################################################################################

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
##############################################


colors <- c("navy","royalblue1","darkgreen","cadetblue2","blue2","mediumseagreen","deeppink","gold","darkorange","firebrick","black","darkcyan","grey","lightpink","springgreen","maroon")


names(colors) <- levels(stint_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

ggplot(stint_phylum, aes(x = type, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+colScale+theme_bw()



###################enviro samples 



env_phylum <- env %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Sample)                                      # Sort data frame alphabetically by phylum



ggplot(env_phylum, aes(x = site, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##mostly proteobacteria and Bacterreoidetes
##actually less diverse at the Phylum level than the stint microbiome

#########################################################recaps


stint_phylum <- recaptures %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

str(stint_phylum)
head(stint_phylum)

unique(stint_phylum$Phylum)

stint_phylum$Phylum<-factor(stint_phylum$Phylum)

stint_phylum$Sample<-factor(stint_phylum$Sample, levels = stint_phylum$Sample[order(stint_phylum$bird_id,stint_phylum$month, stint_phylum$age )])

stint_phylum$Phylum<-factor(stint_phylum$Phylum, level = c("Actinobacteria", "Proteobacteria","Firmicutes",           
                                                           "Bacteria_unclassified", "Tenericutes","Bacteroidetes","Cyanobacteria",        
                                                           "Spirochaetae" ,"Deferribacteres","Fusobacteria"))




colors <- c("navy","royalblue1","darkgreen","cadetblue2","blue2","mediumseagreen","deeppink","gold","darkorange","firebrick")


names(colors) <- levels(stint_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

##plot where recaps are next to each other

ggplot(stint_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per recaptured bird")+
  scale_fill_manual(values = colors)+
facet_wrap(~bird_id, scales = "free")

######################################

stint_class <- recaptures %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by phylum

str(stint_class)
head(stint_class)

unique(stint_class$Class)

stint_class$Class<-factor(stint_class$Class)

stint_class$Sample<-factor(stint_class$Sample, levels = stint_class$Sample[order(stint_class$bird_id, stint_class$age, stint_class$month, stint_class$Class, stint_class$Abundance )])



ggplot(stint_class, aes(x = Sample, y = Abundance, fill = Class)) + 
  
  geom_bar(stat = "identity", position = "fill" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Phyla abundance per sample")+
  facet_wrap(~bird_id, scales = "free")



################################################ ORDINATION ANLYSIS#############################################

##PCA

stint_pcoa <- ordinate(
  physeq = stint, 
  method = "MDS", 
  distance = "unifrac"
)
str(stint_pcoa)

# Plot 
plot_ordination(
  physeq = stint,
  ordination = stint_pcoa,
  color = "site",
  axes = c(1,4),
  
  title = "PCoA of Stint bacterial Communities"
) + 
  
  geom_point(aes(color = site, shape = type), alpha = 0.7, size = 4) +theme_bw()



##########################  NMDS  ###########################################

set.seed(1)
stint_nmds <- ordinate(
  physeq = stint, 
  method = "NMDS", 
  distance = "bray"
)
stint_nmds

cbPalette <- c("red", "darkgreen")
plot_ordination(
  physeq = stint,
  ordination = stint_nmds,
  color = "site")+
  
  
  geom_point(aes(color = site), alpha = 0.7, size = 4)+
  scale_colour_manual(values=cbPalette)+theme_bw()



stint_bray <- phyloseq::distance(stint, method = "bray")
str(stint_bray)

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint))

# Adonis test
adonis(stint_bray ~ site, data = sampledf)


######################################unifrac

set.seed(1)

stint_nmds <- ordinate(
  physeq = stint, 
  method = "NMDS", 
  distance = "unifrac"
)
stint_nmds

plot_ordination(
  physeq = stint,
  ordination = stint_nmds,
  color = "site")+
  
  
  geom_point(aes(color = site), alpha = 0.7, size = 4)+
  scale_colour_manual(values=cbPalette)+theme_bw()



stint_unifrac <- phyloseq::distance(stint, method = "unifrac")
str(stint_unifrac)

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint))

# Adonis test
adonis(stint_unifrac ~ site, data = sampledf)


####################################all samples

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



# Calculate bray curtis distance matrix
stint_bray <- phyloseq::distance(stint, method = "bray")
str(stint_bray)

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint))

# Adonis test
adonis(stint_bray ~ site, data = sampledf)

#This output tells us that our adonis test is significant so we can reject the null hypothesis 
#that our three sites have the same centroid.

# Homogeneity of dispersion test
beta <- betadisper(stint_bray, sampledf$site)
permutest(beta)


##########################################just septermber to look at differences between adults and second years

set.seed(1)

stint_sep_nmds <- ordinate(
  physeq = stint_sep, 
  method = "NMDS", 
  distance = "bray"
)
stint_sep_nmds

plot_ordination(
  physeq = stint_sep,
  ordination = stint_sep_nmds,
  color = "age",
  
  title = "NMDS of Stint bacterial Communities caught at Flinders in September"
) + 
  
  geom_point(aes(color = age), alpha = 0.7, size = 4) +
  theme_bw()+  scale_colour_manual(values=cbPalette)+theme_bw()



# Calculate bray curtis distance matrix
stint_sep_bray <- phyloseq::distance(stint_sep, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint_sep))

# Adonis test
adonis(stint_sep_bray ~ age, data = sampledf)

#################################using unifrac, which takes into account phylogeny

set.seed(1)

stint_sep_nmds <- ordinate(
  physeq = stint_sep, 
  method = "NMDS", 
  distance = "unifrac"
)
stint_sep_nmds

plot_ordination(
  physeq = stint_sep,
  ordination = stint_sep_nmds,
  color = "age",
  
  title = "NMDS of Stint bacterial Communities caught at Flinders in September"
) + 
  
  geom_point(aes(color = age), alpha = 0.7, size = 4) +
  theme_bw()+  scale_colour_manual(values=cbPalette)



# Calculate unifrac distance matrix
stint_sep_unifrac <- phyloseq::distance(stint_sep, method = "unifrac")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint_sep))

# Adonis test
adonis(stint_sep_unifrac ~ age, data = sampledf)

#This output tells us that our adonis test is significant so we can reject the null hypothesis 
#that our three sites have the same centroid.




##########################################jaccard

set.seed(1)

stint_sep_nmds <- ordinate(
  physeq = stint_sep, 
  method = "NMDS", 
  distance = "jaccard"
)
stint_sep_nmds

plot_ordination(
  physeq = stint_sep,
  ordination = stint_sep_nmds,
  color = "age",
  
  title = "NMDS of Stint bacterial Communities caught at Flinders in September"
) + 
  
  geom_point(aes(color = age), alpha = 0.7, size = 4) +
  theme_bw()


############################################

set.seed(1)

stint_flind_nmds <- ordinate(
  physeq = stint_flind, 
  method = "NMDS", 
  distance = "unifrac"
)
stint_flind

plot_ordination(
  physeq = stint_flind,
  ordination = stint_flind_nmds,
  color = "month",
  
  title = "NMDS of Stint bacterial Communities caught at Flinders in September"
) + 
  
  geom_point(aes(color = month), alpha = 0.7, size = 4) +
  theme_bw()



# Calculate bray curtis distance matrix
stint_flind_unifrac <- phyloseq::distance(stint_flind, method = "unifrac")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint_flind))

# Adonis test
adonis(stint_flind_unifrac ~ month, data = sampledf)

# Calculate bray curtis distance matrix
stint_flind_bray <- phyloseq::distance(stint_flind, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(stint_flind))

# Adonis test
adonis(stint_flind_bray ~ month, data = sampledf)

#############################################################################################


############################################################################


################################################recaps#######################################################



set.seed(1)

stint_flind_nmds <- ordinate(
  physeq = stint_flind, 
  method = "NMDS", 
  distance = "unifrac"
)

axes<-data.frame(stint_flind_nmds$points)


Palette <- c("red","darkgreen","darkgrey", "darkorange")


#fig 4

plot_ordination(
  physeq = stint_flind,
  ordination = stint_flind_nmds,
  color = "age",
  shape = "month")+
  
  geom_point(aes(color = age), alpha = 0.7, size = 4) +
  theme_bw()+scale_colour_manual(values=Palette)+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_path(aes(group = bird_id),col = "darkgrey", arrow=arrow(length=unit(0.5,"cm")))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))


str(stint_flind_nmds)



############################################################################################################

########################################################################################################

####alpha diversity

##use unscaled data


plot_richness(stint_all, color = "site", shape = "type")
plot_richness(stint_flind, color = "month", shape = "type")
plot_richness(replicates, color = "site")


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
  r <-  rarefy_even_depth(stint, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
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

alphadiv_rich_shannon <- merge(rich_stats, s, by = "row.names") 
alphadiv_even_inverse <- merge(even_stats, s, by = "row.names") 

ggplot(alphadiv_rich_shannon, aes(x = site, y = mean))+geom_boxplot()+geom_point(aes(col = month))

t.test(alphadiv_rich_shannon$mean ~ alphadiv_rich_shannon$site)

alpha_flind<-subset(alphadiv_rich_shannon, site =="Flinders")
alpha_wtp<-subset(alphadiv_rich_shannon, site =="WTP")

sd(alpha_flind$mean)
sd(alpha_wtp$mean)

model<-lm(mean ~ age+site+month, alphadiv_rich_shannon)

summary(model)

#alphadiv_rich_shannon<-subset(alphadiv_rich_shannon, mean < 5)

ggplot(alphadiv_even_inverse, aes(x = age, y = mean)) +
  geom_boxplot(width = 0.5) +
  geom_point()+
  theme_bw()

ggplot(alphadiv_rich_shannon, aes(x = age, y = mean)) +
  geom_boxplot(width = 0.5) +
  geom_point()+
  theme_bw()

ggplot(alphadiv_rich_shannon, aes(x = site, y = mean)) +
  geom_boxplot(width = 0.5) +
  geom_point()+
  theme_bw()
#+ylim(0.5,5.7)

ggplot(alphadiv_rich_shannon, aes(x = site, y = mean)) +
  geom_jitter(aes(col = type), size = 2, width = 0.1)+
  theme_bw()


ggplot(rich_flind, aes(month, weight))+ geom_boxplot()+geom_jitter(aes(col = age), size = 2, width = 0.2)
ggplot(rich_flind, aes(month, mean))+ geom_boxplot()+geom_jitter(aes(col = age), size = 2, width = 0.2)



##what about adults and second years in Sep?

alphadiv_rich_sep<-subset(alphadiv_rich_shannon, month=="September")

ggplot(alphadiv_rich_sep, aes(x = age, y = mean)) +
  geom_boxplot(width = 0.5) +
  ylab("Observed richness")+xlab("Age")+
  geom_point()+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

t.test(alphadiv_rich_sep$mean ~ alphadiv_rich_sep$age)

alphadiv_rich_sep.adults<-subset(alphadiv_rich_sep, age == "Adult")
alphadiv_rich_sep.sy<-subset(alphadiv_rich_sep, age == "Second year")
sd(alphadiv_rich_sep.adults$mean)
sd(alphadiv_rich_sep.sy$mean)

##but not higher diversity between adults and second years in september



#What about over time at Flinders?

Palette <- c("red","darkgreen","darkgrey", "darkorange")

alphadiv_rich_flind<-subset(alphadiv_rich_shannon, site == "Flinders")

ggplot(alphadiv_rich_flind, aes(x = month, y = mean)) +
  geom_boxplot(width = 0.5) +
  ylab("Observed richness")+xlab("")+
  geom_jitter(aes(col=age),width = 0.05, size = 2)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(aes(group = bird_id, col = age), size = 0.5, linetype = "longdash")+theme(axis.text=element_text(size=14), axis.title=element_text(size=16))+
  scale_color_manual(values = Palette)


library(nlme)



model<-lme(mean ~ age + month, random = ~ 1|bird_id , data = alphadiv_rich_flind)
summary(model)

##mass 

mass.sep<-data.frame(sample_data(stint_sep))

t.test(mass.sep$weight~ mass.sep$age)


################################################


#############################is this reflected in env samples?

plot_richness(env, color = "site")
env_flind
env_wtp
min_lib <- min(sample_sums(env))


# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(env)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(env)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(env)

set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(env, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
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


s <- data.frame(sample_data(env))


alphadiv_env_rich <- merge(rich_stats, s, by = "row.names") 
alphadiv_env_even <- merge(even_stats, s, by = "row.names") 

###plot richness

ggplot(alphadiv_env_rich, aes(x = site, y = mean)) +
  geom_boxplot() 
#+ylim(0.5,5.7)

##plot evenness

ggplot(alphadiv_env_even, aes(x = site, y = mean, color = site)) +
  geom_boxplot() 



##############################fig1 


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

alphadiv_rich_observed <- merge(rich_stats, s, by = "row.names") 
alphadiv_even_inverse <- merge(even_stats, s, by = "row.names") 

cbPalette <- c("black", "darkgrey")

ggplot(alphadiv_rich_observed, aes(x = site, y = mean))+geom_jitter(aes(shape = type, color = site), width = 0.3, size = 2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))+scale_color_manual(values = cbPalette)+
  ylab("Observed richness")+xlab("")+theme(legend.position="none")

t.test(alphadiv_rich_shannon$mean ~ alphadiv_rich_shannon$site)




##################################################OLIGOTYPING CARYNEBACTERIUM


sep.oligo<-read.csv("oligotyping.csv")


sep.oligo<-subset(oligotyping, month == "September")
sep.oligo<-subset(oligotyping, site == "Flinders")
sep.oligo$OTU_other<-sep.oligo$OTU_4+sep.oligo$OTU_5+sep.oligo$OTU_6+sep.oligo$OTU_7+sep.oligo$OTU_9+sep.oligo$OTU_11

sep.oligo<-sep.oligo[,c(2,4,6,7,8,9,27,28,29,42)]

sep.oligo1<-gather(sep.oligo, OTU, proportion, OTU_1:OTU_other, factor_key=TRUE)
str(sep.oligo1)
sep.oligo1$Sample<-as.factor(sep.oligo1$Sample)
sep.oligo1$proportion1<-sep.oligo1$proportion*sep.oligo1$Abundance

sep.oligo1$month<-factor(sep.oligo1$month, level=c("September","January","March"),ordered = T)


sep.oligo1$Sample<-factor(sep.oligo1$Sample, levels = sep.oligo1$Sample[order(sep.oligo1$month, sep.oligo1$age, -sep.oligo1$Abundance)])

ggplot(sep.oligo1, aes(x = Sample, y = proportion1, fill= OTU))+geom_col() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  facet_wrap(~month, scales="free")


###############################################################################

##SourceTracker


# load sample metadata


###below using data sets where all OTUs which were found in the negative control are excluded, because these made up to 3% of OTU sharing

##run pairwise analyses 


##using WTP as a source and then comparing birds from each site to WTP
metadata<- read.csv("metadata_wtp_source.csv",  header=T, row.names=1)

##using Flinders as a source and then comparing birds from each site to Flinders
#metadata<- read.csv("metadata_flind_source.csv",  header=T, row.names=1)

##Comparing two groups of birds to each other
#metadata<- read.csv("metadata_bird_bird_wtp.source.csv",  header=T, row.names=1)
#metadata<- read.csv("metadata_bird_bird_flind.source.csv",  header=T, row.names=1)

##compare two env samples to each other
#metadata<- read.csv("metadata_env_env.csv",  header=T, row.names=1)

metadata<-metadata[order(metadata$Group),]

str(metadata)


# load OTU table
# This 'read.table' command is designed for a 
# QIIME-formatted OTU table.
# namely, the first line begins with a '#' sign
# and actually _is_ a comment; the second line
# begins with a '#' sign but is actually the header

##example code
#otus <- read.table('otus.txt',sep='\t', header=T,row.names=1,check=F,skip=1,comment='')
#otus <- t(as.matrix(otus))

otus<- read.table('count_table_filtered.txt', sep = '\t', header = T, row.names = 1)

otus<-otus[,1:3102]

otus<-as.matrix(otus)
str(otus)
head(otus)



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
source('SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix], rarefaction_depth = 1000)

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

sourcetracker1<-data.frame(results$proportions)

summary(sourcetracker1$soil)
summary(sourcetracker1$Flinders)
summary(sourcetracker1$bird)
sourcetracker1<-sourcetracker1[order(-sourcetracker1$soil),]

write.csv(sourcetracker1,"results_env_env_flind.source_filtered.csv")

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie',include.legend = TRUE)

sourcetracker<-data.frame(results$proportions)
summary(sourcetracker$soil)
summary(sourcetracker$neg)
summary(sourcetracker$Unknown)

write.csv(sourcetracker, "sourcetracker_wtp.csv")
