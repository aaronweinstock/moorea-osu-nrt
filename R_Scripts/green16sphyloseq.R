library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
theme_set(theme_bw())

setwd("/Users/janiskrechter/Desktop/greenjesus/")
getwd()

greenseqtab <- readRDS("/Users/janiskrechter/Desktop/greenjesus/seqtab_final_green.rds")
write.csv(greenseqtab, "/Users/janiskrechter/Desktop/greenjesus/greenseqtab.csv")

greentax <- readRDS("/Users/janiskrechter/Desktop/greenjesus/tax_final_green.rds")
write.csv(greentax, "/Users/janiskrechter/Desktop/greenjesus/greentax.csv")

# First go to your ASV table and transpose it, then load it into phyloseq
OTU_info = greenseqtab[-1,-1]
OTU_infoready <- as.matrix(OTU_info)
rownames(OTU_infoready) <- paste0 ("SAMPLE", 1:nrow(OTU_infoready))
colnames(OTU_infoready) <- paste0 ("ASV", 1:ncol(OTU_infoready))
head(OTU_infoready)
class(OTU_infoready)

otumat = OTU_infoready
head(otumat)

#load in your taxa table, remove seqs or na?"
taxmat_info = greentax[-1,]
taxmat_infoready <- as.matrix(taxmat_info, header = TRUE)
rownames(taxmat_infoready) <- paste0 ("ASV", 1:nrow(taxmat_infoready))
head(taxmat_infoready)
taxmat = taxmat_infoready
taxmat

# Merge
OTU = otu_table(otumat, taxa_are_rows = FALSE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX)

sample_names(physeq)
physeq

# Load in sample data
mappingdata = read.csv("metadata_run1run2_june2.csv", header=FALSE, sep = ",")
sampledata_ready <- as.data.frame(mappingdata)
colnames(sampledata_ready)<-c("SampleID", "CORAL", "Genus", "CT", "POINT_X", 
                              "POINT_Y", "POINT_Z", "SITE", "HABITAT",
                              "SITEID", "ISLANDSIDE", "SEASON", "DHW", 
                              "TIMEPOINT", "ENRICHED")
rownames(sampledata_ready) <- sample_names(physeq)
sampledata = sample_data(sampledata_ready)
tail(sampledata)
class(sampledata)

sample_names()

physeqr = merge_phyloseq(physeq, sampledata)
physeqr


##### multiple sequence alignment

library(ape)
random_tree = rtree(ntaxa(physeqr), rooted=TRUE, tip.label=taxa_names(physeqr))
plot(random_tree)

physeq2 <- merge_phyloseq(random_tree, physeqr)
phy_tree(physeq2)

physeq2

#remove mitochondria and Chloroplasts

physeq_nomito <- subset_taxa(physeq2, (Family!="f__mitochondria"))
physeq_nomito

physeq_nomitochloro <- subset_taxa(physeq_nomito, (Class!="c__Chloroplast"))
physeq_nomitochloro
                                   

#remove singletons #

physeq_nosinglemito <- prune_taxa(taxa_sums(physeq_nomitochloro) > 1, 
                                      physeq_nomitochloro)

physeq_nosinglemito

physeq_nosinglemitoDF = as(tax_table(physeq_nosinglemito), "matrix")


#compare what's left
physeq2
physeq_nosinglemito

##distribution of read counts##
sample_sum_df <- data.frame(sum = sample_sums(physeq_nosinglemito))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

write.csv(sample_sum_df, "readdepth_noprune.csv")

smin <- min(sample_sums(physeq_nosinglemito))
smean <- mean(sample_sums(physeq_nosinglemito))
smax <- max(sample_sums(physeq_nosinglemito))

##which samples have no reads?###
##see data frame from above##

####get rid of samples with reads <3000 and greater that 7500 ###
physeq_nooutlierlow <- prune_samples(sample_sums(physeq_nosinglemito) > 3000, 
                                     physeq_nosinglemito)
physeq_nooutlier <- prune_samples(sample_sums(physeq_nooutlierlow) < 75000, 
                                  physeq_nooutlierlow)

physeq_nosinglemito
physeq_nooutlier

##recheck distribution of reads
##distribution of read counts##
sample_sum_df <- data.frame(sum = sample_sums(physeq_nooutlier))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

smin <- min(sample_sums(physeq_nooutlier))
smean <- mean(sample_sums(physeq_nooutlier))
smax <- max(sample_sums(physeq_nooutlier))

# Rarefy samples to even depth
physeq_rarefied <- rarefy_even_depth(physeq_nooutlier, sample.size = 15000,
                                       rngseed = 711, replace = TRUE, 
                                       trimOTUs = TRUE, verbose = TRUE)
physeq_rarefied






# stacked barplots
#get ready
alpha_barplots <- physeq_rarefied %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum) 

#plot

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


# Plot 1
relativeabundance <- ggplot(alpha_barplots, aes(x = sample_Genus, fill = Phylum)) + 
  facet_grid(HABITAT~.) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    labels = c("ACR", "POC", "POR", "WATER", "SED"), 
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition Corals \n Bacterial Communities by Genus and Habitat")

plot(relativeabundance)




# stacked barplots 2 for endozoicamonas
#get ready
alpha_barplots_endos <- physeq_rarefied %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at Family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Family) 

#plot


# Plot 2
ggplot(alpha_barplots_endos, aes(x = sample_Genus, y = Abundance, fill = Family)) + 
  facet_grid(HABITAT~.) +
  geom_bar(stat = "identity") +
  scale_fill_ordinal()+
  scale_x_discrete(
    labels = c("ACR", "POC", "POR", "WATER", "SED"), 
    drop = FALSE
    ) +
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Family > 5%) \n") +
  ggtitle("Family Composition Corals \n Bacterial Communities by Genus and Habitat")

plot(relativeabundance_family)

#PERMANOVA time?

#NMDS using Bray Curtis

NMDSbray <- ordinate(physeq_rarefied,  method = "NMDS", distance = "bray")
plot_ordination(physeq_rarefied, NMDSbray, color = "Genus",
                shape = "HABITAT") +
  geom_point(size=3)+
  labs(col = "Sample Type", shape = "Habitat")

# PCOA using Wunifrac and Bray Curtis

PCOAwunifrac <- ordinate(physeq_rarefied, method = "PCoA", distance = "wunifrac")
plot_ordination(physeq_rarefied, PCOAbray, color = "Genus",
                shape = "HABITAT") +
  geom_point(size=3)+
  labs(col = "Sample Type", shape = "Habitat")

PCOABray <- ordinate(physeq_rarefied, method = "PCoA", distance = "bray")
plot_ordination(physeq_rarefied, PCOAbray, color = "Genus",
                shape = "HABITAT") +
  geom_point(size=3)+
  labs(col = "Sample Type", shape = "Habitat")

# estimate richness and make a df


richness <- plot_richness(physeq_rarefied, 
              x = "Genus",
              measures = "Shannon")

  


write.csv((estimate_richness(physeq_rarefied, split = TRUE, measures = NULL)), 
          file = "/Users/janiskrechter/Desktop/greenjesus/outputs/estimaterichness.csv")

#add Shannon to mapping data


#heatmap

heatmap <- physeqr_ra %>%
  tax_glom(taxrank = "Phyla") %>%   
  
  plot_heatmap(physeqr_rarephylum2, method = "NMDS", distance = "bray", sample.label = "SampleName")

plot(heatmap(physeqr))


# export to biom for PiCrust

# Extract abundance matrix from the phyloseq object
ASV1 = as(otu_table(physeq_nooutlier), "matrix")
# transpose if necessary
if(taxa_are_rows(physeq_nooutlier)){ASV1 <- t(ASV1)}
# Coerce to data.frame
ASVdf = as.data.frame(ASV1)
write.csv(ASVdf, "ASVdf_forPicrust")

#Extract taxa matrix from the phyloseq object
TAX1 = as(tax_table(physeq_nooutlier), "matrix")
TAXdf = as.data.frame(TAX1)
write.csv(TAXdf, "TAXdf_forPicrust")

library(biomformat)
b <- make_biom(ASVdf, 
               sample_metadata = NULL,
               observation_metadata = TAXdf,
               )



