### Updated October 8 2021
### Updated November 8 2021
### Updated 01262022

### Updated July 17 2021

# Set the working directory
library(dada2)
pathcdat <- "/media/mbb/Sidras_Projects/Chan_paper/primerfastqfiles/"
list.files(pathcdat)
fnFscdat <- sort(list.files(pathcdat, pattern="1.fastq_primer.fastq", full.names = TRUE))
sample.namescdat <- sapply(strsplit(basename(fnFscdat), "_"), `[`, 1)
fnRscdat <- sort(list.files(pathcdat, pattern="2.fastq_primer.fastq", full.names = TRUE))

# To get the quality profile of the forward reads:
plotQualityProfile(fnFscdat[1:2])
plotQualityProfile(fnRscdat[1:2])

# Assigning filenames for the filtered fastq files
filtFscdat <- file.path(pathcdat, "filtered", paste0(sample.namescdat, "_F_filt.fastq"))
names(filtFscdat) <- sample.namescdat
filtRscdat <- file.path(pathcdat, "filtered", paste0(sample.namescdat, "_R_filt.fastq"))
names(filtRscdat) <- sample.namescdat

outcdat <- filterAndTrim(fnFscdat, filtFscdat,
                         maxN=0, maxEE= 3, truncQ = 2, minLen = 50, rm.phix=TRUE,
                         compress=TRUE, multithread=3)

head(outcdat)

# truenaval <- is.na(outcdat) == TRUE
# colnames(truenaval) <- c("readsin", "readsout")
# for (i in 1:length(truenaval)) {
#  # if (truenaval[i] == TRUE) {
#   print(truenaval[i])
#  # }
# }


# Learn the Error rates
errFcdat <- learnErrors(filtFscdat, multithread=TRUE)
## 128716730 total bases in 524831 reads from 2 samples will be used for learning the error rates.
plotErrors(errFcdat, nominalQ=TRUE)

# Dereplication
#derepFscdat <- derepFastq(filtFscdat, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
#names(derepFscdat) <- sample.namescdat
#names(derepRs) <- sample.names

# Applying the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFscdat <- dada(filtFscdat, err=errFcdat, multithread=TRUE)
dadaFscdat[[1]]
# dada-class: object describing DADA2 denoising results
# 38 sequence variants were inferred from 100046 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


# Constructing sequence table
seqtabcdat <- makeSequenceTable(dadaFscdat)
dim(seqtabcdat)
# 137 3743

table(nchar(getSequences(seqtabcdat)))
# 73   74   75   76   77   78   79   80   81   82   86   94   95  100  108  115  116  117  120  137  138  139 
# 2    2    8   13   16    1    2    1    1    5    1    1    1    2    3    2    5    1    2    1    2    1 
# 140  151  177  184  185  187  193  195  206  221  224  227  229  231  232  233  235  237  239  240  241  242 
# 1    1    1    1    2    1    1    2    1    1    2    1    3    1    1    1    2    2    3    2    1   12 
# 243  244  246  247  248  249  250 
# 8    5    3   19    3    9 3580 


#From https://benjjneb.github.io/dada2/tutorial.html#construct-sequence-table
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

# Remove chimeras
seqtab.nochimcdat <- removeBimeraDenovo(seqtabcdat, method="consensus", multithread=TRUE, verbose=TRUE)
# Get lesser bimeras by increasing minFoldParentOverAbundance
# Identified 640 bimeras out of 3743 input sequences.
seqtabnochimnewcdat <- t(seqtab.nochimcdat)
write.table(seqtabnochimnewcdat, file = "/media/mbb/Sidras_Projects/Chan_paper/asvtablechan10082021.txt")  
  
sum(seqtab.nochimcdat)/sum(seqtabcdat)
rowSums(seqtab.nochimcdat)
rowSums(seqtabcdat)

# Track reads thru pipeline
getNcdat <- function(x) sum(getUniques(x))

# According to Dada2 tutorial, make the following modification to the track command: 
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) 
#with getN(dadaFs)
trackcdat <- cbind(outcdat, sapply(dadaFscdat, getNcdat), rowSums(seqtab.nochimcdat))
colnames(trackcdat) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(trackcdat) <- sample.namescdat
head(trackcdat)


# Assign Taxonomy
taxasilvacdat <- assignTaxonomy(seqtab.nochimcdat, "/media/mbb/Sidras_Projects/testfolderurbaniak/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
## Went through!
write.table(taxasilvacdat, file = "/media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/Chantaxafile10022021.txt")
taxa.printcdat <- taxasilvacdat
rownames(taxa.printcdat) <- NULL
head(taxa.printcdat)

# Did not evaluate accuracy with mock data as no mock samples in our data

# Phyloseq
## Need to download phyloseq first
##if (!requireNamespace("BiocManager", quietly = TRUE))
##   install.packages("BiocManager")
##  BiocManager::install("phyloseq")
## Then load library

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

# Steps to create a metadata (samdf file)
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
Sample_Name <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(Sample_Name,1,1)
subject <- substr(Sample_Name,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(sample_name=Sample_Name, Name=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


# Urbaniak Metadata
samdfcdat <- read.table("/media/mbb/Sidras_Projects/Chan_paper/ChanMetadata.txt", header=TRUE)
View(samdfcdat)

# Make phyloseq object
pscdat <- phyloseq(otu_table(seqtab.nochimcdat, taxa_are_rows=FALSE),
                   sample_data(samdfcdat), # originally sample_data
                   tax_table(taxasilvacdat))
pscdat <- prune_samples(sample_names(pscdat) != "Mock", pscdat)
dnacdat <- Biostrings::DNAStringSet(taxa_names(pscdat))
names(dnacdat) <- taxa_names(pscdat)
pscdat <- merge_phyloseq(pscdat, dnacdat)
taxa_names(pscdat) <- paste0("ASV", seq(ntaxa(pscdat)))
pscdat

# Shannon - Simpson: Alpha Diversity
shsicdat <- plot_richness(pscdat, x = "Isolate", measures=c("Shannon", "Simpson"), color = "health_state") + geom_point(size = 2) + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))
#shsicdatpoint <- theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))
#shsicdatfinal <- shsicdatpoint

# Ordination & Bray NMDS: Beta Diversity
ps.propcdat <- transform_sample_counts(pscdat, function(otu) otu/sum(otu))
ord.nmds.braycdat <- ordinate(ps.propcdat, method="NMDS", distance="bray", color="health_state")
braynmdscdat <- plot_ordination(ps.propcdat, ord.nmds.braycdat, color="health_state", title="Bray NMDS") + geom_point(size = 2)

# Abundance Plots
top100cdat <- names(sort(taxa_sums(pscdat), decreasing=TRUE))[1:100]
ps.top100cdat <- transform_sample_counts(pscdat, function(OTU) OTU/sum(OTU))
ps.top100cdat <- prune_taxa(top100cdat, ps.top100cdat)
plot_bar(ps.top100cdat, x = "Sample_Type", fill = "Family") + facet_wrap(~Isolation_source, scales = "free_x")
plot_bar(ps.top100cdat, x = "Sample_Type", fill = "Phylum") + facet_wrap(~Isolation_source, scales = "free_x")


# *** Proportional Abundance Code ***

## Make phyloseq object
pscdat <- phyloseq(otu_table(seqtab.nochimcdat, taxa_are_rows=FALSE),
                   sample_data(samdfcdat), # originally sample_data
                   tax_table(taxasilvacdat))
pscdat <- prune_samples(sample_names(pscdat) != "Mock", pscdat)
dnacdat <- Biostrings::DNAStringSet(taxa_names(pscdat))
names(dnacdat) <- taxa_names(pscdat)
pscdat <- merge_phyloseq(pscdat, dnacdat)
taxa_names(pscdat) <- paste0("ASV", seq(ntaxa(pscdat)))
pscdat

## pscdat file = all sequences
ps.propcdat <- transform_sample_counts(pscdat, function(otu) otu/sum(otu))
ord.nmds.braycdat <- ordinate(ps.propcdat, method="NMDS", distance="bray", color="Sample_Type")
top100cdat <- names(sort(taxa_sums(pscdat), decreasing=TRUE))[1:100]
ps.top100cdat <- transform_sample_counts(pscdat, function(OTU) OTU/sum(OTU))
ps.top100cdat <- prune_taxa(top100cdat, ps.top100cdat)
## ps.top100cdat file = top 100 sequences

## Phylum Level
pscdatp <- tax_glom(ps.top100cdat, "Phylum")
ps0cdatp <- transform_sample_counts(pscdatp, function(x) x / sum(x))
ps1cdatp <- merge_samples(ps0cdatp, "health_state")
ps2cdatp <- transform_sample_counts(ps1cdatp, function(x) x / sum(x))
pcdatp <- plot_bar(ps2cdatp, fill="Phylum")
#plot_bar(ps2cdatp, fill="Phylum")
finalplotcdatp <- pcdatp + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))

## Family Level
pscdatf <- tax_glom(ps.top100cdat, "Family")
ps0cdatf <- transform_sample_counts(pscdatf, function(x) x / sum(x))
ps1cdatf <- merge_samples(ps0cdatf, "Sample_Type")
ps2cdatf <- transform_sample_counts(ps1cdatf, function(x) x / sum(x))
pcdatf <- plot_bar(ps2cdatf, fill="Family")
finalplotcdatf <- pcdatf + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))
#plot_bar(ps2urbggsf, x = "Family", fill = "Sample_Type") + facet_wrap(~Abundance, scales = "free_x")

## Genus Level
pscdatg <- tax_glom(ps.top100cdat, "Genus")
ps0cdatg <- transform_sample_counts(pscdatg, function(x) x / sum(x))
ps1cdatg <- merge_samples(ps0cdatg, "health_state")
ps2cdatg <- transform_sample_counts(ps1cdatg, function(x) x / sum(x))
pcdatg <- plot_bar(ps2cdatg, fill="Genus")
finalplotcdatg <- pcdatg + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

# *** Phylogenetic Tree Code ***

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DECIPHER")

# Load libraries
library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)

## Extract sequences from DADA2 output
sequencescdat <- getSequences(seqtab.nochimcdat)
names(sequencescdat) <- sequencescdat

## Run Sequence Alignment (MSA) using DECIPHER
alignmentcdat <- AlignSeqs(DNAStringSet(sequencescdat), anchor=NA)

## Change sequence alignment output into a phyDat structure
phang.aligncdat <- phyDat(as(alignmentcdat, "matrix"), type="DNA")

## Create distance matrix
dmcdat <- dist.ml(phang.aligncdat)
UPGMAtreecdat <- upgma(dmcdat)
write.tree(UPGMAtreecdat, file = "/media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/09132021_version/upgmachan10202021.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)
                                    

# Microshades Analysis
                                    
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microshades)

# reading in the updated asv table
seqtab2 <- read.table(file = "/media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/11162021_version/chaneditedasvs/asvchannewseqnoch11102021.txt")
seqtab2 <- t(seqtab2)

# reading in the updated taxa file
taxa2 <- read.table(file = "/media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/11162021_version/chaneditedasvs/newtaxachan.txt")

# reading in the metadata tables
metashort <- read.table(file = "/media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/09132021_version/ChanMetadatNS.txt")
metashort <- read.table(file = "/media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/09132021_version/ChanMetadatPBS.txt")
metashort <- read.table(file = "/media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/09132021_version/ChanMetadatNAF.txt")

# creating the phyloseq object and abundance plot at the phylum level
pscdat2 <- phyloseq(otu_table(seqtab2, taxa_are_rows=FALSE),
                   sample_data(metashort), # originally sample_data
                   tax_table(taxa2))
pscdat2 <- prune_samples(sample_names(pscdat2) != "Mock", pscdat2)
dnacdat2 <- Biostrings::DNAStringSet(taxa_names(pscdat2))
names(dnacdat2) <- taxa_names(pscdat2)
pscdat2 <- merge_phyloseq(pscdat2, dnacdat2)
taxa_names(pscdat2) <- paste0("ASV", seq(ntaxa(pscdat2)))
pscdat2
ps.propcdat2 <- transform_sample_counts(pscdat2, function(otu) otu/sum(otu))
ord.nmds.braycdat2 <- ordinate(ps.propcdat2, method="NMDS", distance="bray", color="health_state")
top100cdat2 <- names(sort(taxa_sums(pscdat2), decreasing=TRUE))[1:100]
ps.top100cdat2 <- transform_sample_counts(pscdat2, function(OTU) OTU/sum(OTU))
ps.top100cdat2 <- prune_taxa(top100cdat2, ps.top100cdat2)
pscdatp <- tax_glom(ps.top100cdat2, "Phylum")
ps0cdatp <- transform_sample_counts(pscdatp, function(x) x / sum(x))
ps1cdatp <- merge_samples(ps0cdatp, "health_state")
ps2cdatp <- transform_sample_counts(ps1cdatp, function(x) x / sum(x))
pcdatp <- plot_bar(ps2cdatp, fill="Phylum")
#plot_bar(ps2cdatp, fill="Phylum")
finalplotcdatp <- pcdatp + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))

# Begin Microshades here for Healthy versus Cancer
ps_top100_cdat2 <- merge_samples(ps.top100cdat2,"health_state")

# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
ps_100_prepcdat2 <- prep_mdf(ps_top100_cdat2)

# directory for saving plots/graphs /media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/09132021_version

# Create a color object for the specified data
color_ps_100_cdat2 <- create_color_dfs(ps_100_prepcdat2, group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE, selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))


# Extract
mdf_ps_cdat2 <- color_ps_100_cdat2$mdf
cdf_ps_cdat2 <- color_ps_100_cdat2$cdf

plot_1_c2 <- plot_microshades(mdf_ps_cdat2, cdf_ps_cdat2, group_label = "Phylum Genus")

plot_1_c2 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=18)) +
  theme(axis.text.x = element_text(size= 12))

                                    
# Begin Microshades here for sequences
ps_100_prepcdat <- prep_mdf(ps_top100_cdat)

# directory for saving plots/graphs /media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/09132021_version

# Create a color object for the specified data
color_ps_100_cdat <- create_color_dfs(ps_100_prepcdat, group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE, selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))


# Extract
mdf_ps_cdat <- color_ps_100_cdat$mdf
cdf_ps_cdat <- color_ps_100_cdat$cdf

plot_1_c <- plot_microshades(mdf_ps_cdat, cdf_ps_cdat, group_label = "Phylum Genus")

plot_1_c + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=18)) +
  theme(axis.text.x = element_text(size= 12))

