### Update 01262022
 

### Update July 17 2021

# Set the working directory
library(dada2)
pathudat <- "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/"
list.files(pathudat)
fnFsudat <- sort(list.files(pathudat, pattern=".fastq", full.names = TRUE))
sample.namesudat <- sapply(strsplit(basename(fnFsudat), "_"), `[`, 1)

# To get the quality profile of the forward reads:
plotQualityProfile(fnFsudat[1:2])

# Assigning filenames for the filtered fastq files
filtFsudat <- file.path(pathudat, "filtered", paste0(sample.namesudat, "_F_filt.fastq"))
names(filtFsudat) <- sample.namesudat

outudat <- filterAndTrim(fnFsudat, filtFsudat,
                     maxN=0, maxEE= 3, truncQ=2, minLen = 50, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(outudat)

# Learn the Error rates
errFudat <- learnErrors(filtFsudat, multithread=TRUE)
## 63051027 total bases in 839849 reads from 68 samples will be used for learning the error rates.
plotErrors(errFudat, nominalQ=TRUE)

# Dereplication
#derepFsudat <- derepFastq(filtFsudat, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
#names(derepFsudat) <- sample.namesudat
#names(derepRs) <- sample.names

# Applying the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFsudat <- dada(filtFsudat, err=errFudat, multithread=TRUE)
dadaFsudat[[1]]
# dada-class: object describing DADA2 denoising results
# 317 sequence variants were inferred from 3312 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Constructing sequence table
seqtabudat <- makeSequenceTable(dadaFsudat)
dim(seqtabudat)
# 68 7021

table(nchar(getSequences(seqtabudat)))

#From https://benjjneb.github.io/dada2/tutorial.html#construct-sequence-table
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

# Remove chimeras
seqtab.nochimudat <- removeBimeraDenovo(seqtabudat, method="consensus", multithread=TRUE, verbose=TRUE)
# Get lesser bimeras by increasing minFoldParentOverAbundance
# Identified 78 bimeras out of 7021 input sequences.
write.table(seqtab.nochimudat,file = "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/20211112_version/UrbASVfile.txt")
sum(seqtab.nochimudat)/sum(seqtabudat)
# 0.9938704
rowSums(seqtab.nochimudat)
rowSums(seqtabudat)

# Track reads thru pipeline
getNudat <- function(x) sum(getUniques(x))

# According to Dada2 tutorial, make the following modification to the track command: 
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) 
#with getN(dadaFs)
trackudat <- cbind(outudat, sapply(dadaFsudat, getNudat), rowSums(seqtab.nochimudat))
colnames(trackudat) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(trackudat) <- sample.namesudat
head(trackudat)


# Assign Taxonomy
taxasilvaudat <- assignTaxonomy(seqtab.nochimudat, "/media/mbb/Sidras_Projects/testfolderurbaniak/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
## Went through!
write.table(taxasilvaudat, file = "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/20211112_version/Urbtaxafile.txt")
taxa.printudat <- taxasilvaudat
rownames(taxa.printudat) <- NULL
head(taxa.printudat)

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
samdfudat <- read.table("/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/UrbaniakMetadata.txt", header=TRUE)
View(samdfudat)

# Make phyloseq object
psudat <- phyloseq(otu_table(seqtab.nochimudat, taxa_are_rows=FALSE),
               sample_data(samdfudat), # originally sample_data
               tax_table(taxasilvaudat))
psudat <- prune_samples(sample_names(psudat) != "Mock", psudat)
dnaudat <- Biostrings::DNAStringSet(taxa_names(psudat))
names(dnaudat) <- taxa_names(psudat)
psudat <- merge_phyloseq(psudat, dnaudat)
taxa_names(psudat) <- paste0("ASV", seq(ntaxa(psudat)))
psudat

# Shannon - Simpson: Alpha Diversity
shsiudat <- plot_richness(psudat, x = "Isolation_source", measures=c("Shannon", "Simpson"), color = "Sample_Type") + geom_point(size = 2, position = "jitter")
#shsiudatpoint <- theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))
#shsiudatfinal <- shsiudatpoint

# Ordination & Bray NMDS: Beta Diversity
ps.propudat <- transform_sample_counts(psudat, function(otu) otu/sum(otu))
ord.nmds.brayudat <- ordinate(ps.propudat, method="NMDS", distance="bray", color="Sample_Type")
braynmdsudat <- plot_ordination(ps.propudat, ord.nmds.brayudat, color="Sample_Type", title="Bray NMDS") + geom_point(size = 3)

# Abundance Plots
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE))[1:100]
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)
plot_bar(ps.top100udat, x = "Sample_Type", fill = "Family") + facet_wrap(~Isolation_source, scales = "free_x")
plot_bar(ps.top100udat, x = "Sample_Type", fill = "Phylum") + facet_wrap(~Isolation_source, scales = "free_x")

samdfushort <- read.table(file = "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/20211112_version/UrbMetadata11122021.txt",header = TRUE)
# *** Proportional Abundance Code ***

## Make phyloseq object
psudat <- phyloseq(otu_table(seqtab.nochimudat, taxa_are_rows=FALSE),
                   sample_data(samdfushort), # originally sample_data
                   tax_table(taxasilvaudat))
psudat <- prune_samples(sample_names(psudat) != "Mock", psudat)
dnaudat <- Biostrings::DNAStringSet(taxa_names(psudat))
names(dnaudat) <- taxa_names(psudat)
psudat <- merge_phyloseq(psudat, dnaudat)
taxa_names(psudat) <- paste0("ASV", seq(ntaxa(psudat)))
psudat

## psudat file = all sequences
ps.propudat <- transform_sample_counts(psudat, function(otu) otu/sum(otu))
ord.nmds.brayudat <- ordinate(ps.propudat, method="NMDS", distance="bray", color="Type")
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE))[1:100]
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)
## ps.top100udat file = top 100 sequences

## Phylum Level
psudatp <- tax_glom(ps.top100udat, "Phylum")
ps0udatp <- transform_sample_counts(psudatp, function(x) x / sum(x))
ps1udatp <- merge_samples(ps0udatp, "Sample_Type")
ps2udatp <- transform_sample_counts(ps1udatp, function(x) x / sum(x))
pudatp <- plot_bar(ps2udatp, fill="Phylum")
#plot_bar(ps2udatp, fill="Phylum")
finalplotudatp <- pudatp + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

## Family Level
psudatf <- tax_glom(ps.top100udat, "Family")
ps0udatf <- transform_sample_counts(psudatf, function(x) x / sum(x))
ps1udatf <- merge_samples(ps0udatf, "Sample_Type")
ps2udatf <- transform_sample_counts(ps1udatf, function(x) x / sum(x))
pudatf <- plot_bar(ps2udatf, fill="Family")
finalplotudatf <- pudatf + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))
#plot_bar(ps2urbggsf, x = "Family", fill = "Sample_Type") + facet_wrap(~Abundance, scales = "free_x")

## Genus Level
psudatg <- tax_glom(ps.top100udat, "Genus")
ps0udatg <- transform_sample_counts(psudatg, function(x) x / sum(x))
ps1udatg <- merge_samples(ps0udatg, "Library_Name")
ps2udatg <- transform_sample_counts(ps1udatg, function(x) x / sum(x))
pudatg <- plot_bar(ps2udatg,fill="Genus") + facet_wrap(~Type, scales = "free_x")
finalplotudatg <- pudatg + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))

#### *** To create a proportional abundance genus plot separated by tissue type *** ####
# Read in the ASV table
seqtab.nochimudat <- read.table("~/Downloads/Urbaniak_05252022/UrbASVmodifiedSRRs.txt", header = TRUE) # my ASV table has ASVs as first column
# ASV table should have ASVs in column and sample names that start with SRR in row
View(seqtab.nochimudat)
# seqtab.nochimudat modified in the case where first column has the ASVs 
seqtab.nochimudat1 <- seqtab.nochimudat[,-1]
rownames(seqtab.nochimudat1) <- seqtab.nochimudat[,1] # making ASVs the row names and the SRR sample names are already the column names
seqtab.nochimudat1 <- t(data.frame(seqtab.nochimudat1)) # transpose so that the ASVs are the column names and SRR sample names are row names
View(seqtab.nochimudat1)
# Read in the taxa table
taxasilvaudat <- read.table("~/Downloads/Urbaniak_05252022/Urbtaxafile_R.txt",header = TRUE)
taxasilvaudat <- data.frame(taxasilvaudat)
taxasilvaudat1 <- as.matrix(taxasilvaudat)
# Read in the metadata table
samdfudat <- read.table(file= "~/Downloads/Urbaniak_05252022/urbaniakmeta.txt",header = TRUE)
samdfudat <- data.frame(samdfudat)
## Modify metadata table
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,4] == "cancer"){
    samdfudat[j,4] <- "Malignant"
  }
}
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,4] == "benign"){
    samdfudat[j,4] <- "Benign"
  }
}

samdfudat1 <- data.frame(samdfudat)
View(samdfudat1)

library(phyloseq)
library(ggplot2)
samdfudat12 <- samdfudat1[,-1] # first column is SRR sample names
rownames(samdfudat12) <- samdfudat1[,1] # make the SRR sample names the row names
View(samdfudat12)

#### Prop Abund ####
# breast tumor vs healthy
psurbdat <- phyloseq(otu_table(seqtab.nochimudat1, taxa_are_rows=FALSE),
                      sample_data(samdfudat12), # originally sample_data
                      tax_table(taxasilvaudat1))
dnaurbdat <- Biostrings::DNAStringSet(taxa_names(psurbdat))
names(dnaurbdat) <- taxa_names(psurbdat)
psurbdat <- merge_phyloseq(psurbdat, dnaurbdat)
taxa_names(psurbdat) <- paste0("ASV", seq(ntaxa(psurbdat)))
psurbdat
top100urbdat <- names(sort(taxa_sums(psurbdat), decreasing=TRUE)) [1:100] # plotting top 100 seqs
ps.top100urbdat <- transform_sample_counts(psurbdat, function(OTU) OTU/sum(OTU))
ps.top100urbdat <- prune_taxa(top100urbdat, ps.top100urbdat)
# tax_glom to look at Genus level                                
psurbdatg <- tax_glom(ps.top100urbdat, "Genus")
ps0urbdatg <- transform_sample_counts(psurbdatg, function(x) x / sum(x))
ps1urbdatg <- merge_samples(ps0urbdatg, "samp_num") #yields NAs in columns of samdfchdat within ps1urbdatg
ps2urbdatg <- transform_sample_counts(ps1urbdatg, function(x) x / sum(x))
# plot the proportional abundance with large text         
purbdatg <- plot_bar(ps2urbdatg, fill="Genus") 
purbdatg + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))

# benign vs cancer vs healthy facet wrap plot
psurbdat <- phyloseq(otu_table(seqtab.nochimudat1, taxa_are_rows=FALSE),
                      sample_data(samdfudat12), # originally sample_data
                      tax_table(taxasilvaudat1))
dnaurbdat <- Biostrings::DNAStringSet(taxa_names(psurbdat))
names(dnaurbdat) <- taxa_names(psurbdat)
psurbdat <- merge_phyloseq(psurbdat, dnaurbdat)
taxa_names(psurbdat) <- paste0("ASV", seq(ntaxa(psurbdat)))
psurbdat
top100urbdat <- names(sort(taxa_sums(psurbdat), decreasing=TRUE)) [1:100] # plotting top 100 seqs
ps.top100urbdat <- transform_sample_counts(psurbdat, function(OTU) OTU/sum(OTU))
ps.top100urbdat <- prune_taxa(top100urbdat, ps.top100urbdat)
# tax_glom to look at Genus level                         
psurbdatg <- tax_glom(ps.top100urbdat, "Genus")
ps0urbdatg <- transform_sample_counts(psurbdatg, function(x) x / sum(x))
ps1urbdatg <- merge_samples(ps0urbdatg, "samp_num") #yields NAs in columns of samdfchdat within ps1urbdatg
ps2urbdatg <- transform_sample_counts(ps1urbdatg, function(x) x / sum(x))

# ps2urbdatg sample metadata has NA values for entries in the tissue, samp_type, samp_num, and type columns
# loop through the ps2urbdatg sample metadata and loop through samdfudat metadata
# if row name of ps2urbdatg sample metadata (i.e. BTS##) matches row name of samdfudat metadata
# then add the entry in type column of samdfudat metadata to type column of ps2urbdatg sample metadata
for (i in 1:length(rownames(ps2urbdatg@sam_data))) {
  for (j in 1:length(samdfudat1$tissue)) {
    if (rownames(ps2urbdatg@sam_data)[i] == samdfudat1$samp_num[j]) {
      ps2urbdatg@sam_data$type[i] <- samdfudat1$type[j]
    }
  }
}
# facet wrap by type to group and separate samples by cancer, benign, healthy
purbdatg <- plot_bar(ps2urbdatg, fill="Genus") + facet_wrap(~type, scales = "free_x")
purbdatg + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))
                                    
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
sequencesudat <- getSequences(seqtab.nochimudat)
names(sequencesudat) <- sequencesudat

## Run Sequence Alignment (MSA) using DECIPHER
alignmentudat <- AlignSeqs(DNAStringSet(sequencesudat), anchor=NA)

## Change sequence alignment output into a phyDat structure
phang.alignudat <- phyDat(as(alignmentudat, "matrix"), type="DNA")

## Create distance matrix
dmudat <- dist.ml(phang.alignudat)
UPGMAtreeudat <- upgma(dmudat)
write.tree(UPGMAtreeudat, file = "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/20211112_version/upgmaURB11122021.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)

# Microshades Analysis
                                    
# install.packages("remotes")
# library(remotes)
# remotes::install_github("KarstensLab/microshades")
# remotes::install_github("mikemc/speedyseq")

library(phyloseq)
library(dplyr)
library(ggplot2)
library(microshades)

ps_top100_udat <- merge_samples(ps.top100udat, "Type")
# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
ps100_prepudat <- prep_mdf(ps_top100_udat)

# Create a color object for the specified data
color_ps100udat <- create_color_dfs(ps100_prepudat, group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE, selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))


# Extract
mdf_psudat <- color_ps100udat$mdf
cdf_psudat <- color_ps100udat$cdf

plot_1 <- plot_microshades(mdf_psudat, cdf_psudat, group_label = "Phylum Genus")

plot_1 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=15)) +
  theme(axis.text.x = element_text(size= 10))
