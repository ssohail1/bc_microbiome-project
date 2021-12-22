### Updated October 8 2021
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microshades)

# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
ps_100_prepcdat <- prep_mdf(ps_top100_cdat)

# directory for saving plots/graphs /media/mbb/Sidras_Projects/Chan_paper/VersionControlDADACommands_Chan/09132021_version

# Create a color object for the specified data
color_ps_100_cdat <- create_color_dfs(ps_100_prepcdat, group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE, selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))


# Extract
mdf_ps_cdat <- color_ps_100_cdat$mdf
cdf_ps_cdat <- color_ps_100_cdat$cdf

plot1_c <- plot_microshades(mdf_ps_cdat, cdf_ps_cdat, group_label = "Phylum Genus")

plot1_c + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=18)) +
  theme(axis.text.x = element_text(size= 12))



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
ps1cdatp <- merge_samples(ps0cdatp, "Sample_Type")
ps2cdatp <- transform_sample_counts(ps1cdatp, function(x) x / sum(x))
pcdatp <- plot_bar(ps2cdatp, fill="Phylum")
#plot_bar(ps2cdatp, fill="Phylum")
finalplotcdatp <- pcdatp + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

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
ps1cdatg <- merge_samples(ps0cdatg, "Sample_Type")
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


## Perform Neighbor joining
treeNJcdat <- NJ(dmcdat) # Note, tip order != sequence order

## Internal maximum likelihood
fitcdat = pml(treeNJcdat, data=phang.aligncdat)
### negative edges length changed to 0!

## negative edges length changed to 0!
fitGTRcdat <- update(fitcdat, k=4, inv=0.2)
fitGTRcdat <- optim.pml(fitGTRcdat, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))

# Import into phyloseq
mapcdat <- import_qiime_sample_data("/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/UrbaniakMetadata.txt")

## Assigning species
taxa.pluscdat <- addSpecies(taxasilvacdat, "/media/mbb/Sidras_Projects/testfolderurbaniak/silva_species_assignment_v138.1.fa", verbose=TRUE)
## 602 out of 6943 were assigned to the species level.
## Of which 436 had genera consistent with the input table.Warning message:
##  In UseMethod("depth") :
##  no applicable method for 'depth' applied to an object of class "NULL"

colnames(taxa.pluscdat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
unname(taxa.pluscdat)

## Making phyloseq object
pscdatphy <- phyloseq(otu_table(seqtab.nochimcdat, taxa_are_rows=FALSE), 
                      tax_table(taxa.pluscdat),phy_tree(fitGTRcdat$tree))

## Merge PhyloSeq object with map
pscdatphy <- merge_phyloseq(pscdatphy,mapcdat)
pscdatphy

## Currently, the phylogenetic tree is not rooted Though it is not necessary here, 
## you will need to root the tree if you want to calculate any phylogeny based 
## diversity metrics (like Unifrac)

set.seed(711)
phy_tree(pscdatphy) <- root(phy_tree(pscdatphy), sample(taxa_names(pscdatphy), 1), resolve.root = TRUE)
is.rooted(phy_tree(pscdatphy))
## TRUE

urbaniaktree <- phy_tree(pscdatphy)

#plot(phy_tree(ps))


# For UniFrac Plots

## Unweighted UniFrac
UniFrac(pscdatphy) #Weighted default is False
unweightedcdat <- UniFrac(pscdatphy)
## this works - makes a scatterplot
plot.default(unweightedcdat)

## this doesn't work
plot(unweightedcdat)


## Weighted UniFrac
UniFrac(pscdatphy, weighted = TRUE)
weightedcdat <- UniFrac(pscdatphy, weighted = TRUE)
## this works - makes a scatterplot
plot.default(weightedcdat)

## this doesn't work
plot(weightedcdat)





# *** Additional Code ***


#PCoA 
##first distance matrix
#https://www.rdocumentation.org/packages/hopach/versions/2.32.0/topics/distancematrix

##then pcoa
#https://www.rdocumentation.org/packages/ape/versions/5.4-1/topics/pcoa

#distancematrix(seqtab.nochim, euclid, na.rm=TRUE)
#distmat<-dist(seqtab.nochim, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

#UPDATED 10122020
write.table(Transose_ASV, file = "/media/mbb/Sidras_Projects/Hieken_paper/VersionControl_DadaCommandRevisions/Dada_commandversions/09242020_version/ASV_file.txt", quote = FALSE)
#Writing the ASV file to 09242020_version folder from RStudio 

library(phyloseq)
# data(GlobalPatterns)


#increase font for this: 
#plot_bar(ps.top100, x = "env_biome", fill = "Genus") + facet_wrap(~final_dx, scales = "free_x") #978 873
# End of dada2 Tutorial

# Function to plot the length distribution of reads in seqtab.nochim
plotLengthDistro <- function(st) {
  require(ggplot2)
  tot.svs <- table(nchar(colnames(st)))
  tot.reads <- tapply(colSums(st), nchar(colnames(st)), sum)
  df <- data.frame(Length=as.integer(c(names(tot.svs), names(tot.reads))),
                   Count=c(tot.svs, tot.reads),
                   Type=rep(c("SVs", "Reads"), times=c(length(tot.svs), length(tot.reads))))
  p <- ggplot(data=df, aes(x=Length, y=Count, color=Type)) + geom_point() + facet_wrap(~Type, scales="free_y") + theme_bw() + xlab("Amplicon Length")
  p
}
# Use on your sequence table
plotLengthDistro(seqtab.nochim)
plotLengthDistro(seqtab.nochim) + scale_y_log10()

# Based on this, we'll need to find some way to eliminate the reads / SVs that contain fewer than 50 nts



#First go through this set of commands, so to use the correct ps file as input for making
#the phylum, genus, and family plots.
pscdat <- phyloseq(otu_table(seqtab.nochimcdat, taxa_are_rows=FALSE),
                   sample_data(samdfcdat), # originally sample_data
                   tax_table(taxasilvacdat))
pscdat <- prune_samples(sample_names(pscdat) != "Mock", pscdat)
dnacdat <- Biostrings::DNAStringSet(taxa_names(pscdat))
names(dnacdat) <- taxa_names(pscdat)
pscdat <- merge_phyloseq(pscdat, dnacdat)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(pscdat)))
pscdat
ps.propcdat <- transform_sample_counts(pscdat, function(otu) otu/sum(otu))
ord.nmds.braycdat <- ordinate(ps.propcdat, method="NMDS", distance="bray", color="Sample_Type")
top100cdat <- names(sort(taxa_sums(pscdat), decreasing=TRUE))[1:100]
ps.top100cdat <- transform_sample_counts(pscdat, function(OTU) OTU/sum(OTU))
ps.top100cdat <- prune_taxa(top100cdat, ps.top100cdat)

#Get different results when use ps.top100 vs using ps
ps <- tax_glom(ps.top100, "Phylum")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")

#Have not used ps.top100
ps <- tax_glom(ps, "Family")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Family",)

#Get different results when use ps.top100 vs using ps
ps <- tax_glom(ps, "Genus")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Genus")

#updated Oct 26 2020
# BiocManager::install("hopach")

######Updated 12-17-2020##########
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf), # originally sample_data
               tax_table(taxa_silva))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray", color="Sample_Type")
top700 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:700]
ps.top700 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top700 <- prune_taxa(top700, ps.top700)

#Get different results when use ps.top700 vs using ps
ps <- tax_glom(ps.top700, "Phylum")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")

ps <- tax_glom(ps.top700, "Family")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Family")

#Get different results when use ps.top700 vs using ps
ps <- tax_glom(ps, "Genus")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, x="Sample_Type", fill="Genus") + facet_wrap(~Isolation_source, scales="free_x")

##updated 06292021

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf), # originally sample_data
               tax_table(taxa_silva))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray", color="Sample_Type")
top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)

#Get different results when use ps.top100 vs using ps
ps <- tax_glom(ps.top100, "Phylum")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")

#Have not used ps.top100
ps <- tax_glom(ps, "Family")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Family",)

#Get different results when use ps.top100 vs using ps
ps <- tax_glom(ps, "Genus")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Sample_Type")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Genus")
