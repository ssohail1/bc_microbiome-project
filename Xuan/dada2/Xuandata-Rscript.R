library(dada2)
path <- "~/Sidra_Project/Xuanfastq/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
# Forward and reverse fastq filenames
fnFs <- sort(list.files(path, pattern="_1_primer.fastq", full.names = TRUE)) #change pattern to format of fastq files
fnRs <- sort(list.files(path, pattern="_2_primer.fastq", full.names = TRUE)) #change pattern to format of fastq files
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2]) #plotting quality of forward fastq files
plotQualityProfile(fnRs[1:2]) #plotting quality of reverse fastq files
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_primer_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_primer_filt.fastq"))

#minLen parameter added and truncLen parameter removed to retain more reads
#increasing maxEE = allowing for more errors = less stringent
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE, minLen = 50,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Error Models
errF_ <- learnErrors(filtFs, multithread=TRUE) #error model for forward filtered files
errR_ <- learnErrors(filtRs, multithread=TRUE) #error model for reverse filtered files
plotErrors(errF_, nominalQ=TRUE) #plotting the forward error model
plotErrors(errR_, nominalQ=TRUE) #plotting the reverse error model

#Dereplication
derepFs_ <- derepFastq(filtFs, verbose=TRUE)
derepRs_ <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs_) <- sample.names
names(derepRs_) <- sample.names

#Denoising
dadaFs_ <- dada(derepFs_, err=errF_, multithread=TRUE)
dadaRs_ <- dada(derepRs_, err=errR_, multithread=TRUE)
dadaFs_[[1]]

#Merging the forward and reverse reads
mergers_ <- mergePairs(dadaFs_, derepFs_, dadaRs_, derepRs_, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers_[[1]])

#Making ASV table
seqtab_ <- makeSequenceTable(dadaFs_)
dim(seqtab_)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_)))
seqtab.nochim_ <- removeBimeraDenovo(seqtab_, method="consensus", multithread=TRUE, verbose=TRUE) #ASV table with chimeras removed
dim(seqtab.nochim_)
sum(seqtab.nochim_)/sum(seqtab_)
getN <- function(x) sum(getUniques(x))
track_ <- cbind(out, sapply(dadaFs_, getN), sapply(dadaRs_, getN), sapply(mergers_, getN), rowSums(seqtab.nochim_))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs_, getN) with getN(dadaFs_)
colnames(track_) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_) <- sample.names
head(track_)
#Assigning taxonomy
taxa_ <- assignTaxonomy(seqtab.nochim_, "~/Sidra_Project/silva.nr_v138_1.align", multithread=TRUE) #change "~/Sidra_Project/silva.nr_v138_1.align" to where reference file is located

#Optional step
#optional taxa_ <- addSpecies(taxa_, "~/tax/silva_species_assignment_v128.fa.gz")
taxa.print_ <- taxa_ # Removing sequence rownames for display only
rownames(taxa.print_) <- NULL
head(taxa.print_)
head(taxa_)

samdf_ <- read.table("~/Sidra_Project/Xuanfastq/XuanMetadata.txt", header = TRUE) #change directory to where metadata is located


library(phyloseq) #Need to first install Bioconductor
library(Biostrings)
library(ggplot2)

#Have metadata available so don't need the code for making samdf_

#ps_ <- phyloseq(otu_table(seqtab.nochim_, taxa_are_rows=FALSE),
#               sample_data(samdf_), # originally sample_data
#               tax_table(taxa_))
#ps_ <- prune_samples(sample_names(ps_) != "Mock", ps_)
#dna_ <- Biostrings::DNAStringSet(taxa_names(ps_))
#names(dna_) <- taxa_names(ps_)
#ps_ <- merge_phyloseq(ps_, dna_)
#taxa_names(ps_) <- paste0("ASV", seq(ntaxa(ps_)))
#ps_
#plot_richness(ps, x="environment_.feature.", measures=c("Shannon", "Simpson"), color="environment_.material.")
#ps.prop_ <- transform_sample_counts(ps_, function(otu) otu/sum(otu))
#ord.nmds.bray <- ordinate(ps.prop_, method="NMDS", distance="bray", color="Genus")
#plot_ordination(ps.prop, ord.nmds.bray, color="environment._material.", title="Bray NMDS")
#top20_ <- names(sort(taxa_sums(ps_), decreasing=TRUE))[1:20]
#ps.top20_ <- transform_sample_counts(ps_, function(OTU) OTU/sum(OTU))
#ps.top20_ <- prune_taxa(top20_, ps.top20_)
#plot_bar(ps.top20_, x = "environment_.material.", fill = "Family") #+ facet_wrap(~environment_.material., scales = "free_x")
#plot_bar(ps.top20_, x = "environment_.material.", fill = "Phylum") #+ facet_wrap(~environment_.material., scales = "free_x")
#plot_bar(ps.top20_, x = "environment_.material.", fill = "Genus") #+ facet_wrap(~environment_.material., scales = "free_x")

#Making the phyloseq object
ps_ <- phyloseq(otu_table(seqtab.nochim_, taxa_are_rows=FALSE),
               sample_data(samdf_),
               tax_table(taxa_))
ps_ <- prune_samples(sample_names(ps_) != "Mock", ps_) #removing mock groups from ps_
dna_ <- Biostrings::DNAStringSet(taxa_names(ps_))
names(dna_) <- taxa_names(ps_)
ps_ <- merge_phyloseq(ps_, dna_)
taxa_names(ps_) <- paste0("ASV", seq(ntaxa(ps_)))
ps_
ps.prop_ <- transform_sample_counts(ps_, function(otu) otu/sum(otu)) #getting proportional abundances
ord.nmds.bray <- ordinate(ps.prop_, method="NMDS", distance="bray", color="environment_.material.") #Beta diversity = Bray NMDS plot
top20_ <- names(sort(taxa_sums(ps_), decreasing=TRUE))[1:20]
ps.top20_ <- transform_sample_counts(ps_, function(OTU) OTU/sum(OTU))
ps.top20_ <- prune_taxa(top20_, ps.top20_)

#For visualizing alpha diversity = Shannon Simpson plot                                    
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

                                     
#Generating proportional abundance plots at Phylum, Family, and Genus levels
                                     
#Using the ps_ file from above to generate the Phylum, Family, and Genus plots
ps_ <- tax_glom(ps_, "Phylum")
ps0 <- transform_sample_counts(ps_, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "environment_.material.") #environment_.material. from metadata
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")


ps_ <- phyloseq(otu_table(seqtab.nochim_, taxa_are_rows=FALSE),
               sample_data(samdf_), # originally sample_data
               tax_table(taxa_))
ps_ <- prune_samples(sample_names(ps_) != "Mock", ps_)
dna_ <- Biostrings::DNAStringSet(taxa_names(ps_))
names(dna_) <- taxa_names(ps_)
ps_ <- merge_phyloseq(ps_, dna_)
taxa_names(ps_) <- paste0("ASV", seq(ntaxa(ps_)))
ps_
ps.prop_ <- transform_sample_counts(ps_, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop_, method="NMDS", distance="bray", color="environment_.material.")
top20_ <- names(sort(taxa_sums(ps_), decreasing=TRUE))[1:20]
ps.top20_ <- transform_sample_counts(ps_, function(OTU) OTU/sum(OTU))
ps.top20_ <- prune_taxa(top20_, ps.top20_)

#Using the ps_ file from above to generate the Phylum, Family, and Genus plots                                     
ps_ <- tax_glom(ps_, "Family")
ps0 <- transform_sample_counts(ps_, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "environment_.material.")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Family")


ps_ <- phyloseq(otu_table(seqtab.nochim_, taxa_are_rows=FALSE),
               sample_data(samdf_), # originally sample_data
               tax_table(taxa_))
ps_ <- prune_samples(sample_names(ps_) != "Mock", ps_)
dna_ <- Biostrings::DNAStringSet(taxa_names(ps_))
names(dna_) <- taxa_names(ps_)
ps_ <- merge_phyloseq(ps_, dna_)
taxa_names(ps_) <- paste0("ASV", seq(ntaxa(ps_)))
ps_
ps.prop_ <- transform_sample_counts(ps_, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop_, method="NMDS", distance="bray", color="environment_.material.")
top20_ <- names(sort(taxa_sums(ps_), decreasing=TRUE))[1:20]
ps.top20_ <- transform_sample_counts(ps_, function(OTU) OTU/sum(OTU))
ps.top20_ <- prune_taxa(top20_, ps.top20_)

#Using the ps_ file from above to generate the Phylum, Family, and Genus plots
ps_ <- tax_glom(ps_, "Genus")
ps0 <- transform_sample_counts(ps_, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "environment_.material.")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Genus")

