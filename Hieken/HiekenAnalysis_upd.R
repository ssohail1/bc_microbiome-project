# Updated 09182021
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microshades)

pstop100kdat2 <- merge_samples(ps.top100kdat, "env_biome") # looking at tissues
pstop100kdat2 <- merge_samples(ps.top100kdat, "final_dx") # looking at benign vs malignant
# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
ps_100_prepkdat_2 <- prep_mdf(pstop100kdat2)

# Create a color object for the specified data
color_ps_100kdat_2 <- create_color_dfs(ps_100_prepkdat_2, group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE, selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))


# Extract
mdf_ps_kdat2 <- color_ps_100kdat_2$mdf
cdf_ps_kdat2 <- color_ps_100kdat_2$cdf

plot_1_k2 <- plot_microshades(mdf_ps_kdat2, cdf_ps_kdat2, group_label = "Phylum Genus")

plot_1_k2 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=18)) +
  theme(axis.text.x = element_text(size= 14))


ps_100_prepkdat_ <- prep_mdf(ps_top100kdat_)

# Create a color object for the specified data
color_ps_100kdat_ <- create_color_dfs(ps_100_prepkdat_, group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE, selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))


# Extract
mdf_ps_kdat <- color_ps_100kdat_$mdf
cdf_ps_kdat <- color_ps_100kdat_$cdf

plot_1_k <- plot_microshades(mdf_ps_kdat, cdf_ps_kdat, group_label = "Phylum Genus")

plot_1_k + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=18)) +
  theme(axis.text.x = element_text(size= 14))


# Set the working directory
library(dada2)
pathkdat <- "/media/mbb/Sidras_Projects/Hieken_paper/fastq_files/primer_folder-nogz/primerforwardreads"
list.files(pathkdat)

fnFskdat <- sort(list.files(pathkdat, pattern="_1.fastq_primer.fastq", full.names = TRUE))
samplenameskdat <- sapply(strsplit(basename(fnFskdat), "_"), `[`, 1)

# To get the quality profile of the forward reads:
plotQualityProfile(fnFskdat[1:2])

# Assigning filenames for the filtered fastq.gz files
filtFskdat <- file.path(pathkdat, "filtered", paste0(samplenameskdat, "_F_filt.fastq"))
names(filtFskdat) <- samplenameskdat

outkdat <- filterAndTrim(fnFskdat, filtFskdat,
                     maxN=0, maxEE= 3, truncQ=2, minLen = 50, rm.phix=TRUE,
                     compress=TRUE, multithread=1) # multithread = 2 due to low storage


head(outkdat)

# Learn the Error rates
errFkdat <- learnErrors(filtFskdat, multithread=TRUE)
## 195538305 total bases in 726311 reads from 2 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)

# Applying the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFskdat <- dada(filtFskdat, err=errFkdat, multithread=2)
dadaFskdat[[1]]
## dada-class: object describing DADA2 denoising results
## 77 sequence variants were inferred from 87589 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Constructing sequence table
seqtabkdat <- makeSequenceTable(dadaFskdat)
dim(seqtabkdat)
table(nchar(getSequences(seqtabkdat)))

#From https://benjjneb.github.io/dada2/tutorial.html#construct-sequence-table
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

# Remove chimeras
seqtab.nochimkdat <- removeBimeraDenovo(seqtabkdat, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 2093 bimeras out of 8204 input sequences.
#Get lesser bimeras by increasing minFoldParentOverAbundance
seqtabnochimnewkdat <- t(seqtab.nochimkdat)
write.table(seqtabnochimnewkdat, file = "/media/mbb/Sidras_Projects/Hieken_paper/VersionControl_DadaCommandRevisions/Dada_commandversions/20210823_version/seqtabnochimhieken10082021.txt")

sum(seqtab.nochimkdat)/sum(seqtabkdat)
##[1] 0.9734543  ~ should be close to 1

rowSums(seqtab.nochimkdat)
rowSums(seqtabkdat)

# Track reads thru pipeline
getNkdat <- function(x) sum(getUniques(x))

# According to Dada2 tutorial, make the following modification to the track command: 
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) 
#with getN(dadaFs)
trackkdat <- cbind(outkdat, sapply(dadaFskdat, getNkdat), rowSums(seqtab.nochimkdat))
colnames(trackkdat) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(trackkdat) <- samplenameskdat
head(trackkdat)


# Assign Taxonomy
taxasilvakdat <- assignTaxonomy(seqtab.nochimkdat, "/media/mbb/Sidras_Projects/testfolderurbaniak/silva_nr99_v138.1_train_set.fa", multithread = TRUE)
## Went through!
write.table(taxasilvakdat, file = "/media/mbb/Sidras_Projects/Hieken_paper/VersionControl_DadaCommandRevisions/Dada_commandversions/09242020_version/taxasilvahieken10082021.txt")


taxa.printkdat <- taxasilvakdat
rownames(taxa.printkdat) <- NULL
head(taxa.printkdat)


## taxa <- addSpecies(taxa, "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/silva_species_assignment_v138.fa")

# Did not evaluate accuracy with mock data as no mock samples in our data

# Phyloseq
## Need to download phyloseq first
##if (!requireNamespace("BiocManager", quietly = TRUE))
  ## install.packages("BiocManager")
  ##BiocManager::install("phyloseq")

## Then load library
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")


## Hieken Metadata
samdfkdat <- read.table("~/Desktop/HMetadata_Tissues.txt", header=TRUE)
write.table(samdfkdat, file = "/media/mbb/Sidras_Projects/Hieken_paper/VersionControl_DadaCommandRevisions/Dada_commandversions/20210823_version/samdfhieken.txt")


## Making Phyloseq object
pskdat <- phyloseq(otu_table(seqtab.nochimkdat, taxa_are_rows=FALSE),
               sample_data(samdfkdat), # originally sample_data
               tax_table(taxasilvakdat))
pskdat <- prune_samples(sample_names(pskdat) != "Mock", pskdat)
dnakdat <- Biostrings::DNAStringSet(taxa_names(pskdat))
names(dnakdat) <- taxa_names(pskdat)
pskdat <- merge_phyloseq(pskdat, dnakdat)
taxa_names(pskdat) <- paste0("ASV", seq(ntaxa(pskdat)))
pskdat

## Shannon - Simpson: Alpha Diversity Measure
p1kdat <- plot_richness(pskdat, x = "env_biome", measures=c("Shannon", "Simpson"), color = "final_dx")

## Ordination & Bray NMDS: Beta Diversity Measure
ps.propkdat <- transform_sample_counts(pskdat, function(otu) otu/sum(otu))
ord.nmds.braykdat <- ordinate(ps.propkdat, method="NMDS", distance="bray", color="env_biome")
plot_ordination(ps.propkdat, ord.nmds.braykdat, color="env_biome", title="Bray NMDS")

# ps.top100kdat = top 100 sequences
top100kdat <- names(sort(taxa_sums(pskdat), decreasing=TRUE))[1:100]
ps.top100kdat <- transform_sample_counts(pskdat, function(OTU) OTU/sum(OTU))
ps.top100kdat <- prune_taxa(top100kdat, ps.top100kdat)
plot_bar(ps.top100kdat, x = "env_biome", fill = "Family") + facet_wrap(~final_dx, scales = "free_x")
plot_bar(ps.top100kdat, x = "env_biome", fill = "Phylum") + facet_wrap(~final_dx, scales = "free_x")
#increase font for this: 
plot_bar(ps.top100kdat, x = "env_biome", fill = "Genus") + facet_wrap(~final_dx, scales = "free_x") #978 873
# End of dada2 Tutorial


# Proportional Abundance Code
pskdat <- phyloseq(otu_table(seqtab.nochimkdat, taxa_are_rows=FALSE),
                    sample_data(samdfkdat), # originally sample_data
                    tax_table(taxasilvakdat))
pskdat <- prune_samples(sample_names(pskdat) != "Mock", pskdat)
dnakdat <- Biostrings::DNAStringSet(taxa_names(pskdat))
names(dnakdat) <- taxa_names(pskdat)
pskdat <- merge_phyloseq(pskdat, dnakdat)
taxa_names(pskdat) <- paste0("ASV", seq(ntaxa(pskdat)))
pskdat
ps.propkdat <- transform_sample_counts(pskdat, function(otu) otu/sum(otu))
ord.nmds.braykdat <- ordinate(ps.propkdat, method="NMDS", distance="bray", color="final_dx")
top100kdat <- names(sort(taxa_sums(pskdat), decreasing=TRUE))[1:100]
ps.top100kdat <- transform_sample_counts(pskdat, function(OTU) OTU/sum(OTU))
ps.top100kdat <- prune_taxa(top100kdat, ps.top100kdat)

pskdatp <- tax_glom(ps.top100kdat, "Phylum")
ps0kdatp <- transform_sample_counts(pskdatp, function(x) x / sum(x))
ps1kdatp <- merge_samples(ps0kdatp, "final_dx")
ps2kdatp <- transform_sample_counts(ps1kdatp, function(x) x / sum(x))
pkdatp <- plot_bar(ps2kdatp, fill="Phylum")
#point_measkdatp <- pkdatp + geom_point(size = 3)
finalplotkdatp <- pkdatp + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))


pskdatf <- tax_glom(ps.top100kdat, "Family")
ps0kdatf <- transform_sample_counts(pskdatf, function(x) x / sum(x))
ps1kdatf <- merge_samples(ps0kdatf, "final_dx")
ps2kdatf <- transform_sample_counts(ps1kdatf, function(x) x / sum(x))
pkdatf <- plot_bar(ps2kdatf, fill="Family")
finalplotkdatf <- pkdatf + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))

pskdatg <- tax_glom(ps.top100kdat, "Genus")
ps0kdatg <- transform_sample_counts(pskdatg, function(x) x / sum(x))
ps1kdatg <- merge_samples(ps0kdatg, "final_dx")
ps2kdatg <- transform_sample_counts(ps1kdatg, function(x) x / sum(x))
pkdatg <- plot_bar(ps2kdatg, fill="Genus")
finalplotkdatg <- pkdatg + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title.x = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12))


# Phylogenetic Code

# Load libraries
library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)

# Phylogenetic Tree
## Extract sequences from DADA2 output
sequenceskdat <- getSequences(seqtab.nochimkdat)
names(sequenceskdat)<-sequenceskdat

## Run Sequence Alignment (MSA) using DECIPHER
alignmentkdat <- AlignSeqs(DNAStringSet(sequenceskdat), anchor=NA)

## Change sequence alignment output into a phyDat structure
phang.alignkdat <- phyDat(as(alignmentkdat, "matrix"), type="DNA")

## Create distance matrix
dmkdat <- dist.ml(phang.alignkdat)
UPGMAtreekdat <- upgma(dmkdat)
write.tree(UPGMAtreekdat, file = "/media/mbb/Sidras_Projects/Hieken_paper/VersionControl_DadaCommandRevisions/Dada_commandversions/20210918_version/upgmahieken.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)
