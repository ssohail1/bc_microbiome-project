#### Read and format relevant files ####

# Taxa
taxasilvacdat <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChannewtaxafromR12302021_R.txt",header = TRUE)
taxasilvacdat <- data.frame(taxasilvacdat)
taxasilvacdat1 <- taxasilvacdat[,-1]
rownames(taxasilvacdat1) <- taxasilvacdat[,1]
taxasilvacdat1 <- as.matrix(taxasilvacdat1)


# Metadata
samdfcdatNAF <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChanMetadatNAF_12202021_R.txt",header = TRUE)
samdfcdatNAF <- data.frame(samdfcdatNAF)

samdfcdatPBS <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChanMetadatPBS_12202021_R.txt",header = TRUE)
samdfcdatPBS <- data.frame(samdfcdatPBS)

samdfcdatNS <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChanMetadatNS_12202021_R.txt",header = TRUE)
samdfcdatNS <- data.frame(samdfcdatNS)

samdfcdat <- read.table("~/Documents/Chan10082021/ChanMetadatCombined_krusk_R.txt",header = TRUE)
samdfcdat <- data.frame(samdfcdat)


# ASV table
seqtab.nochimcdat <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChanASV12202021upd_R.txt",header = TRUE)
seqtab.nochimcdat <- data.frame(seqtab.nochimcdat)
seqtab.nochimcdat1 <- seqtab.nochimcdat[,-1]
rownames(seqtab.nochimcdat1) <- seqtab.nochimcdat[,1]

#### Load libraries ####
library(phyloseq)
library(ape)
library(dplyr)
library(ggplot2)
library(vegan)
library(stats)
library(reshape2)


#### Remove PBS and NS samples from ASV table for adonis test for NAF samples ####

for (i in 1:length(colnames(seqtab.nochimcdat1))) {
  for (j in 1: length(rownames(samdfcdatPBS))) {
    if (rownames(samdfcdatPBS)[j] == colnames(seqtab.nochimcdat1)[i]) {
      seqtab.nochimcdat1 <- seqtab.nochimcdat1[,-i]
    }
  }
}
for (i in 1:length(colnames(seqtab.nochimcdat1))) {
  for (j in 1: length(rownames(samdfcdatNS))) {
    if (rownames(samdfcdatNS)[j] == colnames(seqtab.nochimcdat1)[i]) {
      seqtab.nochimcdat1 <- seqtab.nochimcdat1[,-i]
    }
  }
}

#### Remove PBS and NAF samples from ASV table for adonis test for NS samples ####

for (i in 1:length(colnames(seqtab.nochimcdat1))) {
  for (j in 1: length(rownames(samdfcdatPBS))) {
    if (rownames(samdfcdatPBS)[j] == colnames(seqtab.nochimcdat1)[i]) {
      seqtab.nochimcdat1 <- seqtab.nochimcdat1[,-i]
    }
  }
}
for (i in 1:length(colnames(seqtab.nochimcdat1))) {
  for (j in 1: length(rownames(samdfcdatNAF))) {
    if (rownames(samdfcdatNAF)[j] == colnames(seqtab.nochimcdat1)[i]) {
      seqtab.nochimcdat1 <- seqtab.nochimcdat1[,-i]
    }
  }
}

#### Remove NAF and NS samples from ASV table for adonis test for PBS samples ####

for (i in 1:length(colnames(seqtab.nochimcdat1))) {
  for (j in 1: length(rownames(samdfcdatNAF))) {
    if (rownames(samdfcdatNAF)[j] == colnames(seqtab.nochimcdat1)[i]) {
      seqtab.nochimcdat1 <- seqtab.nochimcdat1[,-i]
    }
  }
}
for (i in 1:length(colnames(seqtab.nochimcdat1))) {
  for (j in 1: length(rownames(samdfcdatNS))) {
    if (rownames(samdfcdatNS)[j] == colnames(seqtab.nochimcdat1)[i]) {
      seqtab.nochimcdat1 <- seqtab.nochimcdat1[,-i]
    }
  }
}

#### Phyloseq Analysis ####
seqtab.nochimcdat1 <- t(seqtab.nochimcdat1)
pscdat <- phyloseq(otu_table(seqtab.nochimcdat1, taxa_are_rows=FALSE),
                   sample_data(samdfcdatPBS), # originally sample_data
                   tax_table(taxasilvacdat1))
# 1997 OTUs removed for PBS ; 2049 OTUs removed for NS ; 1927 OTUs removed for NAF
pscdat <- rarefy_even_depth(pscdat,rngseed = 123456789)
pscdat <- prune_samples(sample_names(pscdat) != "Mock", pscdat)
dnacdat <- Biostrings::DNAStringSet(taxa_names(pscdat))
names(dnacdat) <- taxa_names(pscdat)
pscdat <- merge_phyloseq(pscdat, dnacdat)
taxa_names(pscdat) <- paste0("ASV", seq(ntaxa(pscdat)))
pscdat
pssortcdat <- names(sort(taxa_sums(pscdat), decreasing=TRUE))
ps.pssortcdat <- transform_sample_counts(pscdat, function(OTU) OTU/sum(OTU)) # pscdat input to sample_data command downstream

#### Running Adonis test for PBS, NAF, NS samples ####

set.seed(12345)
ps.pssortcdat <- transform_sample_counts(pscdat, function(OTU) OTU/sum(OTU))
asvtab <- data.frame(ps.pssortcdat@otu_table)

vegdibray <- vegdist(asvtab,method = "bray")
sampledf <- data.frame(sample_data(pscdat)) # making a data frame from the sample_data
adonisChan_health <- adonis2(formula=asvtab ~ health_state, data=sampledf, permutations = 1000, method="bray")


asvta <- pcoa(vegdibray)
asvt <- data.frame(asvta$vectors)

## NAF ##
rownames(asvt) == rownames(samdfcdatNAF)
asvt$meta <- samdfcdatNAF$health_state

ggplot(data = asvt,aes(x=Axis.1 ,y= Axis.2)) +
  geom_point(aes(color = meta)) +
  geom_jitter(aes(color = meta)) +
  labs(title = paste("NAF p-value:",adonisChan_health$`Pr(>F)`, sep = ' '))

## PBS ## 
rownames(asvt) == rownames(samdfcdatPBS)
asvt$meta <- samdfcdatPBS$health_state

ggplot(data = asvt,aes(x=Axis.1 ,y= Axis.2)) +
  geom_point(aes(color = meta)) +
  geom_jitter(aes(color = meta)) +
  labs(title = paste("PBS p-value:",adonisChan_health$`Pr(>F)`, sep = ' '))

## NS ## 
rownames(asvt) == rownames(samdfcdatNS)
asvt$meta <- samdfcdatNS$health_state

ggplot(data = asvt,aes(x=Axis.1 ,y= Axis.2)) +
  geom_point(aes(color = meta)) +
  geom_jitter(aes(color = meta)) +
  labs(title = paste("NS p-value:",adonisChan_health$`Pr(>F)`, sep = ' '))
