# Total Sum Scaling

## Hieken
ASVtab <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
ASVtab <- data.frame(ASVtab)
ASVtab2Hiek <- ASVtab[,-1]
# make ASVs the rows
rownames(ASVtab2Hiek) <- ASVtab[,1]
# transpose to have ASVs in the columns
ASVtab2Hiek <- t(ASVtab2Hiek)
seqtab.nochimhkdat12 <- data.frame(ASVtab2Hiek) # original asv table with raw counts
ASVcolsums <- colSums(seqtab.nochimhkdat12)
for (i in 1:length(ASVcolsums)){ 
  for (j in 1:length(seqtab.nochimhkdat12[,1])) { 
    seqtab.nochimhkdat12[j,i] <- round((seqtab.nochimhkdat12[j,i])/(ASVcolsums[i]),digits=10)
  }
}
write.table(seqtab.nochimhkdat12,file = "~/Downloads/HiekenASV12202021upd_percentnormalupd.txt")

## Chan
seqtab.nochimcdat1 <- read.table("~/Downloads/Chan_05242022/ChanASV12202021upd_R.txt",header = TRUE)
seqtab.nochimcdat12 <- data.frame(seqtab.nochimcdat1)
ASVcolsums <- colSums(seqtab.nochimcdat12)
for (i in 1:length(ASVcolsums)){ 
  for (j in 1:length(seqtab.nochimcdat12[,1])) { 
    seqtab.nochimcdat12[j,i] <- round((seqtab.nochimcdat12[j,i])/(ASVcolsums[i]),digits=10)
  }
}
write.table(seqtab.nochimcdat12, "~/Downloads/Chan_05242022/ChanASV-TSS_upd.txt")

## Urbaniak
seqtab.nochimudat <- read.table("~/Downloads/Urbaniak_05252022/UrbASVmodifiedSRRs.txt", header = TRUE)
seqtab.nochimudat1 <- seqtab.nochimudat[,-1]
rownames(seqtab.nochimudat1) <- seqtab.nochimudat[,1]
seqtab.nochimudat1 <- t(data.frame(seqtab.nochimudat1))
seqtab.nochimudat12 <- data.frame(seqtab.nochimudat1)
ASVcolsums <- colSums(seqtab.nochimudat12)
for (i in 1:length(ASVcolsums)){ 
  for (j in 1:length(seqtab.nochimudat12[,1])) { 
    seqtab.nochimudat12[j,i] <- round((seqtab.nochimudat12[j,i])/(ASVcolsums[i]),digits=10)
  }
}
write.table(seqtab.nochimudat12, "~/Downloads/Urbaniak_05252022/UrbASV-TSS_upd.txt")

# Python script for normalizing ASV tables

## Hieken python script
# python percentile_norm.py -i HiekenASV12202021_percentnormal.txt -case Hiekencases.txt -control Hiekencontrols.txt -o Hieken_percentile_norm-proportional.txt

## Chan python script
# python percentile_norm.py -i ChanASVTSS_forpercentnormal.txt -case Chancases.txt -control Chancontrols.txt -o Chan_percentile_norm-proportional.txt

## Urbaniak python script
# python percentile_norm.py -i UrbASVTSS_forpercentnormal.txt -case Urbaniakcases.txt -control Urbaniakcontrols.txt -o Urbaniak_percentile_norm-proportional.txt

### normalized files have the associated taxa of each ASV as column
hiekpernorm <- read.table("~/Downloads/Hieken_percentile_norm-proportional.txt")
urbpernorm <- read.table("~/Downloads/Urbaniak_percentile_norm-proportional.txt")
chanpernorm <- read.table("~/Downloads/Chan_percentile_norm-proportional.txt")

### Before normalization ASV tables with ASVs as columns
casv <- read.table("~/Downloads/Chan_05242022/ChanASV-TSS_upd.txt")
hasv <- read.table(file = "~/Downloads/HiekenASV12202021upd_percentnormalupd.txt")
uasv <- read.table("~/Downloads/Urbaniak_05252022/UrbASV-TSS_upd.txt")
colnames(casv)[1]
colnames(uasv)[1]
colnames(hasv)[1]

### make column names of the pre-normalized ASV tables, the column names of the normalized ASV tables
colnames(chanpernorm) <- colnames(casv)
colnames(urbpernorm) <- colnames(uasv)
colnames(hiekpernorm) <- colnames(hasv)

### transpose so that ASVs are rows and sample names are columns
chanpernorm <- t(chanpernorm)
hiekpernorm <- t(hiekpernorm)
urbpernorm <- t(urbpernorm)

### edit these in text edit to include the word " ASV " before the first sample name
### so can join the files in next step
write.table(chanpernorm,file = "~/Downloads/Chan_05242022/chantsspercnormforR.txt")
write.table(hiekpernorm,file = "~/Downloads/hieken_05242022/hiekentsspercnormforR.txt")
write.table(urbpernorm,file = "~/Downloads/Urbaniak_05252022/urbaniaktsspercnormforR.txt")

# Joining tables

## ASV tables
chanpercnorm <- read.table(file = "~/Downloads/Chan_05242022/chantsspercnormforR.txt", header = TRUE)
hiekpercnorm <- read.table(file = "~/Downloads/hieken_05242022/hiekentsspercnormforR.txt", header = TRUE)
urbpercnorm <- read.table(file = "~/Downloads/Urbaniak_05252022/urbaniaktsspercnormforR.txt", header = TRUE)

HiekUrbASVcombined <- full_join(hiekpercnorm,urbpercnorm,by="ASV")
HiekUrbASVcombined <- data.frame(HiekUrbASVcombined)
HiekUrbChanASVcombined <- full_join(HiekUrbASVcombined,chanpercnorm,by="ASV")
## Modify & format the ASV tables
HiekUrbChanASVcombined1 <- HiekUrbChanASVcombined[,-1]
rownames(HiekUrbChanASVcombined1) <- HiekUrbChanASVcombined[,1]
HiekUrbChanASVcombined1 <- t(HiekUrbChanASVcombined1)
## change the NA values to 0
for (i in 1:length(rownames(HiekUrbChanASVcombined1))){
  ftasv <- is.na(HiekUrbChanASVcombined1[i,])
  for (j in 1:length(ftasv)) {
    if (ftasv[j] == TRUE) {
      HiekUrbChanASVcombined1[i,j] <- 0
    }
  }
}

## Order the ASV table from least to greatest value
newdatafr <- HiekUrbChanASVcombined1[,order(colSums(-HiekUrbChanASVcombined1))]
HiekUrbChanASVcombined1 <- data.frame(newdatafr)

## TAXONOMY tables


## METADATA tables


#### Code #### - need to edit and add comments
#### Read and format input files ####
library(dplyr)
# Urbaniak ASV table
seqtab.nochimudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVmodifiedSRRs.txt", header = TRUE)
seqtab.nochimudat <- data.frame(seqtab.nochimudat)
# Urbaniak Taxa table
taxasilvaudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile_R.txt",header = TRUE)
taxasilvaudat <- data.frame(taxasilvaudat)
# Urbaniak Metadata table
samdfudat <- read.table(file= "~/Downloads/Urbaniak_05252022/urbaniakmeta.txt")
samdfudat <- data.frame(samdfudat)
## Modify Urbaniak metadata table
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,5] == "cancer"){
    samdfudat[j,5] <- "Malignant"
  }
}
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,5] == "benign"){
    samdfudat[j,5] <- "Benign"
  }
}
colnames(samdfudat) <- c("sample_id","sample","tissue","sample_name","sample_type")

# Hieken ASV table
HiekenASVtab <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
HiekenASVtab <- data.frame(HiekenASVtab)
# Hieken Taxa table
Hiekentaxa <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021.txt',header=TRUE)
Hiekentaxa <- data.frame(Hiekentaxa)
colnames(Hiekentaxa)[1] <- "Taxonomy"
# Hieken metadata table
HMeta <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissues.txt")
HMeta <- data.frame(HMeta)
colnames(HMeta) <- c("sample_id","sample","sample_type")

# Chan ASV table
seqtab.nochimcdat <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChanASV12202021upd_R.txt",header = TRUE)
seqtab.nochimcdat <- data.frame(seqtab.nochimcdat)
colnames(seqtab.nochimcdat)[1] <- "NAME"
# Chan taxa table
taxasilvacdat <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChannewtaxafromR12302021_R.txt",header = TRUE)
taxasilvacdat <- data.frame(taxasilvacdat)
colnames(taxasilvacdat)[1] <- "Taxonomy"
# Chan metadata table
samdfcdat <- read.table("~/Documents/Chan10082021/ChanMetadatCombined_.txt",header = TRUE)
samdfcdat <- data.frame(samdfcdat)
colnames(samdfcdat) <- c("sample_id","sample","sample_type")

samdfudat <- samdfudat[-1,]

#### Load libraries ####
library(dplyr)
library(phyloseq)
library(ape)
library(ggplot2)
library(vegan)
# install.packages("ecodist")
library(ecodist)

#### Combine the tables for all three studies ####
# Combine the ASV tables
HiekUrbASVcombined <- full_join(HiekenASVtab,seqtab.nochimudat, by = "NAME")
HiekUrbASVcombined <- data.frame(HiekUrbASVcombined)
HiekUrbChanASVcombined <- full_join(HiekUrbASVcombined,seqtab.nochimcdat,by="NAME")
## Modify & format the ASV tables
HiekUrbChanASVcombined2 <- HiekUrbChanASVcombined[,-1]
rownames(HiekUrbChanASVcombined2) <- HiekUrbChanASVcombined[,1]
HiekUrbChanASVcombined2 <- t(HiekUrbChanASVcombined2)
# Combine the metadata tables
HiekUrbMETADATAcombined <- full_join(HMeta,samdfudat,by = "sample_id")
HiekUrbMETADATAcombined <- data.frame(HiekUrbMETADATAcombined)
HiekUrbChanMETADATAcombined <- full_join(HiekUrbMETADATAcombined,samdfcdat,by="sample_id")
## Modify & format the metadata tables
HiekUrbChanMETADATAcombined1 <- HiekUrbChanMETADATAcombined[,-1]
rownames(HiekUrbChanMETADATAcombined1) <- HiekUrbChanMETADATAcombined[,1]
HiekUrbChanMETADATAcombined1$StudyGroup <- "asv"
for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(HiekUrbChanMETADATAcombined1)[i] == HMeta$sample_id[j]) {
      HiekUrbChanMETADATAcombined1$StudyGroup[i] <- paste("Hiek_", HMeta$sample_type[j], sep = '')
    }
  }
}

for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(HiekUrbChanMETADATAcombined1)[i] == samdfudat$sample_id[j]) {
      HiekUrbChanMETADATAcombined1$StudyGroup[i] <- paste("Urb_", samdfudat$sample_type[j], sep = '')
    }
  }
}

for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(HiekUrbChanMETADATAcombined1)[i] == samdfcdat$sample_id[j]) {
      if (samdfcdat$sample_type[j] == "Healthy_Control_Women") {
        HiekUrbChanMETADATAcombined1$StudyGroup[i] <- "Cha_Healthy"
      }
      if (samdfcdat$sample_type[j] == "Women_with_a_History_of_Breast_Cancer") {
        HiekUrbChanMETADATAcombined1$StudyGroup[i] <- "Cha_Cancer"
      }
    }
  }
}

HiekUrbChanMETADATAcombined1$Tissues <- "asv"
for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(HiekUrbChanMETADATAcombined1)[i] == HMeta$sample_id[j]) {
      HiekUrbChanMETADATAcombined1$Tissues[i] <- paste("Hiek_", HMeta$sample[j], sep = '')
    }
  }
}

for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(HiekUrbChanMETADATAcombined1)[i] == samdfudat$sample_id[j]) {
      HiekUrbChanMETADATAcombined1$Tissues[i] <- paste("Urb_", samdfudat$sample[j], sep = '')
    }
  }
}

for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(HiekUrbChanMETADATAcombined1)[i] == samdfcdat$sample_id[j]) {
        HiekUrbChanMETADATAcombined1$Tissues[i] <- paste("Chan_", samdfcdat$sample[j], sep = '')
    }
  }
}



# Combine the taxa tables
## the combined taxa file is generated from python script editmetadataASVTaxa.py
HiekUrbChanTAXAcombined <- read.table("~/Documents/HiekUrbChan_CombinedTaxa_fin.txt",header = TRUE)
HiekUrbChanTAXAcombined2 <- HiekUrbChanTAXAcombined[,-1]
rownames(HiekUrbChanTAXAcombined2) <- HiekUrbChanTAXAcombined[,1]
HiekUrbChanTAXAcombined2 <- as.matrix(HiekUrbChanTAXAcombined2)

#### Replace the NA values with zeros ####
# Replace NAs with the character zero ("zero") in the metadata table
for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
  ft <- is.na(HiekUrbChanMETADATAcombined1[i,])
  for (j in 1:length(ft)) {
    if (ft[j] == TRUE) {
      HiekUrbChanMETADATAcombined1[i,j] <- "zero"
    }
  }
}
# Replace NAs with the number zero (0) in the ASV table
for (i in 1:length(rownames(HiekUrbChanASVcombined2))){
  ftasv <- is.na(HiekUrbChanASVcombined2[i,])
  for (j in 1:length(ftasv)) {
    if (ftasv[j] == TRUE) {
      HiekUrbChanASVcombined2[i,j] <- 0
    }
  }
}

#### original asv -no modifs beta div ####
D.BCHiek <- as.matrix(vegdist(HiekUrbChanASVcombined2, method="robust.aitchison")) # also have used method="robust.aitchison"
cmdBCurtis <- pcoa(D.BCHiek) #cmdscale(D.BCHiek)
cmdBCurtis <- data.frame(cmdBCurtis[["vectors"]])
cmdBCurtis$meta <- 0
for (i in 1:length(HiekUrbChanMETADATAcombined1[,1])){
  for (j in 1:length(rownames(cmdBCurtis))) {
    if (rownames(cmdBCurtis)[j] == rownames(HiekUrbChanMETADATAcombined1)[i]) {
      cmdBCurtis$meta[j] <- HiekUrbChanMETADATAcombined1[i,9]
    }
  }
}
cmdBCurtis$meta <- 0
for (i in 1:length(HiekUrbChanMETADATAcombined1[,1])){
  for (j in 1:length(rownames(cmdBCurtis))) {
    if (rownames(cmdBCurtis)[j] == rownames(HiekUrbChanMETADATAcombined1)[i]) {
      cmdBCurtis$meta[j] <- HiekUrbChanMETADATAcombined1[i,10]
    }
  }
}

# text(x,y,cmduni$meta,cex=0.8)
ggplot(cmdBCurtis,aes(x=Axis.1,y=Axis.2,color= meta)) +
  #geom_dotplot(y=cmduni[,2],binwidth = 0.004)
  geom_point() +
  stat_ellipse() 


#### Phyloseq Analysis ####
# Rarefaction is not performed
psjoined <- phyloseq(otu_table(HiekUrbChanASVcombined2, taxa_are_rows=FALSE),
                     sample_data(HiekUrbChanMETADATAcombined1), # originally sample_data
                     tax_table(HiekUrbChanTAXAcombined2))
## do not have mock samples in our data so do not need this
# psjoined <- prune_samples(sample_names(psjoined) != "Mock", psjoined)
dnajoined <- Biostrings::DNAStringSet(taxa_names(psjoined))
names(dnajoined) <- taxa_names(psjoined)
psjoined <- merge_phyloseq(psjoined, dnajoined)
taxa_names(psjoined) <- paste0("ASV", seq(ntaxa(psjoined)))
psjoined

## Bray NMDS ordination plot
# ps.propjoined <- transform_sample_counts(psjoined, function(otu) otu/sum(otu))
ord.nmds.brayjoined <- ordinate(psjoined, method="NMDS", distance="bray")
plot_ordination(psjoined, ord.nmds.brayjoined, color="StudyGroup", title="Bray NMDS")

allseqsjoined <- names(sort(taxa_sums(psjoined), decreasing=TRUE)) # [1:700]
psallseqsjoined <- transform_sample_counts(psjoined, function(OTU) OTU/sum(OTU))
psallseqsjoined <- prune_taxa(allseqsjoined, psallseqsjoined)
psjoinedg <- tax_glom(psallseqsjoined, "Phylum")

pscombinejoined <- names(sort(taxa_sums(psjoined), decreasing=TRUE)) # [1:700]
psglomjoined <- transform_sample_counts(psjoined, function(OTU) OTU/sum(OTU))
psglomjoined <- prune_taxa(pscombinejoined, psglomjoined)

## tax_glom function to collapse at the taxa level
psjoinedg <- tax_glom(psglomjoined, "Phylum")
ps0joinedg <- transform_sample_counts(psjoinedg, function(x) x / sum(x))
ps1joinedg <- merge_samples(ps0joinedg, "StudyGroup") #yields NAs in columns of samdfjoined within ps1joinedg
ps2joinedg <- transform_sample_counts(ps1joinedg, function(x) x / sum(x))
pjoinedg <- plot_bar(ps2joinedg, fill= "Phylum")


## Another way to make Bray NMDS or Bray-Curtis plot
ps_braycdat <- phyloseq::distance(psallseqsjoined, method = "bray") # calculating bray curtis distance matrix
pcoabraycurt <- pcoa(ps_braycdat)
# biplot(pcoabraycurt)
pcoabray <- data.frame(pcoabraycurt$vectors)
pcoabray$meta <- "zero"
for (i in 1:length(rownames(pcoabray))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(pcoabray)[i] == HMeta$sample_id[j]) {
      pcoabray$meta[i] <- paste("Hiek_",HMeta$sample_type[j],sep = '')
    }
  }
}
# paste("Hiek_", paste(HMeta$sample_type[j], j, sep = ''),sep = '')
for (i in 1:length(rownames(pcoabray))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(pcoabray)[i] == samdfudat$sample_id[j]) {
      pcoabray$meta[i] <- paste("Urb_", samdfudat$sample_type[j],sep = '')
    }
  }
}

for (i in 1:length(rownames(pcoabray))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(pcoabray)[i] == samdfcdat$sample_id[j]) {
      if (samdfcdat$sample_type[j] == "Healthy_Control_Women") {
        pcoabray$meta[i] <- "Cha_Healthy"
      }
      if (samdfcdat$sample_type[j] == "Women_with_a_History_of_Breast_Cancer") {
        pcoabray$meta[i] <- "Cha_Cancer"
      }
    }
  }
}

pcoabray1 <- data.frame(pcoabray)

# this plot does not label the chan data points for some reason
pcoabray1 %>%
  ggplot(x=pcoabray1[,1], y=pcoabray1[,2],color= pcoabray1$meta) +
  geom_point(x=pcoabray1[,1], y=pcoabray1[,2]) +
  geom_text(x=pcoabray1[,1], y=pcoabray1[,2], label = pcoabray1$meta,nudge_x = 0.09, nudge_y = 0.09, 
            check_overlap = T) +
  xlim(-0.6,0.8) +
  ylim(-0.7,0.8)

x <- pcoabray[100:244,1]
y <- pcoabray[100:244,2]
ggplot(pcoabray[100:244,], aes(x,y)) +
  geom_point() +
  geom_text(label = pcoabray$meta[100:244],nudge_x = 0.09, nudge_y = 0.09, 
            check_overlap = T
  ) +
  xlim(-0.6,0.8) +
  ylim(-0.3,0.8)
rownames(pcoabray) == rownames(HiekUrbChanMETADATAcombined1)

# Define color for each of the 3 iris species
colors <- c("#00AFBB", "#E7B800", "#FC4E07")
colors <- colors[as.character(pcoabray1$meta)]

# Define shapes
shapes <- c(16, 17, 18) 
shapes <- shapes[as.character(pcoabray1$meta)]

# Plot
x <- pcoabray1[,1]
y <- pcoabray1[,2]
plot(x,y, frame = FALSE,
     pch=16)
text(x,y,pcoabray1$meta,cex=0.25)
# legend("topright", legend = unique(pcoabray1$meta),
#        col =  c("#00AFBB", "#E7B800", "#FC4E07"),
#        pch = c(16, 17, 18) )



ps_brayord <- ordinate(psjoined, method="NMDS", distance="bray", color="StudyGroup")
plot_ordination(psjoined,ps_brayord)
sampledf <- data.frame(sample_data(psjoined))
ps_braycdat1 <- data.frame(ps_braycdat)
ggplot(data = sampledf,aes(x=StudyGroup, y=ps_braycdat)) +
  geom_point(aes(color = StudyGroup)) 
#labs(title = paste("PBS p-value",adonisChan_health$aov.tab$`Pr(>F)`[1], sep=" "))



# Pick the relative abundance table
rel_abund_assay <- psallseqsjoined@otu_table

# Calculates Bray-Curtis distances between samples. Because taxa is in
# columns, it is used to compare different samples. We transpose the
# assay to get taxa to columns
bray_curtis_dist <- vegan::vegdist(t(rel_abund_assay), method = "bray")

# PCoA
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)

# All components could be found here: 
# bray_curtis_pcoa$vectors
# But we only need the first two to demonstrate what we can do:
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

# Create a plot
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10)) # makes titles smaller

bray_curtis_plot


#### metamicrobiomeR - for normalizing data ####
# Comparison of bacterial taxa relative abundance up to genus level 
library(devtools)
#install and load package metamicrobiomeR
# install_github("nhanhocu/metamicrobiomeR")
library(metamicrobiomeR) 
#Load other needed packages 
library(knitr)
library(plyr)
# library(dplyr)
library(gdata)
library(gridExtra)
# library(ggplot2)
library(lme4) 
library(lmerTest)
library(mgcv) 
library(meta) 
load("~/Downloads/taxtab.rm7.rda")
load("~/Downloads/taxacom.rm.sex.adjustbfage.rda")
data(taxtab6)
taxacompareMetaAnalysis <- taxa.compare(taxtab=taxtab6.rm[[5]],propmed.rel="gamlss",comvar="gender",adjustvar=c("bf","age.sample"),longitudinal="yes")
psallseqsdf <- data.frame(psallseqsjoined@otu_table)
psallseqsdf$meta <- "zero"
psallseqsdf$type_v <- "zero"
for (i in 1:length(rownames(psallseqsdf))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(psallseqsdf)[i] == HMeta$sample_id[j]) {
      psallseqsdf$meta[i] <- paste("Hiek_",HMeta$sample_type[j],sep = '')
    }
  }
}
# paste("Hiek_", paste(HMeta$sample_type[j], j, sep = ''),sep = '')
for (i in 1:length(rownames(psallseqsdf))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(psallseqsdf)[i] == samdfudat$sample_id[j]) {
      psallseqsdf$meta[i] <- paste("Urb_", samdfudat$sample_type[j],sep = '')
    }
  }
}

for (i in 1:length(rownames(psallseqsdf))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(psallseqsdf)[i] == samdfcdat$sample_id[j]) {
      if (samdfcdat$sample_type[j] == "Healthy_Control_Women") {
        psallseqsdf$meta[i] <- "Cha_Healthy"
      }
      if (samdfcdat$sample_type[j] == "Women_with_a_History_of_Breast_Cancer") {
        psallseqsdf$meta[i] <- "Cha_Cancer"
      }
    }
  }
}

for (i in 1:length(rownames(psallseqsdf))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(psallseqsdf)[i] == HMeta$sample_id[j]) {
      psallseqsdf$type_v[i] <- HMeta$sample_type[j]
    }
  }
}
# paste("Hiek_", paste(HMeta$sample_type[j], j, sep = ''),sep = '')
for (i in 1:length(rownames(psallseqsdf))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(psallseqsdf)[i] == samdfudat$sample_id[j]) {
      psallseqsdf$type_v[i] <- samdfudat$sample_type[j]
    }
  }
}

for (i in 1:length(rownames(psallseqsdf))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(psallseqsdf)[i] == samdfcdat$sample_id[j]) {
      if (samdfcdat$sample_type[j] == "Healthy_Control_Women") {
        psallseqsdf$type_v[i] <- "healthy"
      }
      if (samdfcdat$sample_type[j] == "Women_with_a_History_of_Breast_Cancer") {
        psallseqsdf$type_v[i] <- "Cancer"
      }
    }
  }
}
typeof(psallseqsdf)
str(psallseqsdf)
psallseqsdf1 <- data.frame(psallseqsdf)
# write.table(psallseqsdf1,file = "~/Documents/Chan10082021/metamicrobmetaanalysis.txt")
psmetamicrob <- read.table("~/Documents/Chan10082021/metamicrobmetaanalysis.txt",header = TRUE)
psmetamicrob1 <- data.frame(psmetamicrob)
View(psmetamicrob1)
# did not work:
taxacompareMetaAnalysis <- taxa.compare(taxtab=psmetamicrob1,propmed.rel="gamlss", comvar = "meta", adjustvar = "type_v", longitudinal="no",personid = rownames(psmetamicrob1), pooldata = TRUE)
#phylum
kable(taxcomtab.show(taxcomtab=taxacompareMetaAnalysis,tax.select="none", showvar="genderMale", tax.lev="l2",p.adjust.method="fdr"))
#order
kable(taxcomtab.show(taxcomtab=taxacompareMetaAnalysis,tax.select="none", showvar="genderMale", tax.lev="l4",p.adjust.method="fdr"))
#family
kable(taxcomtab.show(taxcomtab=taxacompareMetaAnalysis,tax.select="none", showvar="genderMale", tax.lev="l5",p.adjust.method="fdr"))
#genus
kable(taxcomtab.show(taxcomtab=taxacompareMetaAnalysis,tax.select="none", showvar="genderMale", tax.lev="l6",p.adjust.method="fdr"))

# taxa.compare
# taxa.compare = edit()

newtaxa.compare <- function (taxtab, propmed.rel = "gamlss", transform = "none", 
          zeroreplace.method = "none", comvar, adjustvar, personid = "personid", 
          longitudinal = "yes", percent.filter = 0.05, relabund.filter = 5e-05, 
          p.adjust.method = "fdr", ...) 
{
  require(gamlss)
  taxdat <- as.data.frame(taxtab)
  taxdat[, comvar] <- gdata::drop.levels(taxdat[, comvar], 
                                         reorder = FALSE)
  if (longitudinal == "yes") {
    taxdat$personid <- as.factor(taxdat[, personid])
  }
  taxlist <- colnames(taxdat)[grep("k__", colnames(taxdat))]
  taxtest <- apply(taxdat[, taxlist], 2, function(x) {
    length(x[!is.na(x) & x > 0])
  })
  taxget <- taxtest[taxtest >= percent.filter * (nrow(taxdat))]
  taxtestm <- apply(taxdat[, taxlist], 2, mean, na.rm = T)
  taxgetm <- taxtestm[taxtestm > relabund.filter]
  taxname <- names(taxget)[names(taxget) %in% names(taxgetm)]
  if (propmed.rel == "gamlss" & transform != "none") {
    stop("gamlss with beta zero-inflated family should only be used for relative abundance without transformation")
  }
  if (transform != "clr" & zeroreplace.method != "none") {
    stop("Zero replacement is only implemented for use with CLR transformation")
  }
  if (transform == "clr" & zeroreplace.method == "none") {
    stop("Zero replacement needs to be done before CLR transformation")
  }
  if (propmed.rel == "lm" & transform == "asin.sqrt") {
    asintransform <- function(p) {
      asin(sqrt(p))
    }
    taxdat[, taxname] <- apply(taxdat[, taxname], 2, asintransform)
  }
  if (propmed.rel == "lm" & transform == "logit") {
    logittransform <- function(p) {
      log(p/(1 - p))
    }
    taxdat[, taxname] <- apply(taxdat[, taxname], 2, logittransform)
  }
  if (propmed.rel == "lm" & transform == "clr") {
    if (zeroreplace.method == "multLN") {
      test0 <- zCompositions::multLN(taxdat[, taxname], 
                                     label = 0, dl = rep(1, length(taxname)))
    }
    if (zeroreplace.method == "multKM") {
      test0 <- zCompositions::multKM(taxdat[, taxname], 
                                     label = 0, dl = rep(1, length(taxname)))
    }
    if (zeroreplace.method == "multRepl") {
      test0 <- zCompositions::multRepl(taxdat[, taxname], 
                                       label = 0, dl = rep(1, length(taxname)))
    }
    if (zeroreplace.method == "lrEM") {
      test0 <- zCompositions::lrEM(taxdat[, taxname], label = 0, 
                                   dl = rep(1, length(taxname)))
    }
    if (zeroreplace.method == "lrDA") {
      test0 <- zCompositions::lrDA(taxdat[, taxname], label = 0, 
                                   dl = rep(1, length(taxname)))
    }
    clrdat <- as.data.frame(compositions::clr(test0))
    taxdat[, taxname] <- clrdat
  }
  require(matrixStats)
  GMPR <- function(comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
    comm[comm < ct.min] <- 0
    if (is.null(colnames(comm))) {
      colnames(comm) <- paste0("S", 1:ncol(comm))
    }
    if (trace) 
      cat("Begin GMPR size factor calculation ...\n")
    comm.no <- numeric(ncol(comm))
    gmpr <- sapply(1:ncol(comm), function(i) {
      if (i%%50 == 0) {
        cat(i, "\n")
      }
      x <- comm[, i]
      pr <- x/comm
      pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
      incl.no <- colSums(!is.na(pr))
      pr.median <- matrixStats::colMedians(pr, na.rm = TRUE)
      comm.no[i] <<- sum(incl.no >= intersect.no)
      if (comm.no[i] > 1) {
        return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
      }
      else {
        return(NA)
      }
    })
    if (sum(is.na(gmpr))) {
      warning(paste0("The following samples\n ", paste(colnames(comm)[is.na(gmpr)], 
                                                       collapse = "\n"), "\ndo not share at least ", 
                     intersect.no, " common taxa with the rest samples! ", 
                     "For these samples, their size factors are set to be NA! \n", 
                     "You may consider removing these samples since they are potentially outliers or negative controls!\n", 
                     "You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n"))
    }
    if (trace) 
      cat("Completed!\n")
    if (trace) 
      cat("Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n")
    names(gmpr) <- names(comm.no) <- colnames(comm)
    attr(gmpr, "NSS") <- comm.no
    return(gmpr)
  }
  if (transform == "gmpr") {
    reltab <- t(taxdat[, taxname])
    gmpr.sf <- GMPR(comm = reltab, intersect.no = 2, ct.min = 0, 
                    trace = TRUE)
    rel.gmpr <- t(t(reltab)/gmpr.sf)
    taxdat[, taxname] <- rel.gmpr
  }
  estisum <- NULL
  for (i in 1:length(taxname)) {
    print(i)
    if (propmed.rel == "lm") {
      if (longitudinal == "yes") {
        fitsum <- try(summary(lme4::glmer(as.formula(paste(taxname[i], 
                                                           paste(c(comvar, adjustvar, "(1|personid)"), 
                                                                 collapse = "+"), sep = "~")), data = taxdat, 
                                          family = gaussian(link = "identity"))))
      }
      if (longitudinal == "no") {
        fitsum <- try(summary(glm(as.formula(paste(taxname[i], 
                                                   paste(c(comvar, adjustvar), collapse = "+"), 
                                                   sep = "~")), data = taxdat, family = "gaussian")))
      }
      if (class(fitsum) == "try-error") {
        cat("Error in model fit, NA introduced.\n")
        fitcoefw <- NULL
        estisum <- plyr::rbind.fill(estisum, fitcoefw)
      }
      if (class(fitsum) != "try-error") {
        if (length(which(rownames(fitsum$coefficients) != 
                         "(Intercept)")) > 1) {
          fitcoef <- as.data.frame(fitsum$coefficients[rownames(fitsum$coefficients) != 
                                                         "(Intercept)", ])
          if (longitudinal == "yes") {
            fitcoef[, "Pr(>|t|)"] <- 1.96 * pnorm(-abs(fitcoef[, 
                                                               "Estimate"]/fitcoef[, "Std. Error"]))
          }
          fitcoef[, "varname"] <- rownames(fitcoef)
          fitcoef[, "id"] <- taxname[i]
          fitcoefw <- reshape(fitcoef, idvar = "id", 
                              timevar = "varname", direction = "wide")
        }
        if (length(which(rownames(fitsum$coefficients) != 
                         "(Intercept)")) == 1) {
          fitcoef <- as.data.frame(matrix(fitsum$coefficients[rownames(fitsum$coefficients) != 
                                                                "(Intercept)", ], ncol = ncol(fitsum$coefficients)))
          rownames(fitcoef) <- rownames(fitsum$coefficients)[rownames(fitsum$coefficients) != 
                                                               "(Intercept)"]
          colnames(fitcoef) <- colnames(fitsum$coefficients)
          if (longitudinal == "yes") {
            fitcoef[, "Pr(>|t|)"] <- 1.96 * pnorm(-abs(fitcoef[, 
                                                               "Estimate"]/fitcoef[, "Std. Error"]))
          }
          fitcoef[, "varname"] <- rownames(fitcoef)
          fitcoef[, "id"] <- taxname[i]
          fitcoefw <- reshape(fitcoef, idvar = "id", 
                              timevar = "varname", direction = "wide")
        }
        if (length(which(rownames(fitsum$coefficients) != 
                         "(Intercept)")) == 0) {
          fitcoefw <- NULL
        }
        estisum <- plyr::rbind.fill(estisum, fitcoefw)
      }
    }
    if (propmed.rel == "gamlss") {
      if (longitudinal == "yes") {
        testdat <- taxdat[, c(taxname[i], comvar, adjustvar, "personid")]
        testdat[, taxname[i]][testdat[, taxname[i]] >= 1] <- 0.9999
        testdat <- stats::na.omit(testdat)
        fitsum <- try(summary(gamlss::gamlss(stats::as.formula(paste(taxname[i], 
                                                                     paste(c(comvar, adjustvar, "random(personid)"), 
                                                                           collapse = "+"), sep = "~")), family = BEZI, 
                                             data = testdat, trace = FALSE), save = TRUE))
      }
      if (longitudinal == "no") {
        testdat <- taxdat[, c(taxname[i], comvar, adjustvar)]
        testdat[, taxname[i]][testdat[, taxname[i]] >= 1] <- 0.9999
        testdat <- stats::na.omit(testdat)
        fitsum <- try(summary(gamlss::gamlss(stats::as.formula(paste(taxname[i], 
                                                                     paste(c(comvar, adjustvar), collapse = "+"), 
                                                                     sep = "~")), family = BEZI, data = testdat, 
                                             trace = FALSE), save = TRUE))
      }
      if (class(fitsum) == "try-error") {
        cat("Error in model fit, NA introduced.\n")
        fitcoefw <- NULL
        estisum <- plyr::rbind.fill(estisum, fitcoefw)
      }
      if (class(fitsum) != "try-error") {
        if (length(which(rownames(fitsum$coef.table) != 
                         "(Intercept)")) > 1) {
          fitcoef <- as.data.frame(fitsum$coef.table[rownames(fitsum$coef.table) != 
                                                       "(Intercept)", ])
          fitcoef[, "varname"] <- rownames(fitcoef)
          fitcoef[, "id"] <- taxname[i]
          fitcoefw <- reshape(fitcoef, idvar = "id", 
                              timevar = "varname", direction = "wide")
        }
        if (length(which(rownames(fitsum$coef.table) != 
                         "(Intercept)")) == 1) {
          fitcoef <- as.data.frame(matrix(fitsum$coef.table[rownames(fitsum$coef.table) != 
                                                              "(Intercept)", ], ncol = ncol(fitsum$coef.table)))
          rownames(fitcoef) <- rownames(fitsum$coef.table)[rownames(fitsum$coef.table) != 
                                                             "(Intercept)"]
          colnames(fitcoef) <- colnames(fitsum$coef.table)
          fitcoef[, "varname"] <- rownames(fitcoef)
          fitcoef[, "id"] <- taxname[i]
          fitcoefw <- reshape(fitcoef, idvar = "id", 
                              timevar = "varname", direction = "wide")
        }
        if (length(which(rownames(fitsum$coef.table) != 
                         "(Intercept)")) == 0) {
          fitcoefw <- NULL
        }
        estisum <- plyr::rbind.fill(estisum, fitcoefw)
      }
    }
  }
  np <- length(colnames(estisum)[grep("Pr(>|t|)", colnames(estisum))])
  if (np > 1) {
    estisum[, sub(".*\\.", "pval.adjust.", colnames(estisum)[grep("Pr(>|t|)", 
                                                                  colnames(estisum))])] <- apply(estisum[, colnames(estisum)[grep("Pr(>|t|)", 
                                                                                                                                  colnames(estisum))]], 2, p.adjust, method = p.adjust.method)
  }
  else {
    estisum[, sub(".*\\.", "pval.adjust.", colnames(estisum)[grep("Pr(>|t|)", 
                                                                  colnames(estisum))])] <- p.adjust(estisum[, colnames(estisum)[grep("Pr(>|t|)", 
                                                                                                                                     colnames(estisum))]], method = p.adjust.method)
    }
    return(estisum)
}

#taxacompareMetaAnalysis <- newtaxa.compare(taxtab=psmetamicrob1,propmed.rel="gamlss", comvar = "meta", adjustvar = "type_v", longitudinal="no", pooldata = TRUE)
for (i in 1:length(psmetamicrob1$type_v)) {
  if (psmetamicrob1$type_v[i] == "healthy") {
    psmetamicrob1$type_v[i] <- as.numeric(0)
  }
}
for (i in 1:length(psmetamicrob1$type_v)) {
  if (psmetamicrob1$type_v[i] == "Cancer") {
    psmetamicrob1$type_v[i] <- as.numeric(2)
  }
}
for (i in 1:length(psmetamicrob1$type_v)) {
  if (psmetamicrob1$type_v[i] == "Benign") {
    psmetamicrob1$type_v[i] <- as.numeric(1)
  }
}
for (i in 1:length(psmetamicrob1$type_v)) {
  if (psmetamicrob1$type_v[i] == "Malignant") {
    psmetamicrob1$type_v[i] <- as.numeric(2)
  }
}

psmetamicrob1 <- psmetamicrob1[-16019] # with the metadata info
psmetamicrob1$type_v <- as.numeric(psmetamicrob1$type_v)
str(psmetamicrob$type_v[1])
# write.table(psmetamicrob1,file="~/Documents/Chan10082021/psmetam.txt")
psmetamicrob <- read.table("~/Documents/Chan10082021/psmetam.txt",header = TRUE)
psmetamicrob$typeva <- "zero"
for (i in 1:length(rownames(psmetamicrob))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(psmetamicrob)[i] == HMeta$sample_id[j]) {
      psmetamicrob$typeva[i] <- paste("Hiek_", HMeta$sample_type[j], sep = '')
    }
  }
}

for (i in 1:length(rownames(psmetamicrob))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(psmetamicrob)[i] == samdfudat$sample_id[j]) {
      psmetamicrob$typeva[i] <- paste("Urb_", samdfudat$sample_type[j], sep = '')
    }
  }
}

for (i in 1:length(rownames(psmetamicrob))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(psmetamicrob)[i] == samdfcdat$sample_id[j]) {
      if (samdfcdat$sample_type[j] == "Healthy_Control_Women") {
        psmetamicrob$typeva[i] <- "Cha_Healthy"
      }
      if (samdfcdat$sample_type[j] == "Women_with_a_History_of_Breast_Cancer") {
        psmetamicrob$typeva[i] <- "Cha_Cancer"
      }
    }
  }
}

for (i in 1:length(rownames(psmetamicrob))) {
  if (psmetamicrob$typeva[i] == "Hiek_Benign") {
    psmetamicrob$typeva[i] <- as.numeric(11)
  }
}
for (i in 1:length(rownames(psmetamicrob))) {
  if (psmetamicrob$typeva[i] == "Hiek_Malignant") {
    psmetamicrob$typeva[i] <- as.numeric(12)
  }
}
for (i in 1:length(rownames(psmetamicrob))) {
  if (psmetamicrob$typeva[i] == "Urb_Benign") {
    psmetamicrob$typeva[i] <- as.numeric(21)
  }
}
for (i in 1:length(rownames(psmetamicrob))) {
  if (psmetamicrob$typeva[i] == "Urb_Malignant") {
    psmetamicrob$typeva[i] <- as.numeric(22)
  }
}
for (i in 1:length(rownames(psmetamicrob))) {
  if (psmetamicrob$typeva[i] == "Urb_healthy") {
    psmetamicrob$typeva[i] <- as.numeric(20)
  }
}
for (i in 1:length(rownames(psmetamicrob))) {
  if (psmetamicrob$typeva[i] == "Cha_Healthy") {
    psmetamicrob$typeva[i] <- as.numeric(30)
  }
}
for (i in 1:length(rownames(psmetamicrob))) {
  if (psmetamicrob$typeva[i] == "Cha_Cancer") {
    psmetamicrob$typeva[i] <- as.numeric(32)
  }
}
seqs <- 1:length(colnames(psmetamicrob))
colnames(psmetamicrob) <- seqs
taxacompareMetaAnalysis <- newtaxa.compare(taxtab=psmetamicrob,propmed.rel="gamlss", comvar = 16019, adjustvar = 16020, longitudinal="no")

#### Percentile Normalization for all three studies - propotional/TSS percentile normalization is in each studies propabund R code ####
# example: Hiekenpropabund_03232022.R or Chanpropabund_05242022.R

#### Chan: Trying the percentile normalization python script GitHub: https://github.com/seangibbons/percentile_normalization#### 
## with sample data
## python: python percentile_norm.py -i baxter_crc_data.txt -case baxter_crc_samples.txt -control baxter_h_samples.txt -o baxter_percentile_norm.txt
## Chan python code: python percentile_norm.py -i ChanASV12202021upd_percentnormal.txt -case Chancases.txt -control Chancontrols.txt -o Chan_percentile_norm.txt
baxterpercentile <- read.table("~/Downloads/baxter_percentile_norm.txt")
newdatf <- data.frame(rep(0,length(rownames(samdfcdat))))
colnames(newdatf) <- c("samples")
for (i in 1:length(rownames(samdfcdat))) {
  if (samdfcdat$sample_type[i] == "Healthy_Control_Women") {
    newdatf$samples[i] <- samdfcdat$sample_id[i]
    #newdatf$type[i] <- paste(samdfcdat$sample[i],"_H",sep='')
  }
}
for (i in 1:length(newdatf$samples)) {
  if (newdatf$samples[i] == 0) {
    newdatf <- newdatf[-i,]
  }
}
for (i in 1:length(newdatf)) {
  if (newdatf[i] == 0) {
    newdatf <- newdatf[-i]
  }
}
# write.table(newdatf,file = "~/Downloads/Chancontrols.txt",row.names = FALSE,col.names = FALSE)
# write.table(newdatf,file = "~/Documents/Chan10082021/Chancontrols.txt",row.names = FALSE)
newdatf <- data.frame(rep(0,length(rownames(samdfcdat))))
colnames(newdatf) <- c("samples")
for (i in 1:length(colnames(seqtab.nochimcdat))) {
  if (samdfcdat$sample_type[i] == "Women_with_a_History_of_Breast_Cancer") {
    newdatf$samples[i] <- samdfcdat$sample_id[i]
    #newdatf$type[i] <- paste(samdfcdat$sample[i],"_C",sep='')
  }
}
for (i in 1:length(newdatf$samples)) {
  if (newdatf$samples[i] == 0) {
    newdatf <- newdatf[-i,]
  }
}
for (i in 1:length(newdatf)) {
  if (newdatf[i] == 0) {
    newdatf <- newdatf[-i]
  }
}
# write.table(newdatf,file = "~/Documents/Chan10082021/Chancases.txt",row.names = FALSE)
# write.table(newdatf,file = "~/Downloads/Chancases.txt",row.names = FALSE,col.names = FALSE)
# files for test python are all in Downloads Folder - upd 05212022
# test python: python percentile_norm.py -i ChanASV12202021upd_percentnormal.txt -case Chancases.txt -control Chancontrols.txt -o Chan_percentile_norm.txt
# ~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChanASV12202021upd_R.txt
# baxterdata <- read.table("~/Downloads/baxter_crc_data.txt")
# already completed the below:
# baxterdata <- read.table("~/Downloads/ChanASV12202021upd_percentnormal.txt")
# baxterdata <- t(baxterdata)
# write.table(baxterdata,file = "~/Downloads/ChanASV12202021upd_percentnormal.txt")

# the ChanASV12202021upd_percnormaltxt.txt file is the ChanASV12202021upd_R.txt file copied from the Documents/Chan10082021 folder to the Downloads folder
baxterdata <- read.table("~/Downloads/ChanASV12202021upd_percnormaltxt.txt")
baxterdata <- t(baxterdata)
View(baxterdata)
one1 <- paste(taxasilvacdat[i,2],taxasilvacdat[i,3],sep = '_')
two2 <- paste(one1,paste(taxasilvacdat[i,4],taxasilvacdat[i,5], sep = '_'), sep = '_')
three3 <- paste(two2, paste(taxasilvacdat[i,6],taxasilvacdat[i,7],sep = '_'), sep = '_')
# baxterdata <- t(seqtab.nochimcdat)

# changing the asv sequences to the asv taxa
for (i in 1:length(colnames(baxterdata))) {
  if (colnames(baxterdata)[i] == taxasilvacdat$Taxonomy[i]) {
    one1 <- paste(taxasilvacdat[i,2],taxasilvacdat[i,3],sep = '_')
    two2 <- paste(one1,paste(taxasilvacdat[i,4],taxasilvacdat[i,5], sep = '_'), sep = '_')
    three3 <- paste(two2, paste(taxasilvacdat[i,6],taxasilvacdat[i,7],sep = '_'), sep = '_')
    colnames(baxterdata)[i] <- three3
  }
}
# write.table(baxterdata,file = "~/Downloads/ChanASV12202021upd_percentnormal.txt")

# storing the control SRR samples in the order that they occur in the asv table - there are 57 samples
newdatf <- data.frame(rep(0,length(rownames(samdfcdat))))
colnames(newdatf) <- "samples"
for (i in 1:length(rownames(baxterdata))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(baxterdata)[i] == samdfcdat$sample_id[j]) {
      if (samdfcdat$sample_type[j] == "Healthy_Control_Women") {
        newdatf$samples[i] <- samdfcdat$sample_id[j]
    #newdatf$type[i] <- paste(samdfcdat$sample[i],"_C",sep='')
      }
    }
  }
}

# storing the cases SRR samples in the order that they occur in the asv table - there are 47 samples
newdatf1 <- data.frame(rep(0,length(rownames(samdfcdat))))
colnames(newdatf1) <- "samples"
for (i in 1:length(rownames(baxterdata))) {
  for (j in 1:length(samdfcdat$sample_id)) {
    if (rownames(baxterdata)[i] == samdfcdat$sample_id[j]) {
      if (samdfcdat$sample_type[j] == "Women_with_a_History_of_Breast_Cancer") {
        newdatf1$samples[i] <- samdfcdat$sample_id[j]
        #newdatf1$type[i] <- paste(samdfcdat$sample[i],"_C",sep='')
      }
    }
  }
}
for (i in 1:length(newdatf$samples)) {
  if (newdatf$samples[i] == 0) {
    newdatf <- newdatf[-i,]
  }
}
for (i in 1:length(newdatf)) {
  if (newdatf[i] == 0) {
    newdatf <- newdatf[-i]
  }
}
for (i in 1:length(newdatf1$samples)) {
  if (newdatf1$samples[i] == 0) {
    newdatf1 <- newdatf1[-i,]
  }
}
for (i in 1:length(newdatf1)) {
  if (newdatf1[i] == 0) {
    newdatf1 <- newdatf1[-i]
  }
}
write.table(newdatf, file = "~/Downloads/Chancontrols.txt",row.names = FALSE, col.names = FALSE)
write.table(newdatf1, file = "~/Downloads/Chancases.txt",row.names = FALSE, col.names = FALSE)

#### Formatting files for input to the percent normalization python script ####
## the ChanASV12202021upd_percentnormal.txt file:
# transpose this file so that the taxa/asvs are the columns and the SRR samples are rows
# Copy the taxa to a separate file
# first change the \t (tabs) to dashes (-)
# then change the \n (new lines) to \t (tabs)
# then copy this modified formatted taxa to the ChanASV12202021upd_percentnormal.txt file
# and remove the " (quotation marks)
# in the ChanASV12202021upd_percentnormal.txt file change the spaces ( ) to \t (tabs)
# go to the first taxa (at the top of the ChanASV12202021upd_percentnormal.txt file) and add a \t (tab)

## the Chancontrols.txt file and the Chancases.txt file
# change the \n (new lines) to \t (tabs)
# remove the " (quotation marks)
# remove the last \t (tab) at the end of the file
# and add a \n (new line)

# View the percent normalization output:
percnormal <- read.table("~/Downloads/Chan_percentile_norm.txt")
View(percnormal)
# sum(percnormal$Bacteria_NA_NA_NA_NA_NA.42)

#### Hieken: Trying the percentile normalization python script GitHub: https://github.com/seangibbons/percentile_normalization#### 
## with sample data
## python: python percentile_norm.py -i baxter_crc_data.txt -case baxter_crc_samples.txt -control baxter_h_samples.txt -o baxter_percentile_norm.txt
## Hieken python: python percentile_norm.py -i HiekenASV12202021upd_percentnormal.txt -case Hiekencases.txt -control Hiekencontrols.txt -o Hieken_percentile_norm.txt

# the seqtabnochimhieken10082021_percnormal.txt file is the seqtabnochimhieken10082021.txt file copied from the Documents/Hieken10082021 folder to the Downloads folder
baxterdata <- read.table("~/Downloads/seqtabnochimhieken10082021_percnormal.txt")
baxterdata <- t(baxterdata)
View(baxterdata)
one1 <- paste(Hiekentaxa[i,2],Hiekentaxa[i,3],sep = '_')
two2 <- paste(one1,paste(Hiekentaxa[i,4],Hiekentaxa[i,5], sep = '_'), sep = '_')
three3 <- paste(two2, paste(Hiekentaxa[i,6],Hiekentaxa[i,7],sep = '_'), sep = '_')
# baxterdata <- t(seqtab.nochimcdat)

# changing the asv sequences to the asv taxa
for (i in 1:length(colnames(baxterdata))) {
  if (colnames(baxterdata)[i] == Hiekentaxa$Taxonomy[i]) {
    one1 <- paste(Hiekentaxa[i,2],Hiekentaxa[i,3],sep = '_')
    two2 <- paste(one1,paste(Hiekentaxa[i,4],Hiekentaxa[i,5], sep = '_'), sep = '_')
    three3 <- paste(two2, paste(Hiekentaxa[i,6],Hiekentaxa[i,7],sep = '_'), sep = '_')
    colnames(baxterdata)[i] <- three3
  }
}
# write.table(baxterdata,file = "~/Downloads/HiekenASV12202021upd_percentnormal.txt")

# storing the control SRR samples in the order that they occur in the asv table - there are 41 samples
newdatf <- data.frame(rep(0,length(rownames(HMeta))))
colnames(newdatf) <- "samples"
for (i in 1:length(rownames(baxterdata))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(baxterdata)[i] == HMeta$sample_id[j]) {
      if (HMeta$sample_type[j] == "Benign") {
        newdatf$samples[i] <- HMeta$sample_id[j]
        #newdatf$type[i] <- paste(HMeta$sample[i],"_C",sep='')
      }
    }
  }
}

# storing the cases SRR samples in the order that they occur in the asv table - there are 57 samples
newdatf1 <- data.frame(rep(0,length(rownames(HMeta))))
colnames(newdatf1) <- "samples"
for (i in 1:length(rownames(baxterdata))) {
  for (j in 1:length(HMeta$sample_id)) {
    if (rownames(baxterdata)[i] == HMeta$sample_id[j]) {
      if (HMeta$sample_type[j] == "Malignant") {
        newdatf1$samples[i] <- HMeta$sample_id[j]
        #newdatf1$type[i] <- paste(HMeta$sample[i],"_C",sep='')
      }
    }
  }
}
for (i in 1:length(newdatf$samples)) {
  if (newdatf$samples[i] == 0) {
    newdatf <- newdatf[-i,]
  }
}
for (i in 1:length(newdatf)) {
  if (newdatf[i] == 0) {
    newdatf <- newdatf[-i]
  }
}
for (i in 1:length(newdatf1$samples)) {
  if (newdatf1$samples[i] == 0) {
    newdatf1 <- newdatf1[-i,]
  }
}
for (i in 1:length(newdatf1)) {
  if (newdatf1[i] == 0) {
    newdatf1 <- newdatf1[-i]
  }
}
write.table(newdatf, file = "~/Downloads/Hiekencontrols.txt",row.names = FALSE, col.names = FALSE)
write.table(newdatf1, file = "~/Downloads/Hiekencases.txt",row.names = FALSE, col.names = FALSE)

#### Formatting files for input to the percent normalization python script ####
## the HiekenASV12202021upd_percentnormal.txt file:
# transpose this file so that the taxa/asvs are the columns and the SRR samples are rows
# Copy the taxa to a separate file
# first change the \t (tabs) to dashes (-)
# then change the \n (new lines) to \t (tabs)
# then copy this modified formatted taxa to the HiekenASV12202021upd_percentnormal.txt file
# and remove the " (quotation marks)
# in the HiekenASV12202021upd_percentnormal.txt file change the spaces ( ) to \t (tabs)
# go to the first taxa (at the top of the HiekenASV12202021upd_percentnormal.txt file) and add a \t (tab)

## the Hiekencontrols.txt file and the Hiekencases.txt file
# change the \n (new lines) to \t (tabs)
# remove the " (quotation marks)
# remove the last \t (tab) at the end of the file
# and add a \n (new line)

# View the percent normalization output:
percnormal <- read.table("~/Downloads/Hieken_percentile_norm.txt")
View(percnormal)
# sum(percnormal$Bacteria_NA_NA_NA_NA_NA.42)


#### Urbaniak: Trying the percentile normalization python script GitHub: https://github.com/seangibbons/percentile_normalization#### 
## with sample data
## python: python percentile_norm.py -i baxter_crc_data.txt -case baxter_crc_samples.txt -control baxter_h_samples.txt -o baxter_percentile_norm.txt
## Urbaniak python: python percentile_norm.py -i UrbaniakASV12202021upd_percentnormal.txt -case Urbaniakcases.txt -control Urbaniakcontrols.txt -o Urbaniak_percentile_norm.txt

# the seqtabnochimUrbaniak10082021_percnormal.txt file is the seqtabnochimUrbaniak10082021.txt file copied from the Documents/Urbaniak10082021 folder to the Downloads folder
baxterdata <- read.table("~/Downloads/UrbASVmodifiedSRRs_percnormal.txt")
baxterdata <- t(baxterdata)
View(baxterdata)
one1 <- paste(taxasilvaudat[i,2],taxasilvaudat[i,3],sep = '_')
two2 <- paste(one1,paste(taxasilvaudat[i,4],taxasilvaudat[i,5], sep = '_'), sep = '_')
three3 <- paste(two2, paste(taxasilvaudat[i,6],taxasilvaudat[i,7],sep = '_'), sep = '_')
# baxterdata <- t(seqtab.nochimcdat)

# changing the asv sequences to the asv taxa
for (i in 1:length(colnames(baxterdata))) {
  if (colnames(baxterdata)[i] == taxasilvaudat$Taxonomy[i]) {
    one1 <- paste(taxasilvaudat[i,2],taxasilvaudat[i,3],sep = '_')
    two2 <- paste(one1,paste(taxasilvaudat[i,4],taxasilvaudat[i,5], sep = '_'), sep = '_')
    three3 <- paste(two2, paste(taxasilvaudat[i,6],taxasilvaudat[i,7],sep = '_'), sep = '_')
    colnames(baxterdata)[i] <- three3
  }
}
# write.table(baxterdata,file = "~/Downloads/UrbaniakASV12202021upd_percentnormal.txt")

# storing the control SRR samples in the order that they occur in the asv table - there are 12 samples
newdatf <- data.frame(rep(0,length(rownames(samdfudat))))
colnames(newdatf) <- "samples"
for (i in 1:length(rownames(baxterdata))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(baxterdata)[i] == samdfudat$sample_id[j]) {
      if (samdfudat$sample_type[j] == "healthy") {
        newdatf$samples[i] <- samdfudat$sample_id[j]
        #newdatf$type[i] <- paste(samdfudat$sample[i],"_C",sep='')
      }
    }
  }
}

# storing the cases SRR samples in the order that they occur in the asv table - there are 30 samples
newdatf1 <- data.frame(rep(0,length(rownames(samdfudat))))
colnames(newdatf1) <- "samples"
for (i in 1:length(rownames(baxterdata))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(baxterdata)[i] == samdfudat$sample_id[j]) {
      if (samdfudat$sample_type[j] == "Malignant") {
        newdatf1$samples[i] <- samdfudat$sample_id[j]
        #newdatf1$type[i] <- paste(samdfudat$sample[i],"_C",sep='')
      }
    }
  }
}
for (i in 1:length(rownames(baxterdata))) {
  for (j in 1:length(samdfudat$sample_id)) {
    if (rownames(baxterdata)[i] == samdfudat$sample_id[j]) {
      if (samdfudat$sample_type[j] == "Benign") {
        newdatf1$samples[i] <- samdfudat$sample_id[j]
        #newdatf1$type[i] <- paste(samdfudat$sample[i],"_C",sep='')
      }
    }
  }
}
for (i in 1:length(newdatf$samples)) {
  if (newdatf$samples[i] == 0) {
    newdatf <- newdatf[-i,]
  }
}
for (i in 1:length(newdatf)) {
  if (newdatf[i] == 0) {
    newdatf <- newdatf[-i]
  }
}
for (i in 1:length(newdatf1$samples)) {
  if (newdatf1$samples[i] == 0) {
    newdatf1 <- newdatf1[-i,]
  }
}
for (i in 1:length(newdatf1)) {
  if (newdatf1[i] == 0) {
    newdatf1 <- newdatf1[-i]
  }
}
write.table(newdatf, file = "~/Downloads/Urbaniakcontrols.txt",row.names = FALSE, col.names = FALSE)
write.table(newdatf1, file = "~/Downloads/Urbaniakcases.txt",row.names = FALSE, col.names = FALSE)

#### Formatting files for input to the percent normalization python script ####
## the UrbaniakASV12202021upd_percentnormal.txt file:
# transpose this file so that the taxa/asvs are the columns and the SRR samples are rows
# Copy the taxa to a separate file
# first change the \t (tabs) to dashes (-)
# then change the \n (new lines) to \t (tabs)
# then copy this modified formatted taxa to the UrbaniakASV12202021upd_percentnormal.txt file
# and remove the " (quotation marks)
# in the UrbaniakASV12202021upd_percentnormal.txt file change the spaces ( ) to \t (tabs)
# go to the first taxa (at the top of the UrbaniakASV12202021upd_percentnormal.txt file) and add a \t (tab)

## the Urbaniakcontrols.txt file and the Urbaniakcases.txt file
# change the \n (new lines) to \t (tabs)
# remove the " (quotation marks)
# remove the last \t (tab) at the end of the file
# and add a \n (new line)

#### Combine TSS Normalized ASV tables ####
# View the percent normalization output:
# percnormal <- read.table("~/Downloads/Urbaniak_percentile_norm.txt")
# View(percnormal)
# sum(percnormal$Bacteria_NA_NA_NA_NA_NA.42)

hiekpernorm <- read.table("~/Downloads/Hieken_percentile_norm-proportional.txt")
urbpernorm <- read.table("~/Downloads/Urbaniak_percentile_norm-proportional.txt")
chanpernorm <- read.table("~/Downloads/Chan_percentile_norm-proportional.txt")

casv <- read.table("~/Downloads/Chan_05242022/ChanASV-TSS_upd.txt")
hasv <- read.table(file = "~/Downloads/HiekenASV12202021upd_percentnormalupd.txt")
uasv <- read.table("~/Downloads/Urbaniak_05252022/UrbASV-TSS_upd.txt")

colnames(casv)[1]
colnames(uasv)[1]
colnames(hasv)[1]

colnames(chanpernorm) <- colnames(casv)
colnames(urbpernorm) <- colnames(uasv)
colnames(hiekpernorm) <- colnames(hasv)

chanpernorm <- t(chanpernorm)
hiekpernorm <- t(hiekpernorm)
urbpernorm <- t(urbpernorm)

# edit these in text edit to include 
# write.table(chanpernorm,file = "~/Downloads/Chan_05242022/chantsspercnormforR.txt")
# write.table(hiekpernorm,file = "~/Downloads/hieken_05242022/hiekentsspercnormforR.txt")
# write.table(urbpernorm,file = "~/Downloads/Urbaniak_05252022/urbaniaktsspercnormforR.txt")

chanpercnorm <- read.table(file = "~/Downloads/Chan_05242022/chantsspercnormforR.txt", header = TRUE)
hiekpercnorm <- read.table(file = "~/Downloads/hieken_05242022/hiekentsspercnormforR.txt", header = TRUE)
urbpercnorm <- read.table(file = "~/Downloads/Urbaniak_05252022/urbaniaktsspercnormforR.txt", header = TRUE)

# # Hieken Taxa table
# Hiekentaxa <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021.txt',header=TRUE)
# Hiekentaxa <- data.frame(Hiekentaxa)
# colnames(Hiekentaxa)[1] <- "Taxonomy"
# # Hieken metadata table
# HMeta <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissues.txt")
# HMeta <- data.frame(HMeta)
# colnames(HMeta) <- c("sample_id","sample","sample_type")
# # Chan taxa table
# taxasilvacdat <- read.table("~/Documents/Chan10082021/ChaneditedASVstaxa/Channewtaxafrom12302021/ChanASVTaxaMetadModified-InputtoRfor-taxaabund/ChannewtaxafromR12302021_R.txt",header = TRUE)
# taxasilvacdat <- data.frame(taxasilvacdat)
# colnames(taxasilvacdat)[1] <- "Taxonomy"
# # Chan metadata table
# samdfcdat <- read.table("~/Documents/Chan10082021/ChanMetadatCombined_.txt",header = TRUE)
# samdfcdat <- data.frame(samdfcdat)
# colnames(samdfcdat) <- c("sample_id","sample","sample_type")
# # Urbaniak Taxa table
# taxasilvaudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile_R.txt",header = TRUE)
# taxasilvaudat <- data.frame(taxasilvaudat)
# # Urbaniak Metadata table
# samdfudat <- read.table(file= "~/Downloads/Urbaniak_05252022/urbaniakmeta.txt")
# samdfudat <- data.frame(samdfudat)
# ## Modify Urbaniak metadata table
# for (j in 1:length(rownames(samdfudat))) {
#   if (samdfudat[j,5] == "cancer"){
#     samdfudat[j,5] <- "Malignant"
#   }
# }
# for (j in 1:length(rownames(samdfudat))) {
#   if (samdfudat[j,5] == "benign"){
#     samdfudat[j,5] <- "Benign"
#   }
# }
# colnames(samdfudat) <- c("sample_id","sample","tissue","sample_name","sample_type")
# 
# 
HiekUrbASVcombined <- full_join(hiekpercnorm,urbpercnorm,by="ASV")
HiekUrbASVcombined <- data.frame(HiekUrbASVcombined)
HiekUrbChanASVcombined <- full_join(HiekUrbASVcombined,chanpercnorm,by="ASV")
## Modify & format the ASV tables
HiekUrbChanASVcombined1 <- HiekUrbChanASVcombined[,-1]
rownames(HiekUrbChanASVcombined1) <- HiekUrbChanASVcombined[,1]
HiekUrbChanASVcombined1 <- t(HiekUrbChanASVcombined1)
# # Combine the metadata tables
# HiekUrbMETADATAcombined <- full_join(HMeta,samdfudat,by = "sample_id")
# HiekUrbMETADATAcombined <- data.frame(HiekUrbMETADATAcombined)
# HiekUrbChanMETADATAcombined <- full_join(HiekUrbMETADATAcombined,samdfcdat,by="sample_id")
# ## Modify & format the metadata tables
# HiekUrbChanMETADATAcombined1 <- HiekUrbChanMETADATAcombined[,-1]
# rownames(HiekUrbChanMETADATAcombined1) <- HiekUrbChanMETADATAcombined[,1]
# HiekUrbChanMETADATAcombined1$StudyGroup <- "asv"
# for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
#   for (j in 1:length(HMeta$sample_id)) {
#     if (rownames(HiekUrbChanMETADATAcombined1)[i] == HMeta$sample_id[j]) {
#       HiekUrbChanMETADATAcombined1$StudyGroup[i] <- paste("Hiek_", HMeta$sample_type[j], sep = '')
#     }
#   }
# }
# 
# for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
#   for (j in 1:length(samdfudat$sample_id)) {
#     if (rownames(HiekUrbChanMETADATAcombined1)[i] == samdfudat$sample_id[j]) {
#       HiekUrbChanMETADATAcombined1$StudyGroup[i] <- paste("Urb_", samdfudat$sample_type[j], sep = '')
#     }
#   }
# }
# 
# for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
#   for (j in 1:length(samdfcdat$sample_id)) {
#     if (rownames(HiekUrbChanMETADATAcombined1)[i] == samdfcdat$sample_id[j]) {
#       if (samdfcdat$sample_type[j] == "Healthy_Control_Women") {
#         HiekUrbChanMETADATAcombined1$StudyGroup[i] <- "Cha_Healthy"
#       }
#       if (samdfcdat$sample_type[j] == "Women_with_a_History_of_Breast_Cancer") {
#         HiekUrbChanMETADATAcombined1$StudyGroup[i] <- "Cha_Cancer"
#       }
#     }
#   }
# }
# # Combine the taxa tables
# ## the combined taxa file is generated from python script editmetadataASVTaxa.py
# HiekUrbChanTAXAcombined <- read.table("~/Documents/HiekUrbChan_CombinedTaxa_fin.txt",header = TRUE)
# HiekUrbChanTAXAcombined1 <- HiekUrbChanTAXAcombined[,-1]
# rownames(HiekUrbChanTAXAcombined1) <- HiekUrbChanTAXAcombined[,1]
# HiekUrbChanTAXAcombined1 <- as.matrix(HiekUrbChanTAXAcombined1)
# 
# #### Replace the NA values with zeros ####
# # Replace NAs with the character zero ("zero") in the metadata table
# for (i in 1:length(rownames(HiekUrbChanMETADATAcombined1))) {
#   ft <- is.na(HiekUrbChanMETADATAcombined1[i,])
#   for (j in 1:length(ft)) {
#     if (ft[j] == TRUE) {
#       HiekUrbChanMETADATAcombined1[i,j] <- "zero"
#     }
#   }
# }
# Replace NAs with the number zero (0) in the ASV table
for (i in 1:length(rownames(HiekUrbChanASVcombined1))){
  ftasv <- is.na(HiekUrbChanASVcombined1[i,])
  for (j in 1:length(ftasv)) {
    if (ftasv[j] == TRUE) {
      HiekUrbChanASVcombined1[i,j] <- 0
    }
  }
}


write.table(t(OrighucASVcombined1),file = "~/Downloads/originalcombinedasvtable.txt",sep = '\t')
#### Prop Abund: TSS Normalized combined studies ####
# can change HiekUrbChanASVcombined1 to OrighucASVcombined1 to do analyses with originalasv combined table
psjoined <- phyloseq(otu_table(HiekUrbChanASVcombined1, taxa_are_rows=FALSE),
                     sample_data(HiekUrbChanMETADATAcombined1), # originally sample_data
                     tax_table(HiekUrbChanTAXAcombined2))
## do not have mock samples in our data so do not need this
# psjoined <- prune_samples(sample_names(psjoined) != "Mock", psjoined)
dnajoined <- Biostrings::DNAStringSet(taxa_names(psjoined))
names(dnajoined) <- taxa_names(psjoined)
psjoined <- merge_phyloseq(psjoined, dnajoined)
taxa_names(psjoined) <- paste0("ASV", seq(ntaxa(psjoined)))
psjoined

## Bray NMDS ordination plot
ps.propjoined <- transform_sample_counts(psallseqsjoined, function(otu) otu/sum(otu))
ord.nmds.brayjoined <- ordinate(ps.propjoined, method="NMDS", distance="bray")
plot_ordination(psjoined, ord.nmds.brayjoined, color="Tissues", title="Bray NMDS")

allseqsjoined <- names(sort(taxa_sums(psjoined), decreasing=TRUE)) [1:800]
psallseqsjoined <- transform_sample_counts(psjoined, function(OTU) OTU/sum(OTU))
psallseqsjoined <- prune_taxa(allseqsjoined, psallseqsjoined)
psjoinedg <- tax_glom(psallseqsjoined, "Phylum")

pscombinejoined <- names(sort(taxa_sums(psjoined), decreasing=TRUE)) [1:800]
psglomjoined <- transform_sample_counts(psjoined, function(OTU) OTU/sum(OTU))
psglomjoined <- prune_taxa(pscombinejoined, psglomjoined)

#### TSS normaliz ASV table cmdscale plot ####
## tax_glom function to collapse at the taxa level
psjoinedg <- tax_glom(psglomjoined, "Genus")
ps0joinedg <- transform_sample_counts(psjoinedg, function(x) x / sum(x))
ps1joinedg <- merge_samples(ps0joinedg, "Tissues") #yields NAs in columns of samdfjoined within ps1joinedg
ps2joinedg <- transform_sample_counts(ps1joinedg, function(x) x / sum(x))
pjoinedg <- plot_bar(ps2joinedg, fill= "Genus")
pjoinedg + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))

library(vegan)
library(ggplot2)
library(ALDEx2)
library(compositions)
library(ape)
#hiekurbchanclr <- clr(HiekUrbChanASVcombined1)
#D.BCHiek <- as.matrix(vegdist(hiekurbchanclr, method="robust.aitchison"))
# factstudygroup <- factor(HiekUrbChanMETADATAcombined1$StudyGroup)
# permanovainfo <- list(HiekUrbChanASVcombined1,D.BCHiek,"bray")
# PERMANOVA(permanovainfo,factstudygroup,nperm = 1000,seed = 1234)

#### Beta diversity with TSS Normalized ASV table ####
D.BCHiek <- as.matrix(vegdist(HiekUrbChanASVcombined1, method="robust.aitchison")) # also used method="bray"
cmdBCurtis <- pcoa(D.BCHiek) #cmdscale(D.BCHiek)
cmdBCurtis <- data.frame(cmdBCurtis[["vectors"]])
cmdBCurtis$meta <- 0
for (i in 1:length(HiekUrbChanMETADATAcombined1[,1])){
  for (j in 1:length(rownames(cmdBCurtis))) {
    if (rownames(cmdBCurtis)[j] == rownames(HiekUrbChanMETADATAcombined1)[i]) {
      cmdBCurtis$meta[j] <- HiekUrbChanMETADATAcombined1[i,9] # StudyGroup with tumor vs healthy samples
    }
  }
}
for (i in 1:length(HiekUrbChanMETADATAcombined1[,1])){
  for (j in 1:length(rownames(cmdBCurtis))) {
    if (rownames(cmdBCurtis)[j] == rownames(HiekUrbChanMETADATAcombined1)[i]) {
      cmdBCurtis$meta[j] <- HiekUrbChanMETADATAcombined1[i,10] # Tissues with breast vs skin vs NAF vs etc... samples
    }
  }
}

# text(x,y,cmduni$meta,cex=0.8)
ggplot(cmdBCurtis,aes(x=Axis.1,y=Axis.2,color= meta)) +
  #geom_dotplot(y=cmduni[,2],binwidth = 0.004)
  geom_point() +
  stat_ellipse() 


#### beta diversity TSS Normalized combined studies ####
library(GUniFrac)
library(MiRKAT)
library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
library(vegan)

sequenceskdat <- getSequences(HiekUrbChanASVcombined1)
names(sequenceskdat) <- sequenceskdat

## Run Sequence Alignment (MSA) using DECIPHER
alignmentkdat <- AlignSeqs(DNAStringSet(sequenceskdat), anchor=NA)

## Change sequence alignment output into a phyDat structure
phang.alignkdat <- phyDat(as(alignmentkdat, "matrix"), type="DNA")

## Create distance matrix
dmkdat <- dist.ml(phang.alignkdat)
UPGMAtreekdat <- upgma(dmkdat)
# write.tree(UPGMAtreekdat, file = "~/Downloads/HiekUrbChanupgmatree.nwk", append = FALSE,
          # digits = 10, tree.names = FALSE)

#### weighted and unweighted UniFrac ####
library(GUniFrac)
library(vegan)
library(ggplot2)
library(ALDEx2)
library(compositions)
library(ape)
library(phyloseq)

hiekurbchantree <- read_tree("~/Downloads/HiekUrbChanupgmatree.nwk")
ASVraref <- data.frame(HiekUrbChanASVcombined1)
ASVtabunifracsHiek <- GUniFrac(ASVraref, hiekurbchantree)$unifracs # , alpha=c(0, 0.5, 1)
D.weightedHiek <- ASVtabunifracsHiek[,,"d_1"]
D.unweightedHiek <- ASVtabunifracsHiek[,,"d_UW"]
# D.BCHiek <- as.matrix(vegdist(ASVtab2Hiek , method="bray"))
KweightedHiek <- D2K(D.weightedHiek)
KunweightedHiek <- D2K(D.unweightedHiek)
# K.BCHiek <- D2K(D.BCHiek)
# MiRKAT(y = ASVtab2Hiek, Ks = KunweightedHiek, out_type = "D", method = "davies")
# ggplot(data= D.unweightedHiek,mapping = HMetaHiek)

UniFrac()

# cmdscale plot
cmdunweunif <- cmdscale(D.unweightedHiek)
cmdweighunif <- cmdscale(D.weightedHiek)
#View(cmdunweunif)
cmduni <- data.frame(cmdunweunif)
cmduniweigh <- data.frame(cmdweighunif)
# x <- cmduni[,1]
# y <- cmduni[,2]
# plot(x,y,asp=1)
# x <- cmduniweigh[,1]
# y <- cmduniweigh[,2]
# plot(x,y,asp=1)
cmduniweigh$meta <- 0
cmduni$meta <- 0

# for unweighted unifrac breast vs skin
# HMetaHiek1 <- read.table(file="~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissuescopy.txt")
for (i in 1:length(HMetaHiek[,1])){
  for (j in 1:length(rownames(cmduni))) {
    if (rownames(cmduni)[j] == rownames(HMetaHiek)[i]) {
      cmduni$meta[j] <- HMetaHiek[i,1]
    }
  }
}
# text(x,y,cmduni$meta,cex=0.8)
ggplot(cmduni,aes(x=X1,y=X2,color= meta)) +
  #geom_dotplot(y=cmduni[,2],binwidth = 0.004)
  geom_point() +
  stat_ellipse() 

# for weighted unifrac breast vs skin
for (i in 1:length(HMetaHiek[,1])){
  for (j in 1:length(rownames(cmduniweigh))) {
    if (rownames(cmduni)[j] == rownames(HMetaHiek)[i]) {
      cmduniweigh$meta[j] <- HMetaHiek[i,1]
    }
  }
}
# text(x,y,cmduniweigh$meta,cex=0.8)
ggplot(cmduniweigh,aes(x=X1,y=X2,color= meta)) +
  #geom_dotplot(y=cmduni[,2],binwidth = 0.004)
  geom_point() +
  stat_ellipse() 


#### Attempting Rarefaction with original raw otu tables ####
# Rarefaction is not performed
psjoined <- phyloseq(otu_table(HiekUrbChanASVcombined1, taxa_are_rows=FALSE),
                     sample_data(HiekUrbChanMETADATAcombined1), # originally sample_data
                     tax_table(HiekUrbChanTAXAcombined1))
# 6431OTUs were removed because they are no longer 
# present in any sample after random subsampling
ASVrarefied <- rarefy_even_depth(psjoined,rngseed = 1274) #set.seed: 123456789
## do not have mock samples in our data so do not need this
# psjoined <- prune_samples(sample_names(psjoined) != "Mock", psjoined)
dnajoined <- Biostrings::DNAStringSet(taxa_names(ASVrarefied))
names(dnajoined) <- taxa_names(ASVrarefied)
ASVrarefied <- merge_phyloseq(ASVrarefied, dnajoined)
taxa_names(ASVrarefied) <- paste0("ASV", seq(ntaxa(ASVrarefied)))
ASVrarefied

## Bray NMDS ordination plot
ps.propjoined <- transform_sample_counts(ASVrarefied, function(otu) otu/sum(otu))
ord.nmds.brayjoined <- ordinate(ps.propjoined, method="NMDS", distance="bray")
plot_ordination(ASVrarefied, ord.nmds.brayjoined, color="StudyGroup", title="Bray NMDS")

allseqsjoined <- names(sort(taxa_sums(ASVrarefied), decreasing=TRUE)) # [1:700]
psallseqsjoined <- transform_sample_counts(ASVrarefied, function(OTU) OTU/sum(OTU))
psallseqsjoined <- prune_taxa(allseqsjoined, psallseqsjoined)
psjoinedg <- tax_glom(psallseqsjoined, "Phylum")

pscombinejoined <- names(sort(taxa_sums(ASVrarefied), decreasing=TRUE)) # [1:700]
psglomjoined <- transform_sample_counts(ASVrarefied, function(OTU) OTU/sum(OTU))
psglomjoined <- prune_taxa(pscombinejoined, psglomjoined)






