#### Input ####

# ASV table is already filtered using prevalence = 10%
ASVtabprev2 <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabmodifHieken03152022__.txt', header = TRUE)
ASVtabprev2 <- data.frame(ASVtabprev2)
ASVtabprev2Hiek <- ASVtabprev2[,-1]
rownames(ASVtabprev2Hiek) <- ASVtabprev2[,1]
ASVtabprev2Hiek <- t(ASVtabprev2Hiek)

# Hieken metadata
HMeta <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissues.txt")
HMeta <- data.frame(HMeta)
HMetaHiek <- HMeta[,-1]
rownames(HMetaHiek) <- HMeta[,1]
colnames(HMetaHiek) <- c("tissue","tumor")

Hiekmet <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/Hmetadata_mirkat.txt',header=TRUE)
HMetaHiek_1 <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/Hmetadata.txt")
Hiekmet$additional <- HMetaHiek_1$V2

# Hieken taxa
hiektaxa <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021.txt',header=TRUE)
hiektaxa <- data.frame(hiektaxa)
hiektaxa1 <- hiektaxa[,-1]
rownames(hiektaxa1) <- hiektaxa[,1]

# load libraries
library(phyloseq)
library(GUniFrac)
library(MicrobiomeStat)
library(ape)

#### TSS ####
ASVcolsums <- colSums(ASVtabprev2Hiek)
for (i in 1:length(ASVcolsums)){ #length(ASVtabprev2Hiek[1,])
  #sumcol <- sum(ASVtabprev2Hiek[,i])
  #print(sumcol)
  for (j in 1:length(ASVtabprev2Hiek[,1])) { 
    print(ASVtabprev2Hiek[j,i])
    ASVtabprev2Hiek[j,i] <- round((ASVtabprev2Hiek[j,i])/(ASVcolsums[i]),digits=10)
  }
}

View(ASVtabprev2Hiek)

# checking if column sums of ASV table is 1
ASVcolsumstaxaabund <- colSums(ASVtabprev2Hiek)

# filtering samples that have a max count less than 0.002 - in this case there were no samples with a max count less than 0.002.
storeasv <- vector()
for (i in 1:length(ASVcolsums)){ #length(ASVtabprev2Hiek[1,])
  #sumcol <- sum(ASVtabprev2Hiek[,i])
  #print(sumcol)
  #for (j in 1:length(ASVtabprev2Hiek[,1])) { 
  if (max(ASVtabprev2Hiek[,i]) <= 0.002) {
    #storeasv <- c(storeasv,colnames(ASVtabprev2Hiek)[i])
    #print(TRUE)
    #print(max(ASVtabprev2Hiek[,i]))
  }
}

# sqrt-transform ASV table
for (i in 1:length(ASVtabprev2Hiek[1,])){
  ASVtabprev2Hiek[,i] <- sqrt(ASVtabprev2Hiek[,i])
}
ASVtabprev2Hiek1 <- t(ASVtabprev2Hiek)

# removing buccal and skin swab, skin tissue, atypia, and DCIS samples
# to compare benign and invasive cancer disease states in breast tissue
storing <- c()

# first remove those samples from the metadata
# and store the sample names to storing
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,1] == "Buccal_Cells") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,1] == "Skin_Swab") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,1] == "Skin_Tissue") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,4] == "Atypia") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,4] == "DCIS") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
# loop through the storing list to remove those sample names from the ASV table
for (i in 1:length(rownames(ASVtabprev2Hiek))) {
  for (j in 1:length(storing)) {
    if (rownames(ASVtabprev2Hiek)[i] == storing[j]){
      ASVtabprev2Hiek <- ASVtabprev2Hiek[-i,]
    }
  }
}
ASVtabprev2Hiek1 <- t(ASVtabprev2Hiek)

# to test with SRR files
linda.obj <- linda(ASVtabprev2Hiek, Hiekmet, formula = '~additional',
                   feature.dat.type = 'proportion',
                   prev.filter = 0.2,
                   p.adj.method = "none", alpha = 0.05)
# volcano plot
linda.plot(linda.obj, 'tumor',
           alpha = 0.05, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)

# to test with ASVs and Taxa
linda.obj1 <- linda(ASVtabprev2Hiek1, Hiekmet, formula = '~additional',
                    feature.dat.type = 'proportion',
                    prev.filter = 0.2,
                    p.adj.method = "none", alpha = 0.05)
# volcano plot
linda.plot(linda.obj1, 'tumor',
           alpha = 0.05, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)


lindaobjpvaladj <- data.frame(linda.obj1[["output"]][["additionalInvCa"]])
lindaobjpvaladj1 <- data.frame(lindaobjpvaladj$padj)
rownames(lindaobjpvaladj1) <- rownames(lindaobjpvaladj)

count <- 0
S <- list()
st <- list()

S <- data.frame("asv",rep(0,length(lindaobjpvaladj1$lindaobjpvaladj.padj)))
S$asv <- rep(0,length(lindaobjpvaladj1$lindaobjpvaladj.padj))
colnames(S) <- c("asv","taxagenus","extra")

for (i in 1:length(lindaobjpvaladj1$lindaobjpvaladj.padj)) {
  if (lindaobjpvaladj1$lindaobjpvaladj.padj[i] <= 0.05){
    S$asv[i] <- rownames(lindaobjpvaladj1)[i]
    count <- count + 1
    st <- c(st,lindaobjpvaladj1$lindaobjpvaladj.padj[i])
  }
}

# run this until all 'asv' rows are removed
for (i in 1:length(lindaobjpvaladj1$lindaobjpvaladj.padj)) {
  if (S$asv[i] == 'asv'){
    S <- S[-i,]
  }
}
S <- data.frame(S)
for (i in 1:length(S$asv)) {
  for (j in 1:length(rownames(hiektaxa1))) {
    if (S$asv[i] == rownames(hiektaxa1)[j]){
      S$taxagenus[i] <- hiektaxa1$Genus[j]
      count <- count + 1
    } 
  }
}

for (i in 1:length(S$asv)) {
  for (j in 1:length(rownames(hiektaxa1))) {
    if (is.na(S$taxagenus[i])) {
      if (S$asv[i] == rownames(hiektaxa1)[j]){
        S$taxagenus[i] <- hiektaxa1$Order[j]
        count <- count + 1
      } 
    }
  }
}

for (j in 1:length(S$asv)) {
  #for (i in 1:length(lindahiek$pvalue)) {
  #for (j in 1:length(S$asv)) {
  if (st[j] <= 0.05){
    S$NA.[j] <- st[j]
  }
}
colnames(S) <- c("ASV","Family-Genus","P-value")
colnames(S) <- c("ASV","Family-Genus","blank", "p-value")
library(phyloseq)
hiektaxa2 <- as.matrix(hiektaxa1)
psudat <- phyloseq(otu_table(ASVtabprev2Hiek, taxa_are_rows=FALSE),
                   sample_data(Hiekmet), # originally sample_data
                   tax_table(hiektaxa2))
dnaudat <- Biostrings::DNAStringSet(taxa_names(psudat))
names(dnaudat) <- taxa_names(psudat)
psudat <- merge_phyloseq(psudat, dnaudat)
# taxa_names(psudat) <- paste0("ASV", seq(ntaxa(psudat)))
psudat
# ps.propudat <- transform_sample_counts(psudat, function(otu) otu/sum(otu))
# ord.nmds.brayudat <- ordinate(ps.propudat, method="NMDS", distance="bray", color="Type")
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE)) 
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)

count <- vector()
library(metagMisc)
propASVdata <- data.frame(psudat@otu_table)
propASVdata <- t(propASVdata)

newdatf <- data.frame(matrix(ncol = length(colnames(propASVdata)), nrow = length(S$ASV)))
newdatf <- t(newdatf)
for (i in 1:length(S[,1])) {
  for (j in 1:length(rownames(propASVdata))) {
    if (S$ASV[i] == rownames(propASVdata)[j]) {
      print(TRUE)
      newdatf[,i] <- propASVdata[j,]
    }
  }
}
rownames(newdatf) <- colnames(propASVdata)
newdatf <- t(newdatf)
newdatf1 <- data.frame(newdatf)
newdatf1$types <- S$`Family-Genus`
# newdatasvl <- list(newdatasv)
# newdatasv1 <- as.data.frame(newdatasvl[[1]])
# newdatasv1$types <- S$`Family-Genus`
# newdatasv1 <- t(data.frame(newdatasv))
# newdatasv1$types <- data.frame(samdfhkdat1$hist)
library(ggplot2)
library(reshape2)

for (i in 1:length(colnames(newdatf1))) {
  for (j in 1:length(storing)) {
    if (colnames(newdatf1)[i] == storing[j]) {
      newdatf1 <- newdatf1[,-i]
    } 
  }
}

colnames(newdatf1)[1:23] <- Hiekmet$additional

colnames(newdatf1)[24] <- "taxa"
newdatf1 <- data.frame(newdatf1)
mndat1 <- melt(newdatf1)
mndat1 <- data.frame(mndat1)

mndat1$variable <- Hiekmet$additional
mndat1 <- data.frame(mndat1)

ggplot(mndat1,aes(x=variable,y=value,fill=taxa)) +
  geom_col() +
  facet_wrap(~taxa)

#### unweighted UniFrac plot ####
# for BBD vs InvCan
ASVtab <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
ASVtab <- data.frame(ASVtab)
ASVtab2Hiek <- ASVtab[,-1]
rownames(ASVtab2Hiek) <- ASVtab[,1]
ASVtab2Hiek <- t(ASVtab2Hiek)
seqtab.nochimhkdat1 <- data.frame(ASVtab2Hiek)

Hiekmet <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/Hmetadata_mirkat.txt',header=TRUE)
HMetaHiek_1 <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/Hmetadata.txt")
Hiekmet$additional <- HMetaHiek_1$V2

storing <- c()
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,1] == "Buccal_Cells") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,1] == "Skin_Swab") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,1] == "Skin_Tissue") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,4] == "Atypia") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(Hiekmet))) {
  if (Hiekmet[i,4] == "DCIS") {
    storing <- c(storing, rownames(Hiekmet)[i])
    Hiekmet <- Hiekmet[-i,]
  }
}
for (i in 1:length(rownames(seqtab.nochimhkdat1))) {
  for (j in 1:length(storing)) {
    if (rownames(seqtab.nochimhkdat1)[i] == storing[j]){
      seqtab.nochimhkdat1 <- seqtab.nochimhkdat1[-i,]
    }
  }
}

hiektaxa <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021.txt',header=TRUE)
hiektaxa <- data.frame(hiektaxa)
hiektaxa1 <- hiektaxa[,-1]
rownames(hiektaxa1) <- hiektaxa[,1]
hiektaxa2 <- as.matrix(hiektaxa1)

psudat <- phyloseq(otu_table(seqtab.nochimhkdat1, taxa_are_rows=FALSE),
                   sample_data(Hiekmet), # originally sample_data
                   tax_table(hiektaxa2))
ASVrarefied <- rarefy_even_depth(psudat,rngseed = 123456789) #set.seed: 123456789

seqtabhk <- data.frame(ASVrarefied@otu_table)
Hiektree <- read.tree(file = "~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/upgmahieken.nwk")
ASVtabunifracsHiek <- GUniFrac(seqtabhk, Hiektree, alpha=c(0, 0.5, 1))$unifracs
D.weightedHiek <- ASVtabunifracsHiek[,,"d_1"]
D.unweightedHiek <- ASVtabunifracsHiek[,,"d_UW"]
# D.BCHiek <- as.matrix(vegdist(ASVtab2Hiek , method="bray"))
KweightedHiek <- D2K(D.weightedHiek)
KunweightedHiek <- D2K(D.unweightedHiek)

cmdunweunif <- cmdscale(D.unweightedHiek)
cmdweighunif <- cmdscale(D.weightedHiek)
#View(cmdunweunif)
cmduni <- data.frame(cmdunweunif)
cmduniweigh <- data.frame(cmdweighunif)
cmduniweigh$meta <- 0
cmduni$meta <- 0

for (i in 1:length(Hiekmet[,1])){
  for (j in 1:length(rownames(cmduni))) {
    if (rownames(cmduni)[j] == rownames(Hiekmet)[i]) {
      print(rownames(cmduni)[j] == rownames(Hiekmet)[i])
      cmduni$meta[j] <- Hiekmet[i,4]
    }
  }
}
ggplot(cmduni,aes(x=X1,y=X2,color= meta)) +
  #geom_dotplot(y=cmduni[,2],binwidth = 0.004)
  geom_point() +
  stat_ellipse() +
  labs(title = "Unweighted UniFrac")

                                         

#### MiRKAT ####
for (i in 1:length(Hiekmet$env_feature)) {
  Hiekmet[i,3] <- as.double(Hiekmet[i,3])
}
Hiekdouble <- Hiekmet[,3]

meerkatsingleHiekunweigUniFrac <- MiRKAT(y= Hiekdouble, Ks = KunweightedHiek, out_type = "D", method = "permutation")
meerkatsingleHiekweighUniFrac <- MiRKAT(y= Hiekdouble, Ks = KweightedHiek, out_type = "D", method = "permutation")
meerkatsingleHiekunweigUniFrac
meerkatsingleHiekweighUniFrac


