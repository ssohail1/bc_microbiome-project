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
#storeasv <- vector()
#for (i in 1:length(ASVcolsums)){ #length(ASVtabprev2Hiek[1,])
  #sumcol <- sum(ASVtabprev2Hiek[,i])
  #print(sumcol)
  #for (j in 1:length(ASVtabprev2Hiek[,1])) { 
  #if (max(ASVtabprev2Hiek[,i]) <= 0.002) {
    #storeasv <- c(storeasv,colnames(ASVtabprev2Hiek)[i])
    #print(TRUE)
    #print(max(ASVtabprev2Hiek[,i]))
  #}
#}

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
library(MiRKAT)
for (i in 1:length(Hiekmet$env_feature)) {
  Hiekmet[i,3] <- as.double(Hiekmet[i,3])
}
Hiekdouble <- Hiekmet[,3]

set.seed(12345)
meerkatsingleHiekunweigUniFrac <- MiRKAT(y= Hiekdouble, Ks = KunweightedHiek, out_type = "D", method = "permutation")
meerkatsingleHiekweighUniFrac <- MiRKAT(y= Hiekdouble, Ks = KweightedHiek, out_type = "D", method = "permutation")
meerkatsingleHiekunweigUniFrac
meerkatsingleHiekweighUniFrac



####### Heatmap - alpha diversity #######
library(phyloseq)
library(reshape2)
library(ggplot2)
library(plyr); library(dplyr)
library(RColorBrewer)
                                         
# load ASV table
ASVtab <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
ASVtab <- read.table('~/Downloads/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
ASVtab <- data.frame(ASVtab)
ASVtab2Hiek <- ASVtab[,-1]
rownames(ASVtab2Hiek) <- ASVtab[,1]
ASVtab2Hiek <- t(ASVtab2Hiek)
seqtab.nochimhkdat1 <- data.frame(ASVtab2Hiek)

# load metadata table
HMetaHiek1 <- read.table(file="~/Downloads/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissuescopy.txt")
HMeta <- data.frame(HMetaHiek1)
HMetaHiek <- HMeta[,-1]
rownames(HMetaHiek) <- HMeta[,1]
colnames(HMetaHiek) <- c("tissue","tumor")
samdfhkdat1 <- data.frame(HMetaHiek)


# load taxa table
hiektaxa <- read.table('~/Downloads/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021.txt',header=TRUE)
hiektaxa <- data.frame(hiektaxa)
hiektaxa1 <- hiektaxa[,-1]
rownames(hiektaxa1) <- hiektaxa[,1]
taxasilvahkdat1 <- as.matrix(hiektaxa1)

# create phyloseq object
pshktisdat <- phyloseq(otu_table(seqtab.nochimhkdat1, taxa_are_rows=FALSE),
                       sample_data(samdfhkdat1), # originally sample_data
                       tax_table(taxasilvahkdat1))
ASVrarefied <- rarefy_even_depth(pshktisdat,rngseed = 12345) #set.seed: 123456789
dnahktisdat <- Biostrings::DNAStringSet(taxa_names(ASVrarefied))
names(dnahktisdat) <- taxa_names(ASVrarefied)
ASVrarefied <- merge_phyloseq(ASVrarefied, dnahktisdat)
taxa_names(ASVrarefied) <- paste0("ASV", seq(ntaxa(ASVrarefied)))
ASVrarefied

# rarefy ASV table
dataASV <- data.frame(ASVrarefied@otu_table)
datatasv <- data.frame(t(dataASV))
colnames(datatasv) <- samdfhkdat1$tissue
newdatf <- data.frame(matrix(ncol = length(colnames(datatasv)), nrow = length(rownames(datatasv))))
count1 <- 0
# add breast tissue samples to newdatf dataframe
for (i in 1:length(colnames(datatasv))) {
  if (colnames(datatasv)[i] == "Breast") {
    newdatf[,i] <- datatasv[,i]
  }
  #colnames(datatasv)[i] <- paste(substr(colnames(datatasv)[i],1,1), i, sep='')
  if (colnames(datatasv)[i] == "Skin_Tissue") {
    count1 <- count1 + 1
  }
}
# remove NA values
for (i in 1:length(colnames(newdatf))) {
  if (is.na(newdatf[1,i]) == TRUE) {
    newdatf <- newdatf[,-i]
  }
}
# all are breast samples so replace breast with "B"
for (i in 1:length(colnames(newdatf))) {
  colnames(newdatf)[i] <- paste("B", i, sep='')
}


# newdatf1 dataframe has the skin_tissue samples
newdatf1 <- data.frame(matrix(ncol = length(colnames(datatasv)), nrow = length(rownames(datatasv))))
# count1 <- 0
for (i in 1:length(colnames(datatasv))) {
  if (colnames(datatasv)[i] == "Skin_Tissue") {
    newdatf1[,i] <- datatasv[,i]
  }
  #colnames(datatasv)[i] <- paste(substr(colnames(datatasv)[i],1,1), i, sep='')
  # if (colnames(datatasv)[i] == "Skin_Tissue") {
  #   count1 <- count1 + 1
  # }
}
# remove NA values
for (i in 1:length(colnames(newdatf1))) {
  if (is.na(newdatf1[1,i]) == TRUE) {
    newdatf1 <- newdatf1[,-i]
  }
}
# all are skin tissue samples so replace skin_tissue with "S"
for (i in 1:length(colnames(newdatf1))) {
  colnames(newdatf1)[i] <- paste("S", i, sep='')
}


# add the skin tissue samples ASV table information to newdatf
## newdatf = all breast samples
## newdatf1 = all skin tissue samples
newdatf[,29:47] <- data.frame(matrix(ncol = count1, nrow = length(rownames(datatasv))))
newdatf[,29:47] <- newdatf1[,1:19]

# add appropriate column names 
for (i in 29:47) {
  colnames(newdatf)[i] <- paste("Skin", sep='')
}
for (i in 1:28) {
  colnames(newdatf)[i] <- paste("Breast", sep='')
}
# add ASV# rownames
for (i in 1:length(rownames(newdatf))) {
  rownames(newdatf)[i] <- paste("ASV",i, sep = '')
}

# modifying dataframes to plot heatmap through base R and ggplot2
newdatf1 <- t(data.frame(newdatf))
meltnewdf <- melt(newdatf1)
ASVdf <- data.frame(meltnewdf)

dataASV1 <- t(datatasv)
ASVmelt <- melt(dataASV1)
ASVdf <- data.frame(ASVmelt)

meltnewdf %>%
  ggplot(aes(Var1, Var2, fill= value)) + 
  geom_tile()

ASVmelt %>%
  ggplot(aes(Var1, Var2, fill= value)) + 
  geom_tile()

adivdist <- 
datatmat <- as.matrix(datatasv)
my_group <- as.numeric(as.factor(substr(rownames(datatmat), 1 , 1)))
colSide <- brewer.pal(9, "Set1")[my_group]
colMain <- colorRampPalette(brewer.pal(8, "Blues"))(47)
colMain1 <- colorRampPalette(brewer.pal(8, "Blues"))(47)
for (i in 1:length(colnames(datatmat))) {
  if (colnames(datatmat)[i] == "Breast") {
    colMain1[i] <- colMain[44]
  }
  if (colnames(datatmat)[i] == "Skin_Tissue") {
    colMain1[i] <- colMain[21]
  }
}
newdat <- as.matrix(newdatf)
heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(x = newdat,ColSideColors = colMain1)
