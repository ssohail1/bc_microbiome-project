# library(devtools)
# install_github("lme4/lme4",dependencies=TRUE)
# install.packages("PERMANOVA)
library(lme4)
library(PERMANOVA)
library(GUniFrac)
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)
comm <- t(throat.otu.tab)
meta.dat <- throat.meta
comm.p <- t(t(comm) / colSums(comm))
linda.obj <- linda(comm.p, meta.dat, formula = '~SmokingStatus+Sex',
                   feature.dat.type = 'proportion',
                   prev.filter = 0.2, is.winsor = TRUE, outlier.pct = 0.03,
                   p.adj.method = "BH", alpha = 0.1
)
linda.plot(linda.obj, c('SmokingStatusSmoker', 'Sexmale'),
           titles = c('Smoke: n v.s. y', 'Sex: female v.s. male'), alpha = 0.1, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)
ASVtab <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
ASVtab <- data.frame(ASVtab)
#ASVtab <- t(ASVtab)
ASVtab2Hiek <- ASVtab[,-1]
rownames(ASVtab2Hiek) <- ASVtab[,1]
ASVtab2Hiek <- t(ASVtab2Hiek)
HMeta <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissues.txt")
HMeta <- data.frame(HMeta)
HMetaHiek <- HMeta[,-1]
rownames(HMetaHiek) <- HMeta[,1]
colnames(HMetaHiek) <- c("tissue","tumor")

# Filter with prevalence = 10% and taxa < 0.2% then square root transform
ASVtabprev <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabwithprevalencepercentsHieken03152022.txt',header = TRUE)
ASVtabprev <- data.frame(ASVtabprev)
#ASVtabprev <- t(ASVtabprev)
hiektaxa <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021.txt',header=TRUE)
hiektaxa <- data.frame(hiektaxa)
hiektaxa1 <- hiektaxa[,-1]
rownames(hiektaxa1) <- hiektaxa[,1]


storeASVs <- vector()
for (i in 1:length(ASVtabprev$percentofzeros)) {
  if (ASVtabprev$percentofzeros[i] >= 90) {
    storeASVs <- c(storeASVs,ASVtabprev[i,])
    ASVtabprev <- ASVtabprev[-i,]
  }
}

#write.table(ASVtabprev,"~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabmodifHieken03152022.txt",row.names = FALSE)
# edit this file written out to folder - edit in Excel to remove countzeros and percentofzeros
ASVtabprev2 <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabmodifHieken03152022__.txt', header = TRUE)
ASVtabprev2 <- data.frame(ASVtabprev2)
ASVtabprev2Hiek <- ASVtabprev2[,-1]
rownames(ASVtabprev2Hiek) <- ASVtabprev2[,1]
ASVtabprev2Hiek <- t(ASVtabprev2Hiek)
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
ASVcolsumstaxaabund <- colSums(ASVtabprev2Hiek)
storeasv <- vector()
for (i in 1:length(ASVcolsums)){ #length(ASVtabprev2Hiek[1,])
  #sumcol <- sum(ASVtabprev2Hiek[,i])
  #print(sumcol)
  #for (j in 1:length(ASVtabprev2Hiek[,1])) { 
  if (max(ASVtabprev2Hiek[,i]) <= 0.002) {
    storeasv <- c(storeasv,colnames(ASVtabprev2Hiek)[i])
    print(TRUE)
    print(max(ASVtabprev2Hiek[,i]))
  }
}

for (i in 1:length(ASVtabprev2Hiek[1,])){
  ASVtabprev2Hiek[,i] <- sqrt(ASVtabprev2Hiek[,i])
}
ASVtabprev2Hiek1 <- t(ASVtabprev2Hiek)
# For SRR files
linda.obj <- linda(ASVtabprev2Hiek, HMetaHiek, formula = '~tumor',
                   feature.dat.type = 'proportion',
                   prev.filter = 0.2,
                   p.adj.method = "none", alpha = 0.05)
linda.plot(linda.obj, 'tumor',
           alpha = 0.05, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)
# For ASVs and Taxa
linda.obj1 <- linda(ASVtabprev2Hiek1, HMetaHiek, formula = '~tumor',
                   feature.dat.type = 'proportion',
                   prev.filter = 0.2,
                   p.adj.method = "none", alpha = 0.05)
linda.plot(linda.obj1, 'tumor',
           alpha = 0.05, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)

lindahiek <- read.table("~/Documents/Hieken10082021/Hiekenlindaobjpvalues.txt",header = TRUE)
lindahiek <- data.frame(lindahiek)
lindaobjupdSRR <- data.frame(linda.obj[["output"]][["tumorMalignant"]]) # p-vals not sig

count <- 0
S <- list()
st <- list()
#linda.obj$output$tumorMalignant$pvalue
S <- data.frame("asv",rep(0,length(lindahiek$pvalue)))
S$asv <- rep(0,length(lindahiek$pvalue))
colnames(S) <- c("asv","taxagenus")

for (i in 1:length(lindahiek$pvalue)) {
  if (lindahiek$pvalue[i] <= 0.05){
    S$asv[i] <- rownames(lindahiek)[i]
    count <- count + 1
    st <- c(st,lindahiek$pvalue[i])
  }
}

for (i in 1:length(lindahiek$pvalue)) {
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
for (i in 1:2) {
  for (j in 1:length(rownames(hiektaxa1))) {
    if (S$asv[i] == rownames(hiektaxa1)[j]){
      S$taxagenus[i] <- hiektaxa1$Family[j]
      count <- count + 1
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
library(phyloseq)
hiektaxa2 <- as.matrix(hiektaxa1)
psudat <- phyloseq(otu_table(ASVtab2Hiek, taxa_are_rows=FALSE),
                   sample_data(HMetaHiek), # originally sample_data
                   tax_table(hiektaxa2))
psudat <- prune_samples(sample_names(psudat) != "Mock", psudat)
dnaudat <- Biostrings::DNAStringSet(taxa_names(psudat))
names(dnaudat) <- taxa_names(psudat)
psudat <- merge_phyloseq(psudat, dnaudat)
taxa_names(psudat) <- paste0("ASV", seq(ntaxa(psudat)))
psudat
# ps.propudat <- transform_sample_counts(psudat, function(otu) otu/sum(otu))
# ord.nmds.brayudat <- ordinate(ps.propudat, method="NMDS", distance="bray", color="Type")
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE)) 
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)
# psudatg <- tax_glom(ps.top100udat, "Genus")
# ps0udatg <- transform_sample_counts(psudatg, function(x) x / sum(x))
# ps1udatg <- merge_samples(ps0udatg, "bencahea") #yields NAs in columns of samdfudat within ps1udatg
# ps2udatg <- transform_sample_counts(ps1udatg, function(x) x / sum(x))
count <- vector()
library(metagMisc)
propASVdata <- data.frame(ps.top100udat@otu_table)
propASVdata <- t(propASVdata)
for (i in 1:length(S$ASV)) {
  for (j in 1:length(rownames(hiektaxa2))) {
    if (S$ASV[i] == rownames(hiektaxa2)[j]) {
      count <- c(count,j)
    }
  }
}
asv <- rep(0,length(count))
for (i in 1:length(count)) {
  asv[i] <- paste("ASV",count[i],sep = "")
}
# asv <- data.frame(asv)
# asv$m <- rep(0,length(count))
# asv$lm <- rep(0,length(count))
# asv$mlm <- rep(0,length(count))
newdatasv <- propASVdata[,count[1:4]]

seqtab.nochimhkdat1 <- data.frame(ASVtab2Hiek)
taxasilvahkdat1 <- as.matrix(hiektaxa1)
# samdfhkdat <- HMetaHiek
samdfhkdat <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/Hmetadata.txt")
samdfhkdat1 <- data.frame(samdfhkdat)
samdfhkdat$asv <- "m"
samdfhkdat1 <- samdfhkdat[,-1]
rownames(samdfhkdat1) <- samdfhkdat[,1]
colnames(samdfhkdat1) <- c("hist","m")
storing <- list()
for (i in 1:length(rownames(samdfhkdat1))) {
  if (samdfhkdat1[i,1] == "Atypia") {
    storing <- c(storing, rownames(samdfhkdat1)[i])
    samdfhkdat1 <- samdfhkdat1[-i,]
  }
}
for (i in 1:length(rownames(samdfhkdat))) {
  if (samdfhkdat1[i,1] == "DCIS") {
    storing <- c(storing, rownames(samdfhkdat1)[i])
    samdfhkdat1 <- samdfhkdat1[-i,]
  }
}
for (i in 1:length(rownames(newdatasv))) {
  for (j in 1:length(storing)) {
    if (rownames(newdatasv)[i] == storing[j]){
      newdatasv <- newdatasv[-i,]
    }
  }
}
for (i in 1:length(rownames(seqtab.nochimhkdat1))) {
  for (j in 1:length(storing)) {
    if (rownames(seqtab.nochimhkdat1)[i] == storing[j]){
      seqtab.nochimhkdat1 <- seqtab.nochimhkdat1[-i,]
    }
  }
}

newdatasv$types <- samdfhkdat1$hist
library(ggplot2)
library(reshape2)
mndat1 <- melt(newdatasv)
# Boxplots
ggplot(mndat1,aes(x=types,y=value)) +
  geom_boxplot() +
  geom_jitter(aes(color = variable), height = 0, width = .2) +
  facet_wrap(~variable)
# Proportional abundance barplots
ggplot(mndat1,aes(x=types,y=value,fill = variable)) +
  geom_col() +
  facet_wrap(~variable)
x1 <- rownames(newdatasv)
# can change ASV121 to ASV122, ASV334, and ASV415
ggplot(newdatasv,aes(x=x1,y=ASV121)) +
  geom_col() +
  facet_wrap(~types)

x1 <- rownames(mndat1)
ggplot(mndat1,aes(x=x1,y=value)) +
  geom_col() +
  facet_wrap(~types)+
  facet_grid(~variable)
# ggplot(newdatasv,aes(x=types,y=ASV121)) +
#   geom_boxplot() +
#   geom_jitter(height = 0, width = .2) #+
#   #facet_wrap(~variable)


pshkdat <- phyloseq(otu_table(seqtab.nochimhkdat1, taxa_are_rows=FALSE),
                   sample_data(samdfhkdat), # originally sample_data
                   tax_table(taxasilvahkdat1))
pshkdat <- prune_samples(sample_names(pshkdat) != "Mock", pshkdat)
dnahkdat <- Biostrings::DNAStringSet(taxa_names(pshkdat))
names(dnahkdat) <- taxa_names(pshkdat)
pshkdat <- merge_phyloseq(pshkdat, dnahkdat)
taxa_names(pshkdat) <- paste0("ASV", seq(ntaxa(pshkdat)))
pshkdat
ps.prophkdat <- transform_sample_counts(pshkdat, function(otu) otu/sum(otu))
ord.nmds.brayhkdat <- ordinate(ps.prophkdat, method="NMDS", distance="bray", color="Type")
top100hkdat <- names(sort(taxa_sums(pshkdat), decreasing=TRUE)) [1:100]
ps.top100hkdat <- transform_sample_counts(pshkdat, function(OTU) OTU/sum(OTU))
ps.top100hkdat <- prune_taxa(top100hkdat, ps.top100hkdat)
pshkdatg <- tax_glom(ps.top100hkdat, "Family")
ps0hkdatg <- transform_sample_counts(pshkdatg, function(x) x / sum(x))
ps1hkdatg <- merge_samples(ps0hkdatg, "hist") #yields NAs in columns of samdfhkdat within ps1hkdatg
ps2hkdatg <- transform_sample_counts(ps1hkdatg, function(x) x / sum(x))
phkdatg <- plot_bar(ps2hkdatg, fill="Family") #+ facet_wrap(~bencahea, scales = "free_x")
finalplothkdatg <- phkdatg + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))
finalplothkdatg

#colnames(samdfhkdat) <- c("Run","type")
pshkdat <- phyloseq(otu_table(seqtab.nochimhkdat1, taxa_are_rows=FALSE),
                   sample_data(samdfhkdat), # originally sample_data
                   tax_table(taxasilvahkdat1))
pshkdat <- prune_samples(sample_names(pshkdat) != "Mock", pshkdat)
dnahkdat <- Biostrings::DNAStringSet(taxa_names(pshkdat))
names(dnahkdat) <- taxa_names(pshkdat)
pshkdat <- merge_phyloseq(pshkdat, dnahkdat)
taxa_names(pshkdat) <- paste0("ASV", seq(ntaxa(pshkdat)))
pshkdat
ps.prophkdat <- transform_sample_counts(pshkdat, function(otu) otu/sum(otu))
ord.nmds.brayhkdat <- ordinate(ps.prophkdat, method="NMDS", distance="bray", color="tumor")
top100hkdat <- names(sort(taxa_sums(pshkdat), decreasing=TRUE)) [1:100]
ps.top100hkdat <- transform_sample_counts(pshkdat, function(OTU) OTU/sum(OTU))
ps.top100hkdat <- prune_taxa(top100hkdat, ps.top100hkdat)
pshkdatg <- tax_glom(ps.top100hkdat, "Genus")
ps0hkdatg <- transform_sample_counts(pshkdatg, function(x) x / sum(x))
ps1hkdatg <- merge_samples(ps0hkdatg, "histology_cat") #yields NAs in columns of samdfhkdat within ps1hkdatg
ps2hkdatg <- transform_sample_counts(ps0hkdatg, function(x) x / sum(x))
pudatg <- plot_bar(ps2hkdatg, x= "histology_cat", fill="Genus") #+ facet_wrap(~bencahea, scales = "free_x")
finalplotudatg <- pudatg + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))
finalplotudatg
pshhkdata <- psmelt(ps2hkdatg)
View(pshhkdata)



library(lme4)
library(ape)
library(MiRKAT)
library(GUniFrac)
library(vegan)
hbbdinv <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/Hmetadata.txt")
Hmetgroup <- HMeta[,3]
Hiektree <- read.tree(file = "~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/upgmahieken.nwk")
ASVtabunifracsHiek_perma <- GUniFrac(ASVtabprev2Hiek, Hiektree, alpha=c(0, 0.5, 1))$unifracs
PERMANOVA(Distance = ASVtabunifracsHiek_perma,group = Hmetgroup,seed = "1234",PostHoc = "none")

#ASVtabsub_rff <- Rarefy(ASVtabsub,1024) #$ASVtabsub_rff
# listnums <- ASVtabsub_rff[1]
# listnums <- listnums[1][["otu.tab.rff"]]
# ASVtab <- as.numeric(as.matrix(ASVtab))
ASVtabunifracsHiek <- GUniFrac(ASVtab2Hiek, Hiektree, alpha=c(0, 0.5, 1))$unifracs
D.weightedHiek <- ASVtabunifracsHiek[,,"d_1"]
D.unweightedHiek <- ASVtabunifracsHiek[,,"d_UW"]
D.BCHiek <- as.matrix(vegdist(ASVtab2Hiek , method="bray"))
KweightedHiek <- D2K(D.weightedHiek)
KunweightedHiek <- D2K(D.unweightedHiek)
K.BCHiek <- D2K(D.BCHiek)
# MiRKAT(y = ASVtab2Hiek, Ks = KunweightedHiek, out_type = "D", method = "davies")
# ggplot(data= D.unweightedHiek,mapping = HMetaHiek)

# cmdscale plot
cmdunweunif <- cmdscale(D.unweightedHiek)
View(cmdunweunif)
cmduni <- data.frame(cmdunweunif)
x <- cmduni[,1]
y <- cmduni[,2]
plot(x,y,asp=1)
cmduni$meta <- 0
for (i in 1:length(hbbdinv[,1])){
  #for (j in 1:length(rownames(cmduni))) {
    if (rownames(cmduni)[i] == hbbdinv[i,1]) {
      cmduni$meta[i] <- hbbdinv[i,2]
    }
  }
text(x,y,cmduni$meta,cex=0.8)

# MiRKAT code
## Formatting metadata for input to MiRKAT
Hiekmet <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/Hmetadata_mirkat.txt',header=TRUE)
Hiekmet <- data.frame(Hiekmet)

for (i in 1:42) {
  Hiekmet[i,3] <- as.double(Hiekmet[i,3])
}
Hiekdouble <- Hiekmet[,3]

# works!
KslistHiek <- list(KweightedHiek,KunweightedHiek,K.BCHiek)
meerkatsingleHiekBC <- MiRKAT(y= Hiekdouble, Ks = K.BCHiek, out_type = "D", method = "permutation") # can be Kunweighted or K.BC
meerkatsingleHiekunweigUniFrac <- MiRKAT(y= Hiekdouble, Ks = KunweightedHiek, out_type = "D", method = "permutation")
meerkatsingleHiekweighUniFrac <- MiRKAT(y= Hiekdouble, Ks = KweightedHiek, out_type = "D", method = "permutation")
meerkatmultipleHiek <- MiRKAT(y= Hiekdouble,Ks = KslistHiek, out_type = "D", nperm = 9999, method = "permutation")





# install.packages("edgeR")
library(edgeR)
edgeRUsersGuide()
DGEList(counts = ASVtab2Hiek, samples = HMeta)
