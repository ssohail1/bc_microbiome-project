ASVtab <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
ASVtab <- data.frame(ASVtab)
#ASVtab <- t(ASVtab)
ASVtab2Hiek <- ASVtab[,-1]
rownames(ASVtab2Hiek) <- ASVtab[,1]
ASVtab2Hiek <- t(ASVtab2Hiek)
HMeta <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissues.txt")
HMeta <- data.frame(t(HMeta))

# Filter with prevalence = 10% and taxa < 0.2% then square root transform
ASVtabprev <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabwithprevalencepercentsHieken03152022.txt',header = TRUE)
ASVtabprev <- data.frame(ASVtabprev)
#ASVtabprev <- t(ASVtabprev)

for (i in 58:length(ASVtabprev$percentofzeros)) {
  if (ASVtabprev$percentofzeros[i] >= 90) {
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

library(ape)
library(MiRKAT)
library(GUniFrac)
library(vegan)
Hiektree <- read.tree(file = "~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/upgmahieken.nwk")
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
