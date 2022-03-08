

# Appending the SRR file that has the corresponding HS or BTS metadata associated to m
# Use m to filter out those SRR files from ASVtab
ASVtab <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVtransposed.txt',header = FALSE)
ASVtab <- data.frame(ASVtab)
Urbmetatoremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbMetaremovefiles.txt',header=TRUE)
Urbremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtoremove.txt', header=FALSE)
m <- rep(0,27)
for (i in 1:27) {
 # print(i)
  for (j in 1:68) {
    if (Urbremove[i,1] == Urbmetatoremove[j,2]) {
      m[i] <- Urbmetatoremove[j,1]
    }
  }
}

for (i in 1:length(m)) {
  if (m[i] == 0) {
    m <- m[-i]
  } else {
    print(FALSE)
  }

}

m <- data.frame(m)
count <- 0
for (i in 2:27) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab[1,i]) {
        ASVtab <- ASVtab[,-i]
      } else{
      count <- count+1
      }
  }
}

for (i in 27:60) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab[1,i]) {
      ASVtab <- ASVtab[,-i]
    } else{
      count <- count+1
    }
  }
}

for (i in 2:50) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab[1,i]) {
      ASVtab <- ASVtab[,-i]
    } else{
      count <- count+1
    }
  }
}

for (i in 2:43) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab[1,i]) {
      ASVtab <- ASVtab[,-i]
    } else{
      count <- count+1
    }
  }
}

write.table(ASVtab,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVmodifiedSRRs.txt", row.names = FALSE, col.names = FALSE)
ASVtab1 <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVmodifiedSRRs.txt", header = TRUE)
ASVtab1 <- data.frame(ASVtab1)
#ASVtab <- t(ASVtab)
ASVtab2 <- ASVtab1[,-1]
rownames(ASVtab2) <- ASVtab1[,1]
ASVtab2 <- t(ASVtab2)
# ASVtab <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVtransposed.txt',header = TRUE)
# ASVtab <- data.frame(ASVtab)
# #ASVtab <- t(ASVtab)
# ASVtab3 <- ASVtab[,-1]
# rownames(ASVtab3) <- ASVtab[,1]
# ASVtab3 <- t(ASVtab3)
# use compositions for clr
for (i in 1:6943) {
  ASVtab2[,i] <- ASVtab2[,i] + 0.50
  
}

# install.packages('compositions')
library(ape)
library(MiRKAT)
library(GUniFrac)
library(vegan)
library(compositions)
library(cluster)
library(ALDEx2)
clrASVtab2 <- clr(ASVtab2)
D.euclid <- as.matrix(vegdist(clrASVtab2 , method="euclidean"))
kmeansclusEuclid <- pam(D.euclid, k=2)
clusplot(kmeansclusEuclid,shade=TRUE)

# for (i in range(2,69)) {
#   ASVtab[,i] <- as.numeric(ASVtab[,i])
# }
#View(ASVtab[,2])

#if (!require("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")

#BiocManager::install("ALDEx2")
library(ALDEx2)
#data(ASVtab)
#data(selex)
#subset for efficiency
#selex <- selex[1201:1600,]
#conds <- c(rep("NS", 7), rep("S", 7))
str(ASVtab[,2])
conditions <- c(rep('BT',25),rep('H',5),rep('BT',1),rep('H',4),rep('BT',1),rep('H',3),rep('BT',3))
x1 <- aldex.clr(ASVtab1[2:43], conditions, mc.samples = 128, denom = "all", verbose = FALSE)

#ASVtabsub <- ASVtab2[1:42, 1:6943]
# ASVtabsub <- t(ASVtab2)
# for (i in 1:42) {
#   for (j in 1:6943) {
#   ASVtabsub[j,i] <- as.numeric(ASVtabsub[j,i])
#   }
# }
# clrASVtab2 <- data.frame(t(clrASVtab2))
# row.names(clrASVtab2) <- NULL
# #colnames(clrASVtab2) <- NULL
# for (i in 1:42) {
#   for (j in 1:6943) {
#     clrASVtab2[j,i] <- as.numeric(clrASVtab2[j,i])
#   }
# }


x <- aldex(ASVtab1[2:43], conditions, mc.samples=128, denom="all", test="t", effect=TRUE, include.sample.summary=TRUE, verbose=FALSE, iterate=FALSE )
#x1 <- aldex.clr(ASVtabsub, conditions, mc.samples = 128, denom = "all", verbose = FALSE)
write.table(x, '~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaldexObject03052022.txt')
#x1 <- data.frame(x1)
write.table(x1, '/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/02152022_version/UrbaldexClrObjASVs.csv')
#getDirichletInstances(.object)
monteCarloDirInstances <- numMCInstances(x1)
#monteCarloDirInstances1 <- getDirichletInstances(x1)
write.table(monteCarloDirInstances, '/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/02152022_version/UrbdirichMCwithAldObj.txt')
#write.table(monteCarloDirInstances1, '/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/02152022_version/UrbdirichMCwithAldClrObj.txt')
library(ape)
library(MiRKAT)
library(GUniFrac)
library(vegan)
Urbtree <- read.tree(file = "~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/upgmaURB11122021.nwk")
ASVtabsub_rff <- Rarefy(ASVtabsub,1024) #$ASVtabsub_rff
listnums <- ASVtabsub_rff[1]
listnums <- listnums[1][["otu.tab.rff"]]
ASVtab <- as.numeric(as.matrix(ASVtab))
ASVtabunifracs <- GUniFrac(ASVtab2, Urbtree, alpha=c(0, 0.5, 1))$unifracs
D.weighted <- ASVtabunifracs[,,"d_1"]
D.unweighted <- ASVtabunifracs[,,"d_UW"]
D.BC <- as.matrix(vegdist(ASVtab2 , method="bray"))
# for constructing kernel matrices for MiRKAT
Kweighted <- D2K(D.weighted)
Kunweighted <- D2K(D.unweighted)
K.BC <- D2K(D.BC)
Urbmet <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMetadata--.txt',header=TRUE)
Urbmet <- as.vector(Urbmet)

# Error message
# Error in model.frame.default(formula = y ~ X1 - 1, drop.unused.levels = TRUE) : 
#   invalid type (list) for variable 'y'
MiRKAT(y = Urbmet, Ks = c(Kweighted,Kunweighted,K.BC), out_type = "D", method = "davies")
 
# install.packages("cluster")
library(cluster)
kmeansclusBray <- pam(D.BC, k=2) #bray-curtis
clusplot(kmeansclusBray,shade=TRUE)
kmeansclusUniWeigh <- pam(D.weighted, k=2) #Unifrac weighted
clusplot(kmeansclusUniWeigh,shade=TRUE)

#yields error message
# Error in princomp.default(x, scores = TRUE, cor = ncol(x) > 2) : 
#   cannot use 'cor = TRUE' with a constant variable
kmeansclusUniUnweigh <- pam(D.unweighted, k=2) #Unifrac unweighted
clusplot(kmeansclusUniUnweigh,shade=TRUE)



#install.packages("stats")
library(stats)
kmeansalgorithm <- kmeans(x=D.BC,centers=2)

#For UniFrac Plots
##Unweighted UniFrac
UniFrac(pskdatphy) #Weighted default is False
unweightedkdat <- UniFrac(pskdatphy)
weightedunifracplot <- UniFrac(pskdatphy, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
plot(weightedunifracplot)
# library(GUniFrac)
# unifracs <- GUniFrac(seqtab.nochimkdat, psktree, alpha = c(0, 0.5, 1))
# weightedUni <- unifracs[, , "d_0"]
# unweightedUni <- unifracs[, , "d_UW"]
##Weighted UniFrac
UniFrac(pskdatphy, weighted = TRUE)
weightedkdat <- UniFrac(pskdatphy, weighted = TRUE)
plot(weightedkdat, type = "p", pch)

# library(dada2)
# library(DECIPHER)
# library(phangorn)
# library(ggplot2)
# library(phyloseq)
#ASVtab <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVtransposed_.txt',header = FALSE)
#ASVtab <- data.frame(ASVtab)
# for (i in range(2,69)) {
#   ASVtab[,i] <- as.numeric(ASVtab[,i])
# }
# asvt <- t(ASVtab)
# ## Extract sequences from DADA2 output
# sequencesudat <- getSequences(asvt)
# names(sequencesudat) <- sequencesudat
# 
# ## Run Sequence Alignment (MSA) using DECIPHER
# alignmentudat <- AlignSeqs(DNAStringSet(sequencesudat), anchor=NA)
# 
# ## Change sequence alignment output into a phyDat structure
# phang.alignudat <- phyDat(as(alignmentudat, "matrix"), type="DNA")
# 
# ## Create distance matrix
# dmudat <- dist.ml(phang.alignudat)
# UPGMAtreeudat <- upgma(dmudat)
# write.tree(UPGMAtreeudat, file = "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/VersionControl_DadaCommandRevisions/20211112_version/upgmaURB11122021.nwk", append = FALSE,
#            digits = 10, tree.names = FALSE)

