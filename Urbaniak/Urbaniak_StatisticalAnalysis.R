
# Formatting ASV table for clr() function
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
ASVtab2 <- ASVtab1[,-1]
rownames(ASVtab2) <- ASVtab1[,1]
ASVtab2 <- t(ASVtab2)
# use compositions for clr
# for (i in 1:6943) {
#   ASVtab2[,i] <- ASVtab2[,i] + 0.50
#   
# }


# ALDEx2 code
# install.packages('compositions')
library(ape)
library(MiRKAT)
library(GUniFrac)
library(vegan)
library(compositions)
library(cluster)
library(ALDEx2)
library(tibble)
library(dplyr)

ASVtab12 <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVmodifiedSRRs.txt", header = FALSE)
ASVtab12 <- data.frame(ASVtab12)
Urbmetatoremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbMetaremovefiles.txt',header=TRUE)
Urbremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtoremovebenign.txt', header=FALSE)

m <- rep(0,length(Urbremove[,1]))
for (i in 1:length(Urbremove[,1])) {
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
for (i in 2:14) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab12[1,i]) {
      ASVtab12 <- ASVtab12[,-i]
    } else{
      count <- count+1
    }
  }
}

for (i in 14:26) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab12[1,i]) {
      ASVtab12 <- ASVtab12[,-i]
    } else{
      count <- count+1
    }
  }
}

for (i in 26:33) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab12[1,i]) {
      ASVtab12 <- ASVtab12[,-i]
    } else{
      count <- count+1
    }
  }
}

for (i in 2:33) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == ASVtab12[1,i]) {
      ASVtab12 <- ASVtab12[,-i]
    } else{
      count <- count+1
    }
  }
}

#write.table(ASVtab12,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVremovebenign.txt", row.names = FALSE, col.names = FALSE)
ASVtab12 <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVremovebenign.txt", header = TRUE)
ASVtab12 <- data.frame(ASVtab12)
ASVtab211 <- ASVtab12[,-1]
rownames(ASVtab211) <- ASVtab12[,1]
conditions <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
#x <- aldex(ASVtab12[2:43], conditions, mc.samples=128, denom="all", test="t", effect=TRUE, include.sample.summary=TRUE, verbose=FALSE, iterate=FALSE)

x1 <- aldex.clr(ASVtab211, conditions, mc.samples = 128, denom = "all", verbose = FALSE)
x1ttest <- aldex.ttest(x1)
x1effect <- aldex.effect(x1,include.sample.summary = TRUE,CI = TRUE,verbose = FALSE)
x1aldex_out <- data.frame(x1ttest, x1effect)
#write.table(x1aldex_out,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/aldexoutput03172022_allsamps03172022.txt")#, row.names = FALSE)
x1alddata <- rownames_to_column(x1aldex_out,"ASVs") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than ttest welch
  dplyr::select(ASVs, we.eBH, wi.eBH, effect, overlap) %>%
  data.frame()

Urbtaxa <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile.txt",header = TRUE)
Urbtaxa <- data.frame(Urbtaxa)

count <- vector()
x1alddata$Taxa <- 0
for (i in 1:length(x1alddata$ASVs)) {
  # nm <- x1alddata$samp_num[i]
  for (j in 1:length(Urbtaxa[,1])) {
    if (Urbtaxa[j,1] == x1alddata$ASVs[i]){
      nmj <- as.character(Urbtaxa[j,2:7])
      if (is.na(nmj[length(nmj)]) == TRUE) {
        x1alddata$Taxa[i] <- nmj[length(nmj)-1]
        count <- c(count,length(nmj)-1)
      } else{
        x1alddata$Taxa[i] <- nmj[length(nmj)]
        count <- c(count,length(nmj))
      }
    }
  }
}


# x1alddata$ASVs <- 0
# for (i in 1:length(x1alddata[,1])) {
#   n <- x1alddata$samp_num[i]
#   x1alddata$ASVs[i] <- ASVtab12[n,1]
# }

# count tells us whether genus - 6 or family - 5 was added
print(count)

# x1alddata1 <- as.data.frame(x1alddata)
# for (i in 1:7) {
#   x1alddata[,i] <- data.frame(x1alddata[,i])
# }

# Has the ASVs in the aldex output
write(x1alddata,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/x1aldexoutput03182022_ver0318.txt")

xeffect <- aldex.effect(x1,include.sample.summary = FALSE,CI = TRUE,verbose = FALSE)
aldex_out <- data.frame(x1ttest, xeffect)
xalddata <- rownames_to_column(aldex_out,"samp_num") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than ttest welch
  dplyr::select(samp_num, we.eBH, wi.eBH, effect, overlap) %>%
  data.frame()
View(xalddata)
xalddata$ASVs <- 0
for (i in 1:length(xalddata[,1])) {
  n <- xalddata$samp_num[i]
  xalddata$ASVs[i] <- ASVtab12[n,1]
}

# Has the ASVs in the aldex output
write.table(xalddata,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/aldexoutput03172022_ver0317.txt") #, row.names = FALSE)

# K-means Clustering code
# library(ape)
# library(MiRKAT)
# library(GUniFrac)
# library(vegan)

ASVtab21 <- ASVtab12[,-1]
rownames(ASVtab21) <- ASVtab12[,1]
ASVtab21 <- t(ASVtab21)

Urbtree <- read.tree(file = "~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/upgmaURB11122021.nwk")
# ASVtabsub_rff <- Rarefy(ASVtabsub,1024) #$ASVtabsub_rff
# listnums <- ASVtabsub_rff[1]
# listnums <- listnums[1][["otu.tab.rff"]]
# ASVtab <- as.numeric(as.matrix(ASVtab))
ASVtabunifracs <- GUniFrac(ASVtab21, Urbtree, alpha=c(0, 0.5, 1))$unifracs
D.weighted <- ASVtabunifracs[,,"d_1"]
D.unweighted <- ASVtabunifracs[,,"d_UW"]
D.BC <- as.matrix(vegdist(ASVtab21 , method="bray"))
# for constructing kernel matrices for MiRKAT
Kweighted <- D2K(D.weighted)
Kunweighted <- D2K(D.unweighted)
K.BC <- D2K(D.BC)

library(stats)
# install.packages("cluster")
# PAM uses k-means medoids approach
pamclusBray <- pam(D.BC, k = 2) #bray-curtis
pamBray <- data.frame(pamclusBray$clustering)
pamBray$samptype <- 0
for (i in 1:42) {
  for (j in 1:68) {
    if (row.names(pamBray)[i] == Urbmetatoremove[j,1]) {
      pamBray$samptype[i] <- Urbmetatoremove[j,2]
    }
  }
}
clusplot(pamclusBray,shade=TRUE)
# kmeans uses k-means clustering approach
kmeansclusBray <- kmeans(D.BC,centers = 2) #bray-curtis
kmeansBray <- data.frame(kmeansclusBray$cluster)
kmeansBray$samptype <- 0
for (i in 1:42) {
  for (j in 1:68) {
    if (row.names(kmeansBray)[i] == Urbmetatoremove[j,1]) {
      kmeansBray$samptype[i] <- Urbmetatoremove[j,2]
    }
  }
}
clusplot(kmeansclusBray,shade=TRUE)
kmeansclusUniWeigh <- pam(D.weighted, k=2) #Unifrac weighted
clusplot(kmeansclusUniWeigh,shade=TRUE)

#yields error message
# Error in princomp.default(x, scores = TRUE, cor = ncol(x) > 2) : 
#   cannot use 'cor = TRUE' with a constant variable
kmeansclusUniUnweigh <- pam(D.unweighted, k=2) #Unifrac unweighted
clusplot(kmeansclusUniUnweigh,shade=TRUE)
# Error message: Error in princomp.default(x, scores = TRUE, cor = ncol(x) > 2) : 
# cannot use 'cor = TRUE' with a constant variable

#install.packages("stats")
library(stats)
kmeansalgorithm <- kmeans(x=D.BC,centers=2)

# MiRKAT code
## Formatting metadata for input to MiRKAT
Urbmet <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMetadata--.txt',header=FALSE)
Urbmet <- data.frame(Urbmet)

count <- 0
for (i in 2:27) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == Urbmet[i,1]) {
      Urbmet <- Urbmet[-i,]
    } else{
      count <- count+1
    }
  }
}

for (i in 27:60) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == Urbmet[i,1]) {
      Urbmet <- Urbmet[-i,]
    } else{
      count <- count+1
    }
  }
}

for (i in 2:50) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == Urbmet[i,1]) {
      Urbmet <- Urbmet[-i,]
    } else{
      count <- count+1
    }
  }
}

for (i in 2:43) {
  for (j in 1:length(m[,1])) {
    #if (var1 != NULL) {
    if (m[j,1] == Urbmet[i,1]) {
      Urbmet <- Urbmet[-i,]
    } else{
      count <- count+1
    }
  }
}

write.table(Urbmet,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMet--.txt",row.names = FALSE,col.names = FALSE)
# remove the NN header from the UrbaniakMet--.txt file and change spaces to tabs
Urbmet <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMet--.txt",header = TRUE)
Urbmet <- data.frame(Urbmet)
for (i in 1:42) {
  Urbmet[i,2] <- as.double(Urbmet[i,2])
}
Urbdouble <- Urbmet[,2]
# Error message
# Error in model.frame.default(formula = y ~ X1 - 1, drop.unused.levels = TRUE) : 
#   invalid type (list) for variable 'y'

# works!
meerkatsingleBC <- MiRKAT(y= Urbdouble, Ks = K.BC, out_type = "D", method = "permutation") # can be Kunweighted or K.BC
meerkatsingleunweigUniFrac <- MiRKAT(y= Urbdouble, Ks = Kunweighted, out_type = "D", method = "permutation")
meerkatsingleweighUniFrac <- MiRKAT(y= Urbdouble, Ks = Kweighted, out_type = "D", method = "permutation")
Kslist <- list(Kweighted,Kunweighted,K.BC)
meerkatmultiple <- MiRKAT(y= Urbdouble,Ks = Kslist, out_type = "D", nperm = 9999, method = "permutation")


library(cluster)
for (i in 1:6943) {
  ASVtab21[,i] <- ASVtab21[,i] + 0.50
  
}
rownames(ASVtab21)
rownames(ASVtab21) <- c("BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "H", "H", "H", "H", "H", "BT", "H", "H", "H", "H", "BT", "H", "H", "H", "BT", "BT", "BT")
clrASVtab21 <- clr(ASVtab21)
D.euclid <- as.matrix(vegdist(clrASVtab21 , method="euclidean"))
kmeansclusEuclid <- pam(D.euclid, k=2)
clusplot(kmeansclusEuclid,shade=TRUE)
clusplot(kmeansclusEuclid,shade=TRUE,labels = 2,add=TRUE)

# euclidclus <- data.frame(kmeansclusEuclid$clustering)
# euclidclus$samptype <- 0
# for (i in 1:32) {
#   for (j in 1:68) {
#     if (row.names(euclidclus)[i] == Urbmetatoremove[j,1]) {
#       euclidclus$samptype[i] <- Urbmetatoremove[j,2]
#     }
#   }
# }



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


# Example data
# library(MiRKAT)
# library(GUniFrac)
# data(throat.tree)
# data(throat.otu.tab)
# data(throat.meta)
# attach(throat.meta)
# #3.2 Prepare the data
# set.seed(123)
# Male = (Sex == "Male")**2
# Smoker =(SmokingStatus == "Smoker") **2  # a "double" object with 0s and 1s and nothing else
# anti = (AntibioticUsePast3Months_TimeFromAntibioticUsage != "None")^2
# cova = cbind(Male, anti)

# do not use daisy - it assigns clusters wrong - all samples were 1 except for 1 sample
# D.eucliddaisy <- daisy(x = clrASVtab2, metric = "euclidean")
# kmeansclusEucliddaisy <- pam(D.eucliddaisy, k=2)
# clusplot(kmeansclusEucliddaisy,shade=TRUE)
# daisyclus <- data.frame(kmeansclusEucliddaisy$clustering)
# daisyclus$samptype <- 0
# for (i in 1:42) {
#   for (j in 1:68) {
#     if (row.names(daisyclus)[i] == Urbmetatoremove[j,1]) {
#       daisyclus$samptype[i] <- Urbmetatoremove[j,2]
#     }
#   }
# }
# 
# m0 <- matrix(NA, 4, 0)
# rownames(m0)
# 
# m2 <- cbind(1, 1:4)
# colnames(m2, do.NULL = FALSE)
# colnames(m2) <- c("x","Y")
# rownames(m2) <- rownames(m2, do.NULL = FALSE, prefix = "Obs.")
# m2

