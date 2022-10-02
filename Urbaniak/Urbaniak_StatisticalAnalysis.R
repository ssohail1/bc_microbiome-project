#### updated code - need to clean up and add comments ####
# ALDEx2
library(ape)
library(MiRKAT)
library(GUniFrac)
library(vegan)
library(compositions)
library(cluster)
library(ALDEx2)
library(tibble)
library(dplyr)
library(phyloseq)
# ASV table has 32 samples - ~12 benign samples removed and 26 samples removed for quality control - based on paper
ASVtab12 <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVremovebenign.txt", header = TRUE)
ASVtab12 <- data.frame(ASVtab12)
ASVtab211 <- ASVtab12[,-1]
rownames(ASVtab211) <- ASVtab12[,1]
conditions <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))

# running aldex2 functions
x1 <- aldex.clr(ASVtab211, conditions, mc.samples = 128, denom = "all", verbose = FALSE)
x1ttest <- aldex.ttest(x1)
x1effect <- aldex.effect(x1,include.sample.summary = TRUE,CI = TRUE,verbose = FALSE)
x1aldex_out <- data.frame(x1ttest, x1effect)
#write.table(x1aldex_out,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/aldexoutput03172022_allsamps03172022.txt")#, row.names = FALSE)
x1alddata <- rownames_to_column(x1aldex_out,"ASVs") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than ttest welch
  dplyr::select(ASVs, we.eBH, wi.eBH, effect, overlap) %>%
  data.frame()

# adding in taxa info to aldex dataframe
Urbtaxa <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile.txt",header = TRUE)
Urbtaxa <- data.frame(Urbtaxa)
View(Urbtaxa)

count <- vector()
x1alddata$Taxa <- 0
for (i in 1:length(x1alddata$ASVs)) {
  for (j in 1:length(Urbtaxa[,1])) {
    if (Urbtaxa[j,1] == x1alddata$ASVs[i]){
      # if the ASVs in Urbtaxa matches the ASVs in aldex dataframe
      # then add the Kingdom to Genus to nmj
      nmj <- as.character(Urbtaxa[j,2:7])
      # if the last column in nmj - the Genus level is NA
      # then add the second to last column - the Family level in the 
      # taxa column of the aldex dataframe and add the index j of where the ASV was
      # found in Urbtaxa to count variable 
      if (is.na(nmj[length(nmj)]) == TRUE) {
        x1alddata$Taxa[i] <- nmj[length(nmj)-1]
        count <- c(count,j)
      # if the last column in nmj - the Genus level is not NA
      # then add the Genus to the taxa column of the aldex dataframe
      # and add the index j of where the ASV was found in Urbtaxa to count variable 
      } else{
        x1alddata$Taxa[i] <- nmj[length(nmj)]
        count <- c(count,j)
      }
    }
  }
}
View(x1alddata)

# Urbremben = the benign samples that need to be removed from samdfudat
# Urbmetatoremove = has all the samples and the associated sample type metadata
Urbremben <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbremoveben.txt',header=TRUE)
Urbmetatoremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbmetremove42samp.txt',header=TRUE)
samdfudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMet--.txt",header = TRUE)
samdfudat <- data.frame(samdfudat)

# healthy = labels of which samples are healthy
# ca = labels of which samples are canc.
healthy <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofhea.txt",header = TRUE)
healthy <- data.frame(healthy)
ben <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofben.txt",header = TRUE)
ben <- data.frame(ben)
ca <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofca.txt",header = TRUE)
ca <- data.frame(ca)

# new column of where the sample type will be stored
samdfudat$Types <- 0
for (i in 1:length(Urbmetatoremove[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    # if the sample name in samdfudat matches sample name in Urbmetatoremove
    # then add the sample type info (BTS vs HS) to the Types column for index j in samdfudat
    if (rownames(samdfudat)[j] == Urbmetatoremove[i,1]){
      samdfudat$Types[j] <- Urbmetatoremove[i,2]
    }
  }
}
# After adding that info for all 42 samples in samdfudat
# then can remove the specific ~10 benign samples
for (i in 1:length(Urbremben[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat$Types[j] == Urbremben[i,1]){
      samdfudat <- samdfudat[-j,]
    }
  }
}

# adding the sample type label info (healthy vs. canc.) to bencahea column in samdfudat
# for healthy
for (i in 1:length(healthy[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == healthy[i,1]){
      samdfudat$bencahea[j] <- healthy[i,2]
    }
  }
}
# for canc.
for (i in 1:length(ca[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == ca[i,1]){
      samdfudat$bencahea[j] <- ca[i,2]
    }
  }
}
# replacing 1 in Sample_Type column in samdfudat with
# associated info
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,2] == 1){
    samdfudat$Sample_Type[j] <- "Breast_Tumor"
  }
}
# replacing 0 in Sample_Type column in samdfudat with
# associated info
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,2] == 0){
    samdfudat$Sample_Type[j] <- "Healthy"
  }
}
# running phyloseq analyses to get proportional abundances of the samples
# in the ASV table
samdfudat <- data.frame(samdfudat)
View(samdfudat)
seqtab.nochimudat1 <- t(ASVtab211)
library(phyloseq)
psudat <- phyloseq(otu_table(seqtab.nochimudat1, taxa_are_rows=FALSE),
                   sample_data(samdfudat), # originally sample_data
                   tax_table(taxasilvaudat1))
psudat <- prune_samples(sample_names(psudat) != "Mock", psudat)
dnaudat <- Biostrings::DNAStringSet(taxa_names(psudat))
names(dnaudat) <- taxa_names(psudat)
psudat <- merge_phyloseq(psudat, dnaudat)
taxa_names(psudat) <- paste0("ASV", seq(ntaxa(psudat)))
psudat
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE)) 
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)

# propASVdata has the proportional abundance otu table
propASVdata <- data.frame(ps.top100udat@otu_table)
# propASVdata <- t(propASVdata)
# colnames(propASVdata) <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
# 1704/(sum(ASVtab211[,1]))
# 878/(sum(ASVtab211[,2]))
# 1098/(sum(ASVtab211[,3]))

# count has the indices of where the taxa in Urbtaxa matched with the aldex output
# to get ASV names in format of ASV24 ASV100 etc. so that it matches the ASV format
# in the otu table
asvs <- rep(0,length(count))
for (i in 1:length(count)) {
  asvs[i] <- paste("ASV",count[i],sep = "")
}

count <- vector()
for (i in 1:length(asvs)) {
  for (j in 1:length(colnames(propASVdata))) {
    if (asvs[i] == colnames(propASVdata)[j]) {
      count <- c(count,j)
      # boxplot(propASVdata[j,])
    }
  }
}

x1alddata1 <- rownames_to_column(x1aldex_out,"ASVs") %>%
  dplyr::select(ASVs) %>%
  data.frame()
x1alddata1 <- data.frame(x1alddata1[1:32,1])
x1alddata1$zero <- rep(0,32)
rownames(x1alddata1) <- rownames(propASVdata)

# to get the graphs with ASV labels
newdatf <- cbind(x1alddata1,propASVdata[,count[1:20]])
View(newdatf)
newdatf$tissuetyp <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
newdatf <- newdatf[,-1]
newdatf <- newdatf[,-1]
library(reshape2)
mndat <- melt(newdatf)
ggplot(mndat,aes(x=tissuetyp,y=value)) +
  geom_boxplot() +
  geom_jitter(aes(color = variable), height = 0, width = .2) +
  facet_wrap(~variable)

# to get that graphs with taxa labels
newdatftaxa <- cbind(x1alddata1,propASVdata[,count[1:20]])
newdatftaxa <- newdatftaxa[,-1]
newdatftaxa <- newdatftaxa[,-1]
colnames(newdatftaxa) <- x1alddata$Taxa
newdatftaxa$tissuetyp <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
mndattaxa <- melt(newdatftaxa)
ggplot(mndattaxa,aes(x=tissuetyp,y=value)) +
  geom_boxplot() + # added aes(fill = variable)
  #geom_jitter(height = 0, width = .2) +
  scale_y_continuous(breaks = round(seq(min(mndattaxa$value), max(mndattaxa$value), by = 0.0005),1)) +
  facet_wrap(~variable)


#### ALDEx2 ####
library(ape)
library(MiRKAT)
library(GUniFrac)
library(vegan)
library(compositions)
library(cluster)
library(ALDEx2)
library(tibble)
library(dplyr)
library(phyloseq)
# ASV table has 32 samples - ~12 benign samples removed and 26 samples removed for quality control - based on paper
ASVtab12 <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVremovebenign.txt", header = TRUE)
ASVtab12 <- data.frame(ASVtab12)
ASVtab211 <- ASVtab12[,-1]
rownames(ASVtab211) <- ASVtab12[,1]
conditions <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))

# running aldex2 functions
x1 <- aldex.clr(ASVtab211, conditions, mc.samples = 128, denom = "all", verbose = FALSE)
x1ttest <- aldex.ttest(x1)
x1effect <- aldex.effect(x1,include.sample.summary = TRUE,CI = TRUE,verbose = FALSE)
x1aldex_out <- data.frame(x1ttest, x1effect)
#write.table(x1aldex_out,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/aldexoutput03172022_allsamps03172022.txt")#, row.names = FALSE)
x1alddata <- rownames_to_column(x1aldex_out,"ASVs") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than ttest welch
  dplyr::select(ASVs, we.eBH, wi.eBH, effect, overlap) %>%
  data.frame()

# adding in taxa info to aldex dataframe
Urbtaxa <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile.txt",header = TRUE)
Urbtaxa <- data.frame(Urbtaxa)

count <- vector()
x1alddata$Taxa <- 0
for (i in 1:length(x1alddata$ASVs)) {
  for (j in 1:length(Urbtaxa[,1])) {
    if (Urbtaxa[j,1] == x1alddata$ASVs[i]){
      # if the ASVs in Urbtaxa matches the ASVs in aldex dataframe
      # then add the Kingdom to Genus to nmj
      nmj <- as.character(Urbtaxa[j,2:7])
      # if the last column in nmj - the Genus level is NA
      # then add the second to last column - the Family level in the 
      # taxa column of the aldex dataframe and add the index j of where the ASV was
      # found in Urbtaxa to count variable 
      if (is.na(nmj[length(nmj)]) == TRUE) {
        x1alddata$Taxa[i] <- nmj[length(nmj)-1]
        count <- c(count,j)
      # if the last column in nmj - the Genus level is not NA
      # then add the Genus to the taxa column of the aldex dataframe
      # and add the index j of where the ASV was found in Urbtaxa to count variable 
      } else{
        x1alddata$Taxa[i] <- nmj[length(nmj)]
        count <- c(count,j)
      }
    }
  }
}

# Urbremben = the benign samples that need to be removed from samdfudat
# Urbmetatoremove = has all the samples and the associated sample type metadata
Urbremben <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbremoveben.txt',header=TRUE)
Urbmetatoremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbmetremove42samp.txt',header=TRUE)
samdfudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMet--.txt",header = TRUE)
samdfudat <- data.frame(samdfudat)

# healthy = labels of which samples are healthy
# ca = labels of which samples are canc.
healthy <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofhea.txt",header = TRUE)
healthy <- data.frame(healthy)
ben <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofben.txt",header = TRUE)
ben <- data.frame(ben)
ca <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofca.txt",header = TRUE)
ca <- data.frame(ca)

# new column of where the sample type will be stored
samdfudat$Types <- 0
for (i in 1:length(Urbmetatoremove[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    # if the sample name in samdfudat matches sample name in Urbmetatoremove
    # then add the sample type info (BTS vs HS) to the Types column for index j in samdfudat
    if (rownames(samdfudat)[j] == Urbmetatoremove[i,1]){
      samdfudat$Types[j] <- Urbmetatoremove[i,2]
    }
  }
}
# After adding that info for all 42 samples in samdfudat
# then can remove the specific ~10 benign samples
for (i in 1:length(Urbremben[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat$Types[j] == Urbremben[i,1]){
      samdfudat <- samdfudat[-j,]
    }
  }
}

# adding the sample type label info (healthy vs. canc.) to bencahea column in samdfudat
# for healthy
for (i in 1:length(healthy[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == healthy[i,1]){
      samdfudat$bencahea[j] <- healthy[i,2]
    }
  }
}
# for canc.
for (i in 1:length(ca[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == ca[i,1]){
      samdfudat$bencahea[j] <- ca[i,2]
    }
  }
}
# replacing 1 in Sample_Type column in samdfudat with
# associated info
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,2] == 1){
    samdfudat$Sample_Type[j] <- "Breast_Tumor"
  }
}
# replacing 0 in Sample_Type column in samdfudat with
# associated info
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,2] == 0){
    samdfudat$Sample_Type[j] <- "Healthy"
  }
}
# running phyloseq analyses to get proportional abundances of the samples
# in the ASV table
samdfudat <- data.frame(samdfudat)
View(samdfudat)
seqtab.nochimudat1 <- t(ASVtab211)
library(phyloseq)
psudat <- phyloseq(otu_table(seqtab.nochimudat1, taxa_are_rows=FALSE),
                   sample_data(samdfudat), # originally sample_data
                   tax_table(taxasilvaudat1))
psudat <- prune_samples(sample_names(psudat) != "Mock", psudat)
dnaudat <- Biostrings::DNAStringSet(taxa_names(psudat))
names(dnaudat) <- taxa_names(psudat)
psudat <- merge_phyloseq(psudat, dnaudat)
taxa_names(psudat) <- paste0("ASV", seq(ntaxa(psudat)))
psudat
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE)) 
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)

# propASVdata has the proportional abundance otu table
propASVdata <- data.frame(ps.top100udat@otu_table)
# propASVdata <- t(propASVdata)
# colnames(propASVdata) <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
# 1704/(sum(ASVtab211[,1]))
# 878/(sum(ASVtab211[,2]))
# 1098/(sum(ASVtab211[,3]))

# count has the indices of where the taxa in Urbtaxa matched with the aldex output
# to get ASV names in format of ASV24 ASV100 etc. so that it matches the ASV format
# in the otu table
asvs <- rep(0,length(count))
for (i in 1:length(count)) {
  asvs[i] <- paste("ASV",count[i],sep = "")
}

count <- vector()
for (i in 1:length(asvs)) {
  for (j in 1:length(colnames(propASVdata))) {
    if (asvs[i] == colnames(propASVdata)[j]) {
      count <- c(count,j)
      # boxplot(propASVdata[j,])
    }
  }
}

x1alddata1 <- rownames_to_column(x1aldex_out,"ASVs") %>%
  dplyr::select(ASVs) %>%
  data.frame()
x1alddata1 <- data.frame(x1alddata1[1:32,1])
x1alddata1$zero <- rep(0,32)
rownames(x1alddata1) <- rownames(propASVdata)

# to get the graphs with ASV labels
newdatf <- cbind(x1alddata1,propASVdata[,count[1:20]])
View(newdatf)
newdatf$tissuetyp <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
newdatf <- newdatf[,-1]
newdatf <- newdatf[,-1]
library(reshape2)
mndat <- melt(newdatf)
ggplot(mndat,aes(x=tissuetyp,y=value)) +
  geom_boxplot() +
  geom_jitter(aes(color = variable), height = 0, width = .2) +
  facet_wrap(~variable)

# to get that graphs with taxa labels
newdatftaxa <- cbind(x1alddata1,propASVdata[,count[1:20]])
newdatftaxa <- newdatftaxa[,-1]
newdatftaxa <- newdatftaxa[,-1]
colnames(newdatftaxa) <- x1alddata$Taxa
newdatftaxa$tissuetyp <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
mndattaxa <- melt(newdatftaxa)
ggplot(mndattaxa,aes(x=tissuetyp,y=value)) +
  geom_boxplot() +
  geom_jitter(aes(color = variable), height = 0, width = .2) +
  facet_wrap(~variable)


# K-means clustering 
library(stats)
library(cluster)
library(GUniFrac)
library(ape)
library(vegan)
ASVtab21 <- ASVtab12[,-1]
rownames(ASVtab21) <- ASVtab12[,1]
ASVtab21 <- t(ASVtab21)

# reading in Urbaniak tree file
Urbtree <- read.tree(file = "~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/upgmaURB11122021.nwk")

# Getting UniFrac distances through GUniFrac with ASV table and tree file as input
rownames(ASVtab21) <- c("BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "H", "H", "H", "H", "H", "BT", "H", "H", "H", "H", "BT", "H", "H", "H", "BT", "BT", "BT")
ASVtabunifracs <- GUniFrac(ASVtab21, Urbtree, alpha=c(0, 0.5, 1))$unifracs
D.weighted <- ASVtabunifracs[,,"d_1"]
D.unweighted <- ASVtabunifracs[,,"d_UW"]
D.BC <- as.matrix(vegdist(ASVtab21 , method="bray"))

# PAM uses k-means medoids approach (from cluster package)
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
clusplot(pamclusBray,shade=TRUE,main = paste("Bray-Curtis K-means Cluster"),labels = 2,add=TRUE)

kmeansclusUniWeigh <- pam(D.weighted, k=2) #Unifrac weighted
clusplot(kmeansclusUniWeigh,shade=TRUE)
clusplot(kmeansclusUniWeigh,shade=TRUE,labels = 2,add=TRUE)

kmeansclusUniUnweigh <- pam(D.unweighted, k=2) #Unifrac unweighted
clusplot(kmeansclusUniUnweigh,shade=TRUE)
clusplot(kmeansclusUniUnweigh,shade=TRUE,labels = 2,add=TRUE)

# Euclidean distance k-means clustering
# CLR-transforming ASV table
# adding a 0.5 uniform prior to CLR-transformation
library(cluster)
library(compositions)
library(vegan)

ASVtab21 <- ASVtab12[,-1]
rownames(ASVtab21) <- ASVtab12[,1]
ASVtab21 <- t(ASVtab21)

for (i in 1:6943) {
  ASVtab21[,i] <- ASVtab21[,i] + 0.50
  
}

rownames(ASVtab21)
rownames(ASVtab21) <- c("BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "H", "H", "H", "H", "H", "BT", "H", "H", "H", "H", "BT", "H", "H", "H", "BT", "BT", "BT")
clrASVtab21 <- clr(ASVtab21)
D.euclid <- as.matrix(vegdist(clrASVtab21 , method="euclidean"))
kmeansclusEuclid <- pam(D.euclid, k=2)
clusplot(kmeansclusEuclid,shade=TRUE)
clusplot(kmeansclusEuclid,shade=TRUE,main = paste("Euclid K-means Cluster"),labels = 2,add=TRUE,)



