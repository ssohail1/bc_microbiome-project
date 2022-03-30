# Prop - Abund
# BiocManager::install(version='3.14')
# BiocManager::install("microbiome")
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
library(metagMisc); packageVersion("metagMisc")

export_ps <- phyloseq_to_df(ps)
library(microbiome)
# data(peerj32)
# p <- boxplot_abundance(peerj32$phyloseq, x='time', y='Akkermansia',
#                        line='subject')

seqtab.nochimudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVmodifiedSRRs.txt", header = TRUE)
seqtab.nochimudat <- data.frame(seqtab.nochimudat)
seqtab.nochimudat1 <- seqtab.nochimudat[,-1]
rownames(seqtab.nochimudat1) <- seqtab.nochimudat[,1]
seqtab.nochimudat1 <- t(seqtab.nochimudat1)

taxasilvaudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile_R.txt",header = TRUE)
taxasilvaudat <- data.frame(taxasilvaudat)
taxasilvaudat1 <- taxasilvaudat[,-1]
rownames(taxasilvaudat1) <- taxasilvaudat[,1]
taxasilvaudat1 <- as.matrix(taxasilvaudat1)

Urbmetatoremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbmetremove42samp.txt',header=TRUE)
samdfudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMet--.txt",header = TRUE)
samdfudat <- data.frame(samdfudat)
healthy <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofhea.txt",header = TRUE)
healthy <- data.frame(healthy)
ben <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofben.txt",header = TRUE)
ben <- data.frame(ben)
ca <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofca.txt",header = TRUE)
ca <- data.frame(ca)
samdfudat$Types <- 0
for (i in 1:length(Urbmetatoremove[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
  if (rownames(samdfudat)[j] == Urbmetatoremove[i,1]){
    samdfudat$Types[j] <- Urbmetatoremove[i,2]
  }
  }
}
samdfudat$bencahea <- 0
for (i in 1:length(ben[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == ben[i,1]){
      samdfudat$bencahea[j] <- ben[i,2]
    }
  }
}
for (i in 1:length(healthy[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == healthy[i,1]){
      samdfudat$bencahea[j] <- healthy[i,2]
    }
  }
}
for (i in 1:length(ca[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == ca[i,1]){
      samdfudat$bencahea[j] <- ca[i,2]
    }
  }
}
for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,2] == 1){
      samdfudat$Sample_Type[j] <- "Breast_Tumor"
  }
}
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,2] == 0){
    samdfudat$Sample_Type[j] <- "Healthy"
  }
}
samdfudat <- data.frame(samdfudat)
View(samdfudat)

#samdfudat$Types <- c(rep('BT',25),rep('H',5),rep('BT',1),rep('H',4),rep('BT',1),rep('H',3),rep('BT',3))
# install.packages("remotes")
# remotes::install_github("KarstensLab/microshades")
# install.packages("cowplot")
# remotes::install_github("mikemc/speedyseq")
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microshades)
library(speedyseq)

psudat <- phyloseq(otu_table(seqtab.nochimudat1, taxa_are_rows=FALSE),
                   sample_data(samdfudat), # originally sample_data
                   tax_table(taxasilvaudat1))
psudat <- prune_samples(sample_names(psudat) != "Mock", psudat)
dnaudat <- Biostrings::DNAStringSet(taxa_names(psudat))
names(dnaudat) <- taxa_names(psudat)
psudat <- merge_phyloseq(psudat, dnaudat)
taxa_names(psudat) <- paste0("ASV", seq(ntaxa(psudat)))
psudat
ps.propudat <- transform_sample_counts(psudat, function(otu) otu/sum(otu))
ord.nmds.brayudat <- ordinate(ps.propudat, method="NMDS", distance="bray", color="Type")
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE)) [1:100]
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)
psudatg <- tax_glom(ps.top100udat, "Genus")
ps0udatg <- transform_sample_counts(psudatg, function(x) x / sum(x))
#ps1udatg <- merge_samples(ps0udatg, "Types") yields NAs in columns of samdfudat within ps1udatg
ps2udatg <- transform_sample_counts(ps0udatg, function(x) x / sum(x))
pudatg <- plot_bar(ps2udatg, x= "Types", fill="Genus") + facet_wrap(~bencahea, scales = "free_x")
finalplotudatg <- pudatg + theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))
finalplotudatg

boxplot_abundance(d=ps.top100udat, x="Types", y='Bacillus')

ps_top100_udat <- merge_samples(ps.top100udat, "Types")
# Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
ps100_prepudat <- prep_mdf(ps_top100_udat)

# Create a color object for the specified data
color_ps100udat <- create_color_dfs(ps100_prepudat, group_level = "Phylum", subgroup_level = "Genus", cvd = TRUE, selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes'))

# Extract
mdf_psudat <- color_ps100udat$mdf
cdf_psudat <- color_ps100udat$cdf

plot_1 <- plot_microshades(mdf_psudat, cdf_psudat, group_label = "Phylum Genus")

plot_1 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=15)) +
  theme(axis.text.x = element_text(size= 10))


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

#write.table(ASVtab,"~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVmodifiedSRRs.txt", row.names = FALSE, col.names = FALSE)
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
# BiocManager::install("ALDEx2")
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
        count <- c(count,j)
      } else{
        x1alddata$Taxa[i] <- nmj[length(nmj)]
        count <- c(count,j)
        #count <- c(count,length(nmj))
      }
    }
  }
}

Urbremben <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbremoveben.txt',header=TRUE)
Urbmetatoremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbmetremove42samp.txt',header=TRUE)
samdfudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMet--.txt",header = TRUE)
samdfudat <- data.frame(samdfudat)
healthy <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofhea.txt",header = TRUE)
healthy <- data.frame(healthy)
ben <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofben.txt",header = TRUE)
ben <- data.frame(ben)
ca <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofca.txt",header = TRUE)
ca <- data.frame(ca)
samdfudat$Types <- 0
for (i in 1:length(Urbmetatoremove[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (rownames(samdfudat)[j] == Urbmetatoremove[i,1]){
      samdfudat$Types[j] <- Urbmetatoremove[i,2]
    }
  }
}
for (i in 1:length(Urbremben[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat$Types[j] == Urbremben[i,1]){
      samdfudat <- samdfudat[-j,]
    }
  }
}


for (i in 1:length(healthy[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == healthy[i,1]){
      samdfudat$bencahea[j] <- healthy[i,2]
    }
  }
}
for (i in 1:length(ca[,1])) {
  for (j in 1:length(rownames(samdfudat))) {
    if (samdfudat[j,3] == ca[i,1]){
      samdfudat$bencahea[j] <- ca[i,2]
    }
  }
}
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,2] == 1){
    samdfudat$Sample_Type[j] <- "Breast_Tumor"
  }
}
for (j in 1:length(rownames(samdfudat))) {
  if (samdfudat[j,2] == 0){
    samdfudat$Sample_Type[j] <- "Healthy"
  }
}
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
# ps.propudat <- transform_sample_counts(psudat, function(otu) otu/sum(otu))
# ord.nmds.brayudat <- ordinate(ps.propudat, method="NMDS", distance="bray", color="Type")
top100udat <- names(sort(taxa_sums(psudat), decreasing=TRUE)) 
ps.top100udat <- transform_sample_counts(psudat, function(OTU) OTU/sum(OTU))
ps.top100udat <- prune_taxa(top100udat, ps.top100udat)
# psudatg <- tax_glom(ps.top100udat, "Genus")
# ps0udatg <- transform_sample_counts(psudatg, function(x) x / sum(x))
# ps1udatg <- merge_samples(ps0udatg, "bencahea") #yields NAs in columns of samdfudat within ps1udatg
# ps2udatg <- transform_sample_counts(ps1udatg, function(x) x / sum(x))

library(metagMisc)
propASVdata <- data.frame(ps.top100udat@otu_table)
propASVdata <- t(propASVdata)
colnames(propASVdata) <- c(rep("BT",15),rep("H",5),rep("BT",1),rep("H",4),rep("BT",1),rep("H",3),rep("BT",3))
# 1704/(sum(ASVtab211[,1]))
# 878/(sum(ASVtab211[,2]))
# 1098/(sum(ASVtab211[,3]))
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
  #filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than ttest welch
  dplyr::select(ASVs) %>%
  data.frame()
x1alddata1 <- data.frame(x1alddata1[1:32,1])
x1alddata1$zero <- rep(0,32)
# colnames(x1alddata1) <- "ASV"
rownames(x1alddata1) <- rownames(propASVdata)
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

# for (i in 1:length(colnames(newdatf))-1) {
#   # pdf(file= "samples1_.pdf", onefile = TRUE)
#   print(ggplot(newdatf,aes(x=tissuetyp,y=paste("ASV",count[i]))) +
#     geom_boxplot()) +
#     labs(x="",y="ASV")
# }
#asvs <- data.frame(asvs)
dataps2udatg <- psmelt(ps.top100udat)
export_ps <- phyloseq_to_df(ps.top100udat)
# dataps2udatg1 <- dataps2udatg[,-1]
# dataps2udatg1 <- t(dataps2udatg1)
# colnames(dataps2udatg1) <- dataps2udatg[,1]
# dataps2udatg1 <- t(dataps2udatg1)
#psgraph <- rep(1,length(x1alddata$Taxa))
ps2graphslist <- list()
# plots <- list()  # new empty list
# for (i in 1:6) {
#   p1 = qplot(1:10, rnorm(10), main = i)
#   plots[[i]] <- p1  # add each plot into plot list
# }
# multiplot(plotlist = plots)
library(ggplot2)
for (i in 1:length(x1alddata$Taxa)){
 # pdf(file= "samples11.pdf", onefile = TRUE)
  ps2data <- column_to_rownames(dataps2udatg,"OTUs") %>%
    filter(Genus == x1alddata$Taxa[i])  %>% # here we chose the wilcoxon output rather than ttest welch
    filter(OTU == asvs[i])
    dplyr::select(OTUs,Abundance, Sample, Genus, OTU) %>%
    data.frame()
 # ps2data %>%
  ggplot(data=ps2data,aes(x= Sample, y= Abundance, fill= Genus)) +
  geom_boxplot()
}
  #ps2graphslist[[i]] <- p1

  # ps2graphslist[i] <- p
  #ggsave("arrange.pdf", arrangeGrob(grobs = l), device = "pdf",limitsize = FALSE)
  
 # ggsave("sigUrbgraphboxplot.pdf",limitsize = FALSE)
  # pdf("plots.pdf")
  # for (i in 1:20) {
  #   print(ps2graphslist[[i]])
  # }
  # dev.off()
#}


# ps2data <- rownames_to_column(dataps2udatg,"OTUs") %>%
#   filter(Genus == "Bacillus")  %>% # here we chose the wilcoxon output rather than ttest welch
#   dplyr::select(OTUs, Abundance, Sample_Type, Genus) %>%
#   data.frame()
# 
# ps2data %>%
#   ggplot( aes(x=Sample_Type, y=Abundance, fill=Genus)) +
#   geom_boxplot() +

# for (i in 1:length(ps2data$Sample_Type)){
#   if (ps2data$Sample_Type[i] == "Breast_Tumor") {
#     ps2data$Sample_Type[i] <- 2
#   } else {
#     if (ps2data$Sample_Type[i] == "Healthy") {
#       ps2data$Sample_Type[i] <- 1
#   }
#   }
# }
# #for (i in 1:length(ps2data$Sample_Type)){
# ps2data$Sample_Type <- as.numeric(ps2data$Sample_Type)
# #}
# for (i in 1:length(ps2data$Sample_Type)){
#   ps2data$Abundance[i] <- as.numeric(ps2data$Abundance[i])
# }
# boxplot(x=ps2data$Sample_Type,y=ps2data$Abundance)


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
library(compositions)
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

