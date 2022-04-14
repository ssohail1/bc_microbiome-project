# # Urbaniak Taxa table
# Urbaniaktaxa <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile.txt",header = TRUE)
# Urbaniaktaxa <- data.frame(Urbaniaktaxa)
# Urbaniaktaxa1 <- Urbaniaktaxa[,-1]
# rownames(Urbaniaktaxa1) <- Urbaniaktaxa[,1]
# Urbaniak ASV table
seqtab.nochimudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbASVmodifiedSRRs.txt", header = TRUE)
seqtab.nochimudat <- data.frame(seqtab.nochimudat)
# seqtab.nochimudat1 <- seqtab.nochimudat[,-1]
# rownames(seqtab.nochimudat1) <- seqtab.nochimudat[,1]
# seqtab.nochimudat1 <- t(seqtab.nochimudat1)
# Urbaniak taxa table
taxasilvaudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbtaxafile_R.txt",header = TRUE)
taxasilvaudat <- data.frame(taxasilvaudat)
# taxasilvaudat1 <- taxasilvaudat[,-1]
# rownames(taxasilvaudat1) <- taxasilvaudat[,1]
# taxasilvaudat1 <- as.matrix(taxasilvaudat1)
# Urbaniak metadata table
# Urbmetatoremove <- read.table('~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/Urbmetremove42samp.txt',header=TRUE)
# samdfudat <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/UrbaniakMet--.txt",header = TRUE)
# samdfudat <- data.frame(samdfudat)
# healthy <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofhea.txt",header = TRUE)
# healthy <- data.frame(healthy)
# ben <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofben.txt",header = TRUE)
# ben <- data.frame(ben)
# ca <- read.table("~/Documents/Urbaniak_02152022_version/MicrobAnalystInputs/listofca.txt",header = TRUE)
# ca <- data.frame(ca)
samdfudat <- read.table(file= "~/Downloads/urbaniakmeta.txt")
samdfudat <- data.frame(samdfudat)
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

# Hieken ASV table
HiekenASVtab <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/seqtabnochimhieken10082021.txt',header = TRUE)
HiekenASVtab <- data.frame(HiekenASVtab)
# HiekenASVtab2 <- HiekenASVtab[,-1]
# rownames(HiekenASVtab2) <- HiekenASVtab[,1]
# Hieken Taxa table
Hiekentaxa <- read.table('~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/taxasilvahieken10082021.txt',header=TRUE)
Hiekentaxa <- data.frame(Hiekentaxa)
# Hiekentaxa1 <- Hiekentaxa[,-1]
# rownames(Hiekentaxa1) <- Hiekentaxa$TAXONOMY
# Hieken metadata table
HMeta <- read.table("~/Documents/Hieken10082021/MicrobiomeAnalyst_Inputs/HMetadata_Tissues.txt")
HMeta <- data.frame(HMeta)
# HMetaHiek <- HMeta[,-1]
# rownames(HMetaHiek) <- HMeta[,1]
colnames(HMeta) <- c("sample_id","sample","sample_type")
colnames(samdfudat) <- c("sample_id","sample","tissue","sample_name","sample_type")

library(dplyr)
library(phyloseq)
HiekUrbASVcombined <- full_join(HiekenASVtab,seqtab.nochimudat, by = "NAME")
HiekUrbASVcombined <- data.frame(HiekUrbASVcombined)
HiekUrbASVcombined1 <- HiekUrbASVcombined[,-1]
rownames(HiekUrbASVcombined1) <- HiekUrbASVcombined[,1]
HiekUrbASVcombined1 <- t(HiekUrbASVcombined1)

HiekUrbTAXAcombined <- full_join(Hiekentaxa,taxasilvaudat, by = "Taxonomy")
HiekUrbTAXAcombined <- data.frame(HiekUrbTAXAcombined)
HiekUrbTAXAcombined1 <- HiekUrbTAXAcombined[,-1]
rownames(HiekUrbTAXAcombined1) <- HiekUrbTAXAcombined[,1]
HiekUrbTAXAcombined1 <- as.matrix(HiekUrbTAXAcombined1)

HiekUrbMETADATAcombined <- full_join(HMeta,samdfudat,by = "sample_id")
HiekUrbMETADATAcombined <- data.frame(HiekUrbMETADATAcombined)
HiekUrbMETADATAcombined1 <- HiekUrbMETADATAcombined[,-1]
rownames(HiekUrbMETADATAcombined1) <- HiekUrbMETADATAcombined[,1]

psudat <- phyloseq(otu_table(HiekUrbASVcombined1, taxa_are_rows=FALSE),
                   sample_data(HiekUrbMETADATAcombined1), # originally sample_data
                   tax_table(HiekUrbTAXAcombined1))






match <- list()
count1 <- 0
for (i in 1:length(rownames(Urbaniaktaxa1))) {
  for (j in 1:length(rownames(Hiekentaxa1))) {
    if (rownames(Urbaniaktaxa1)[i] == rownames(Hiekentaxa1)[j]) {
      match <- c(match, rownames(Urbaniaktaxa1)[i])
    } else {
      count1 <- count1 + 1
    }
  }
}


hiekenGreengenes <- read.table("/media/mbb/Sidras_Projects/Hieken_paper/HiekenGreengenes-joininginR.txt", header = TRUE)
hiekenSilva <- read.table("/media/mbb/Sidras_Projects/Hieken_paper/HiekenSilva-joininginR.txt", header = TRUE)
urbaniakGreengenes <- read.table("/media/mbb/Sidras_Projects/Urbaniak_paper/UrbaniakGreengenes-joininginR.txt", header = TRUE)
urbaniakSilva <- read.table("/media/mbb/Sidras_Projects/Urbaniak_paper/UrbaniakSilva-joininginR.txt", header = TRUE)

library(dplyr)
file <- inner_join(urbaniakGreengenes, urbaniakSilva, by="Taxa") %>%
  group_by(Taxa)
#this one looks like the correct file
file_one <- full_join(urbaniakGreengenes, urbaniakSilva)
file_two <- full_join(hiekenGreengenes, hiekenSilva)
write.table(file_one, file = "/media/mbb/Sidras_Projects/Urbaniak_paper/urbaniaktaxajoin.txt")
write.table(file_two, file = "/media/mbb/Sidras_Projects/Hieken_paper/hiekentaxajoin.txt")

# join formats - info

# inner_join(): includes all rows in x and y.
# left_join(): includes all rows in x.
# right_join(): includes all rows in y.
# full_join(): includes all rows in x or y.

# if don't specify common column name then will perform a natural join, 
# using all variables in common across x and y. 
# A message lists the variables so that you can check they're correct 
# suppress the message by supplying by explicitly.
# To join by different variables on x and y, use a named vector. 
# For example, by = c("a" = "b") will match x$a to y$b.
# To join by multiple variables, use a vector with length > 1. 
# For example, by = c("a", "b") will match x$a to y$a and x$b to y$b. 
# Use a named vector to match different variables in x and y. 
# For example, by = c("a" = "b", "c" = "d") will match x$a to y$b and x$c to y$d.
# To perform a cross-join, generating all combinations of x and y, use by = character().
inner_join(x,y,by="common column name")



