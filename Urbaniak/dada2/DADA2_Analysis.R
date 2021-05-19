### NOTE: the Urbaniak data has merged forward and reverse reads so only using the 
### forward read notation as there are no reverse reads

# Set the working directory
library(dada2)
setwd("/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/") #change to directory where files are located
dir()
path <- "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/"
list.files(path)
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# To get the quality profile of the forward reads:
plotQualityProfile(fnFs[1:2])

# Assigning filenames for the filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
names(filtFs) <- sample.names

#minLen = 50 so that there are at least 50 nucleotides in all samples for the assignTaxonomy command
#truncLen is omitted as the merged reads are pre-filtered
out <- filterAndTrim(fnFs, filtFs,
                     maxN=0, maxEE= 3, truncQ=2, minLen = 50, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


head(out)

# Learn the Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
##63051027 total bases in 839849 reads from 68 samples will be used for learning the error rates.

#Plot errors to visualize easily
plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Applying the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#Sample 1 - 28252 reads in 3312 unique sequences.
#Sample 2 - 24050 reads in 1738 unique sequences.
#Sample 3 - 5140 reads in 1302 unique sequences.
#Sample 4 - 5114 reads in 1298 unique sequences.
#Sample 5 - 6931 reads in 1625 unique sequences.
#Sample 6 - 5866 reads in 1334 unique sequences.
#Sample 7 - 36562 reads in 5019 unique sequences.
#Sample 8 - 2498 reads in 706 unique sequences.
#Sample 9 - 12784 reads in 1786 unique sequences.
#Sample 10 - 11652 reads in 2223 unique sequences.
#Sample 11 - 10662 reads in 1955 unique sequences.
#Sample 12 - 10320 reads in 2138 unique sequences.
#Sample 13 - 7452 reads in 1578 unique sequences.
#Sample 14 - 2270 reads in 681 unique sequences.
#Sample 15 - 3332 reads in 717 unique sequences.
#Sample 16 - 4392 reads in 1052 unique sequences.
#Sample 17 - 6137 reads in 1357 unique sequences.
#Sample 18 - 5947 reads in 1420 unique sequences.
#Sample 19 - 5206 reads in 1117 unique sequences.
#Sample 20 - 4660 reads in 1171 unique sequences.
#Sample 21 - 8366 reads in 1675 unique sequences.
#Sample 22 - 3534 reads in 777 unique sequences.
#Sample 23 - 4232 reads in 997 unique sequences.
#Sample 24 - 5016 reads in 1290 unique sequences.
#Sample 25 - 5271 reads in 1253 unique sequences.
#Sample 26 - 4607 reads in 1026 unique sequences.
#Sample 27 - 4799 reads in 1153 unique sequences.
#Sample 28 - 23086 reads in 1918 unique sequences.
#Sample 29 - 9518 reads in 1670 unique sequences.
#Sample 30 - 20969 reads in 3269 unique sequences.
#Sample 31 - 7239 reads in 2065 unique sequences.
#Sample 32 - 6800 reads in 1431 unique sequences.
#Sample 33 - 10394 reads in 2322 unique sequences.
#Sample 34 - 13393 reads in 1670 unique sequences.
#Sample 35 - 16193 reads in 2297 unique sequences.
#Sample 36 - 49182 reads in 4123 unique sequences.
#Sample 37 - 8002 reads in 1532 unique sequences.
#Sample 38 - 5497 reads in 1459 unique sequences.
#Sample 39 - 21982 reads in 2853 unique sequences.
#Sample 40 - 3239 reads in 810 unique sequences.
#Sample 41 - 69984 reads in 12888 unique sequences.
#Sample 42 - 42815 reads in 8432 unique sequences.
#Sample 43 - 22717 reads in 5272 unique sequences.
#Sample 44 - 23675 reads in 4746 unique sequences.
#Sample 45 - 29291 reads in 6211 unique sequences.
#Sample 46 - 3092 reads in 879 unique sequences.
#Sample 47 - 4841 reads in 1433 unique sequences.
#Sample 48 - 8226 reads in 2231 unique sequences.
#Sample 49 - 14872 reads in 3463 unique sequences.
#Sample 50 - 12757 reads in 2300 unique sequences.
#Sample 51 - 15070 reads in 2795 unique sequences.
#Sample 52 - 10976 reads in 2350 unique sequences.
#Sample 53 - 6798 reads in 1464 unique sequences.
#Sample 54 - 11316 reads in 2263 unique sequences.
#Sample 55 - 6099 reads in 1300 unique sequences.
#Sample 56 - 5404 reads in 1305 unique sequences.
#Sample 57 - 17585 reads in 2193 unique sequences.
#Sample 58 - 4174 reads in 1065 unique sequences.
#Sample 59 - 5955 reads in 1427 unique sequences.
#Sample 60 - 7160 reads in 1703 unique sequences.
#Sample 61 - 5734 reads in 1648 unique sequences.
#Sample 62 - 5076 reads in 1239 unique sequences.
#Sample 63 - 26723 reads in 4933 unique sequences.
#Sample 64 - 20099 reads in 4889 unique sequences.
#Sample 65 - 17398 reads in 3830 unique sequences.
#Sample 66 - 6175 reads in 1482 unique sequences.
#Sample 67 - 2534 reads in 726 unique sequences.
#Sample 68 - 2757 reads in 806 unique sequences.

dadaFs[[1]]
#dada-class: object describing DADA2 denoising results
#317 sequence variants were inferred from 3312 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Constructing sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
##[1]   68 7021

table(nchar(getSequences(seqtab)))
#51   53   54   55   56   57   58   63   65   66   67   68   69   70   71   72   73   74   75 
#2    6    1    7    6    2    1    2    1    1    1    6   11   14  117  430  621 1090 2116 
#76   77   78   79   80   81   82   83   84   85   86   87   88   89   90   91   92   93   94 
#405 1094  296  418   95   61   92   24    8   11   13    6   15    1    1    6    1    3    5 
#95   96  100  101  102 
#14   14    1    1    1 


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 78 bimeras out of 7021 input sequences.
#minFoldParentOverAbundance not required here

sum(seqtab.nochim)/sum(seqtab)
##[1] 0.9938704  ~ should be close to 1

rowSums(seqtab.nochim)
rowSums(seqtab)

# Track reads thru pipeline
getN <- function(x) sum(getUniques(x))

# According to Dada2 tutorial, make the following modification to the track command: 
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) 
#with getN(dadaFs)
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)


# Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/silva.nr_v138.align", multithread = TRUE)
## Went through!

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# Did not evaluate accuracy with mock data as no mock samples in our data

# Phyloseq
## Need to download phyloseq first
##if (!requireNamespace("BiocManager", quietly = TRUE))
##   install.packages("BiocManager")
##  BiocManager::install("phyloseq")
## Then load library
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
# Steps to create a metadata (samdf file)
# Already have a metadata - so can skip this
theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
Sample_Name <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(Sample_Name,1,1)
subject <- substr(Sample_Name,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(sample_name=Sample_Name, Name=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# Already have metadata so can read that in
samdf <- read.table("/media/mbb/Sidras_Projects/Urbaniak_paper/fastqfiles/UrbaniakMetadata.txt", header=TRUE)
View(samdf)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf), # originally sample_data
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
#Shannon-Simpson Alpha-diversity plot
plot_richness(ps, x = "Isolation_source", measures=c("Shannon", "Simpson"), color = "Sample_Type")
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray", color="Sample_Type")
#Bray NMDS Beta-diversity plot
plot_ordination(ps.prop, ord.nmds.bray, color="Sample_Type", title="Bray NMDS")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
#Can have the fill = be at the Genus level as well
plot_bar(ps.top20, x = "Sample_Type", fill = "Family") + facet_wrap(~Isolation_source, scales = "free_x")
plot_bar(ps.top20, x = "Sample_Type", fill = "Phylum") + facet_wrap(~Isolation_source, scales = "free_x")

# End of dada2 Tutorial
