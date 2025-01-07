# bc_microbiome-project
## Sidra 
## Project - Analysis of Breast Cancer Microbial Datasets
#### Background
The microbiome is a collection of bacteria that reside in tissues where in cancer it has been implicated in a variety of tissues, and specifically distinct microbial communities have been linked to different types of cancers [[Hieken et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4971513/)]. The association of microbial communities to human health has prompted researchers to investigate the microbial ecology of various tissues leading to a surplus of microbial data. This in turn has led to the development of various software and tools that can analyze microbial data and disentangle the biological variation from amplicon sequencing errors [[Callahan et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)]. Here, I analyze data from _The Microbiome of Aseptically Collected Human Breast Tissue in Benign and Malignant Disease_ study by [Hieken et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4971513/),  _The Microbiota of Breast Tissue and Its Association with Breast Cancer_ study by [Urbaniak et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4968547/), and _Characterization of the microbiome of nipple aspirate fluid of breast cancer survivors_ study by [Chan et al.](https://www.nature.com/articles/srep28061) by running them through DADA2 [[Callahan et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)]. The data can be downloaded from the [SRA](https://www.ncbi.nlm.nih.gov/sra) using the following accession numbers:  
Hieken et al. PRJNA335375  
Urbaniak et al. SRP076038  
Chan et al. PRJNA314877

#### Packages and software to install before running
R: Bioconductor, dada2, phyloseq, Biostrings, ggplot2, DECIPHER, phangorn, hopach
### Files in repo
- Chan folder:
  - ChanRscript_upd.R: R script for [DADA2](https://benjjneb.github.io/dada2/tutorial.html), phylogenetic, and microshades analyses
  - graphs: proportional abundance plots from phyloseq and MicrobiomeAnalyst
- Hieken folder:
  - HiekenAnalysis_upd.R: R script for [DADA2](https://benjjneb.github.io/dada2/tutorial.html), phylogenetic, and microshades analyses
  - graphs: proportional abundance plots from phyloseq and MicrobiomeAnalyst
- Urbaniak folder:
  - UrbaniakAnalysis_upd.R: R script for [DADA2](https://benjjneb.github.io/dada2/tutorial.html), phylogenetic, and microshades analyses
  - graphs: proportional abundance plots from phyloseq

### How to use
Clone repository into personal directory using this command,  
```
git clone https://github.com/ssohail1/bc_microbiome-project.git
```


To move into bc_microbiome-project directory use `cd`,  
```
cd bc_microbiome-project
```
