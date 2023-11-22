# BeeeDNA
Data Analysis Scripts for Sickel et al. 'BEE-quest of the nest: A novel method for eDNA-based, nonlethal detection of cavity-nesting hymenopterans and other arthropods' (in press)

## Study design (short version)
105 vacated nest tubes of cavity nesting Hymenoptera were selected, covering nests of different sizes, parasitised and otherwise mixed species nests.
We compared two sample types, obtained from those nests (fecal pellets and swab heads). 
DNA was extracted and two different COI fragments were amplified, referred to as COI.short (~ 205 bp; primers fwhF2&fwhR2n ([Elbrecht et al., 2019](https://peerj.com/articles/7745/))) and COI.long (~ 313 bp; primers mlCOIintF&HCO2198 ([Folmer et al. 1994](https://pubmed.ncbi.nlm.nih.gov/7881515/); [Leray et al. 2013](https://frontiersinzoology.biomedcentral.com/articles/10.1186/1742-9994-10-34))).
We followed a dual-tagging strategy adapted from [Herbold et al., 2015](https://www.frontiersin.org/articles/10.3389/fmicb.2015.00731/full) and [Elbrecht & Steinke, 2019](https://onlinelibrary.wiley.com/doi/full/10.1111/fwb.13220) (upload tag file), samples were amplified in duplicate.


## Third party tools
+ [BCdatabaser](https://github.com/molbiodiv/bcdatabaser) 
+ [MIDORI UNIQ GB 241](http://www.reference-midori.info/download.php#)
+ [VSEARCH](https://github.com/torognes/vsearch)
+ [cutadapt](https://github.com/marcelm/cutadapt)
+ R (see used packages further down)

## Overview
The workflow contains of two bash scripts and one config file. The config file contains information about used sample-tags, primer sequences, reference databases, etc.

## Main steps
### Preparation of raw data ('_prepRawDataCOI.sh')
1. Sort by fragment size (short & long), via primer sequences, remove duplicates between long & short fragment
2. Demultiplexing and trimming of primer (+tag) sequences

### Further processing ('_processing_COI.sh')
1. Merge read pairs
2. Quality filtering
3. Dereplication
4. Denoising
5. Chimera checking
6. Make community table
7. Taxonomic classification

## Data analysis in R
see file R_script_Analysis_BeeeDNA.R

### 0. load required packages
+ tidyverse (v.1.3.2)
+ phyloseq (v.1.38.0)
+ rstatix (v.0.7.0)
+ PMCMRplus (v.1.9.3) 
+ broom (v.1.0.1)
+ parsnip (v.1.0.1) 
+ yardstick (v.1.1.0)
+ cowplot  (v.1.1.1)
+ EnvStats (v.2.7.0)
+ ggsignif (v.0.6.3)
### 1. Data import
+ import data containing all samples, both fragments
+ clean up data, e.g. converting & ordering categorical variables
### 2. Analysis of DNA quantity & quality, short&long fragment
#### DNA quantity
+ determine standard metrics such as mean, standard deviation, minimum, maximum of DNA yields, for all samples, but also for feces and swabs separately, as well as different nest sizes
+ test for normal distribution using Shapiro-Wilk test (not normally distributed)
+ test for differences between sample types & nest sizes regarding DNA yields using Kruskal-Wallis test, posthoc-testing: Kruskal-Nemenyi’s All-Pairs Rank 
+ visualisation via violin and scatter plot
Comparison Test
#### DNA quality
+ DNA quality = visible band in agarose gel electrophoresis after 1st PCR
+ analysed for short and long fragment separately
+ determined number of succsessfully amplified samples as well as numbers of cases where duplicate samples were in agreement (all samples, sample types, nest sizes)
### 3. Import and preparation of sequencing results, short fragment only
+ import community table (ASV-based) & taxonomy table, merge with sample data (already imported) to *phyloseq*-object
+ subset to Metazoa, analyse laboratory controls (remove contaminant taxa, positive control taxon (*Apis mellifera*) and control samples)
+ agglomerate to species / genus / family
+ subset to Arthropoda / Hymenoptera & Diptera / Hymenoptera, remove unclassified species
### 4. Analysis of sequencing results
#### 4.1. Experiment 1: Detection of Hymenoptera
+ determine richness estimates, count Hymenoptera detections, count detections of mixed species
+ visualise Hymenoptera detection rates as bar plot
+ compare Hymenoptera detection rates between sample types, nest sizes, etc. using Fisher's exact test
+ analyse generalised linear models (categorical variables converted to numerical), to assess whether explanatory variables, DNA yield and/or electrophoresis can be used to predict Hymenoptera detection rates
#### 4.2. Experiment 2: Mixed species detection
+ similar to 4.1., but counting cases where number of species detected > 1 (Hymenoptera & Diptera)
+ no visualisation
#### 4.3. Validation: Validation via morphological data
+ import and clean up morphological data
+ compare detected taxa across all samples / different taxonomic levels
+ sample-wise comparison
+ visualise as bar plots
