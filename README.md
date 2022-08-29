# BeeeDNA
Data Analysis Scripts for Sickel et al. 'eDNA allows non-lethal detection of cavity nesting wild bees and their parasitoids'

## Study design (short version)
105 vacated nest tubes of cavity nesting Hymenoptera collected as part the F.R.A.N.Z.-project https://www.franz-projekt.de/website/english-summary were selected, covering nests of different sizes, parasitised and otherwise mixed species nests.
We compared two sample types, obtained from those nests (fecal pellets and swab heads). 
DNA was extracted and two different COI fragments were amplified, referred to as COI.short (~ 302 bp; primers fwhF2&fwhR2n (Elbrecht et al., 2019)) and COI.long (~ 413 bp; primers mlCOIintF&HCO2198 (Folmer et al. 1994; Leray et al. 2013)).
We followed a dual-tagging strategy adapted from Herbold et al., 2015 and Elbrecht&Steinke, 2018 (upload tag file), samples were amplified in duplicate.


## Third party tools
+ COI via bc-databaser 
+ MIDORI
+ VSEARCH & cutadapt
+ R (packages: dplyr, vegan, phyloseq)

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
