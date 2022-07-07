---
title: "Calling Chromosomal Copy Number Changes"
author: "Madyson Colton"
date: "22/04/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# RNAseqCNV 

Implementing the RNAseqCNV caller into the SAHMRI bioinformatics pipeline following variant calling. This tool can be used to determine aneuploidy in samples missing karyotype data.


## Installing on the HPC

R is required to run RNAseqCNV. An virtual environment was created to install the latest version of R and related packages.

```{bash Renv, eval=F, echo=T}
 mamba create --name Renv r-base r-essentials
 mamba activate Renv
 mamba install -c bioconda bioconductor-biocinstaller

```

To install RNAseqCNV and dependencies, run the installCNV.R script. The RNAseqCNV package is stored on GitHub, so installation of devtools is required for installation. 

```{bash devtools, eval=F, echo=T}
mamba install -c conda-forge r-devtools

```

```{bash install command, eval=F, echo=T}
Rscript installCNV.R

```


```{r installCNV.R script, eval=F, echo=T}

library(devtools)

devtools::install_github(repo = "honzee/RNAseqCNV")

```


Some other dependent packages may have to be installed via conda. For me, ggpubr was causing the installation to fail.

```{bash ggpubr, eval=F, echo=T}
mamba install -c conda-forge r-ggpubr

```

## Locate and Tidy Data

To determine CNVs, the tool requires both sequencing counts and SNV frequencies.

The counts data should have only 2 columns - the ensembl IDs and the read counts - and no header.
The SNV information can be provided as a .csv file, or as a VCF (or vcf.gz) file following the GATK pipeline. This file should contain headers, and RNAseqCNV will require the chromosome number, start position, read depth and mutate allele frequency.


The RNAseqCNV wrapper then requires the following:

* Config file
  + Directory to store the output files
  + Directory containing the counts files
  + Directory containing the vcf files
  
* Metadata file, containing 3 columns with no header
  + Sample ID
  + Path to that sample's counts file
  + Path to that sample's vcf file
  

### Config file

An example of the config file is given below. The variable names must not be changed. The out_dir path should be the same as the outDir for the variant calling pipeline, and should end in a "/"; the count_dir and snv_dir paths should not have a "/" at the end. 

Due to the nature of RNAseqCNV, a config file will have to be generated for each sample in order for the output to be saved in the same output directory as the variant caller results.

```{bash config file, eval=F, echo=T}
out_dir = "~/local_NF_ALL_runs/variants/outDir/"
count_dir = "~/local_NF_ALL_runs/variants/outDir"
snv_dir = "~/local_NF_ALL_runs/variants/outDir"


```


### Count data

The counts data need to be converted from a ".rData" object to a ".csv" file with two columns. Following the variant calling pipeline, the output of featureCounts can be used to extract this information.

### Metadata file

For the variant calling pipeline, a sample sheet is provided. This sample sheet can be used to generate the metadata file.

#### Sample Sheet

```{r samplesheet example, echo = FALSE, warning=FALSE}
sampleSheet <- read.csv("Test/samplesheet.example.csv")
knitr::kable(sampleSheet, caption = "An example of the samplesheet")

```


#### Metadata Sheet


```{r metadata example, echo = FALSE}
metaTable <- read.csv("metaData.csv")
knitr::kable(metaTable, caption = "An example of the metadata sheet")

```


### Running RNAseqCNV

To run the tool, the path to the config and metadata files should be provided and specified within the RNAseqCNV_wrapper() function. The format of the variant file should also be specified ("vcf") and the genome version (typically "hg19" in our pipeline). Here's the basic requirements:

```{r RNAseqCNV.R, eval=F, echo=T}

library(RNAseqCNV)

configPath <- "~/local_NF_ALL_runs/variants/cnv_config"
metaPath <- paste0("~/local_NF_ALL_runs/variants/cnv_metadata", date, ".csv")


RNAseqCNV_wrapper(config = configPath,
                  metadata = metaPath,
                  snv_format = "vcf",
                  genome_version = "hg19",
                  arm_lvl = FALSE)

```


# CNV Caller for Pipeline

All the steps required to run RNAseqCNv on a new sample have been included in a script to follow the variant calling pipeline. The full script is shown below.

The script creates one metadata table based on the variant calling samplesheet. The only input field required to be changed is the date in the format YYYYMMDD; this is the suffix for the samplesheet and will form the suffix for the metadata sheet.

A new RNAseqCNV directory is then created within each individual sample output directory. The counts are extracted for each sample and saved in this output directory. A config file is also created for each sample and saved in this directory. 

Finally, RNAseqCNV is ran for each sample.

```{r cnv_caller.R, message=FALSE, warning=FALSE, eval = FALSE}

# # # # # # # # # # # # # # # # # # # # # # # # #
#            Calling Chromosomal CNV            #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Date for file name suffix in format YYYYMMDD
# This is the only field to edit each time
date <- "220503"

# Load sample sheet file (make sure this is your own working dir)
sampleSheet <- read.csv(paste0("~/local_NF_ALL_runs/variants/samplesheet_", date, ".csv"))

# Load libraries
library(RNAseqCNV)
library(dplyr) 

# Creates an RNAseqCNV output directory for each sample
# Creates counts sheet and config file for each sample
# Runs RNAseqCNV for each sample

for(i in 1:nrow(sampleSheet)){
  samplePath <- sampleSheet[i,]
  
  # Create RNAseqCNV output directory
  dir.create(paste0("outDir/",samplePath$group, "/", samplePath$filename,"/RNAseqCNV"), 
             recursive = TRUE)
  
### Metadata ###
  
  # Generate path to vcf file
  metaData <- samplePath %>% 
    mutate(varNames = paste0(samplePath$group, "/", samplePath$filename, "/Variants/", samplePath$filename, ".raw.vcf.gz")) %>%
    mutate(countNames = paste0(samplePath$group, "/", samplePath$filename, "/RNAseqCNV/sampleCounts.csv"))
  
  # Generate metadata table
  metaData <- metaData %>% 
    select(c("filename", "countNames", "varNames"))
  
  # Save table in working directory
  write.table(metaData, 
              paste0("outDir/",samplePath$group, "/", samplePath$filename,"/RNAseqCNV/cnv_metadata.csv"), 
              col.names = FALSE, 
              sep=",", 
              row.names = FALSE, 
              quote = FALSE, 
              eol = "\r")
  
  print(paste0("Metadata table generated for ", samplePath$filename))
  
### Counts ###
  
  # Use sample sheet to find featureCounts output 
  sampleCounts <- get(load(paste0("outDir/",samplePath$group, "/", samplePath$filename, "/featureCounts/counts.Rdata")))
 
  # Extract gene counts
  sampleCounts <- data.frame(sampleCounts$counts) 
  sampleCounts$genes <- rownames(sampleCounts)
  sampleCounts <- sampleCounts[,c(2,1)]
  
  # Save as csv file
  write.table(sampleCounts, 
              paste0("outDir/",samplePath$group, "/", samplePath$filename,"/RNAseqCNV/sampleCounts.csv"), 
              col.names = FALSE, 
              sep=",", 
              row.names = FALSE, 
              quote = FALSE, 
              eol = "\r")
  
  print(paste0("Counts generated for ", samplePath$filename))

### Configs ###
  
  # Create a config file for each sample in the RNAseqCNV directory
  outPath <- paste0("outDir/",samplePath$group, "/", samplePath$filename,"/RNAseqCNV/")
  
  cat(paste0("out_dir = \"~/local_NF_ALL_runs/variants/", outPath, "\""),
      file=paste0(outPath, "config"),sep="\n")
  
  cat(paste0("count_dir = \"outDir\""),
      file=paste0(outPath, "config"),append=TRUE, sep="\n")
  
  cat(paste0("snv_dir = \"outDir\""),
      file=paste0(outPath, "config"),append=TRUE, sep="\n") 
  
  print(paste0("Config generated for ", samplePath$filename))
  
### RNAseqCNV ###

  configPath <- paste0("outDir/",samplePath$group, "/", samplePath$filename,"/RNAseqCNV/config")
  metaPath <- paste0("outDir/",samplePath$group, "/", samplePath$filename,"/RNAseqCNV/cnv_metadata.csv")
  
  RNAseqCNV_wrapper(config = configPath,
                  metadata = metaPath,
                  snv_format = "vcf",
                  genome_version = "hg19",
                  arm_lvl = TRUE)
  
  print(paste0("RNAseqCNV complete for ", samplePath$filename))
}

print("RNAseqCNV complete.")

```



