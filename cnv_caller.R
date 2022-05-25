# # # # # # # # # # # # # # # # # # # # # # # # #
#           Calling Chromosomal CNVs            #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Date for file name suffix in format YYYYMMDD
# This is the only field to edit
date <- "220525"

# Load sample sheet file
sampleSheet <- read.csv(paste0("~/local_NF_ALL_runs/variants/samplesheet_", date, ".csv"))

# Load libraries
library(RNAseqCNV)
library(dplyr) 

# Creates an RNAseqCNV directory for each sample
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
