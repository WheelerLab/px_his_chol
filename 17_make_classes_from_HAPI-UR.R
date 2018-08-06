#makes .classes file for RFMix from HAPI-UR
library(dplyr)
library(data.table)
args <- commandArgs(trailingOnly=T)
phind_file_name <- args[1] #name of .phind file that HAPI-UR outputted
pop_file_name <- args[2] #name of pop file with three columns: FID, IID, and pop
test_pop <- args[3] #code of population to be tested

phind_file_name <- "ancestry_pipeline/HCHS_chr22/phase_chr22.phind"
pop_file_name <- "ancestry_pipeline/HCHS_chr22/ordered_pops.txt"
test_pop <- "HCHS"

#parse .phind file
phind <- fread(phind_file_name, header = F)
phind$V2 <- NULL
phind$V3 <- NULL
phind$FID <- gsub("\\:.*", "", phind$V1)
phind$IID <- gsub(".*:\\s*|_.*", "", phind$V1)
phind$V1 <- NULL

#add pops
pop_file <- fread(pop_file_name)
