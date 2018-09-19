#MAKE THIS MORE READABLE
#converts local ancestry output from MOSAIC to a more human-readable format

library(data.table)
library(dplyr)
library(MOSAIC)
library(tibble)
library(reshape2)
args <- commandArgs(trailingOnly = T)
input_file_name <- args[1] #localanc.RData produced by MOSAIC
output_file_name <- args[2] #output prefix

load(input_file_name)
df <- reshape2::melt(localanc[[1]])
  #Var1 is each of the 3 possible ancestries
  #Var2 is haplotype
  #Var3 is location (bp)
g_loc <- tibble::rownames_to_column(as.data.frame(g.loc), "Var3")
g_loc$Var3 <- as.integer(g_loc$Var3)
loc_anc <- dplyr::left_join(df, g_loc)
loc_anc$Var3 <- NULL
colnames(loc_anc) <- c("ancestry", "haplotype_num", "prob", "bp")
phind <- fread("/home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/HCHS_chr22.phind", header = F)
phind <- phind %>% dplyr::select(V1)
colnames(phind) <- "haplotype"
phind <- tibble::rownames_to_column(phind, "haplotype_num")
phind$haplotype_num <- as.integer(phind$haplotype_num)
loc_anc <- left_join(loc_anc, phind)
loc_anc$haplotype_num <- NULL
snpfile <- fread("/home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/snpfile.22")
snpfile <- snpfile %>% dplyr::select(V3, V4)
colnames(snpfile) <- c("cM", "bp")
loc_anc <- left_join(loc_anc, snpfile)
pdf(paste(output_file_name,".pdf"))
plot_localanc(chrnos, g.loc, localanc)
dev.off()
fwrite(loc_anc, paste(output_file_name, ".csv"), sep = ",", row.names = F, col.names = T, quote = F, na = "NA")
#Okay but how do I tell which ancestry is which
