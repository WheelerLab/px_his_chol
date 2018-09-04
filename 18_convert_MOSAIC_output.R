#converts local ancestry output from MOSAIC to a more human-readable format
load("/home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/MOSAIC_RESULTS/localanc_HCHS_3way_1-100_22-22_684_60_0.99_100.RData")
library(data.table)
library(dplyr)
library(MOSAIC)
library(tibble)
library(reshape2)
df <- reshape2::melt(localanc[[1]])
  #Var1 is each of the 3 possible ancestries
  #Var2 is haplotype
  #Var3 is location
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
#pdf("/home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/MOSAIC_PLOTS/plot_localanc.pdf")
#plot_localanc(chrnos, g.loc, localanc)
#dev.off()
fwrite(loc_anc, "/home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/loc_anc.csv", sep = ",", row.names = F, col.names = T, quote = F, na = "NA")
#Okay but how do I tell which ancestry is which
