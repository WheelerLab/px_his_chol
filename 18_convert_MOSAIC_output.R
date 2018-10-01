#converts local ancestry output from MOSAIC to a more human-readable format
#example input - Rscript 18_convert_MOSAIC_output.R --input_file_name localanc_HCHS_3way_1-12135_1-1_24754_60_0.99_100.RData --phind_file_name /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/HCHS_chr1.phind --snpfile_file_name /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/snpfile.1 --output_file_name local_anc_1_HCHS
#surprisingly this doesn't take horribly long to produce a 31m line file at the end

library(argparse)
library(data.table)
library(dplyr)
library(MOSAIC)
library(tibble)
library(reshape2)

parser <- ArgumentParser()
parser$add_argument("--input_file_name",help="localanc.RData produced by MOSAIC")
parser$add_argument("--phind_file_name",help="MOSAIC input admixed.phind")
parser$add_argument("--snpfile_file_name", help="MOSAIC input snpfile")
parser$add_argument("--output_file_name", help="Output prefix. Default = 'local_anc'.", default = "local_anc")
parser$add_argument("--probability_threshold", help = "Threshold to cut off ancestry probability. Default = 0.8.", default = 0.8)
args <- parser$parse_args()

input_file_name <- args$input_file_name
phind_file_name <- args$phind_file_name
snpfile_file_name <- args$snpfile_file_name
output_file_name <- args$output_file_name
probability_threshold <- as.integer(args$probability_threshold)

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
loc_anc <- subset(loc_anc, prob > probability_threshold) #keep ancestry probabilities > 0.8
phind <- fread(phind_file_name, header = F)
phind <- phind %>% dplyr::select(V1)
colnames(phind) <- "haplotype"
phind <- tibble::rownames_to_column(phind, "haplotype_num")
phind$haplotype_num <- as.integer(phind$haplotype_num)
loc_anc <- left_join(loc_anc, phind)
loc_anc$haplotype_num <- NULL
snpfile <- fread(snpfile_file_name)
snpfile <- snpfile %>% dplyr::select(V3, V4)
colnames(snpfile) <- c("cM", "bp")
loc_anc <- left_join(loc_anc, snpfile)
anc_codes <- fread("/home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/ancestry_codes.txt", header = T)
loc_anc <- left_join(loc_anc, anc_codes, by = "ancestry")
loc_anc$ancestry <- NULL
loc_anc$cM <- NULL
#pdf(paste(output_file_name, ".pdf", sep = ""))
#plot_localanc(chrnos, g.loc, localanc) #error here?
#dev.off()
output_file_name <- paste(output_file_name, ".csv", sep = "")
fwrite(loc_anc, output_file_name, sep = ",", row.names = F, col.names = T, quote = F, na = "NA")
