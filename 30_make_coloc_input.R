library(data.table)
library(dplyr)
library(R.utils)

"%&%" = function(a,b) paste(a,b,sep="")
phenos <- c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")

frq <- fread("/home/angela/px_his_chol/MESA_compare/HIS_rm_missnps.frq")
frq <- frq %>% select(SNP, MAF)
meqtl <- fread("/home/lauren/files_for_revisions_plosgen/meqtl_results/MESA/HIS_Nk_10_PFs_chr20pcs_3.meqtl.cis.2018-04-04.txt.gz")
meqtl$se <- meqtl$beta / meqtl$statistic
meqtl$n_samples <- 352
meQTL_for_COLOC <- left_join(meqtl, frq, by = c("snps" = "SNP"))
meQTL_for_COLOC <- meQTL_for_COLOC %>% select(gene, snps, MAF, pvalue, beta, se)
colnames(meQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
meQTL_for_COLOC <- meQTL_for_COLOC[complete.cases(meQTL_for_COLOC),]

GEMMA_result <- fread("/home/angela/px_his_chol/GEMMA/output/CHOL_rank_sig_snps_p-0.txt", header = T)
GEMMA_result$chr_pos <- paste(gsub("chr", "", GEMMA_result$chr), GEMMA_result$ps, sep = ":")
GEMMA_for_COLOC <- GEMMA_result %>% select(rs, beta, se, af)
GEMMA_for_COLOC$sample_size <- 11103
colnames(GEMMA_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")

snps_in_both <- intersect(GEMMA_for_COLOC$panel_variant_id, meQTL_for_COLOC$variant_id)
GEMMA_for_COLOC <- subset(GEMMA_for_COLOC, panel_variant_id %in% snps_in_both)
meQTL_for_COLOC <- subset(meQTL_for_COLOC, variant_id %in% snps_in_both)

fwrite(meQTL_for_COLOC, "/home/angela/px_his_chol/COLOC/HIS_for_COLOC_chr20.txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
gzip("/home/angela/px_his_chol/COLOC/HIS_for_COLOC_chr20.txt", destname = "/home/angela/px_his_chol/COLOC/HIS_for_COLOC_chr10.txt.gz")
fwrite(GEMMA_for_COLOC, "/home/angela/px_his_chol/COLOC/CHOL_rank_chr20_for_COLOC.txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
gzip("/home/angela/px_his_chol/COLOC/CHOL_rank_chr20_for_COLOC.txt", destname = "/home/angela/px_his_chol/COLOC/CHOL_rank_chr20_for_COLOC.txt.gz")

