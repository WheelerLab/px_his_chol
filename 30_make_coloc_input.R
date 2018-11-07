library(data.table)
library(dplyr)
library(R.utils)

"%&%" = function(a,b) paste(a,b,sep="")
phenos <- c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")
chrs <- c(1:22)
pops <- c("AFA", "CAU", "HIS", "AFHI", "ALL") #do combined pops later
pops_sample_size <- c(233, 352, 578, 585, 1163) #R doesn't have dicts so we're doing it a slgihtly more ratchet way

for(pop in 1:length(pops)){ #read in pop's .frq file for MAF
  frq <- fread("/home/angela/px_his_chol/MESA_compare/" %&% pops[pop] %&% ".frq")
  frq <- frq %>% select(SNP, MAF)

  for(pheno in phenos){ #read in GEMMA output file
    GEMMA_result <- fread("/home/angela/px_his_chol/GEMMA/output/" %&% pheno %&% "_sig_snps_p-0.txt", header = T)
    GEMMA_result$chr_pos <- paste(gsub("chr", "", GEMMA_result$chr), GEMMA_result$ps, sep = ":")
    GEMMA_for_COLOC <- GEMMA_result %>% select(rs, beta, se, af) #subset to COLOC input
    GEMMA_for_COLOC$sample_size <- 11103
    colnames(GEMMA_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")
    GEMMA_for_COLOC <- GEMMA_for_COLOC[complete.cases(GEMMA_for_COLOC),] #COLOC does not like missing values
    
    for(chr in chrs){ #yes triple loops are ratchet
      system("zcat -f /home/lauren/files_for_revisions_plosgen/meqtl_results/MESA/" %&% pops[pop] %&% "_Nk_10_PFs_chr" %&% chr %&% "pcs_3.meqtl.cis.* > /home/angela/px_his_chol/COLOC/COLOC_input/meQTL_input.txt") #fread doesn't seem to like wildcards so we're gonna do this the ugly way
      meqtl <- fread("/home/angela/px_his_chol/COLOC/COLOC_input/meQTL_input.txt") #read in matrix eQTL results
      meqtl$se <- meqtl$beta / meqtl$statistic #make your own standard error since it's not in the meQTL output
      meqtl$n_samples <- pops_sample_size[pop]
      meQTL_for_COLOC <- left_join(meqtl, frq, by = c("snps" = "SNP")) #add freq to COLOC input
      meQTL_for_COLOC <- meQTL_for_COLOC %>% select(gene, snps, MAF, pvalue, beta, se) #subset to COLOC input
      colnames(meQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
      meQTL_for_COLOC <- meQTL_for_COLOC[complete.cases(meQTL_for_COLOC),]

      snps_in_both <- intersect(GEMMA_for_COLOC$panel_variant_id, meQTL_for_COLOC$variant_id) #is there a better way to do this? Probably. Do I feel like figuring it out? Nah.
      GEMMA_for_COLOC_chr <- subset(GEMMA_for_COLOC, panel_variant_id %in% snps_in_both)
      meQTL_for_COLOC_chr <- subset(meQTL_for_COLOC, variant_id %in% snps_in_both)

      fwrite(meQTL_for_COLOC_chr, "/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_" %&% pops[pop] %&% "_chr" %&% chr %&% "_" %&% pheno %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
      gzip("/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_" %&% pops[pop] %&% "_chr" %&% chr %&% "_" %&% pheno %&% ".txt", destname = "/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_" %&% pops[pop] %&% "_chr" %&% chr %&% "_" %&% pheno %&% ".txt.gz") #script may only take .gz values so can't hurt to be too careful
      fwrite(GEMMA_for_COLOC_chr, "/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_" %&% pops[pop] %&% "_chr" %&% chr %&% "_" %&% pheno %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
      gzip("/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_" %&% pops[pop] %&% "_chr" %&% chr %&% "_" %&% pheno %&% ".txt", "/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_" %&% pops[pop] %&% "_chr" %&% chr %&% "_" %&% pheno %&% ".txt.gz")
      print("Completed with " %&% pops[pop] %&% ", chr" %&% chr %&% ", for " %&% pheno %&% ".")
    }
  }
}
