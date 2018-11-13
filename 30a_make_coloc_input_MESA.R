library(data.table)
library(dplyr)
library(R.utils)

"%&%" = function(a,b) paste(a,b,sep="")
phenos <- c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")
chrs <- c(1:22)
pops <- c("AFA", "CAU", "HIS", "AFHI", "ALL") #do combined pops later
pops_sample_size <- c(233, 578, 352, 585, 1163) #R doesn't have dicts so we're doing it a slgihtly more ratchet way
sig_gene_SNPs <- fread("/home/angela/px_his_chol/COLOC/COLOC_input/sig_gene_SNPs_84.txt", header = F) #so we don't run all the SNPs b/c it takes forever
sig_gene_SNPs <- sig_gene_SNPs$V1

for(pop in 1:length(pops)){ #read in pop's .frq file for MAF
  frq <- fread("/home/angela/px_his_chol/MESA_compare/" %&% pops[pop] %&% ".frq")
  frq <- frq %>% dplyr::select(SNP, MAF)

  for(pheno in phenos){ #read in GEMMA output file
    GEMMA_result <- fread("/home/angela/px_his_chol/GEMMA/output/" %&% pheno %&% "_sig_snps_p-0.txt", header = T)
    GEMMA_result$chr_pos <- paste(gsub("chr", "", GEMMA_result$chr), GEMMA_result$ps, sep = ":")
    GEMMA_for_COLOC <- GEMMA_result %>% dplyr::select(rs, beta, se, af) #subset to COLOC input
    GEMMA_for_COLOC$sample_size <- 11103
    colnames(GEMMA_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")
    GEMMA_for_COLOC <- GEMMA_for_COLOC[complete.cases(GEMMA_for_COLOC),] #COLOC does not like missing values
    #GWAS_write <- data.frame(panel_variant_id = character(), effect_size = numeric(), standard_error = numeric(), frequency = numeric(), sample_size = numeric(), stringsAsFactors = F) 
    GWAS_write <- GEMMA_for_COLOC
    eQTL_write <- data.frame(gene_id = character(), variant_id = character(), maf = numeric(), pval_nominal = numeric(), slope = numeric(), slope_se = numeric(), stringsAsFactors = F) 
    
    for(chr in chrs){ #yes triple loops are ratchet
      system("zcat -f /home/lauren/files_for_revisions_plosgen/meqtl_results/MESA/" %&% pops[pop] %&% "_Nk_10_PFs_chr" %&% chr %&% "pcs_3.meqtl.cis.* > /home/angela/px_his_chol/COLOC/COLOC_input/meQTL_input.txt") #fread doesn't seem to like wildcards so we're gonna do this the ugly way
      meqtl <- fread("/home/angela/px_his_chol/COLOC/COLOC_input/meQTL_input.txt", nThread = 40) #read in matrix eQTL results
      meqtl$se <- meqtl$beta / meqtl$statistic #make your own standard error since it's not in the meQTL output
      meqtl$n_samples <- pops_sample_size[pop]
      meQTL_for_COLOC <- left_join(meqtl, frq, by = c("snps" = "SNP")) #add freq to COLOC input
      meQTL_for_COLOC <- meQTL_for_COLOC %>% dplyr::select(gene, snps, MAF, pvalue, beta, se) #subset to COLOC input
      colnames(meQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
      meQTL_for_COLOC <- meQTL_for_COLOC[complete.cases(meQTL_for_COLOC),]

      #GWAS_write <- rbind(GWAS_write, GEMMA_for_COLOC_chr)
      eQTL_write <- rbind(eQTL_write, meQTL_for_COLOC)
    }
    
    snps_in_both <- intersect(GWAS_write$panel_variant_id, eQTL_write$variant_id) #is there a better way to do this? Probably. Do I feel like figuring it out? Nah.
    snps_in_all <- intersect(snps_in_both, sig_gene_SNPs)
    GWAS_write <- subset(GWAS_write, panel_variant_id %in% snps_in_all)
    eQTL_write <- subset(eQTL_write, variant_id %in% snps_in_all)
    GWAS_write <- GWAS_write[order(GWAS_write$gene_id),]
    eQTL_write <- eQTL_write[order(eQTL_write$gene_id),]
    
    fwrite(eQTL_write, "/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
    gzip("/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", destname = "/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz") #script may only take .gz values so can't hurt to be too careful
    fwrite(GWAS_write, "/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
    gzip("/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", "/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz")
    print("Completed with " %&% pops[pop] %&% ", for " %&% pheno %&% ".")
  }
}

