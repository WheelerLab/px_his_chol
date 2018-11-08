library(data.table)
library(dplyr)
library(R.utils)

"%&%" = function(a,b) paste(a,b,sep="")
phenos <- c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")
chrs <- c(1:22)
pops <- fread("/home/angela/px_his_chol/MESA_compare/GTEx_WB/SNPs_n_samples/sig_gene_in_tiss/Fst_num_genes_n_samples.csv")
pops_sample_size <- pops$Sample_size
pops <- pops$Tiss
pops <- gsub("TW_", "", pops) #GTEx eQTL files don't have prefix or suffix
pops <- gsub("_0.5.db", "", pops)
sig_gene_SNPs <- fread("/home/angela/px_his_chol/COLOC/GTEx_SNPs_in_sig_genes.txt", header = F) #so we don't run all the SNPs b/c it takes forever
sig_gene_SNPs <- sig_gene_SNPs$V1

for(pop in 1:length(pops)){ #read in pop's .frq file for MAF
  frq <- fread("/home/angela/px_his_chol/COLOC/GTEx_WHLBLD.frq") #we're just gonna use WHLBLD dosages for everyone
  bim <- fread("/home/angela/px_his_chol/MESA_compare/GTEx_WB/GTEx_WHLBLD.bim") #for adding cpos
  bim$cpos <- paste(bim$V1, bim$V4, sep = "_")
  bim <- bim %>% select(V2, cpos)
  colnames(bim) <- c("SNP", "cpos")
  frq <- left_join(frq, bim, by = "SNP")
  frq <- frq %>% select(SNP, MAF, cpos)
  
  #DOWNLOADED V6 FROM GTEx WEBSITE
  meqtl <- fread("/home/angela/px_his_chol/COLOC/GTEx_V6_eqtl/GTEx_Analysis_v6p_all-associations/" %&% pops[pop] %&% "_Analysis.v6p.all_snpgene_pairs.txt.gz", nThread = 30) #read in matrix eQTL results
  
  #for testing
  #meqtl <- meqtl[1:1000000,]
  
  #meqtl$se <- meqtl$beta / meqtl$statistic #make your own standard error since it's not in the meQTL output
  meqtl$n_samples <- pops_sample_size[pop]
  
  #why can't you just be in normal cpos format
  meqtl$cpos <- gsub("^([^_]*_[^_]*)_.*$", "\\1", meqtl$variant_id) #https://stackoverflow.com/questions/7449564/regex-return-all-before-the-second-occurrence
  
  meQTL_for_COLOC <- left_join(meqtl, frq, by = "cpos") #add freq to COLOC input
  meQTL_for_COLOC <- meQTL_for_COLOC %>% select(gene_id, SNP, MAF, pval_nominal, slope, slope_se) #subset to COLOC input
  colnames(meQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
  meQTL_for_COLOC <- meQTL_for_COLOC[complete.cases(meQTL_for_COLOC),]
  
  for(pheno in phenos){ #read in GEMMA output file
    GEMMA_result <- fread("/home/angela/px_his_chol/GEMMA/output/" %&% pheno %&% "_sig_snps_p-0.txt", header = T)
    GEMMA_result$chr <- gsub("chr", "", GEMMA_result$chr)
    GEMMA_result$cpos <- paste(GEMMA_result$chr, GEMMA_result$ps, sep = "_")
    GEMMA_for_COLOC <- GEMMA_result %>% select(rs, beta, se, af, cpos) #subset to COLOC input
    GEMMA_for_COLOC$sample_size <- 11103
    colnames(GEMMA_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size", "variant_id")
    GEMMA_for_COLOC <- GEMMA_for_COLOC[complete.cases(GEMMA_for_COLOC),] #COLOC does not like missing values
    
    snps_in_both <- intersect(GEMMA_for_COLOC$panel_variant_id, meQTL_for_COLOC$variant_id) #is there a better way to do this? Probably. Do I feel like figuring it out? Nah.
    snps_in_all <- intersect(snps_in_both, sig_gene_SNPs)
    GEMMA_for_COLOC_chr <- subset(GEMMA_for_COLOC, panel_variant_id %in% snps_in_all)
    meQTL_for_COLOC_chr <- subset(meQTL_for_COLOC, variant_id %in% snps_in_all)
    
    fwrite(meQTL_for_COLOC_chr, "/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_TW_" %&% pops[pop] %&% "_0.5.db_" %&% pheno %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
    gzip("/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_TW_" %&% pops[pop] %&% "_0.5.db_" %&% pheno %&% ".txt", destname = "/home/angela/px_his_chol/COLOC/COLOC_input/eQTL_TW_" %&% pops[pop] %&% "_0.5.db_" %&% pheno %&% ".txt.gz") #script may only take .gz values so can't hurt to be too careful
    fwrite(GEMMA_for_COLOC_chr, "/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_TW_" %&% pops[pop] %&% "_0.5.db_" %&% pheno %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
    gzip("/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_TW_" %&% pops[pop] %&% "_0.5.db_" %&% pheno %&% ".txt", "/home/angela/px_his_chol/COLOC/COLOC_input/GWAS_TW_" %&% pops[pop] %&% "_0.5.db_" %&% pheno %&% ".txt.gz")
    print("Completed with " %&% pops[pop] %&% ", for " %&% pheno %&% ".")
  }
}
