#backward elimination of significant genes to determine which ones are independent
library(data.table)
library(MASS)
library(dplyr)

setwd("/home/angela/px_his_chol/COLOC/backward_elimination/")
"%&%" = function(a,b) paste(a,b,sep="")
sig_gene_HCHS <- fread("/home/angela/px_his_chol/MESA_compare/GTEx_WB/sig_gene_HCHS_no_minus.csv")
pheno <- fread('/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt', header = T)
pheno <- pheno %>% dplyr::select(FID, IID, CHOL_rank, HDL_rank, TRIG_rank, LDL_rank)
phenos <- c("CHOL", "HDL", "TRIG", "LDL")
tissues <- fread("/home/angela/px_yri_chol/PrediXcan/database_tissues.txt", header = F)
tissues <- tissues$V1 #get list of tissues
tissues <- c("AFA.db", "AFHI.db", "ALL.db", "CAU.db", "HIS.db", tissues) #add MESA tissues
back_elim <- pheno[, 2] #extract only IID
sig_gene_HCHS$gene_tiss_combos <- paste(sig_gene_HCHS$tissue, sig_gene_HCHS$gene, sep = "_")
gene_tiss_combos <- sig_gene_HCHS$gene_tiss_combos
gene_tiss_combos <- c("IID", gene_tiss_combos)
sig_gene_HCHS$gene_tiss_combos <- NULL

#get table of all predicted expressions
#for(tiss_name in tissues){
#  print("Running analyses on " %&% tiss_name %&% ".")
#  tiss <- fread('/home/angela/px_his_chol/PrediXcan/' %&% tiss_name %&% '_predicted_expression.txt', header = F, nThread = 30) #load in predicted expression file
#  tiss$V1 <- NULL #remove FID
#  tiss_header <- paste(tiss_name, "_", gsub("\\..*", "", c(tiss[1,])), sep = "") #remove all after . and add tiss name
#  tiss <- tiss[2:nrow(tiss),] #remove header line
#  colnames(tiss) <- tiss_header #assign header line
#  colnames(tiss)[1] <- "IID" #rename beginning
#  tiss_gene_in_tiss <- c("IID", colnames(tiss)) #add IID to list
#  tiss_gene_in_tiss <- intersect(tiss_gene_in_tiss, gene_tiss_combos)
#  for_back_elim <- tiss %>% select(tiss_gene_in_tiss) #select IIDs and intersection of tiss-gene combos
#  back_elim <- left_join(back_elim, for_back_elim, by = "IID") #add to list to back-elim
#}

#fwrite(back_elim, "/home/angela/px_his_chol/COLOC/backward_elimination/all_predicted_expression.csv", row.names = F, col.names = T, sep = ",", na = "NA", quote = F)
#save.image(file = "back_elim.RData") #this likes to crash b/c it's too big
back_elim <- fread("/home/angela/px_his_chol/COLOC/backward_elimination/all_predicted_expression.csv")
back_elim <- left_join(back_elim, pheno, by = "IID")

#make clusters of genes
gene_clusters <- sig_gene_HCHS %>% dplyr::select(pheno, chr, gene, tissue)
gene_clusters <- gene_clusters[order(pheno, chr),] 
gene_clusters$first_occur <- paste(gene_clusters$pheno, gene_clusters$chr, gene_clusters$gene, sep = "_") #instead of checking across three columns for unique, just check one
gene_clusters <- gene_clusters[match(unique(gene_clusters$gene), gene_clusters$gene),]
gene_clusters$first_occur <- NULL
gene_clusters <- aggregate(data = gene_clusters, gene ~ pheno + chr, FUN=paste) #https://stackoverflow.com/questions/38125125/combine-rows-which-have-same-value-in-two-columns-r
gene_clusters$pheno_chr <- paste(gene_clusters$pheno, gene_clusters$chr, sep = "_") #instead of checking across three columns for unique, just check one
gene_clusters$pheno <- NULL
gene_clusters$chr <- NULL

for(pheno_chr in gene_clusters$pheno_chr){
  #pheno_chr <- (gene_clusters$pheno_chr)[1]
  pheno_name <- c(strsplit(pheno_chr, "_"))
  chr <- pheno_name[[1]][2]
  pheno_name <- pheno_name[[1]][1] #why is this so complicated
  pheno_name_rank <- paste(pheno_name, "rank", sep = "_")
  
  genes_to_test <- subset(gene_clusters, pheno_chr == paste(pheno_name, chr, sep = "_"))
  genes_to_test <- genes_to_test[1, 1][[1]] #why is this data structure so weird
  tissues_to_test <- subset(sig_gene_HCHS, gene %in% genes_to_test & chr == as.numeric(chr) & pheno == pheno_name)
  tissue_gene <- paste(tissues_to_test$tissue, tissues_to_test$gene, sep = "_")
  tissue_gene <- c("IID", tissue_gene)
  
  #there was probably a MUCH less convoluted way to do this but here we are
  
  print("Started making models for " %&% pheno_name %&% ", chr" %&% chr %&% ".")
  tiss_gene_to_keep <- c(pheno_name_rank, intersect(colnames(back_elim), tissue_gene))
  back_elim_cluster <- back_elim %>% dplyr::select(tiss_gene_to_keep)
  back_elim_cluster <- back_elim_cluster[complete.cases(back_elim_cluster),]
  back_elim_cluster$IID <- NULL
  predictor_genes <- colnames(back_elim_cluster)[2:length(colnames(back_elim_cluster))] #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
  fmla <- paste(pheno_name_rank,  " ~ ", paste(predictor_genes, collapse= "+"))
  fmla <- gsub(" - ", "_", fmla)
  #fmla <- as.formula(gsub(" - ", "_", paste(pheno_name_rank,  " ~ ", paste(predictor_genes, collapse= "+"))))
  
  fwrite(back_elim_cluster, pheno_chr %&% "_all_tiss_gene_before_back_elim.csv", row.names = F, col.names = T, sep = ",", na = NA, quote = F)
  all_tiss_gene <- lm(fmla, data = back_elim_cluster) 
  saveRDS(all_tiss_gene, file = pheno_chr %&% "_all_tiss_gene.rds")
  print("Finished making full model for " %&% pheno_name_rank %&% ".")
  
  #using stepAIC, is there a better option somewhere?
  back_elim_complete <- stepAIC(all_tiss_gene, direction = "backward", trace = FALSE) #perform backward analysis on full model
  saveRDS(back_elim_complete, file = pheno_chr %&% "_back_elim.rds")
  print("Finished making backward-eliminated model for " %&% pheno_name_rank %&% ".")
}


