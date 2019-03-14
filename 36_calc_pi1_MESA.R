library(data.table)
library(dplyr)
library(qvalue)
"%&%" = function(a,b) paste(a,b,sep="")
setwd("/home/angela/px_his_chol/MESA_compare/")

disc_thres <- c("1", "0.05", "0.01", "0.005", "0.001", "0.0005", "0.0001") #discovery threshold
rep_pop <- c("MESA_HIS", "MESA_CAU")
all_PX <- fread("three_pops_filtered.csv")

for(h in c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")){
  pi1_P_matrix <- matrix(NA, nrow = length(disc_thres), ncol = length(rep_pop)) #intiialize to hold pi1 vals
  rownames(pi1_P_matrix) <- disc_thres
  colnames(pi1_P_matrix) <- c("MESA_HIS", "MESA_CAU")
  P_cols <- c("MESA_HIS_P", "MESA_CAU_P")
  for (i in 1:length(disc_thres)){
    thres <- as.numeric(disc_thres[i])
    for (j in 1:length(P_cols)){ #triple loops may or may not be ratchet
      pops <- all_PX %>% dplyr::select(c("HCHS_P", "gene", "model", "pheno", P_cols[j]))
      pops <- pops[!duplicated(pops), ]
      pops_filtered <- subset(pops, HCHS_P <= thres & pheno == h)
      pops_filtered <- pops_filtered %>% dplyr::select(P_cols[j])
      colnames(pops_filtered) <- "pop"
      pops_filtered <- na.omit(pops_filtered)
      if(nrow(pops_filtered) > 0){
        tryCatch({
          qobj <- qvalue(pops_filtered$pop)
          png("hists/" %&% P_cols[j] %&% "_" %&% disc_thres[i] %&% "_" %&% h %&% ".png")
          hist(pops_filtered$pop, main = P_cols[j] %&% ", HCHS P < " %&% disc_thres[i] %&% ", " %&% h)
          dev.off()
          pi1 <- signif(1 - qobj$pi0, 4)
          pi1_P_matrix[i, j] <- pi1
        }, error=function(e){})
      }
    }
  }
  write.table(pi1_P_matrix, h %&% "_pi1_matrix.txt", quote = F, sep = "\t")
}

#all phenos
pi1_P_matrix <- matrix(NA, nrow = length(disc_thres), ncol = length(rep_pop)) #intiialize to hold pi1 vals
rownames(pi1_P_matrix) <- disc_thres
colnames(pi1_P_matrix) <- c("MESA_HIS", "MESA_CAU")
P_cols <- c("MESA_HIS_P", "MESA_CAU_P")
for (i in 1:length(disc_thres)){
  thres <- as.numeric(disc_thres[i])
  for (j in 1:length(P_cols)){ #triple loops may or may not be ratchet
    pops <- all_PX %>% dplyr::select(c("HCHS_P", "gene", "model", "pheno", P_cols[j]))
    pops <- pops[!duplicated(pops), ]
    pops_filtered <- subset(pops, HCHS_P <= thres)
    pops_filtered <- pops_filtered %>% dplyr::select(P_cols[j])
    colnames(pops_filtered) <- "pop"
    pops_filtered <- na.omit(pops_filtered)
    if(nrow(pops_filtered) > 0){
      tryCatch({
        qobj <- qvalue(pops_filtered$pop)
        png("hists/" %&% P_cols[j] %&% "_" %&% disc_thres[i] %&% ".png")
        hist(pops_filtered$pop, main = P_cols[j] %&% ", HCHS P < " %&% disc_thres[i])
        dev.off()
        pi1 <- signif(1 - qobj$pi0, 4)
        pi1_P_matrix[i, j] <- pi1
      }, error=function(e){})
    }
  }
}
write.table(pi1_P_matrix, "pi1_matrix.txt", quote = F, sep = "\t")
