#Calculates pi1 across phenotypes b/w HCHS and MESA HIS
#Code from https://github.com/WheelerLab/DivPop/blob/master/matrix_eQTL/09_pi1_allpops.R
library(qvalue)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
"%&%" = function(a,b) paste(a,b,sep="")

#open.dir <- "Z:/" #for mount
open.dir <- "/home/angela/" #for direct

source(open.dir %&% "px_yri_chol/GWAS/qqman.r")
pheno_list <- c("CHOL", "HDL", "TRIG", "LDL")
pheno_list2 <- c("CHOL", "HDL", "TRIG", "LDL")
nphenos <- length(pheno_list)
nphenos2 <- length(pheno_list2)
pi1_matrix <- matrix(NA, nrow = nphenos, ncol = nphenos2)
rownames(pi1_matrix) <- c("CHOL", "HDL", "TRIG", "LDL")
colnames(pi1_matrix) <- c("CHOL", "HDL", "TRIG", "LDL")

#genes
for(i in 1:length(pheno_list)){ #discovery, HCHS MultiXcan (no related)
  refinfile <- fread(open.dir %&% "MultiXcan/HCHS_wo_relateds/" %&% pheno_list[i] %&% "_HCHS.txt")
  refinfile <- refinfile[complete.cases(refinfile$pvalue),]
  png(open.dir %&% "px_his_chol/pi1/" %&% pheno_list[i] %&% "_HCHS_hist.png")
  hist(refinfile$pvalue, main = pheno_list[i] %&% " HCHS",  xlim = c(0, 1), breaks = 20)
  dev.off()
  #png(open.dir %&% "px_his_chol/pi1/" %&% pheno_list[i] %&% "_HCHS_qq.png")
  #qq(refinfile$pvalue, plot = T, main = pheno_list[i] %&% " HCHS")
  #dev.off()
  for(j in 1:length(pheno_list2)){ #replication, MESA HIS MultiXcan
    testinfile <- fread(open.dir %&% "MultiXcan/" %&% pheno_list2[j] %&% "_MESA_HIS.txt")
    png(open.dir %&% "px_his_chol/pi1/" %&% pheno_list2[j] %&% "_MESA_hist.png")
    hist(testinfile$pvalue, main = pheno_list2[j] %&% " MESA",  xlim = c(0, 1), breaks = 20)
    dev.off()
    #png(open.dir %&% "px_his_chol/pi1/" %&% pheno_list2[j] %&% "_MESA_qq.png")
    #qq(testinfile$pvalue, plot = T, main = pheno_list2[j] %&% " MESA")
    #dev.off()
    refinfile$FDR <- p.adjust(refinfile$pvalue)
    testinfile$FDR <- p.adjust(testinfile$pvalue)
    pheno1 <- refinfile
    pheno2 <- testinfile
    pheno1fdr05 <- dplyr::filter(pheno1, FDR < 0.05)
    fdr <- dim(pheno1fdr05)
    pheno2tested <- inner_join(pheno1fdr05, pheno2, by = "gene")
    over <- dim(pheno2tested)
    pheno2pvals <- pheno2tested$pvalue.y
    pheno2pvals <- pheno2pvals[complete.cases(pheno2pvals)]
    png(open.dir %&% "px_his_chol/pi1/" %&% pheno_list[i] %&% "_HCHS_" %&% pheno_list2[j] %&% "_MESA.png")
    hist(pheno2pvals, main = pheno_list[i] %&% " HCHS, " %&% pheno_list2[j] %&% " MESA", xlim = c(0, 1), breaks = 20)
    dev.off()
    tryCatch({ #if qval calculation doesn't work, skip to next one
      qobjpheno2 <- qvalue(p = pheno2pvals)
      hist(qobjpheno2)
      pi1 <- 1 - qobjpheno2$pi0
      pi12<-signif(pi1,4)
      pi1_matrix[i,j] <- pi12
      print(pi1_matrix)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

write.table(pi1_matrix, open.dir %&% "px_his_chol/pi1/pi1values_MultiXcan_HCHS_wo_r_MESA_HIS.txt", quote = F, sep = "\t")
melted_pi1 <- melt(pi1_matrix, na.rm = TRUE)
melted_pi1 <- melted_pi1[nrow(melted_pi1):1,]
ggplot(data = melted_pi1, aes(x=Var2, y=Var1, fill=value)) + geom_tile(color = "white") + scale_fill_gradient2(low = "white", high = "black", mid = "dark green", midpoint = 0.5, limit = c(0,1), name = "Pi1") + ggtitle("Pi1 for MultiXcan HCHS vs. MESA HIS") + ylab("HCHS, w/o relateds\n(discovery)") + xlab("MESA_HIS\n(replication)") + scale_y_discrete(limits = rev(levels(melted_pi1$Var2)))
