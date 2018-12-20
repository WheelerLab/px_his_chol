source("/home/angela/px_yri_chol/GWAS/qqman.r")
library(data.table)
library(dplyr)
BP_Chrome <- fread("/home/angela/px_yri_chol/PrediXcan/BP_Chrome.txt")

#HCHS
HCHS_CHOL <- fread("/home/angela/px_his_chol/GEMMA/output/CHOL_rank.assoc.txt")
HCHS_CHOL$V12 <- as.numeric(HCHS_CHOL$V12)
HCHS_CHOL <- HCHS_CHOL[complete.cases(HCHS_CHOL),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_qq.png", width = 3, height = 3, units = "in", res = 300)
par(mfrow=c(2,2))
qq(HCHS_CHOL$V12, main = "HCHS, CHOL")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_CHOL_qq.pdf")
#qq(HCHS_CHOL$V12, main = "HCHS, CHOL")
#dev.off()

'
HCHS_CHOL <- HCHS_CHOL %>% dplyr::select(V2, V1, V3, V12)
colnames(HCHS_CHOL) <- c("SNP", "CHR", "BP", "P")
HCHS_CHOL$SNP <- as.factor(HCHS_CHOL$SNP)
HCHS_CHOL$CHR <- as.integer(gsub("chr", "", HCHS_CHOL$CHR))
HCHS_CHOL$BP <- as.integer(HCHS_CHOL$BP)
HCHS_CHOL$P <- as.double(HCHS_CHOL$P)
HCHS_CHOL <- HCHS_CHOL[complete.cases(HCHS_CHOL),]
HCHS_CHOL <- HCHS_CHOL[order(HCHS_CHOL$CHR, HCHS_CHOL$BP),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_CHOL_man.png")
manhattan(HCHS_CHOL, main = "HCHS, CHOL")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_CHOL_man.pdf")
manhattan(HCHS_CHOL, main = "HCHS, CHOL")
dev.off()
'

HCHS_HDL <- fread("/home/angela/px_his_chol/GEMMA/output/HDL_rank.assoc.txt")
HCHS_HDL$V12 <- as.numeric(HCHS_HDL$V12)
HCHS_HDL <- HCHS_HDL[complete.cases(HCHS_HDL),]
#png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_HDL_qq.png")
qq(HCHS_HDL$V12, main = "HCHS, HDL")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_HDL_qq.pdf")
#qq(HCHS_HDL$V12, main = "HCHS, HDL")
#dev.off()

HCHS_TRIG <- fread("/home/angela/px_his_chol/GEMMA/output/TRIG_rank.assoc.txt")
HCHS_TRIG$V12 <- as.numeric(HCHS_TRIG$V12)
HCHS_TRIG <- HCHS_TRIG[complete.cases(HCHS_TRIG),]
#png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_TRIG_qq.png")
qq(HCHS_TRIG$V12, main = "HCHS, TRIG")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_TRIG_qq.pdf")
#qq(HCHS_TRIG$V12, main = "HCHS, TRIG")
#dev.off()

HCHS_LDL <- fread("/home/angela/px_his_chol/GEMMA/output/LDL_rank.assoc.txt")
HCHS_LDL$V12 <- as.numeric(HCHS_LDL$V12)
HCHS_LDL <- HCHS_LDL[complete.cases(HCHS_LDL),]
#png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_LDL_qq.png")
qq(HCHS_LDL$V12, main = "HCHS, LDL")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/HCHS_LDL_qq.pdf")
#qq(HCHS_LDL$V12, main = "HCHS, LDL")
dev.off()

all_PX <- fread("/home/angela/px_his_chol/GEMMA/pred_exp/output/all_PX_results.csv")
HCHS_CHOL_PX <- subset(all_PX, pheno == "CHOL")
HCHS_CHOL_PX$start_bp <- NULL
HCHS_CHOL_PX <- HCHS_CHOL_PX[complete.cases(HCHS_CHOL_PX),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_qq_PX.png", width = 3, height = 3, units = "in", res = 300)
par(mfrow=c(2,2))
#qq(HCHS_CHOL_PX$p_wald, main = "HCHS, CHOL")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_CHOL_qq.pdf")
#qq(HCHS_CHOL_PX$p_wald, main = "HCHS, CHOL")
#dev.off()

HCHS_CHOL_PX <- subset(all_PX, pheno == "CHOL")
HCHS_CHOL_PX$start_bp <- NULL
HCHS_CHOL_PX <- HCHS_CHOL_PX[complete.cases(HCHS_CHOL_PX),]
#png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_CHOL_qq.png")
qq(HCHS_CHOL_PX$p_wald, main = "HCHS, CHOL")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_CHOL_qq.pdf")
#qq(HCHS_CHOL_PX$p_wald, main = "HCHS, CHOL")
#dev.off()

HCHS_HDL_PX <- subset(all_PX, pheno == "HDL")
HCHS_HDL_PX$start_bp <- NULL
HCHS_HDL_PX <- HCHS_HDL_PX[complete.cases(HCHS_HDL_PX),]
#png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_HDL_qq.png")
qq(HCHS_HDL_PX$p_wald, main = "HCHS, HDL")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_HDL_qq.pdf")
#qq(HCHS_HDL_PX$p_wald, main = "HCHS, HDL")
#dev.off()

HCHS_TRIG_PX <- subset(all_PX, pheno == "TRIG")
HCHS_TRIG_PX$start_bp <- NULL
HCHS_TRIG_PX <- HCHS_TRIG_PX[complete.cases(HCHS_TRIG_PX),]
#png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_TRIG_qq.png")
qq(HCHS_TRIG_PX$p_wald, main = "HCHS, TRIG")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_TRIG_qq.pdf")
#qq(HCHS_CHOL_PX$p_wald, main = "HCHS, TRIG")
#dev.off()

HCHS_LDL_PX <- subset(all_PX, pheno == "LDL")
HCHS_LDL_PX$start_bp <- NULL
HCHS_LDL_PX <- HCHS_LDL_PX[complete.cases(HCHS_LDL_PX),]
#png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_LDL_qq.png")
qq(HCHS_LDL_PX$p_wald, main = "HCHS, LDL")
#dev.off()
#pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/HCHS_LDL_qq.pdf")
#qq(HCHS_CHOL_PX$p_wald, main = "HCHS, LDL")
dev.off()

'
#MESA
MESA_CHOL <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/output/CHOL_rank.assoc.txt")
MESA_CHOL$p_wald <- as.numeric(MESA_CHOL$p_wald)
MESA_CHOL <- MESA_CHOL[complete.cases(MESA_CHOL),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_CHOL_qq.png")
qq(MESA_CHOL$p_wald, main = "MESA HIS, CHOL")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_CHOL_qq.pdf")
qq(MESA_CHOL$p_wald, main = "MESA HIS, CHOL")
dev.off()

MESA_HDL <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/output/HDL_rank.assoc.txt")
MESA_HDL$p_wald <- as.numeric(MESA_HDL$p_wald)
MESA_HDL <- MESA_HDL[complete.cases(MESA_HDL),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_HDL_qq.png")
qq(MESA_HDL$p_wald, main = "MESA HIS, HDL")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_HDL_qq.pdf")
qq(MESA_HDL$p_wald, main = "MESA HIS, HDL")
dev.off()

MESA_TRIG <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/output/TRIG_rank.assoc.txt")
MESA_TRIG$p_wald <- as.numeric(MESA_TRIG$p_wald)
MESA_TRIG <- MESA_TRIG[complete.cases(MESA_TRIG),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_TRIG_qq.png")
qq(MESA_TRIG$p_wald, main = "MESA HIS, TRIG")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_TRIG_qq.pdf")
qq(MESA_TRIG$p_wald, main = "MESA HIS, TRIG")
dev.off()

MESA_LDL <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/output/LDL_rank.assoc.txt")
MESA_LDL$p_wald <- as.numeric(MESA_LDL$p_wald)
MESA_LDL <- MESA_LDL[complete.cases(MESA_LDL),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_LDL_qq.png")
qq(MESA_LDL$p_wald, main = "MESA HIS, LDL")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/SNP/MESA_LDL_qq.pdf")
qq(MESA_LDL$p_wald, main = "MESA HIS, LDL")
dev.off()

MESA_CHOL_PX <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/pred_exp/output/CHOL_rank_all.csv")
MESA_CHOL_PX$start_bp <- NULL
MESA_CHOL_PX <- MESA_CHOL_PX[complete.cases(MESA_CHOL_PX),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_CHOL_qq.png")
qq(MESA_CHOL_PX$p_wald, main = "MESA, CHOL")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_CHOL_qq.pdf")
qq(MESA_CHOL_PX$p_wald, main = "MESA, CHOL")
dev.off()

MESA_HDL_PX <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/pred_exp/output/HDL_rank_all.csv")
MESA_HDL_PX$start_bp <- NULL
MESA_HDL_PX <- MESA_HDL_PX[complete.cases(MESA_HDL_PX),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_HDL_qq.png")
qq(MESA_HDL_PX$p_wald, main = "MESA, HDL")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_HDL_qq.pdf")
qq(MESA_CHOL_PX$p_wald, main = "MESA, HDL")
dev.off()

MESA_TRIG_PX <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/pred_exp/output/TRIG_rank_all.csv")
MESA_TRIG_PX$start_bp <- NULL
MESA_TRIG_PX <- MESA_TRIG_PX[complete.cases(MESA_TRIG_PX),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_TRIG_qq.png")
qq(MESA_TRIG_PX$p_wald, main = "MESA, TRIG")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_TRIG_qq.pdf")
qq(MESA_CHOL_PX$p_wald, main = "MESA, TRIG")
dev.off()

MESA_LDL_PX <- fread("/home/angela/px_his_chol/MESAreplication/GEMMA/pred_exp/output/LDL_rank_all.csv")
MESA_LDL_PX$start_bp <- NULL
MESA_LDL_PX <- MESA_LDL_PX[complete.cases(MESA_LDL_PX),]
png("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_LDL_qq.png")
qq(MESA_LDL_PX$p_wald, main = "MESA, LDL")
dev.off()
pdf("/home/angela/px_his_chol/Manuscript_figures/qqman/PX/MESA_LDL_qq.pdf")
qq(MESA_CHOL_PX$p_wald, main = "MESA, LDL")
dev.off()
'