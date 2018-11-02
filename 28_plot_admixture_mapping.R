library(data.table)
library(ggplot2)
HDL_16 <- fread("/home/angela/px_his_chol/GEMMA/output/HDL_rank_sig_snps_p-0.txt")
#gene <- c(110230418, 110251661)
#CHOL_rank_chr1_gene <- subset(CHOL_rank_chr1, (gene[1] - 1000000) <= ps & ps <= (gene[2] + 1000000))
#CHOL_rank_chr1_gene$location <- ""
#for (row in 1:nrow(CHOL_rank_chr1_gene)) {
#  if(as.numeric(CHOL_rank_chr1_gene[row, 3]) > gene[1] & as.numeric(CHOL_rank_chr1_gene[row, 3]) < gene[2]) {
#    CHOL_rank_chr1_gene[row, 15] <- "in gene"
#  }else{
#    CHOL_rank_chr1_gene[row, 15] <- "not in gene"
#  }
#}
#CHOL_rank_chr1_gene$location <- as.factor(CHOL_rank_chr1_gene$location)

HDL_NAT_16 <- fread("/home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/output/chr16_HDL_rank_NAT.assoc.txt")
HDL_IBS_16 <- fread("/home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/output/chr16_HDL_rank_IBS.assoc.txt")
HDL_YRI_16 <- fread("/home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/output/chr16_HDL_rank_YRI.assoc.txt")

#entire chr
ggplot() + 
  #geom_point(data = HDL_16, aes(x = ps, y = -log10(p_wald))) +
  geom_line(data = HDL_NAT_16, aes(x = ps, y = -log10(p_wald), color = "NAT")) +
  geom_line(data = HDL_IBS_16, aes(x = ps, y = -log10(p_wald), color = "IBS")) +
  geom_line(data = HDL_YRI_16, aes(x = ps, y = -log10(p_wald), color = "YRI")) + 
  coord_cartesian(xlim = c(min(HDL_IBS_16$ps), max(HDL_IBS_16$ps)), ylim = c(0, 2.5)) + 
  #coord_cartesian(xlim = c(min(HDL_IBS_16$ps), max(HDL_IBS_16$ps)), ylim = c(0, 10)) + 
  ggtitle("Admixture mapping, HDL") +
  xlab("bp (chr. 16)") + 
  ylab("-log10(p)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + 
  guides(color = guide_legend(title = "Ancestry"), fill = guide_legend(title = "Location")) + 
  scale_color_brewer(type = "div", palette = "Dark2") + 
  scale_fill_brewer(palette = "Paired", direction = -1)
