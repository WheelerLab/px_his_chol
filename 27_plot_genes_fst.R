#Plots various models' sample size, significant gene associations, and Fst against each other
library(ggrepel)
library(ggplot2)

#MESA
genes_fst <- fread("/home/angela/px_his_chol/MESA_compare/genes_Fst.csv", header = T)
genes_fst_noGTEx <- genes_fst[1:5,]
fig <- ggplot(genes_fst_noGTEx, aes(x = Model_Sample_Size, y = Significant_gene_associations, label = Model)) + 
  #geom_raster() + 
  geom_point(aes(fill = log10(Sig_gene_Fst)), shape = 21, size = 10) + #, size = Model_Sample_Size), shape = 21) + 
  coord_cartesian() + 
  geom_label_repel(aes(label = Model), box.padding = 0.35, point.padding = 0.15, segment.color = NA, hjust = 0.5, direction = "x") +
  xlab("Model sample size") + 
  ylab("Significant gene associations") +
  labs(fill = expression(paste("log10(", F[ST], ")"))) + 
  theme_bw(12) +
  scale_fill_continuous(type = "viridis") + 
  scale_x_continuous(limits = c(100, 1250)) + 
  scale_y_continuous(limits = c(8.5, 18.5)) + 
  theme(text = element_text(size = 15))
print(fig)

#GTEx
Fst_num_genes_n_samples <- read_csv("/home/angela/px_his_chol/MESA_compare/GTEx_WB/SNPs_n_samples/Fst_num_genes_n_samples.csv")
Fst_num_genes_n_samples <- Fst_num_genes_n_samples[complete.cases(Fst_num_genes_n_samples),]
fig <- ggplot(Fst_num_genes_n_samples, aes(x = Sample_size, y = Sig_gene_associations, label = Tiss_abb)) + 
  #geom_raster() + 
  geom_point(aes(fill = log10(Mean_fst)), shape = 21, size = 10) + #, size = Model_Sample_Size), shape = 21) + 
  #coord_cartesian() + 
  #geom_label_repel(aes(label = Tiss_abb), box.padding = 0.35, point.padding = 0.15, segment.color = NA, hjust = 0.5, direction = "x") +
  xlab("Model sample size") + 
  ylab("Significant gene associations") +
  labs(fill = expression(paste("log10(", F[ST], ")"))) + 
  theme_bw(12) +
  scale_fill_continuous(type = "viridis") + 
  scale_x_continuous(limits = c(50, 400)) + 
  scale_y_continuous(limits = c(1, 14)) + 
  theme(text = element_text(size = 15))
print(fig)
