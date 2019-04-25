#Make scree plot of KING PCs to justify 5
#WHY AREN'T EIGENVALUES STORED IN A FILE UGH RUN AGAIN
library(ggplot2)

eigen <- c(4081.32, 2665.49, 1324.21, 1194.32, 1170.90, 1132.73, 1090.39, 1045.27, 988.82, 946.72, 932.58, 920.77, 913.25, 865.50, 851.32, 835.64, 820.70, 781.40, 775.50, 763.24)
eigen_percent <- eigen/sum(eigen)
eigen_percent <- as.data.frame(eigen_percent)
eigen_percent$index <- as.numeric(rownames(eigen_percent))
colnames(eigen_percent) <- c("percent", "index")

#plot(eigen_percent$eigen_percent, eigen_percent$index, ylab = "Percent of variance", xlab("PC #"))
  #why aren't you working

scree <- ggplot() + 
  geom_point(data = eigen_percent, aes(x = index, y = percent)) + 
  theme_bw() + 
  labs(x = "PC #", y = "Percent of variance explained") + 
  theme(text = element_text(size = 15))
scree

pdf("/home/angela/px_his_chol/Manuscript_figures/S2Fig.pdf", width = 4, height = 4)
print(scree)
dev.off()

tiff("/home/angela/px_his_chol/Manuscript_figures/S2Fig.tiff", width = 15.24, height = 7.62, units = 'cm', res = 300, compression = 'lzw')
print(scree)
dev.off()
