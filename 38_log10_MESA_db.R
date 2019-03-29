#compare correlations in the MESA model associations b/w HCHS and MESA HIS, MESA CAU, GLGC
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
"%&%" = function(a,b) paste(a,b,sep="")
setwd("/home/angela/px_his_chol/MESA_compare/MESA_db/")

pops <- c("AFA", "CAU", "HIS", "AFHI", "ALL")
MESA_db <- fread("MESA_db_assoc.csv")
MESA_db <- MESA_db %>% dplyr::select(tissue, HCHS_P, pheno, MESA_HIS_P, MESA_CAU_P, GLGC_P)
MESA_db_log10_P <- MESA_db %>% dplyr::select(tissue, HCHS_P, pheno)
MESA_db_log10_P$HCHS_log_P <- -log10(MESA_db$HCHS_P)
MESA_db_log10_P$MESA_CAU_log_P <- -log10(MESA_db$MESA_CAU_P)
MESA_db_log10_P$MESA_HIS_log_P <- -log10(MESA_db$MESA_HIS_P)
MESA_db_log10_P$GLGC_log_P <- -log10(MESA_db$GLGC_P)
thres <- c(1, 0.05, 0.01, 0.005, 0.001)#, 0.0005, 0.0001)
rep_pops <- c("MESA_HIS", "MESA_CAU", "GLGC")

#USING PLOT_GRID
#log_P_results
pops_plots <- list()
for(k in 1:length(pops)){
  cor_mat <- matrix(nrow = 0, ncol = length(rep_pops)) 
  colnames(cor_mat) <- rep_pops
  db <- subset(MESA_db_log10_P, tissue == pops[k])
  for(i in thres){
    thres_df <- subset(db, HCHS_P < as.numeric(i))
    thres_df <- thres_df %>% dplyr::select(HCHS_log_P, MESA_HIS_log_P, MESA_CAU_log_P, GLGC_log_P)
    if(nrow(thres_df[complete.cases(thres_df)]) > 0){
      cor_df <- as.data.frame(cor(thres_df, use = "complete.obs", method = "pearson"))
      cor_df$HCHS_log_P <- NULL
      cor_df <- cor_df[1,]
      rownames(cor_df) <- as.character(i)
      cor_mat <- rbind(cor_mat, cor_df)
    }
  }
  cor_mat$thres <- rownames(cor_mat)
  cor_df <- as.data.frame(cor_mat)
  fwrite(cor_df, pops[k] %&% "_cor_log_P.csv", row.names = F, quote = F)
  
  #plotting
  if(k == 1){
    cor_df$thres <- as.numeric(cor_df$thres)
    cor_df <- cor_df[order(as.numeric(cor_df$thres)),]
    cor_df$thres <- as.factor(cor_df$thres)
    melted <- melt(cor_df, id.vars = "thres")
    pop_plot <- ggplot(data = melted, aes(x = thres, y = value, group = variable)) + 
      geom_line(aes(color = variable)) +
      scale_x_discrete(limits = rev(levels(melted$thres))) +
      ylim(-0.5, 1) + 
      scale_color_brewer(palette = "Set1", labels = c("MESA HIS", "MESA CAU", "GLGC")) + 
      labs(title = pops[k], color = "Pop.", y = "Correlation of -log10(P)", x = "") + 
      theme_bw() + 
      theme(text = element_text(size = 15), legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))
    pops_plots[[k]] <- pop_plot
  }else if(k == 3){
    cor_df$thres <- as.numeric(cor_df$thres)
    cor_df <- cor_df[order(as.numeric(cor_df$thres)),]
    cor_df$thres <- as.factor(cor_df$thres)
    melted <- melt(cor_df, id.vars = "thres")
    pop_plot <- ggplot(data = melted, aes(x = thres, y = value, group = variable)) + 
      geom_line(aes(color = variable)) +
      scale_x_discrete(limits = rev(levels(melted$thres))) + 
      ylim(-0.5, 1) +  
      scale_color_brewer(palette = "Set1", labels = c("MESA HIS", "MESA CAU", "GLGC")) + 
      labs(title = pops[k], color = "Pop.", y = "Correlation of -log10(P)", x = "HCHS P < #") + 
      theme_bw() + 
      theme(text = element_text(size = 15), legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    pops_plots[[k]] <- pop_plot
  }else if(k == 5){
    cor_df$thres <- as.numeric(cor_df$thres)
    cor_df <- cor_df[order(as.numeric(cor_df$thres)),]
    cor_df$thres <- as.factor(cor_df$thres)
    melted <- melt(cor_df, id.vars = "thres")
    pop_plot <- ggplot(data = melted, aes(x = thres, y = value, group = variable)) + 
      geom_line(aes(color = variable)) +
      scale_x_discrete(limits = rev(levels(melted$thres))) +
      ylim(-0.5, 1) + 
      scale_color_brewer(palette = "Set1", labels = c("MESA HIS", "MESA CAU", "GLGC")) + 
      labs(title = pops[k], color = "Pop.", y = "Correlation of -log10(P)", x = "") + 
      theme_bw() + 
      theme(text = element_text(size = 15), legend.position = "right", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    pops_plots[[k]] <- pop_plot
  }else{
    cor_df$thres <- as.numeric(cor_df$thres)
    cor_df <- cor_df[order(as.numeric(cor_df$thres)),]
    cor_df$thres <- as.factor(cor_df$thres)
    melted <- melt(cor_df, id.vars = "thres")
    pop_plot <- ggplot(data = melted, aes(x = thres, y = value, group = variable)) + 
      geom_line(aes(color = variable)) +
      scale_x_discrete(limits = rev(levels(melted$thres))) + 
      ylim(-0.5, 1) + 
      scale_color_brewer(palette = "Set1", labels = c("MESA HIS", "MESA CAU", "GLGC")) + 
      labs(title = pops[k], color = "Pop.", y = "Correlation of -log10(P)", x = "") + 
      theme_bw() + 
      theme(text = element_text(size = 15), legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    pops_plots[[k]] <- pop_plot
  }
}

#PLOTTING
p1 <- pops_plots[[1]]
p2 <- pops_plots[[2]]
p3 <- pops_plots[[3]]
p4 <- pops_plots[[4]]
p5 <- pops_plots[[5]]

plot_grid(p1, p2, p3, p4, p5, align = "h", nrow = 1, rel_heights = c(1, 1, 1, 1, 1), rel_widths = c(1.25, 1, 1, 1, 1.85)) #rel widths depend on the print height



