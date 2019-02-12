#NOT FOR HCHS BUT CLOSE ENOUGH REPO
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyr)
setwd("/home/angela/loc_anc_plots/")
FB <- fread("1000G_AFA_chr22.rfmix.2.ForwardBackward.txt", header = F) #Forward-Backward: Have *.ForwardBackward.txt suffix. Not output unless --forward-backward option included. Same format as Viterbi, except instead  of one column per haplotype, each haplotype gets one column per ancestry.  The value is the posterior probability of that ancestry at that SNP in that  haplotype. For example, if there are 2 reference populations, the first row  will be: <probability of ancestry 1 at SNP 1 in admixed haplotype 1>,  <probability of ancestry 2 at SNP 1 in admixed haplotype 1>, <probability of  ancestry 1 at SNP 1 in admixed haplotype 2>, etc. If doing phase correction option, these posterior probabilities are calculated given the rephrasing calculated for that iteration.
phind <- fread("1000G_AFA_chr22.phind", header = F)
snp_locations <- fread("1000G_AFA_chr22.snp_locations")
colnames(snp_locations) <- "snp_locations"
'
hap_anc <- c(crossing(phind$V1, c("IBS", "NAT", "YRI")))
colnames(FB) <- paste(hap_anc$`phind$V1`, hap_anc$`c("IBS", "NAT", "YRI")`, sep = "_")
'
hap_anc <- c(crossing(phind$V1, c("CEU", "YRI")))
colnames(FB) <- paste(hap_anc$`phind$V1`, hap_anc$`c("CEU", "YRI")`, sep = "_")

p <- list() #https://stackoverflow.com/questions/9315611/grid-of-multiple-ggplot2-plots-which-have-been-made-in-a-for-loop

#take out 10 random people
#inds <- seq(1, ncol(FB), by = 6) #each person is 6 cols
inds <- seq(829, ncol(FB), by = 4)
for(iteration in c(1:10)){
  ind <- sample(inds, 1)
  ind_prob <- FB[, ind:(ind + 3)]
  '
  ind_prob$IBS <- ind_prob[, 1] + ind_prob[, 4]
  ind_prob$NAT <- ind_prob[, 2] + ind_prob[, 5]
  ind_prob$YRI <- ind_prob[, 3] + ind_prob[, 6]
  ind_prob <- ind_prob %>% dplyr::select(IBS, NAT, YRI)
  '
  ind_prob$CEU <- ind_prob[, 1] + ind_prob[, 3]
  ind_prob$YRI <- ind_prob[, 2] + ind_prob[, 4]
  ind_prob <- ind_prob %>% dplyr::select(CEU, YRI)
  
  ind_prob <- cbind(snp_locations, ind_prob)
  ind_prob <- ind_prob[!duplicated(ind_prob$snp_locations),]
  ind_melt <- melt(ind_prob, id = "snp_locations")
  colnames(ind_melt) <- c("snp_locations", "Ancestry", "prob_of_ancestry")
  if(iteration < 10){
    p[[iteration]] <- ggplot() + 
      geom_bar(data = ind_melt, aes(y = prob_of_ancestry, x = snp_locations, col = Ancestry, fill = Ancestry), stat = "identity") + 
      scale_fill_brewer(palette = "Set1") + 
      scale_color_brewer(palette = "Set1") + 
      theme_bw() +
      theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") #+
      #facet_grid(snp_locations ~ .)#, 
        #plot.margin = unit(c(0, .5, -1.5, 0)))
  }else{
    p[[iteration]] <- ggplot() + 
      geom_bar(data = ind_melt, aes(y = prob_of_ancestry, x = snp_locations, col = Ancestry, fill = Ancestry), stat = "identity") + 
      scale_fill_brewer(palette = "Set1") + 
      scale_color_brewer(palette = "Set1") + 
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "bottom",
            plot.margin = unit(c(0, 0, 0, 0), "cm")) +  
      xlab("Position on chr. 22 (cM)") + 
      ylab("Probability of ancestry")
  }
}

g_legend <- function(a.gplot){ #https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(p[[10]])

p3 <- grid.arrange(grobs = p, mylegend, nrow = 10)
p3
