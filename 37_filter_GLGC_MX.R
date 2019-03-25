#filter GLGC MetaXcan results to include only MESA or those w/ green flags
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
setwd("/home/angela/MetaXcan-master/GLGC/")

abbv <- fread("/home/angela/tiss_abbreviations.txt")
abbv2 <- fread("/home/angela/tiss_abbreviations_no_TW.txt")
colnames(abbv2)[1] <- "tissu"
abbv <- left_join(abbv, abbv2, by = "tiss")

results <- data.frame(gene = character(), gene_name = character(), effect_size = numeric(), pvalue = numeric(), pheno = character(), model = character())
for(j in c("CHOL", "HDL", "TRIG", "LDL")){
  for(i in 1:nrow(abbv)){
    model <- abbv[i, 1]
    if(i <= 5){ #MESA
      MX_results <- fread(model %&% "_" %&% j %&% ".csv")
    }else{ #GTEx
      MX_results <- fread(model %&% "/" %&% j %&% "_rank.csv")
    }
    MX_results <- MX_results %>% dplyr::select(gene, gene_name, effect_size, pvalue)
    MX_results$pheno <- j %&% "_rank"
    MX_results$model <- abbv[i, 3]
    results <- rbind(results, MX_results)
  }
}
colnames(results)[2] <- "genename"

#filter by flags
flags <- fread("/home/angela/px_his_chol/GEMMA/pred_exp/output/sig_thres/flags.txt")
flags$gene <- NULL
results <- left_join(results, flags, by = c("model", "genename"))
results <- results %>% mutate(flag = if_else(is.na(flag), "MESA", flag))
fwrite(results, "all_MX_results.csv", quote = F, sep = ",", na = "NA", nThread = 10)
pass_results <- subset(results, flag %in% c("MESA", "green"))
fwrite(pass_results, "pass_MX_results.csv", quote = F, sep = ",", na = "NA", nThread = 10)
