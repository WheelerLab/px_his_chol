#subsets input for GEMMA based on individuals present in analysis
#for use only in HCHS

library(argparse)
library(data.table)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--ind_file_name", help="Individual file produced by 19_loc_anc.py", required = T)
parser$add_argument("--output_suffix", help = "Suffix for output.")
args <- parser$parse_args()

ind_file_name <- args$ind_file_name
if(is.null(args$output_suffix)){
  output_suffix <- ""
}else{
  output_suffix <- args$output_suffix
}

load("/home/angela/px_his_chol/LMM-OPS/run_LMM-OPS_chr22.RData") #KINGmat is stored in here
ind <- fread(ind_file_name, header = F, sep = ":")

#make subsetted relatedness matrix
KINGmat_v2 <- KINGmat[ind$V2,, drop = FALSE]
KINGmat_v2 <- KINGmat_v2[,ind$V2, drop = FALSE]
relatedness_file_name <- paste("relatedness_", output_suffix, ".txt", sep = "")
fwrite(as.data.frame(KINGmat_v2), relatedness_file_name, sep = "\t", row.names = F, col.names = F, quote = F, na = "NA")

#make subsetted phenotype file
pheno <- fread("/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt", header = T)
pheno <- left_join(ind, pheno, by = c("V2" = "IID"))
pheno <- pheno %>% select(CHOL_rank, HDL_rank, TRIG_rank, LDL_rank)
pheno_file_name <- paste("pheno_", output_suffix, ".txt", sep = "")
fwrite(pheno, pheno_file_name, sep = "\t", row.names = F, col.names = F, quote = F, na = "NA")

#make subsetted covariates file
region <- fread("/home/angela/px_his_chol/editedPheno/14_regions.txt", header = F)
region[region == "Unknown"] <- NA
region$V3 <- as.integer(as.factor(region$V3))
pcs <- fread("/home/angela/px_his_chol/KING/kingpc.ped")
pcs <- pcs %>% select(V2, V7, V8, V9, V10, V11)
covariates <- left_join(pcs, region, by = "V2")
covariates <- left_join(ind, covariates, by = "V2")
covariates$V1.x <- NULL
covariates$V2 <- NULL
covariates_file_name <- paste("covariates_", output_suffix, ".txt", sep = "")
fwrite(covariates, covariates_file_name, sep = "\t", row.names = F, col.names = F, quote = F, na = "NA")

