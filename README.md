# px_his_chol
Scripts used in studying the lipid traits of the Hispanic Community Health Study and the Multi-Ethnic Study of Atherosclerosis, Hispanic

01_PrediXcan_dosages_to_GEMMA_BIMBAM.py: converts PrediXcan-style dosages (from UMich_vcf2px.py or Sanger_vcf2px.py) to GEMMA BIMBAM style

02_PrediXcan_dosages_to_GEMMA_anno.py: takes information from PrediXcan-style dosages to make SNP annotation files for GEMMA

03_pred_exp_to_GEMMA_input.R: converts predicted expression (PrediXcan output) to GEMMA "genotype" input

04_make_KING_GRM.R: makes GEMMA-style GRM using KING output

05_sig_genes_GEMMA_pred_exp.R: finds most significant genes (both FDR < 0.05 and p < bonferroni correction) from GEMMA predicted expression LMM output

06_sig_SNP_GEMMA.py: finds significant SNPs in a series of GEMMA files

07_1000_to_PCA.py: extract YRI, CEU, and PEL from 1000G to use as anchors in Hispanic PCA

08_compare_MESA.py: extract SNP weights from .db files

09_make_bfile_to_RFMix: converts from PLINK --recode beagle output to RFMix input (HAS SEVERE MEMORY ISSUES)

10_chromosome_painting_for_R.py: make "chromosome painting" input (for R's ggplot) from LAMP output

11_filter_Native_American_WGS.sh: Convert and filter individual vcf files to PLINK format

12_pi1.R: Calculates pi1 across phenotypes b/w HCHS and MESA HIS

13_make_KING.R: Runs KING and converts output to both matrix and GCTA format

14_make_eigenvalue_plot.R: make principal component eigenvalue plot as seen in PAGE fig. 4B

15_ancestry_pipeline.sh: personal version of https://github.com/armartin/ancestry_pipeline
