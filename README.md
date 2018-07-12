# px_his_chol
Scripts used in studying the lipid traits of the Hispanic Community Health Study and the Multi-Ethnic Study of Atherosclerosis

01_PrediXcan_dosages_to_GEMMA_BIMBAM.py: converts PrediXcan-style dosages (from UMich_vcf2px.py or Sanger_vcf2px.py) to GEMMA BIMBAM style

02_PrediXcan_dosages_to_GEMMA_anno.py: takes information from PrediXcan-style dosages to make SNP annotation files for GEMMA

03_pred_exp_to_GEMMA_input.R: converts predicted expression (PrediXcan output) to GEMMA "genotype" input

04_make_KING_GRM.R: makes GEMMA-style GRM using KING output

05_sig_genes_GEMMA_pred_exp.R: finds most significant genes (both FDR < 0.05 and p < bonferroni correction) from GEMMA predicted expression LMM output

06_sig_SNP_GEMMA.py: finds significant SNPs in a series of GEMMA files

07_1000_to_PCA.py: extract YRI, CEU, and PEL from 1000G to use as anchors in Hispanic PCA

08_compare_MESA.py: extract SNP weights from .db files

09_make_bfile_to_RFMix: converts from PLINK --recode beagle output to RFMix input (HAS SEVERE MEMORY ISSUES)

10_chromosome_painting_for_R.py: make "chromsome painting" input (for R's ggplot) from LAMP output

11_filter_SGDP.sh: #Pipeline to merge and convert all downloaded genotypes from SGDP
