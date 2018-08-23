# px_his_chol
Scripts used in studying the lipid traits of the Hispanic Community Health Study and the Multi-Ethnic Study of Atherosclerosis, Hispanic

ASHG abstract: Plasma lipid levels are risk factors for cardiovascular disease, a leading cause of death worldwide. While many studies have been conducted on lipid genetics, they mainly comprise individuals of European ancestry and thus their transferability to diverse populations is unclear. We performed genome- and transcriptome-wide association studies of four lipid traits in the Hispanic Community Health Study (HCHS) cohort (n = 11,103), with origins in Mexico, Cuba, Puerto Rico, Central America, the Dominican Republic, and South America. We tested our findings for replication in the Hispanic population in the Multi-Ethnic Study of Atherosclerosis (MESA) (n = 822) and compared our results to larger, predominantly European ancestry meta-analyses. In both our GWAS and TWAS, we used a linear mixed model to control for relatedness and included the first five genotypic principal components and geographic region as covariates. In our GWAS, five previously-implicated SNPs reached significance (P < 5 x 10-8). After predicting gene expression levels with PrediXcan software using multi-tissue models built in the Genotype-Tissue Expression Project (GTEx) and multi-ethnic monocyte models built in MESA, we tested genes for association with four lipid phenotypes (total cholesterol, HDL, LDL, triglycerides). This revealed 255 significant gene-phenotype associations (FDR < 0.05) with 84 unique significant genes, many of which occurred across multiple phenotypes, tissues and MESA populations. Of these significant genes, 36 were previously implicated within the GWAS catalog, such as CETP, PSRC1, and DOCK7, and 20/36 replicated in MESA and 30/36 replicated in the Genetic Lipid Global Consortium (GLGC). We found 48 of the significant genes are novel for any lipid association, including TSNAXIP1, which associated with HDL in eight tissues, and C19orf52, which associated with total cholesterol, HDL, and LDL in two tissues. Of the 58 novel gene-phenotype associations found significant, 19 replicated independently in MESA and 42 replicated independently in GLGC at P < 0.05, with 8 associations replicating in both. The largest and most diverse MESA expression prediction model, including African Americans, Caucasians, and Hispanics, had more significant genes (18, n = 1,163) compared to the Caucasian model (8, n = 578), indicating that to fully characterize the impact of genetic variation between populations, larger studies in non-European ancestry populations are needed.

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

16_make_qq.R: makes qq and aggregate gene lists from GEMMA/PrediXcan output

17_make_classes_from_HAPI-UR: makes .classes file for RFMix from HAPI-UR (use in ancestry pipeline)
