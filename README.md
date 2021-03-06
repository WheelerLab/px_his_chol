# px_his_chol
Scripts used in studying the lipid traits of the Hispanic Community Health Study and the Multi-Ethnic Study of Atherosclerosis, Hispanic, with our manuscript avaialble in [bioRxiv](https://www.biorxiv.org/content/early/2018/12/28/507905). This repository is mainly for documentation. For more user-friendly, better commented, and better explained scripts, please refer to [Ad_PX_pipe](https://github.com/aandaleon/Ad_PX_pipe).

* ~~01_PrediXcan_dosages_to_GEMMA_BIMBAM.py: defunct, see 02_PrediXcan_dosages_to_GEMMA.py~~

* 02_PrediXcan_dosages_to_GEMMA_anno.py: takes information from PrediXcan-style dosages to make SNP annotation and BIMBAM files for GEMMA

* 03_pred_exp_to_GEMMA_input.R: converts predicted expression (PrediXcan output) to GEMMA "genotype" input

* ~~04_make_KING_GRM.R: defunct, use 13_make_KING.R~~

* ~~05_sig_genes_GEMMA_pred_exp.R: finds most significant genes (both FDR < 0.05 and p < bonferroni correction) from GEMMA predicted expression LMM output~~

* 06_sig_SNP_GEMMA.py: finds significant SNPs in a series of GEMMA files

* 07_1000_to_PCA.py: extract YRI, CEU, and PEL from 1000G to use as anchors in Hispanic PCA

* 08_compare_MESA.py: extract SNP weights from .db files

* ~~09_make_bfile_to_RFMix: converts from PLINK --recode beagle output to RFMix input (HAS SEVERE MEMORY ISSUES)~~

* ~~10_chromosome_painting_for_R.py: make "chromosome painting" input (for R's ggplot) from LAMP output~~

* ~~11_filter_Native_American_WGS.sh: Convert and filter individual vcf files to PLINK format~~

* ~~12_pi1.R: Calculates pi1 across phenotypes b/w HCHS and MESA HIS~~

* 13_make_KING.R: Runs KING and converts output to both matrix and GCTA format

* 14_make_eigenvalue_plot.R: make principal component eigenvalue plot as seen in PAGE fig. 4B

* ~~15_ancestry_pipeline.sh: personal version of https://github.com/armartin/ancestry_pipeline; phase PLINK files with HAPI-UR to be used in RFMix~~

* 16_make_qq.R: makes qq and aggregate gene lists from GEMMA/PrediXcan output

* ~~17_make_classes_from_HAPI-UR: makes .classes file for RFMix from HAPI-UR (use in ancestry pipeline)~~

* ~~18_convert_MOSAIC_output.R: converts local ancestry output from MOSAIC to a more human-readable format. Output to be used in 19_loc_anc.py.~~

* ~~19_MOSAIC_loc_anc.py: "imputes" local ancestry between markers to use local ancestry as a dosage in GEMMA.~~

* ~~20_GEMMA_wrapper_local_anc_cov.py: wrapper for SNP-by-SNP level GEMMA~~

* ~~21_make_GEMMA_input.R: subsets input for GEMMA based on individuals present in analysis~~

* ~~22_GEMMA_wrapper_admixture_mapping.py: uses local ancestry as a dosage to run an ancestry-by-ancestry level admixture mapping analysis~~

* ~~23_split_haplotype_file.py: MOSAIC input with all haplotypes is too large, so this program chunks the 24k haplotypes into 12 equal chunks for an easier time on the memory~~ 

* 24_BIMBAM_to_PX_dosages.py: converts from BIMBAM format (comma delimited) to PrediXcan dosage format

* ~~25_RFMix_loc_anc.py: translates HAPI-UR and RFMix output into GEMMA-style input dosages to be used in 22_GEMMA_wrapper_admixture_mapping.py~~

* ~~26_local_ancestry_wrapper.py: Goes from the RFMix local ancestry estimates (25), PrediXcan-style dosages (2), and HAPI-UR phasings (15) to full admixture mapping in GEMMA (21, 22)~~

* ~~27_plot_genes_fst.R: plots various models' sample size, significant gene associations, and Fst against each other~~

* ~~28_plot_admixture_mapping.R: makes Manhattan plot with ancestry mixture overlays by ancestry~~

* 29_backward_elimination_sig_genes.R: backward elimination of significant genes by phenotype and chromosome to determine which ones are independent at a locus
 
* 30a/b_make_coloc_input_MESA/GTEx.R: converts from GEMMA results and matrix eQTL format or GTEx download format to input for run_COLOC.py script - https://github.com/hakyimlab/summary-gwas-imputation/wiki/Running-Coloc

* 31_make_COLOC_table.R: process COLOC output and join to novel significant gene associations

* 32_plot_pred_exp.R: make scatterplot of predicted expression vs. phenotype
 
* 33_get_genename_R2.py: pull gene names and R2 from MESA db files through sqlite3 library
 
* 34_make_pred_exp_result_table.R: extracts both all and significant (P < 9.654e-6) PrediXcan results; adds gene name and starting bp
 
Strikethrough: depreciated and not used for final analysis

Abstract: Plasma lipid levels are risk factors for cardiovascular disease, a leading cause of death worldwide. While many studies have been conducted on lipid genetics, they mainly comprise individuals of European ancestry and thus their transferability to diverse populations is unclear. We performed genome-wide (GWAS) and imputed transcriptome-wide association studies of four lipid traits in the Hispanic Community Health Study/Study of Latinos cohort (n = 11,103), tested the findings for replication in the Hispanic population in the Multi-Ethnic Study of Atherosclerosis (MESA HIS, n = 1,297), and compared the results to the larger, predominantly European ancestry meta-analysis by the Global Lipids Genetic Consortium (GLGC, n = 196,475). GWAS revealed both known and potentially novel SNPs. In the imputed transcriptome-wide association study in multiple tissues and ethnicities, we found 101 significant gene-phenotype associations (P < 9.65 x 10-6) with 36 unique significant genes, many of which occurred across multiple phenotypes, tissues, and multi-ethnic populations. We found 21 of the 36 significant genes are novel for any lipid association. Of these, 15/21 replicated in MESA HIS, 14/21 replicated in GLGC, and 10/21 associations replicated in both cohorts (P < 0.05). We identified genes that associate independently from nearby genes and colocalize with expression quantitative trait loci (eQTLs), indicating a possible mechanism of gene regulation in lipid level variation. We also investigated prediction in multi-ethnic versus European-only models. To fully characterize the genetic architecture of lipid traits in diverse populations, larger studies in non-European ancestry populations are needed.
