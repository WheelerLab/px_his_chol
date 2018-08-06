#Personal version of https://github.com/armartin/ancestry_pipeline

#RUNNING THE SAMPLE TOY DATA
cd /home/angela/px_his_chol/ancestry_pipeline/toy_data/
#for i in {1..22};
for i in {22..22};
do 
  #Note: the original .bed/.bim/.fam already has combined reference/test
  #1. Split into chromosomes
  plink --bfile ACB_example --chr ${i} --make-bed --out ACB_example_chr${i};
  plink --bfile ACB_example_chr${i} --geno 0.01 --maf 0.01 --make-bed --out ACB_example_chr${i};
  #2. Check data is aligned to reference
  shapeit -check --input-ref /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.hap.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.legend.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3.sample -B ACB_example_chr${i} --input-map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --output-log ACB_example_chr${i}.mendel
  #3. Phasing the dataset using the reference panel of haplotypes (long step)
  shapeit --input-ref /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.hap.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.legend.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3.sample -B ACB_example_chr${i} --duohmm --input-map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --exclude-snp ACB_example_chr${i}.mendel.snp.strand.exclude --output-max ACB_example_chr${i}.haps.gz ACB_example_chr${i}.sample
  #4. Convert to RFMix input
  python ../shapeit2rfmix.py --shapeit_hap_ref ACB_example_chr${i}.haps.gz --shapeit_hap_admixed ACB_example_chr${i}.haps.gz --shapeit_sample_ref ACB_example_chr${i}.sample --shapeit_sample_admixed ACB_example_chr${i}.sample --ref_keep ACB_example.ref --admixed_keep ACB_example.notref --chr ${i} --genetic_map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --out ACB_example
done

#5. Fix classes file for RFMix
python ../classes.py --ref CEU_example.keep,YRI_example.keep --sample CEU_YRI_ACB.sample --out CEU_YRI_ACB.classes
  #MUST KEEP CONSISTENT WITH INDVIDUAL ORDER IN THE ALLELES FILE

#6. Run RFMix (very very long step)
cd /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/
#for i in {1..22};
for i in {22..22}
do
  python RunRFMix.py -e 2 -w 0.2 --num-threads 40 --use-reference-panels-in-EM --forward-backward PopPhased /home/angela/px_his_chol/ancestry_pipeline/toy_data/ACB_example_chr${i}.alleles /home/angela/px_his_chol/ancestry_pipeline/toy_data/CEU_YRI_ACB.classes /home/angela/px_his_chol/ancestry_pipeline/toy_data/ACB_example_chr${i}.snp_locations -o /home/angela/px_his_chol/ancestry_pipeline/toy_data/CEU_YRI_ACB_chr${i}.rfmix
done

#7. Collapse RFMix output into TRACTS-compatible bed files (needs all chrs. to be complete first)
python ../collapse_ancestry.py --rfmix CEU_YRI_ACB_chr22.rfmix.2.Viterbi.txt --snp_locations ACB_example_chr22.snp_locations --fbk CEU_YRI_ACB_chr22.rfmix.2.ForwardBackward.txt --fbk_threshold 0.9 --ind HG02481 --ind_info CEU_YRI_ACB.sample --pop_labels EUR,AFR --chrX --out HG02481

###CONTINUE HERE

#RUNNING ACTUAL DATA (chr 22)
#when running rest of the data, do merging/filtering first, as well as adding cM to .bim
cd /home/angela/px_his_chol/ancestry_pipeline/HCHS_chr22/
for i in {22..22};
do 
  #1. Phasing the dataset using the reference panel of haplotypes (long step)
  ##Large (>12k) cohort, so use HAPI-UR instead of SHAPEIT
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p merged_chr${i}_filtered_ordered -w 64 -o phase_chr${i}
  
  #2. Make additional files for RFMix input
  awk '{print $3}' phase_chr${i}.phsnp > phase_chr${i}.snp_locations
  Rscript make_classes_from_HAPI-UR.R phase_chr22.phind ordered_pops.txt HCHS



done

#The RFMix portion needs to be changed b/c of 3-way admixture (rather than 2)
