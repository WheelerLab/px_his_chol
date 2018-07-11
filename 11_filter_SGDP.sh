#Pipeline to merge and convert all downloaded genotypes from SGDP
cd /home/angela/px_his_chol/SGDP_filtered/

#1. Extract one chromosome at a time
  #Maybe do all at the same time when pipeline is done
bcftools view /home/wheelerlab3/Data/SGDP/LP6005441-DNA_A04.annotated.nh.vcf.gz --threads 70 --regions 22 > test_chr22.vcf.gz

#2. Filter by FL = 9
  #"filter level in range 0-9 or no value (non-integer: N,?) with zero being least reliable; to threshold at FL=n, use all levels n-9"
bcftools filter -i 'FORMAT/FL=="9"' test_chr22.vcf.gz --threads 70 > test_chr22_FL9.vcf

#3. Create index to be used in annotation
bcftools index -f --tbi test_chr22_FL9.vcf.gz > test_chr22_FL9.vcf.gz.tbi

#4. Add rsid, REF, and ALT alleles
bcftools annotate -a anno/homo_sapiens.vcf.gz -c CHROM,POS,ID,REF,ALT --threads 70 test_chr22_FL9.vcf.gz > test_chr22_FL9_anno.vcf.gz

#5. Print known biallelic sites only (ID column is not ".")
bcftools view --known -m2 -M2 -v snps --threads 70 -o test_chr22_FL9_anno_known.vcf.gz test_chr22_FL9_anno.vcf.gz

#ID list
declare -a arr=("LP6005441-DNA_A04" "LP6005441-DNA_A12" "LP6005441-DNA_B04" "LP6005441-DNA_B12" "LP6005441-DNA_E10" "LP6005441-DNA_F10" "LP6005441-DNA_G06" "LP6005441-DNA_G07" "LP6005441-DNA_H06" "LP6005441-DNA_H07" "LP6005443-DNA_A12" "LP6005443-DNA_E11" "LP6005443-DNA_F11" "LP6005443-DNA_G11" "LP6005443-DNA_H11" "LP6005519-DNA_D01" "LP6005519-DNA_G02" "LP6005677-DNA_D01" "LP6005677-DNA_E01" "LP6005677-DNA_F01")
cd /home/angela/px_his_chol/SGDP_filtered/
for ind in "${arr[@]}"
do
  bcftools filter --threads 10 -i 'FORMAT/FL=="9"' --output-type z --output $ind.FL9 /home/wheelerlab3/Data/SGDP/$ind.annotated.nh.vcf.gz 
  bcftools index --threads 10 -f --tbi -o $ind.FL9.vcf.gz $ind.FL9.vcf.gz
  bcftools annotate --threads 10 -a anno/homo_sapiens.vcf.gz -c CHROM,POS,ID,REF,ALT --output-type z --output $ind.FL9.anno $ind.FL9.vcf.gz
  bcftools view --threads 10 --known -m2 -M2 -v snps --output-type z -o $ind.FL9.anno.known $ind.FL9.anno
done

# merge all samples together
#bcftools merge [OPTIONS] A.vcf.gz B.vcf.gz […] 
