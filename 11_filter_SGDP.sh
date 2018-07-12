#Pipeline to merge and convert all downloaded genotypes from SGDP
cd /home/angela/px_his_chol/SGDP_filtered/ 
declare -a arr=("LP6005441-DNA_A04" "LP6005441-DNA_A12" "LP6005441-DNA_B04" "LP6005441-DNA_B12" "LP6005441-DNA_E10" "LP6005441-DNA_F10" "LP6005441-DNA_G06" "LP6005441-DNA_G07" "LP6005441-DNA_H06" "LP6005441-DNA_H07" "LP6005443-DNA_A12" "LP6005443-DNA_E11" "LP6005443-DNA_F11" "LP6005443-DNA_G11" "LP6005443-DNA_H11" "LP6005519-DNA_D01" "LP6005519-DNA_G02" "LP6005677-DNA_D01" "LP6005677-DNA_E01" "LP6005677-DNA_F01") 

for ind in "${arr[@]}"
do
  echo "Starting filtering on" $ind
  #1. Filter by FL = 9
    #"filter level in range 0-9 or no value (non-integer: N,?) with zero being least reliable; to threshold at FL=n, use all levels n-9"
  bcftools filter --threads 40 -i 'FORMAT/FL=="9"' --output $ind.FL9.vcf /home/wheelerlab3/Data/SGDP/$ind.annotated.nh.vcf.gz
  bgzip $ind.FL9.vcf
  echo "Done filtering on" $ind
  #2. Create index to be used in annotation
  bcftools index --threads 40 -f --tbi -o $ind.FL9.vcf.tbi $ind.FL9.vcf.gz
  echo "Done indexing on" $ind
  #3. Add rsid, REF, and ALT alleles
  bcftools annotate --threads 40 -a anno/homo_sapiens.vcf.gz -c CHROM,POS,ID,REF,ALT --output $ind.FL9.anno.vcf $ind.FL9.vcf.gz
  echo "Done annotating on" $ind
  #4. Print known biallelic sites only (ID column is not ".")
  bcftools view --threads 40 --known -m2 -M2 -v snps -o $ind.FL9.anno.known.vcf $ind.FL9.anno.vcf
  bgzip $ind.FL9.anno.vcf
  bgzip $ind.FL9.anno.known.vcf
  echo "Done viewing on" $ind
done

#merge all samples together
#bcftools merge [OPTIONS] A.vcf.gz B.vcf.gz […] 
