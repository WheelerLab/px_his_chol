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
  bcftools index --threads 40 -f --tbi -o $ind.FL9.vcf.tbi $ind.FL9.vcf.gz
  echo "Done indexing on" $ind
done

#2. Merge samples
bcftools merge --threads 40 --merge all LP*.FL9.vcf.gz > SGDP.FL9.vcf
bgzip SGDP.FL9.vcf
bcftools index --threads 40 -f --tbi SGDP.FL9.vcf.gz > SGDP.FL9.vcf.tbi
echo "Samples have been merged"

#3. Add rsid, REF, and ALT
bcftools annotate --threads 40 -a /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423.vcf.gz -c CHROM,POS,ID SGDP.FL9.vcf.gz > SGDP.FL9.rsid.vcf
  #add rsid
bgzip SGDP.FL9.rsid.vcf
bcftools index --threads 40 -f --tbi SGDP.FL9.rsid.vcf.gz > SGDP.FL9.rsid.vcf.tbi
echo "rsids have been added"
bcftools +fixref SGDP.FL9.rsid.vcf.gz -Ob -o SGDP.FL9.match.vcf --threads 40 -- -d -f /home/angela/human_g1k_v37.fasta -i /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423.vcf.gz
  #add REF/ALT
bcftools sort SGDP.FL9.match.vcf -Ob -o SGDP.FL9.match.sorted.vcf
echo "Samples have been sorted"

#4. Filter to SNPs only
bcftools view --threads 40 --known -m2 -M2 -v snps SGDP.FL9.match.sorted.vcf > SGDP.FL9.match.sorted.known.vcf
bgzip SGDP.FL9.match.sorted.vcf
bgzip SGDP.FL9.match.sorted.known.vcf
bcftools index --threads 40 -f --tbi SGDP.FL9.match.sorted.known.vcf.gz > SGDP.FL9.match.sorted.known.vcf.tbi
echo "Samples have been filtered to SNPs only"
