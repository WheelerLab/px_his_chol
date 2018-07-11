#1. Extract one chromosome at a time
bcftools view /home/wheelerlab3/Data/SGDP/LP60f.gz --threads 70 --regions 1 > test_chr1.vcf.gz

#2. Filter by FL = 9
  #"filter level in range 0-9 or no value (non-integer: N,?) with zero being least reliable; to threshold at FL=n, use all levels n-9"
bcftools filter -i 'FORMAT/FL=="9"' test_chr1.vcf.gz --threads 70 > test_chr1_FL9.vcf

#3. Make annotation files from 1000G phase 3
python /home/angela/px_his_chol/SGDP_filtered/convert_SGDP_to_PLINK.py
