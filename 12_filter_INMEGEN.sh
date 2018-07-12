
#1. Filter by passing SNPs
bcftools filter -i 'FORMAT/FT=="PASS"' /home/wheelerlab3/Data/INMEGEN/12G/MAY1.vcf.gz > MAY1.PASS.vcf
  #leaves 3.8m SNPs
bgzip MAY1.PASS.vcf
bcftools index --threads 40 -f --tbi MAY1.PASS.vcf.gz > MAY1.PASS.vcf.tbi
  
#2. Why is the world is the rsid in info and not in the ID column (and why are there SNPs from like db88 how vintage)
bcftools annotate --threads 40 -a /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423.vcf.gz -c CHROM,POS,ID,ALT MAY1.PASS.vcf.gz > MAY1.PASS.anno.vcf

