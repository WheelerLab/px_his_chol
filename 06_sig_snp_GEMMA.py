#finds significant SNPs in a series of GEMMA files
pheno_list = ['CHOL_rank', 'HDL_rank', 'TRIG_rank', 'LDL_rank']

for pheno in pheno_list:
  sig_snp = open("/home/angela/px_his_chol/GEMMA/output/" + pheno + "_sig_snps.txt", "w")
  sig_snp.write("chr\trs\tps\tn_miss\tallele1\tallele0\taf\tbeta\tse\tl_remle\tl_mle\tp_wald\tp_lrt\tp_score\n")
  for chrom in range(1, 23):
    for line in open("/home/angela/px_his_chol/GEMMA/output/" + pheno + "_chr" + str(chrom) + ".assoc.txt", "r"):
      arr = line.strip().split()
      (chr, rs, ps, n_miss, allele1, allele0, af, beta, se, l_remle, l_mle, p_wald, p_lrt, p_score) = arr[0:14]
      if line.startswith("chr"):
        continue
      elif float(p_wald) < 5e-8: #only sig. SNPs
        is_sig_snp = (chr + "\t" + rs + "\t" + ps + "\t" + n_miss + "\t" + allele1 + "\t" + allele0 + "\t" + af + "\t" + beta + "\t" + se + "\t" + l_remle + "\t" + l_mle + "\t" + p_wald + "\t" + p_lrt + "\t" + p_score + "\n")
        sig_snp.write(is_sig_snp)
  sig_snp.close()
