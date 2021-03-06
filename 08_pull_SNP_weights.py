#pulls SNPs from sig. genes in MESA models
import sqlite3
import pandas
pops = ["AFA", "AFHI", "ALL", "CAU", "HIS"]
input_genenames = ["ATG4C", "ATP13A1", "ATXN1L", "DOCK7", "PSRC1", "USP1", "CETP", "DUS2L", "EPOR", "GFOD2", "LPL", "NLRC5", "NUTF2", "PLA2G15", "TMEM205", "ATG4C", "BRE", "CCDC121", "DOCK7", "EIF2B4", "LPL", "NRBP1", "SNX17", "USP1", "ZNF513", "ATP13A1", "HMGCR", "PSRC1"]
input_genenames = set(input_genenames) #only unique names
db_path = "/home/lauren/files_for_revisions_plosgen/en_v7/dbs/"
output_path = "/home/angela/px_his_chol/MESA_compare/"

for pop in pops:
  for input_genename in input_genenames:
    gene_output = []
    #find gene from gene_name
    for line in open("/home/angela/px_yri_chol/PrediXcan/ChrENGene_forRenaming.txt"):
      arr = line.strip().split()
      (CHR, gene, genename) = arr[0:3]
      genename = genename.replace('"', '')
      if input_genename == genename:
        break
      
    #look for a gene and its SNPs in db
    conn = sqlite3.connect(db_path + pop + "_imputed_10_peer_3_filtered.db")
    c = conn.cursor()
    c.execute('select * from weights where gene = ' + gene + ';')
    for row in c:
      gene_output.append([row[0], row[1], row[5]])
    conn.close()
    gene_output_df = pandas.DataFrame(gene_output, columns = ['gene', 'rs', 'weight'])
    gene_output_df.to_csv(output_path + pop + "_" + input_genename + "_sig_gene_weights.txt", index = False)
