#a subset of UMich_vcf2_px that translate dosage w/ cpos to dosage w/ rsid
import gzip
chrs = range(1, 23)

for chro in chrs:
    c = str(chro)
    
    ##make dictionary: keys -> positions values -> rsids
    posdict = {}
    snpfile = "/home/wheelerlab3/Data/imputation_reference/1000g/1000GP_Phase3_chr"+c+".legend.gz"
    for line in gzip.open(snpfile):
        if(line.startswith('id')):
            continue
        arr = line.strip().split()
        if(arr[0].count(":") == 3): #exclude CNV
            (rs, pos, a0, a1) = arr[0].split(":")
            cpos = str(c) + ":" + str(pos)
            if(rs.startswith("rs")):
                posdict[cpos] = rs

    #read in cpos file and replace column 2 w/ rsid
    indosagefile = "dosages/chr" + c + ".txt.gz"
    outdosagefile = gzip.open("dosages_rsid/chr" + c + ".txt.gz","wb")
    for line in gzip.open(indosagefile):
        arr = line.strip().split()
        cpos = arr[1]
        if cpos in posdict:
            arr[1] = posdict[arr[1]]
            outdosage_SNP = " ".join(arr) + "\n"
            outdosagefile.write(outdosage_SNP)
    print("Completed with chr. " + c + ".")






