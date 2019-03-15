# ACE_GWAS

These are scripts for GenR. GenR contains files that are imputed vcf files. The stages of GWAS analyses include:

1. Convert vcf to Plink format files
2. Quality control
3. GWAS


## Step 1: VCF to Plink

Let's next convert it to plink binary files keeping only the files we need

```{bash}
for i in {1..22}; do ./plink --vcf ./Yourdirectory/Yourfilenamechr${i}.dose.vcf.gz --make-bed --out ./Yourdirectory/Yourfilename_chr${i}  --const-fid 0; done
```

Let's do some QC. Here, we are removing SNPs with 0.01 < ALT_Freq < 0.99. This will ensure that all our SNPs have a minor allele freq > 0.01. We also want to remove poorly imputed SNPs. We do this by removing SNPs with Rsq < 0.3. Finally, we create a file to that can be passed on to the --extract command in Plink. We will do this in R

```{R}
for (i in 1:22){
  a = read.table(paste0("chr", i, ".info"), header = T)
  a$Rsq = as.numeric(as.character(a$Rsq))
  b = subset(a, Rsq > 0.3)
  b = subset(b, ALT_Frq > 0.01 & ALT_Frq < 0.99)
  write.table(b[,1], file = paste0("chr", i, "extract.txt"), row.names = F, col.names = T, quote = F)
}

```

Next, let's combine all the files, and recode the SNP IDs. To recode the SNP IDs, you need to download the VCF file that matches your build (GRCh37), and then manipulate it to get the file you need to pass that onto plink. We will restrict it to only the common variants.

```{bash}
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz
zgrep -v "^##" common_all_20170710.vcf.gz | cut -f1-3 > fileforrecoding.txt
awk '{print $1":"$2"\t"$3}' < fileforrecoding.txt > plinkrecodingfile.txt

```
Merge the files, update SNP name, and do QC

```{bash}
for i in {1..22}; do ./plink --bfile ./Yourdirectory/Yourfilename_chr${i} --extract ./Yourdirectory/chr${i}exclude.txt --make-bed --maf 0.05 --geno 0.05 --hwe 0.000001 --out ./Yourdirectory/Yourfilename2_chr${i}; done
      
./plink --bfile ./Yourdirectory/Yourfilename2_chr22 --merge-list ./Yourdirector/mergelist.txt --make-bed -biallelic-only --out ./Yourdirectory/Yourfilename3

for i in {1..22}; do ./plink --bfile ./Yourdirectory/Yourfilename2_chr${i} --exclude ./Yourdirectory/Yourfilename3.missnp --make-bed --out ./Yourdirectory/Yourfilename2_v2_chr${i}; done

./plink --bfile ./Yourdirectory/Yourfilename2_v2_chr22 --merge-list 1Mv1mergelist2.txt --make-bed -biallelic-only --out ./Yourdirectory/Yourfilename3

./plink --bfile ./Yourdirectory/Yourfilename3 --maf 0.01 --update-name ~/SFARI/liftOverPlink/plinkrecodingfile.txt --hwe 0.000001 --geno 0.05
```

## Step 3
Finally, update the Fam file, as this gets messed up in the whole process. To do this, open the fam file of the imputed merged version, and the fam file of the non-imputed version in R.

```{R}
library(data.table)
library(tidyr)
setwd()

fileimputed = fread("Yourfilename3 .fam")

fileimputed = fileimputed %>% separate(V2, into = c('FID', 'IID'), sep = 6)
setnames(filenonimputed, 2, "IID")

merged = merge(fileimputed, filenonimputed, by = "IID")

merged$oldIID = paste0(merged$FID, merged$IID)

setnames(merged, "IID", "NewIID")
setnames(merged, "V1.x", "oldFID")
setnames(merged, "V1.y", "newFID")

write.table(merged[,c("oldFID", "oldIID", "newFID", "NewIID")], file = "1Mv1updatenames.txt", row.names = F, col.names = F, quote = F)
write.table(filenonimputed[,1:4], file = "updateparents.txt", row.names = F, col.names = F, quote = F)
write.table(filenonimputed[,c(1,2,5)], file = "updatesex.txt", row.names = F, col.names = F, quote = F)
write.table(filenonimputed[,c(1,2,6)], file = "updatepheno.txt", row.names = F, col.names = F, quote = F)
```
