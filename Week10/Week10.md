Undergraduate Coding Club - Kelley Lab 2024
================
Week 10
================

This week we will be discucssing principal component analysis (PCA) as a method to better understanding genomic similarities and differences between individuals. I will work through my PCA plotting script [script link]

## Preparation Steps  

Before plotting our PCA and better understanding how to modify differnet aspects of our plots we need to actually produce the data that goes into the PCA. Although you will not be running this (I have pre-prepared the files needed for the tutorial I have put some example scripts below so you can get a feel for what this entails.  


First we will take our raw concatenated VCF file, which contains high quality genotype calls from across the genome and filter out SNPs with missing data.

```
#!/bin/bash
#SBATCH --job-name=filt
#SBATCH --time=0-12:00:00
#SBATCH --partition=256x44
#SBATCH --output=geno_filt.%j.out # output file
#SBATCH --error=geno_filt.%j.err # error file
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB 

module load htslib/1.17 sambamba picard mosdepth samtools bcftools

mkdir -p /hb/scratch/rdekayne/genotype_temp_concat_filt

#no missing filter for SNPs
bcftools filter -Oz PmexWGS160_mindepth7_minqual30.vcf.gz -e 'INFO/MAC < 1 | AN < 320' | bcftools view -m2 -M2 -v snps | bcftools sort -Oz --temp-dir /hb/scratch/rdekayne/genotype_temp_concat_filt -o PmexWGS160_mindepth7_minqual30_mac1_AN320.vcf.gz

touch genotype_stats.done
```

In the script above we use `-e 'INFO/MAC < 1 | AN < 320'` to specify that we want to exclude `-e` SNPs with a minor allele count (i.e. the count of less requent allele calls) of < 1 i.e. invariant sites and SNPs where the count of alleles `AN` is less than 320 - because we have 160 individual calls our diploid allele count is 160*2=360 meaning we will exclude any sites with missing data.  

Next we use `bcftools view -m2 -M2` to retain only sites with a minimum `-m2` of 2 alleles and a maximum `-M2` of two alleles - meaning we now have only biallelic SNPs.  

Finally, we sort the output using `bcftools sort` specifying a temporary directory that we already made in the first line of the script above using `--temp-dir /hb/scratch/rdekayne/genotype_temp_concat_filt`.   

Great! Now we have our filtered VCF file - this filtered file will have many fewer SNPs than the full file due to our stringent filtering. In our case we retain almost 5.5 million SNPs which sounds like a lot but remember our genome contains ~800 million sites!  

Next we will actually use PLINK to make our PCA.  

