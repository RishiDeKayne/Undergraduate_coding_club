Undergraduate Coding Club - Kelley Lab 2024
================
Week 10
================

This week we will be discucssing principal component analysis (PCA) as a method to better understanding genomic similarities and differences between individuals. The data will be real data from my project carrying out sequencing on 160 _Poecilia mexicana_.

![pca](images/pca.png) 

## Preparation Steps  

Before producing and plotting our PCA and learning to modify differnet aspects of our plots we need to actually produce the data that goes into the PCA. Although you will not be running this (I have pre-prepared the files needed for the tutorial) I have put some example scripts below so you can get a feel for what this entails.  


First we will take our raw concatenated VCF file, which contains high quality genotype calls from across the genome and filter out SNPs with missing data:  

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

Finally, we sort the output using `bcftools sort` specifying a temporary directory that we already made in the first line of the script above using `--temp-dir /hb/scratch/rdekayne/genotype_temp_concat_filt`. For more details on using BCFtools check out the [BCFtools manual](https://samtools.github.io/bcftools/bcftools.html)

Great! Now we have our filtered VCF file - this filtered file will have many fewer SNPs than the full file due to our stringent filtering. In our case we retain almost 5.5 million SNPs which sounds like a lot but remember our genome contains ~800 million sites!  

Next we will actually use PLINK to make our PCA:  

```
#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --time=0-12:00:00
#SBATCH --partition=128x24
#SBATCH --output=pca.%j.out # output file
#SBATCH --error=pca.%j.err # error file
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB 

module load plink

cd /hb/home/rdekayne/PmexWGS160/06_PCA

plink --vcf PmexWGS160_mindepth7_minqual30_mac1_AN320.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out 160_unfilt_out

touch /hb/home/rdekayne/PmexWGS160/06_PCA/160_pca.done
```

This script moves into my PCA directory using `cd`. Then I run PLINK specifying my VCF file (`--vcf`). Since PLINK is made for human chromosomes I tell PLINK to `--allow-extra-chr` so that it is not expecting only human named chromosomes and I specify the fomat in which I want my variant IDs i.e. how SNPs are coded in my output `--set-missing-var-ids @:#`. `--make-bed` means there will be a `.bed` file as output which I often use for other downstream applications and finally the most important part - `--pca` produces our PCA files. We also specify a prefix with PLINK which will be placed before all output files and we do that with `--out`. Look at the [PLINK manual](https://www.cog-genomics.org/plink/) to learn more about the 100s of things PLINK can do.  

Our key PLINK output files are:  
- `160_unfilt_out.eigenval` which contains the eigenvalues from our analysis (these will be converted later to % variance explained by each PC axis).
- `160_unfilt_out.eigenvec` which contains the eigenvectors from our analysis (there are what we will use to actually plot where individuals sit in PC space relative to one another.

I have also made a `.csv` file with background information on all 160 individuals which you can find [HERE](https://github.com/RishiDeKayne/Undergraduate_coding_club/blob/main/Week10/Samples_selected_background_shape_col.csv).  

Now it's time to head to [the R script](). If you are following along after today's session you will need to download the R script along with the background, `.eigenvec`, and `.eigenval` files. Make a folder on your machine and place all of these files there. Open the R script using R studio and change the working directory e.g. with `setwd` to the directory that contains all of the important files. You will also probably need to update the filepath of the background info file. Also be mindful you may have to install and load some of the libraries I use in the tutorial before you are able to run through the script.  

A fantastic resource written by friends and former colleagues of mine Mark Ravinet and Joana Meier contains example scripts and explanations for a bunch of adaptation and speciation genomics analyses. It can be found [HERE](https://speciationgenomics.github.io/pca/) - I have linked to the PCA page which you'll notice is very similar to my R script but take some time to check out the other information there - this is a super valuable resource and I still find myself coming back to this we page frequently!

