Undergraduate Coding Club - Kelley Lab 2024
================
Week 6
================

This is an explanation of how to go from running scripts to execute commands on a single file to submitting array jobs. The following scripts are just examples and will require modifying to run your own jobs on the cluster!

## Single file jobs

At this point you have likely submitted a single-file job for something like FastQC as follows:

```
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --partition=128x24
#SBATCH --output=QC.out # output file
#SBATCH --error=QC.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=4GB # Memory limit of 4GB

#this line loads the modules
module load fastqc/0.12.1

#this line runs fastqc on one thread (-t 1) outputs the resulting .html
#file to a specific directory (-o /hb/home/rdekayne/PmexWGS160/01_QC/fastqc_reports/)
fastqc -t 1 -o /hb/home/USER/fastqc_reports/ /hb/groups/BLAH/INDIV1_R1_001.fastq.gz

#and then 'touches' this file so we know the script has been run through
#(NOTE: this file will be made even if the job fails and is only an indication
#that the job no longer sits in the queue or is currently being run)
touch INDIV1_R1_001.fastq.gz.done
```

## Multi-file analysis in 1 job

This approach, rather than running a single file, will run through a list of individual names. This can be particularly useful if we already have these names in a text file (and if not we can always make that list). In this case the list of individuals are stored one per line in `sample_list.txt`. Note, this is one job that will carry out all fastqc analysis in series (one after another within the same job) - to see how you can run many jobs in parallel (at once) see below.

```
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --partition=128x24
#SBATCH --output=QC.out # output file
#SBATCH --error=QC.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=4GB # Memory limit of 4GB

#this line loads the modules
module load fastqc/0.12.1

#this line takes our sample list and runs/reads through it line by line in a loop
#the bit between DO and DONE is the bit being looped.

#Sample name can be specified in the rest of the script using ${sample}

#Since fastqc needs to be run on the R1 and R2 file for each individual here we can run both
#one after another for each sample ID (stored as ${sample})

cat sample_list.txt | while read sample
do
fastqc -t 1 -o /hb/home/USER/fastqc_reports/ /hb/groups/BLAH/${sample}_R1_001.fastq.gz
fastqc -t 1 -o /hb/home/USER/fastqc_reports/ /hb/groups/BLAH/${sample}_R2_001.fastq.gz
done

#and again we have a 'touch' statement indicating the loop has finished and the script has run through to the end
touch all_samples.done

```

## Multi-file analysis in an array

This is the  most powerful approach (and probably most broadly applicable). Again we will use a file to split our job up but

Rather than running a single file, will run through a list of individual names. This can be particularly useful if we already have these names in a text file (and if not we can always make that list). In this case the list of individuals are stored one per line in `sample_list.txt`.

```
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --time=0-12:00:00 # Wall clock time limit in Days-Hours:min:seconds
#SBATCH --partition=128x24
#SBATCH --output=QC.out # output file
#SBATCH --error=QC.err # error file
#SBATCH --ntasks=1 # Run 1 job
#SBATCH --ntasks-per-node=1 # One task per computer
#SBATCH --cpus-per-task=2 # 2 CPUs per job
#SBATCH --mem=4GB # Memory limit of 4GB

#this line loads the modules
module load fastqc/0.12.1
