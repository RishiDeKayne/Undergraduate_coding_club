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

#This line loads the modules
module load fastqc/0.12.1

#This line runs fastqc on one thread (-t 1) outputs the resulting .html
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

#This line loads the modules
module load fastqc/0.12.1

#Then we use cat to read through our sample list line by line in a loop
#assign whatever is in that line to the variable `sample`
#and then execute the bit between DO and DONE in a loop.

#The sample name variable can be specified in the rest of the script using ${sample}

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

This is the  most powerful approach (and probably most broadly applicable). Again we will use a file to split our job up by individual but this time the different indicuals will be run in parallel rather than series. This means we can spread the jobs across many nodes rather than using a single job on a single node. For jobs that take a lot of time this is important - 100x1h jobs can finish in 1h since jobs can be distributed across 100 nodes rather than taking 100h on a single node. NOTE: the lines after the module load are the key part here where we extract the correct value of `ind` and `indiv_name` from our sample file using the array number or job ID. 

Array jobs like this also have an importnat distinction in that they are submitted differently to regular jobs. Rather than `sbatch my_job.sh` which you would use for the jobs above you need to specify the range of value for the array. This correspons to the number of items in your list file, which could be file names, contigs, or in our case individuals. So here we would have as many rows as we have individuals. Say way have 10 individuals we need to tell the cluster we will run 10 jobs ranting from job number 1 to job number 10. We do this using the following command: `sbatch --array=1-10 my_script.sh`. If the jobs are computationally intensive it is key to specify that not all jobs are to be submitted/run at once. We can do this using `%`. Say we have 10 jobs as before but only want 2 to be run at any one time (to avoid irritating other cluster users) we can run this as `sbatch --array=1-10%2 my_script.sh`.

```
#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --time=0-12:00:00
#SBATCH --partition=128x24
#SBATCH --output=mapping.%j.out # output file
#SBATCH --error=mapping.%j.err # error file
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16GB 

#This line loads the modules
module load fastqc/0.12.1

#Here we produce two variables
#The first, `ind` will take the number of the array job we specify
#If we say array 1-10 this `ind` value will be every value 1-10
ind=${SLURM_ARRAY_TASK_ID}

#Then we produce an individual name variable `indiv_name`
#by extracting a line from our sample_list.txt file
#Sinc we use `ind` the line being extracted corresponds to the array number
#i.e. for an array of 1-10 the indiv_name will range from the 1st row/line to the 10th
indiv_name=$(cat sample_list.txt | sed -n ${ind}p)

#Then we run fastqc as before for the R1 and R2 reads for that individual
#We also add a touch command at the end to make a file indiciating that run finished
fastqc -t 1 -o /hb/home/USER/fastqc_reports/ /hb/groups/BLAH/${indiv_name}_R1_001.fastq.gz && touch ${indiv_name}_R1.done
fastqc -t 1 -o /hb/home/USER/fastqc_reports/ /hb/groups/BLAH/${indiv_name}_R2_001.fastq.gz && touch ${indiv_name}_R2.done
```


