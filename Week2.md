Undergraduate Coding Club - Kelley Lab 2024
================
Week 2 - Hummingbird tips
================

This week we will be discussing individual workflows for scientific projects. Additionally, we will cover some things to consider when working on the Hummingbird cluster.

## Before running jobs
- *Slurm Basics*: Familiarize yourself with Slurm concepts such as partitions (queues), jobs, job scripts, nodes, and resources (CPU, memory, GPUs).
- *Read Documentation*: This applies to both Hummingbird itself ([HERE](https://hummingbird.ucsc.edu/) you can find info about Hummingbird) and the queue system Slurm ([HERE](https://www.carc.usc.edu/user-information/user-guides/hpc-basics/slurm-cheatsheet) you can find a cheat-sheet for Slurm) which has comprehensive documentation available online.
- *Project File and Folder Organization*:
  - E.g. `data`, `analysis`, `output`
    - `data` should have raw data and should remain untouched/unaltered (could contain `.fastq`, `.bam`, `.vcf` files e.g.)
    - `analysis` is where the bulk of your work should be carried out taking files from Data and doing something to them
    - `output` is where we can send/copy useful output files, figures, summary stats which can then be shared and manipulated without the risk of altering the raw data or analysis
  - Within analysis you can then have folders for each part e.g. `01_QC`, `02_Mapping`, `03_Genotyping`.
  - Within those folders you can then have different folders for different parameter combinations e.g. `FST_50kb_windows`, `FST_100kb_windows`

## While running jobs
- *Job Submission*: Always submit jobs using `sbatch` rather than running them interactively (`srun`), especially for longer tasks. This ensures better resource management and job tracking.
- *Software Environment*: Load required modules (`module load`) for your software packages and tools at the beginning of your job scripts. Ensure your environment is set up correctly.
    - Conda
    - Already installed modules
- *Job Scripts*: Write clear and well-commented job submission scripts (`.sh` files). Include necessary directives (`#SBATCH`) for resource requests, job name, output/error logs, etc.
- *Resource Requests*: Specify resource requirements accurately in your job scripts (`#SBATCH --mem`, `#SBATCH --cpus-per-task`) based on the needs of your analysis tools. You can also limit the number of jobs run at the same time even within a bigger batch array using `#SBATCH --array=[1-18]%4` - in this example the `%4` specifies the number of jobs running at once.
- *Partition Selection*: Choose the appropriate partition (`-p` or `--partition`) for your jobs. Different partitions may have different resource limits or priorities.
- *Job Dependencies*: Use job dependencies (`--dependency`) when necessary to ensure jobs run in sequence or based on completion of other jobs.
- *Error Handling*: Redirect standard output and standard error to log files (`#SBATCH --output`, `#SBATCH --error`) to capture any errors or debugging information.
- *Check Job Status*: Use `squeue` to check the status of your jobs (`squeue -u your_username`). This helps you monitor progress and identify issues.
- *Monitor Resource Usage*: Use `sacct` or `sstat` to monitor resource usage of your jobs after they have finished. This helps optimize resource requests for future jobs.

## After running jobs
- *Clean Up*: Remove temporary files and directories after job completion to avoid cluttering the file system and using unnecessary disk space.
- Learn from Others: Collaborate with eachother and seek advice from more experienced Hb users. 
  - Improving coding skills
  - Learning new/more efficient tools/ways of doing tasks
  - Learning new software/programs
- *Stay Updated with Hb*: Keep yourself updated with cluster announcements, maintenance schedules, and any changes in cluster policies or procedures.
- *Back up*: Stay on top of knowing which files are ‘key’ files that need backing up
