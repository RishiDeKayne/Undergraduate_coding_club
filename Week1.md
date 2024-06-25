

File types
.txt
Unformatted file
Problems arise when certain tools need input to be structured in a certain way
Differences in txt files based on your operating system
windows to UNIX (dos2unix / unix2dos ) text file converter, converts text files with dos or mac line breaks to unix line breaks and vice versa, does take double quotation marks and turns them into 2 single quotations which can be a really important issue when using dos2unix to create a linux .txt from a windows generated .txt file 
.sh
Means it is a shell script 
Job script made to be submitted 
Look up (for example) google coding best practices - makes it easier for collaborative work

.fastq
Files we get from the sequencing companies that include DNA sequences and the quality score that associates with them
Sequencing is really accurate in the middle of the read but the ends are more prone to issues so the fastq files help you decide if the quality is good
Informative output of the distribution of the quality score
Helps you figure out what could have caused issues; quality score is really consistent, and then in a certain portion of the read, it drops dramatically → could be for multiple reasons, ex is someone bumped the sequencing machine
They come in specific format
.fasta/.fna
All fasta files start with a greater than symbol ‘>’
Header which is the name of the read (single line description) and then followed by the read itself (sequence data)
Typically the file format used for genome assembly 
Reference used to map files against
Text-based format represent nucleotide sequences/peptide sequences in single-letter codes
.sam/.bam
SequenceAlignmentMap/BinaryAlignmentMap
Sam is an uncompressed version of the bam file
Represents aligned sequences
Has a header which describes the entire file and then followed by the alignment which contains all the read information
Can compress it to different extents → takes computer a while to compress & decompress
.vcf
Variant Call format
Probably the most useful file
They are the files you deal with once you call snips
Row —> different snips
Columns –> individual
“./.” → can’t call a variant from an individual 
Other variants we call are going to be 0 and 1 
In really complex data you might even see 0, 1, and 2
Each number represents a different variant/allele
"0/0"-both reference alleles, "0/1" 1 reference allele, 1 SNP
Easiest to filter (may filter them to exclude missing data)

.GFF/.GTF
Nightmare format
Used to represent the annotations within genomes
Have a nested structure which makes them really complicated
Can be broken very easily when you try to manipulate them
Can match a segment of it with SNPs to understand whats going on in there
Counts matrices
"count matrices" or "expression matrices" contain a measure of gene expression for every gene in every sample/individual
Rows → 1 per gene
Column → 1 per individual
How many read sequenced (RNA) correlate with each gene for individual
How many times do you have a gene match to individual
Are things similar or different to each other? 

A lot of these file formats have a standard format they follow, so it is good to look at the manual before approaching the task

Mapping:
Mapping takes the short read that we cut up and compares them to the reference

Zipping files/tar
Zipping files is simply compressing them
Can do this by converting them to some type of binary format
Some formats let us compress things to different extents 
Ex: bam files
Problem is that it takes a long time to compress them and a longer time to uncompress them
Problem with compress files is that you can’t just read them, have to use something to read them or else your screen will look like a mess
Gzip 
Zcat → pipe it into less (otherwise it will be chaotic) 


HELP TIPS:
When in hummingbird/terminal/all that → after writing the first letter or first few letters, you can press tab and it will fill in the rest of the file name
When in file:
Arrow key to go up and down
If the file is really long/big, can use the space key to go to the next page
In less, you can do -S and it will organize the file and allows for horizontal scrolling
GT:PL:DP:GQ → Genotype : ? : depth : genotype quality → ex: 0/0:0,54,127:18:63

SPECIAL WARNING: DO NOT USE rm -r* at any point
-it will delete every file and directory at the level it is ran at and everything below directory wise. The -r option means delete even directories. By itself rm without the option -r does not delete directories and kicks back an error if you try to use it to remove a directory. That is why the rm -r is so dangerous, the addition of the * is a variable that means everything before*after is included in the group.  Ex. head star* will display every file that starts with star no matter what characters follow with the head command.  A * without any additional modifiers includes everything.  So again to reiterate a rm -r* selects EVERYTHING to be DELETED
