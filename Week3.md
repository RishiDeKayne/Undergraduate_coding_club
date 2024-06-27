Undergraduate Coding Club - Kelley Lab 2024
================
Week 3
================

This week we will be taking a look at some different tips, tricks, and features of bioinformatics code you may come across or need to write yourself. 
These tools can all be used in conjunction one another so I have structured this workshop to go through examples of how to use each tool separately before comining them in some more complicated challenges at the end.

## Data Files

First we are going to make an example phenotype file we will use during this workshop. 

```
echo 'Individual_Number Spawning_Depth Standard_Length Sex Colour' > sample_info.txt
echo 'BL_100 4 35 M Blue' >> sample_info.txt
echo 'BL_101 5 32 F Blue' >> sample_info.txt
echo 'BL_102 4 36 M Blue' >> sample_info.txt
echo 'BL_400 4 34 M Blue' >> sample_info.txt
echo 'BL_432 10 41 F Blue' >> sample_info.txt
echo 'NR_700 24 18 F Red' >> sample_info.txt
echo 'NR_703 31 16 F Red' >> sample_info.txt
echo 'NR_704 33 19 F Red' >> sample_info.txt
echo 'NR_870 32 17 M Red' >> sample_info.txt
echo 'NR_871 31 18 F Red' >> sample_info.txt
```

Now convert spaces to tabs so that we can use awk easily (you can specify the delimiter but to keep it simple here we will keep it at tabs)
```
sed 's/ /\t/g' sample_info.txt > sample_info_tabs.txt
```

## Piping (| and >)
Piping (|) allows you to redirect the output of one command as input to another command, while redirection (>) is used to redirect output to a file.  

Count number of organisms
```
wc -l sample_info_tabs.txt
```

Pipe output to cut to select only the organism names
```
cut -f 1 sample_info_tabs.txt
```

Redirect output to a new file
```
cut -f 1 sample_info_tabs.txt > individuals.txt
```

Now we can make a version of this file for individuals with a spawning depth < 20m
We use `awk` (which we will cover in more detail next week) to filter based on column 2 (indicated by `$2`) and then return all columns (`{print $0}`)  
```
awk '$2 < 20 {print $0}' sample_info_tabs.txt
```

The pipe symbol allows us to also further manipulate this output.
For example we can sort this output baed on the value of column 2 (using `sort -k 2`) from deepest to shallowest.
```
awk '$2 < 20 {print $0}' sample_info_tabs.txt | sort -k 2
```

## Variables
Variables are used to store values for later use in scripts. They can hold user input, filenames, or any other data needed in your script.  

Assign a variable
```
organism="Organism1"
```

Use the variable in a command
```
grep "$organism" phenotypes.txt
```

## Control Structures and Conditional Statements (if, elif, else, fi)
Control structures allow you to make decisions based on conditions within your script.  

```
#!/bin/bash

organism="Organism3"
trait_value=$(grep "$organism" phenotypes.txt | cut -f 2)

if [ $trait_value -gt 10 ]; then
    echo "$organism has a trait value greater than 10."
elif [ $trait_value -eq 10 ]; then
    echo "$organism has a trait value equal to 10."
else
    echo "$organism has a trait value less than 10."
fi
```
