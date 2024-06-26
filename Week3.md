Undergraduate Coding Club - Kelley Lab 2024
================
Week 3
================

This week we will be taking a look at some different tips, tricks, and features of bioinformatics code you may come across or need to write yourself. 
These tools can all be used in conjunction one another so I have structured this workshop to go through examples of how to use each tool separately before comining them in some more complicated challenges at the end.

## Data Files

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
sed -i 's/ /\t/g' sample_info.txt
```
## Piping (| and >)


## Variables

```
test
```

## Control Structures and Conditional Statements (if, elif, else, fi)

