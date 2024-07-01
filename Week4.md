Undergraduate Coding Club - Kelley Lab 2024
================
Week 4
================

This week we will be 

## Data Files

As with last week we are going to make an example phenotype file we will use during this workshop. 

```
#Individual_Number Spawning_Depth Standard_Length Sex Colour
echo 'BL_100 4 35 M Blue' > sample_info.txt
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

## Part 1: AWK

AWK is a domain-specific language designed for file and text processing. For bioinformatics work we will mostly use it for filtering and summarising large text files just like `sed` and `grep` (which we will cover below).  

This week we will start by showing how you can print the whole contents of a file with AWK
```
awk '{print}' sample_info.txt
```

With AWK we can also do this a different way specifying `$0` which means 'all columns'
```
awk '{print $0}' sample_info.txt
```

We can also print only specific columns, for example if we just want the name of the individual (column 1 `Individual_Number`) and the sex of the fish (column 4 `Sex`)
```
awk '{print $1, $4}' sample_info.txt
```

Here's where things get useful for summarising - let's say we want to calculate the average for a given column - Spawning depth (column 2). We can also do this with AWK
```
awk '{sum += $2} END {print "Average Spawning_Depth:", sum/NR}' sample_info.txt
```

Another useful way we use AWK is to filter columns based on some condition (we covered conditional statements in more detail in week 3 HERE)
If we only want to return female fish we can do the following:
```
awk '$4 == "F" {print}' sample_info.txt
```

As we covered in more detail last week you can imagine this can be useful for subsetting sample lists. If we only wanted to run some code for Female fish in our dataset we could prepare an input sample file with something like this. 
```
awk '$4 == "F" {print}' sample_info.txt | cut -f 1 > female_sample_list.txt
```

But wait!! This doesn't only extract the individual names even though we specified the first column with `cut -f1` WHY?
Note: last week we first converted spaces to tabs - this week we did not so we need to specify the delimiter for `cut` to work. Let's do that again correctly and overwrite our original file.
```
awk '$4 == "F" {print}' sample_info.txt | cut -d ' ' -f 1 > female_sample_list.txt
```

Last week we learned to use `|` to join together commands including things like counting output rows with `wc -l` after a command.  
But we can do things like count the output directly with AWK too. Here we will count the number of red fish in our set.
```
awk '$5 == "Red" {count++} END {print "Number of Red samples:", count}' sample_info.txt
```

Let's say we went to send a collaborator some output and they have asked for it to be in a `.csv` format ready to import to R or python. We can use awk to subset our data, in this case to extract specific columns, as well as specifying what delimiter to put between them. You can see this could be another way to replace our spaces with tabs or commas if we wanted.
```
awk -v OFS=',' '{print $1, $2, $3}' sample_info.txt > data_for_collaborator.txt
cat data_for_collaborator.txt
```

We can also use AWK to specify a range of rows we want to extract. If we wanted only rows 2-4 we can run
```
awk 'NR==2,NR==4 {print}' sample_info.txt
```

As well as `count` we can also `sum` a column. For example the sum of the spawning depth values can be summed as follows
```
awk '{sum += $2} END {print "Total Spawning_Depth:", sum}' sample_info.txt
```

Finally, we can do a kind of find and replace, simlar to that we will cover with `sed` below using AWK. We found out that our blue fish were actually wronly tagged and in fact they were green. We can replace blue with green in column 5 as follows. If we wanted we could also asign this to a new data frame.
```
awk '{$5 = ($5 == "Blue" ? "Green" : $5)} 1' sample_info.txt
```




