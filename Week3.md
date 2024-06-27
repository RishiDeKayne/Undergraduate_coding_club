Undergraduate Coding Club - Kelley Lab 2024
================
Week 3
================

This week we will be taking a look at some different tips, tricks, and features of bioinformatics code you may come across or need to write yourself. 
These tools can all be used in conjunction one another so I have structured this workshop to go through examples of how to use each tool separately before comining them in some more complicated challenges at the end.

## Data Files

First we are going to make an example phenotype file we will use during this workshop. 

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
If we only wanted specific columns we can specify those within the print statement with `$1` for column 1 etc.
```
awk '$2 < 20 {print $0}' sample_info_tabs.txt
```

The pipe symbol allows us to also further manipulate this output.
For example we can sort this output based on the value of column 2 (using `sort -k 2`) from deepest to shallowest.
```
awk '$2 < 20 {print $0}' sample_info_tabs.txt | sort -k 2
```

We can keep adding commands using pipe. 
Firstly we might want to identify unique values in our sorted set. We can do this with `uniq`.
```
awk '$2 < 20 {print $2}' sample_info_tabs.txt | sort | uniq
```

You can see this way we have only three values, 10, 5, and 4. But what if we were looking at unique gene names and had 10,000+?
We can use `wc -l` to count the number of lines in this output without even having to make an intermediate file. 
```
awk '$2 < 20 {print $2}' sample_info_tabs.txt | sort | uniq | wc -l
```

## Variables
Variables are used to store values for later use in scripts. They can hold user input, filenames, or any other data needed in your script.  

Assign a variable
```
individual="NR_703"
```

Use the variable in a command. Here we use `grep` which allows us to use regular expressions to extract information from files.
In it's simplest form below grep is searching for the contents of our variable within the `.txt` file provided.
```
grep "$individual" sample_info_tabs.txt
```

We can also do this for a whole group. This time we want to extract all blue fish
```
color='Blue'
grep "$color" sample_info_tabs.txt
```

As before we can build on this to get other features. This time we will look for the number of red fish in the sample set
```
color='Red'
grep "$color" sample_info_tabs.txt | wc -l
```

We can also use variables to store a specific file name - this is something you might do if you are carrying out an action on a specific file type.
For example you may decide to carry out some filtering on a `.vcf` file - you want to make a filtering script once and then you can use it for any vcf file by specifying your desired input and output file names as variables in the script. Here we will produce subsetted versions of our `sample_info_tabs.txt` file based on the color of the fish to extract only the color and the standard length.

```
color1='Blue'
color2='Red'
input_file="sample_info_tabs.txt"
output_file_blue="Blue_SL.txt"
output_file_red="Red_SL.txt"

grep "$color1" $input_file | awk '{print $1,$3}' > $output_file_blue && grep "$color2" $input_file | awk '{print $1,$3}' > $output_file_red

#look in each file with cat
cat $output_file_blue
cat $output_file_red

```

## Control Structures and Conditional Statements (if, elif, else, fi)

In bash scripting, `if`, `elif`, `else`, and `fi` are control structures used for conditional execution of commands. They allow you to control the flow of your script based on conditions, enabling you to make decisions and execute different commands or actions accordingly.  

An `if` statement starts a block of code that will be executed only if a specified condition is true.

```
if [ condition ]; then
    # Commands to execute if condition is true
fi
```

- `[ condition ]` is a conditional expression enclosed within square brackets. 
- We can specify a bunch of different conditions including (-eq, -ne, -lt, -gt, -le, -ge), file checks (-f, -d, -r, -w, -x), and logical operations (&&, ||). In your own time look up what some of these mean and think about how you could use them in a script
- `then` marks the beginning of the block of code to execute **if** the condition is true.
- `fi` (reverse of `if`) marks the end of the if block.

An `elif` statement allows you to check additional conditions if the preceding if condition (or any previous elif condition) was false.

```
if [ condition1 ]; then
    # Commands to execute if condition1 is true
elif [ condition2 ]; then
    # Commands to execute if condition2 is true
else
    # Commands to execute if all conditions are false
fi
```

- You can have multiple `elif` blocks to check multiple conditions sequentially.
- The `else` block (if it is present) will execute if none of the preceding conditions (`if` or `elif`) were true.
- As before `fi` ends the entire conditional structure.

Here is an example using a few of the things we have learned today:
```
echo $individual
trait_value=$(grep "$individual" $input_file | cut -f 2)

if [ $trait_value -gt 10 ]; then
    echo "$individual has spawning depth greater than 10."
elif [ $trait_value -eq 10 ]; then
    echo "$individual has a spawning depth equal to 10."
else
    echo "$individual has a spawning depth less than 10."
fi

#verify the output
echo "$individual has a spawning depth of: " && grep "$individual" $input_file | cut -f 2
```
