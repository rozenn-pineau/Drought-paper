# Drought-paper
Suite of notes and scripts for the drought project.






# Ancestry mapping using ancestry hmm

## Processing results from ancestry_hmm output

Ancestry_hmm gives as an output one file per sample, with one line per position and three columns, corresponding to the probability for each site to come from ancestry 1, both or ancestry 2. The three columns add up to 1. 
We choose a threshold of 0.9 to assign ancestry at each site using the following script. 


### Step (1) : gather all samples in one file

```
#!/bin/bash

folder_name=$(basename "$PWD")
output_file="${folder_name}_values.txt"

#initialize dataset with chrom and pos
awk 'BEGIN { OFS="\t" } {print $1,$2}' P11_Ag_11_T.posterior > $output_file

#loop through each individual and paste information to previous version of the file
for file in *.posterior; do

    header_name=$(basename "$file" .posterior)
    
    awk -v header="$header_name" 'NR == 1 { print header; }
    
    NR > 1 {
    
    result = "NA";
    if ($3 > 0.9) result = "0";
    else if ($4 > 0.9) result = "1";
    else if ($5 > 0.9) result = "1";
    else if ($6 > 0.9) result = "2";
    else if ($7 > 0.9) result = "1";
    else if ($8 > 0.9) result = "2";
    print result;
    }' $file > tmp
    
    paste $output_file tmp > tmp2
    
    mv tmp2 $output_file

done

#check: file is 786262 lines and 284 columns (chromosome, position and 282 individuals). 

```

### Visualize : plot ancestry by chromosome in R

