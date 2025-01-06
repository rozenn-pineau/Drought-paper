# Drought-paper
Suite of notes and scripts for the drought project.






### Ancestry mapping using ancestry hmm

## Processing results from ancestry_hmm output

# Step (1) : gather all samples in one file
Threshold for probability of 0.9

```
#!/bin/bash
#initialize dataset
folder_name=$(basename "$PWD")
output_file="${folder_name}_values.txt"

for file in *.posterior; do

    header_name=$(basename "$file" .posterior)
    
    awk -v header="$header_name" 'NR == 1 { print header; }
    
    NR > 1 {
    
    result = "NA";
    if ($3 > 0.9) result = "0";
    else if ($4 > 0.9) result = "1";
    else if ($5 > 0.9) result = "2";
    print result;
    }' $file > tmp
    
    paste $output_file tmp > tmp2
    
    mv tmp2 $output_file

done

```


