#!/bin/bash

# Define the input file and output file
input_file="two_pulse_flexible_prop_2_values.txt"
output_file="mean_ancestry_1Mb_windows.txt"

# Initialize the output file
echo -e "Chromosome\tStart_Position\tEnd_Position\tIndividual_ID\tMean_Ancestry" > $output_file

# Sort the input file by chromosome and position
sorted_file="sorted_${input_file}"
sort -k1,1 -k2,2n $input_file > $sorted_file

# Process each 1Mb window using awk
awk '
BEGIN {
    # Define window size and initialize variables
    window_size = 1000000
    start_position = -1
}

{
    # Read chromosome and position
    chrom = $1
    pos = $2
    
    # Check if the position is within the current window
    if (start_position == -1 || pos >= end_position || chrom != current_chrom) {
        # If not, finalize the current window and start a new window
        if (start_position != -1) {
            # Calculate mean ancestry for each individual in the current window
            for (i = 3; i <= NF; i++) {
                mean = sum[i] / count
                print current_chrom "\t" start_position "\t" end_position "\t" (i - 2) "\t" mean >> "'$output_file'"
            }
        }
        
        # Initialize new window
        current_chrom = chrom
        start_position = pos
        end_position = start_position + window_size
        
        # Initialize summation and count for new window
        for (i = 3; i <= NF; i++) {
            sum[i] = 0
        }
        count = 0
    }
    
    # Accumulate sums and count for the current window
    for (i = 3; i <= NF; i++) {
        sum[i] += $i
    }
    count += 1
}
END {
    # Finalize the last window
    if (start_position != -1) {
        for (i = 3; i <= NF; i++) {
            mean = sum[i] / count
            print current_chrom "\t" start_position "\t" end_position "\t" (i - 2) "\t" mean >> "'$output_file'"
        }
    }
}
' $sorted_file

echo "Mean ancestry calculation in 1Mb windows complete. Results saved to $output_file."
