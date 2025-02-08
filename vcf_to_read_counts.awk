#!/bin/awk -f

#This script reads in a vcf and outputs the number of alleles of each ancestry. 2,0 would be homz reference / 1,1 --> hets / 0,2 hets alternative
#two columns per sample 

BEGIN { FS=OFS="\t" }

/^#CHROM/ {
    # Print header for output
#    printf "CHROM\tPOS"
#    for (i = 10; i <= NF; i++) {
#        printf "\t%s_ref\t%s_alt", $i, $i
#    }
#    print ""
    next
}

# Skip other header lines
/^#/ {
    next
}

# Process genotype data
{
    # Output the position (CHROM, POS)
#    printf "%s\t%s", $1, $2

    # Process genotype columns (starting from the 10th field)
    for (i = 10; i <= NF; i++) {
        split($i, fields, ":")
        split(fields[1], genotype, /[\/|]/)

        # Count ref and alt alleles
        ref_count = 0
        alt_count = 0

        for (j = 1; j <= length(genotype); j++) {
            if (genotype[j] == "0") {
                ref_count++
            } else if (genotype[j] == "1") {
                alt_count++
            }
        }

        # Output the counts for this sample
        printf "%d\t%d\t", ref_count, alt_count
    }

    print ""
}
