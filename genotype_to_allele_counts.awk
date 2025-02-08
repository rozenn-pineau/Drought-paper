#!/bin/awk -f

#this script reads in a vcf file and outputs a txt file with the chromosome, position, number of allele a and number of allele A
BEGIN { FS=OFS="\t" }

# Skip headers
/^#/ {
    next
}

# Process genotype data
{
    # Initialize total counts for the site
    total_ref_count = 0
    total_alt_count = 0

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

        # Sum counts for this site
        total_ref_count += ref_count
        total_alt_count += alt_count
    }

    # Output the position (CHROM, POS) and total counts
    printf "%s\t%s\t%d\t%d\n", $1, $2, total_ref_count, total_alt_count
}
