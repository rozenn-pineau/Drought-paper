#!/bin/bash
#SBATCH --job-name=vcf_plink
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=20:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50GB

cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/

# Input and output files
input_file="two_pulse_flexible_prop_2_values.txt"
output_file="two_pulse_flexible_prop_2_values.vcf"

# Create the VCF header
cat <<EOL > $output_file
##fileformat=VCFv4.2
##source=CustomScript
##INFO=<ID=.,Number=1,Type=String,Description="Custom info">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOL
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$(head -1 $input_file | cut -f3- | tr -s ' ' '\t')" >> $output_file


# Transform the input data into VCF format
tail -n +2 $input_file | while read -r line; do
    # Parse the fields
    chrom=$(echo "$line" | awk '{print $1}')
    pos=$(echo "$line" | awk '{print $2}')
    genotypes=$(echo "$line" | cut -f3-)

    # Convert genotypes to VCF GT format
    formatted_genotypes=$(echo "$genotypes" | awk '
    {
        for (i = 1; i <= NF; i++) {
            if ($i == 0) $i = "0|0";      # Homozygous reference
            else if ($i == 1) $i = "1|0";      # Heterozygote
            else if ($i == 2) $i = "1|1"; # Homozygous alternate
            else $i = "./.";              # Missing data or unrecognized value
        }
        print $0
    }' | sed 's/ /\t/g')

    # Placeholder values for REF, ALT, ID, QUAL, FILTER, INFO
    ref="A" # Adjust this to your data's REF allele
    alt="T" # Adjust this to your data's ALT allele
    id="." # No variant ID provided
    qual="."
    filter="PASS"
    info="."
    format="GT"

    # Print the VCF record
    echo -e "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$formatted_genotypes" >> $output_file
done

