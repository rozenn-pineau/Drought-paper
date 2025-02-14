# Drought-paper
Suite of notes and scripts for the drought project.


# Outline

[Drought dataset : ancestry mapping using ancestry hmm](#Drought-dataset---ancestry-mapping-using-ancestry-hmm)

[Processing results from ancestry_hmm output](#Processing-results-from-ancestry_hmm-output)

[GWAS on ancestry calls](#GWAS-on-ancestry-calls)

[Drought selection experiment trajectory analyses](#Drought-selection-experiment-trajectory-analyses)

[Herbarium dataset: ancestry mapping using ancestry hmm](#Herbarium-dataset---ancestry-mapping-using-ancestry-hmm)

        [Calling ancestry on drought-informative sites only](#Calling-ancestry-on-drought-informative-sites-only)
    
        [Calling ancestry on more sites](#Calling-ancestry-on-more-sites)



# Drought dataset - ancestry mapping using ancestry hmm


[Ancestry_hmm](https://github.com/russcd/Ancestry_HMM) : tool to infer ancestry at input positions in the genome (Corbett-Detig, R. and Nielsen, R., 2017.)


### Step (1) : define var rudis versus and var tuberculatus pure ancestry individuals
To define the "ancestry panel" and not lose any of the drought data, we used an additional dataset, variants from common garden experiments used in Kreiner et al, 2022 that were realigned on a newer version of the reference genome. We ran structure with k=2, and kept individuals with pure ancestry with a threshold of 0.00001. We identified 44 pure var. rudis and 21 pure var. tuberculatus samples.

```
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

my_bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/commongarden_allfiltsnps_193_hap2_numericChr_filt.bed

for K in 2;

do admixture --cv $my_bed $K ;

done
```

### Step (2) : calculate Fst on ancestry sites with MAF <0.05
We filtered the pure ancestry variant file (50,811,811 mutations) for minimum allele frequency of 0.05. To focus the analysis on ancestry-informative sites, we kept the sites with Fst values situated in the 75th percentile and above (2,089,620 variants, Fst calculated using VCFtools version 0.1.16).


```
#calculate fst per site with magf > 0.05
module load vcftools

VCF=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/ancestry.vcf.recode.vcf
vcftools --vcf ${VCF} --maf 0.05 --recode --stdout > fst/ancestry_maf.vcf

vcftools --vcf fst/ancestry_maf.vcf \
--weir-fst-pop var_rudis_samp.list \
--weir-fst-pop var_tub_samp.list \
--out fst/ancestry_maf

```

Rscript to choose Fst threshold : [fst_on_ancestry.Rmd](https://github.com/rozenn-pineau/Drought-paper/blob/main/fst_on_ancestry.Rmd)

### Step (3) : keep the variants common to the ancestry variant file and the drought variant file
The ancestry was filtered for the high Fst sites. This filtered file was used to keep the intersection between the ancestry variant file and the drought variant file (bcftools isec).
At this step, we have two variant files with the same variants for the two populations.
Now, we need to generate the rho information between each site.


### Step (4) : calculate rho between each site
Based on an estimation of ld for (nearly) the whole genome, we calculated, chromosome per chromosome, the linkage value between two consecutive sites. To do this, we fit a monotonic spline on the cumulative distribution of the mean recombination rates. 

Important note: we had to exclude the sites that were outside of the region defined by the recombination map (the monotonic spline does not extrapolate outside of boundaries). Then, we kept track of which sites to filter out of the variant files. 

Rscript to calculate rho between sites : [calculate_ldhat_between_sites.Rmd](https://github.com/rozenn-pineau/Drought-paper/blob/main/calculate_ldhat_between_sites.Rmd).



### Step (5) : extract allele counts from the ancestry variant file

[genotype_to_allele_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/genotype_to_allele_counts.awk)

exanple output :
22 22 for 22 individuals, all heterozygotes
00 44 for 22 individuals, all homs alternative

### Step (6) : get genotypes for drought panel

2,0 for homs reference
1,1 for hets
0,2 for homs alternative

[vcf_to_read_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/vcf_to_read_counts.awk)

### Step (7) : put the file together

Paste columns together to make the full file, following instructions on ancestry_hmm github page. 







# Processing results from ancestry_hmm output

Ancestry_hmm gives as an output one file per sample, that has the following header :
```
chrom	position	2,0,0	1,1,0	1,0,1	0,2,0	0,1,1	0,0,2
```
"2,0,0" column has the probability that the loci was homozygote for var. tuberculatus before hybridization, and stayed homozygote for var. tuberculatus at pulse 2 and 3 (0 chromosomes from var. rudis). GT = 0

"1,1,0" column has the probability that the loci was heterozygote for var. tuberculatus before hybridization, and became heterozygote for var. rudis at pulse 2 (1 chromosome from var. rudis). GT = 1

"1,0,1" column has the probability that the loci was heterozygote for var. tuberculatus before hybridization, and became heterozygote for var. rudis at pulse 3 (1 chromosome from var. rudis). GT = 1

"0,2,0" column has the probability that the loci was homozygote for var. rudis before hybridization, and stayed homozygote for var. rudis (0 chromosome from var. tuberculatus). GT = 2

"0,1,1" column has the probability that the loci was homozygote for var. tuberculatus before hybridization, and became gain one var. rudis loci at each pulse (1 chromosome from var. rudis). GT = 2

"0,0,2" column has the probability that the loci was homozygote for var. rudis before hybridization, and became homzygote for var. rudis at pulse 3 (1 chromosome from var. rudis). GT = 2


The probabilitiees in the columns add up to 1. 

We choose a threshold of 0.8 to assign ancestry at each site using the following script. 


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
    if ($3 > 0.8) result = "0";
    else if ($4 > 0.8) result = "1";
    else if ($5 > 0.8) result = "1";
    else if ($6 > 0.8) result = "2";
    else if ($7 > 0.8) result = "2";
    else if ($8 > 0.8) result = "2";
    print result;
    }' $file > tmp
    
    paste $output_file tmp > tmp2
    
    mv tmp2 $output_file

done

#check: file is 786262 lines and 284 columns (chromosome, position and 282 individuals).

```

### Visualize : plot ancestry by chromosome in R

We used the following script [plot_ancestry.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/plot_ancestry.R) to visualize ancestry per chromosome for each smaple in the dataset, ordered by longitude. 

To make sure independent runs gave consistent ancestry calls, I compared the output from two runs that had different starting proportions (-.3 and -.37 versus  -.5 and -.17) : [compare_ancestry_hmm_runs.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/compare_ancestry_hmm_runs.R). I tested 1000 random sites on two chromosomes and dit not identify differences in ancestry calls between models, suggesting the robustness of the results. 


# GWAS on ancestry calls

## (1) Using a "manual" logistic model with and without covariate

We ran logistic regression model in R based on the ancestry calls and using PC1 as a covariate using the following R scripts : without covariate - [gwas_ancestry_logistic_v2.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/gwas_ancestry_logistic_v2.R) and with covariate - [gwas_ancestry_logistic_covariate_v2.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/gwas_ancestry_logistic_covariate_v2.R).
These scripts require the sample order vector ("samp_order_noout.txt"), the phenotypes ("drought_phenos.txt") and the PC information ("PCA_output_noout.txt").


## (2) using GEMMA with and without relatedness matrix

### Prepping the data for GWAS

To use the ancestry calls with gemma, I first had to reformat it to bimbam format :

```

#!/bin/bash


cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs
file=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs/two_pulse_flexible_prop_2_values_cutoff_08_v2.txt
output=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs/two_pulse_flexible_prop_2_values_cutoff_08_v2.gemma
#remove first line
head -n +1 $file > samp.order

tail -n -786261 $file > tmp

#remove outlier columns
awk -F'\t' '{for(i=1;i<=NF;i++) {if($i == "P16_Nat_1_T") printf("Column %d is outlier #1\n", i-1)}}' samp.order
awk -F'\t' '{for(i=1;i<=NF;i++) {if($i == "P12_Nat_14_T") printf("Column %d is outlier #2\n", i-1)}}' samp.order
#columns 29 and 104
cut -d$'\t' --complement -f 29,104 tmp > tmp1

#add fake ref to alt allele columns in positions 3 and 4
awk 'BEGIN { FS = OFS = "\t" } {
    for (i = 1; i <= NF; i++) {
        if (i == 3) printf "A\tT\t"; # Add "A" and "T" before column 3
        printf "%s%s", $i, (i < NF ? OFS : ""); # Print original column with tab separator
    }
    print ""; # Add newline after each row
}' tmp1 > tmp2

#replace first tab by ":", replace every remaining tab by commas, and delete column 2
awk -F'\t' '{$1 = $1 ":" $2; print}' tmp2 > tmp3
awk 'BEGIN { FS=" "; OFS="," } {$1=$1; print}' tmp3 > tmp4
cut -d ',' --complement -f2 tmp4 > $output


rm -I tmp*

```

Running GEMMA (we allowed up to 20% missing data per sample for each SNP):

```
#!/bin/bash
#SBATCH --job-name=gemma
#SBATCH --output=gemma.out
#SBATCH --error=gemma.err
#SBATCH --time=25:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100G   # memory per cpu-core

cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/gemma_gwas/two_pulse_flexible_prop_2

anc_file=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs/two_pulse_flexible_prop_2_values_cutoff_08_v2.gemma
phenos=/scratch/midway3/rozennpineau/drought/ancestry_hmm/manual_gwas/phenos_ordered_noout.txt
geno_relatedness=/scratch/midway3/rozennpineau/drought/ancestry_hmm/gemma_gwas/relatedness_matrix_geno_ordered.txt

module load gemma
gemma -g $anc_file -p $phenos -miss 0.1 -lm #gwas not corrected for population structure
#-miss 0.1 is to accept up to 10% missing data per sample per SNP, will be interpolated from existing data

#gemma gwas corrected for pop structure
##with genotype-based relatedness matrix
gemma -g $anc_file -p $phenos -k $geno_relatedness -km 1 -miss 0.1 -lmm -o geno_corrected_gemma_gwas

##with anc covariate, first generating relatedness matrix (-gk)

gemma -g $anc_file -p $phenos -gk 1 -o ancestryrelatedness #generate relatedness
gemma -g $anc_file  -p $phenos -k output/ancestryrelatedness.cXX.txt -miss 0.1 -lmm -o ancestry_corrected_gemma_gwas

```

For the genotype-based relatedness matrix, a bit of sample reordering was necessary and done in R :

```
#goal: to re order the relatedness matrix built from genotypes such that it matches the ancestry call order
rm(list= ls())

#load relatedness matrix 
mat <- read.table("/project/kreiner/drought/gwas_episode3/output/merged_numericChr_nooutliers_maf01.cXX.txt", sep = "\t", header = F)#280 x 280

#load order of samples in rm
order_mat <- read.table("/project/kreiner/drought/gwas_episode3/prep_files/merged_numericChr_nooutliers_maf01.fam", sep = " ", header=F)

#load order of samples in ancestry calls
order_goal <- read.table("/scratch/midway3/rozennpineau/drought/ancestry_hmm/manual_gwas/samp_order_noout.txt", sep = "\t", header = F)

#order matrix based on order of samples in order_goal
#template : vector1[order(match(vector1,vector2))] 
new_order <- order(match(order_mat[,1], order_goal[,1]))
new_mat <- mat[new_order, new_order]

#export
write.table(new_mat, "/scratch/midway3/rozennpineau/drought/ancestry_hmm/gemma_gwas/relatedness_matrix_geno_ordered.txt", sep = "\t", col.names = F, row.names = F)
```
### Analyzing associated sites in R and correcting p-values

We analyzed the runs from different probability cutoffs (0.5 to 0.9) and different NA content tolerance in GEMMA command line (0.05 -default-, 0.1 and 0.2). 
We compared QQplots and Manhanttan plots from GWAS with and without covariate (ancestry and genotype-based relatedness matrices). 
---> work in progress R script to be updated : [plot_gwas_ancestry.R)](https://github.com/rozenn-pineau/Drought-paper/blob/main/plot_gwas_ancestry.R)

We obtained a promising GWAS when looking at ancestry-corrected gwas.

### Clumping based on LD
Next step is to clump together the sites that may be in LD. We do this using plink.

First, we export the **inflated** and **FDR corrected** p-values from the association file (see above script). We will use these inflated p values to clump sites based on LD (and significance).




(1) make vcf based on table
```
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

echo "VCF file created: $output_file"


#fill in the ID field by combining chrom and pos
awk 'NR <= 5 {print; next} {OFS="\t"; $3 = $1 ":" $2; print}' two_pulse_flexible_prop_2_values.vcf > two_pulse_flexible_prop_2_values_ID.vcf
```

(2) create plink family files (.ped and .map)
```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2
#make plink family files and create ID based on position information
plink --vcf two_pulse_flexible_prop_2_values_ID.vcf --out two_pulse_flexible_prop_2_values --allow-extra-chr --recode --double-id 
```

(3) clump sites using plink

```
#add ID column to association file
sed -e 's/rs/ID/g' ancestry_corrected_inflated_gemma_gwas.assoc.txt > ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt        

#run plink clump
plink --file two_pulse_flexible_prop_2_values --clump ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt --clump-p1 0.05 --clump-field FDR --clump-kb 100 --out two_pulse_flexible_prop_2_clumped_100kb --allow-no-sex --allow-extra-chr --clump-snp-field ID


#bfile is 0.9 cutoff for ancestry calls
#association file from gemma gwas with --miss 0.2
```

Output from plink --clump :


```
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --clump ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt
  --clump-field FDR
  --clump-kb 100
  --clump-p1 0.05
  --clump-snp-field ID
  --file two_pulse_flexible_prop_2_values
  --out two_pulse_flexible_prop_2_clumped_100kb

192953 MB RAM detected; reserving 96476 MB for main workspace.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (786261 variants, 282 people).
--file: two_pulse_flexible_prop_2_clumped_100kb-temporary.bed +
two_pulse_flexible_prop_2_clumped_100kb-temporary.bim +
two_pulse_flexible_prop_2_clumped_100kb-temporary.fam written.
786261 variants loaded from .bim file.
282 people (0 males, 0 females, 282 ambiguous) loaded from .fam.
Ambiguous sex IDs written to two_pulse_flexible_prop_2_clumped_100kb.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 282 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.889586.
786261 variants and 282 people pass filters and QC.
Note: No phenotypes present.
--clump: 36 clumps formed from 893 top variants.
Results written to two_pulse_flexible_prop_2_clumped_100kb.clumped


```


Prep the files to analyze allelic trajectories :

(1) make bed file to filter ancestry calls vcf file

```
#!/bin/bash

file=two_pulse_flexible_prop_2_clumped_100kb.clumped

awk -F ' ' '{ print $1, $4, $5}' $file > temp.bed #extract cols 1 and 4 for position information, 5 for p value

awk '{print ($2 + 1) " " $4 }' temp.bed > pos.bed #remove 1 from position in file

paste temp.bed pos.bed | awk -v OFS='\t' '{print $1, $2, $4, $3}' > full.bed

tail -n +2 full.bed > full2.bed #remove first line

head -n 36 full2.bed > ancestry_gwas_filtered_sites.bed

```

(2) filter the vcf file based on the bed file

```
module load vcftools
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2
vcftools --vcf two_pulse_flexible_prop_2_values_ID.vcf --bed ancestry_gwas_filtered_sites.bed --out two_pulse_flexible_prop_2_values_ID_filtered --recode
```

output 

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf two_pulse_flexible_prop_2_values_ID.vcf
        --out two_pulse_flexible_prop_2_values_ID_filtered.vcf
        --recode
        --bed ancestry_gwas_filtered_sites.bed

After filtering, kept 282 out of 282 Individuals
Outputting VCF file...
        Read 36 BED file entries.
After filtering, kept 35 out of a possible 786261 Sites
Run Time = 4.00 seconds
````
Why did we lose one site here ?


Getting the ancestry calls (in the form of genotypes) from the filtered vcf for more downstream analyses :
```
bgzip -f two_pulse_flexible_prop_2_values_ID_filtered.recode.vcf

tabix -f two_pulse_flexible_prop_2_values_ID_filtered.recode.vcf.gz

two_pulse_flexible_prop_2]$ bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' two_pulse_flexible_prop_2_values_ID_filtered.recode.vcf.gz > two_pulse_flexible_prop_2_values_ID_filtered_GT.txt
```


## Checking : Does ancestry predict response to drought ?

Our expectation is that var. rudis ancestry is better adapted to drought than var. tuberculatus. 
Do we see this in our ancestry calls ? This is also simply a way to make sure that our pipeline was coded correctly. 



# Drought selection experiment trajectory analyses

## upload scripts for trajectory analyses

### Compare to randomized distribution
We calculate the allele frequency change at each site using the following script : [calculate_anc_calls_fq.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/calculate_anc_calls_fq.R).



# Herbarium dataset - ancestry mapping using ancestry hmm

## Calling ancestry on drought-informative sites only


### Step (1) : filter herbarium vcf files for drought adaptive sites

<ins>(1) merge * herbarium * vcf for each chromosome together in one file </ins>

```
cd /cds3/kreiner/herbarium
bcftools concat herb108_193_2_Scaffold_10_filteredsnps.vcf.gz herb108_193_2_Scaffold_16_filteredsnps.vcf.gz herb108_193_2_Scaffold_6_filteredsnps.vcf.gz herb108_193_2_Scaffold_11_filteredsnps.vcf.gz herb108_193_2_Scaffold_1_filteredsnps.vcf.gz herb108_193_2_Scaffold_7_filteredsnps.vcf.gz herb108_193_2_Scaffold_12_filteredsnps.vcf.gz herb108_193_2_Scaffold_2_filteredsnps.vcf.gz herb108_193_2_Scaffold_8_filteredsnps.vcf.gz herb108_193_2_Scaffold_13_filteredsnps.vcf.gz herb108_193_2_Scaffold_3_filteredsnps.vcf.gz herb108_193_2_Scaffold_9_filteredsnps.vcf.gz herb108_193_2_Scaffold_14_filteredsnps.vcf.gz herb108_193_2_Scaffold_4_filteredsnps.vcf.gz herb108_193_2_Scaffold_15_filteredsnps.vcf.gz herb108_193_2_Scaffold_5_filteredsnps.vcf.gz > merged.vcf

```

<ins>(2) change to numeric chromosome names</ins>

```
#numeric chr names
cat merged.vcf | \sed s/^Scaffold_//g  > merged_numericChr.vcf
```

<ins>(3) make bed file for the 35 FDR clumped significant drought sites AND the bed file from the inflated FDR significant sites from the GWAS, * before clumping * (893 sites)</ins>

```
awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2}' two_pulse_flexible_prop_2_values_ID_filtered_GT.txt > two_pulse_flexible_prop_2_values_ID_filtered_GT.bed

#the other bed file was made from the R script and sent to the cluster
#rozennpineau@midway3.rcc.uchicago.edu:/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/FDR_non_clumped_significant_sites.bed
```

<ins>(4) extract lines from vcf based on bed file</ins>

```
cd /cds3/kreiner/herbarium/
module load vcftools

#35 FDR clumped sites
vcf=merged_numericChr.vcf
bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/two_pulse_flexible_prop_2_values_ID_filtered_GT.bed

vcftools --vcf $vcf --positions $bed --recode --stdout > /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/herb_35FDR_clumped.vcf

#893 FDR non clumped sites
bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/FDR_non_clumped_significant_sites.bed
vcftools --vcf $vcf --positions $bed --recode --stdout > /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/herb_893FDR_non_clumped.vcf

```

I found 1 site in common between the 35 FDR clumped significant sites, and 92 with the 893 before clumping FDR sites.

I will work with the set of 92 and clump based on LD using plink, as previously done for the drought vcf. 

<ins>(5) clump</ins>


The association file has more sites than the filtered vcf, filter the association file for the specific sites only:

```
#!/bin/bash

cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium

#make bed from vcf file
vcf=herb_893FDR_non_clumped.vcf 
bed=herb_893FDR_non_clumped_plus1.bed

awk 'BEGIN {OFS="\t"} 
     !/^#/ {print $1, $2, $2+1}' "$vcf" > "$bed" #position in vcf is shifted by one

# Input files
bed=herb_893FDR_non_clumped_plus1.bed
assoc=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt
out=FDR_non_clumped_significant_sites_filtered.assoc.txt

# Define column positions for chrom and pos in the TXT file (1-based index)
CHROM_COL=13  
POS_COL=14    


# Convert column positions to AWK's 1-based indexing
awk -v chrom_col="$CHROM_COL" -v pos_col="$POS_COL" '
    NR==FNR {sites[$1][$3]; next}  
    { 
        chrom = $chrom_col; 
        pos = $pos_col;
        if (chrom in sites) 
            for (s in sites[chrom]) 
                if (pos >= s && pos <= s) 
                    print $0;
    }' "$bed" "$assoc" > "$out"

echo "Filtering complete. Results saved in $out"

```

The association file has 92 sites. However, they are not in the same order as the vcf file. Let's reorder them both. 


```
#Sort the vcf file :
grep "^#" herb_893FDR_non_clumped.vcf > herb_893FDR_non_clumped_sorted.vcf
grep -v "^#" herb_893FDR_non_clumped.vcf | sort -k1,1V -k2,2g >> herb_893FDR_non_clumped_sorted.vcf

#check
grep -v "#" herb_893FDR_non_clumped_sorted.vcf | awk '{print $1,$2}' #looks good

#positions in vcf are off by 1, fix :

grep "^#" herb_893FDR_non_clumped_sorted.vcf > herb_893FDR_non_clumped_sorted_POSfixed.vcf
grep -v "^#" herb_893FDR_non_clumped_sorted.vcf | awk  '{OFS="\t"; $2 = $2+1; print}' >> herb_893FDR_non_clumped_sorted_POSfixed.vcf

#check
grep -v "#" herb_893FDR_non_clumped_sorted_POSfixed.vcf | awk '{print $1,$2}' #looks good

#add the ID field in the vcf
grep "^#" herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed.vcf > header
grep -v "^#" herb_893FDR_non_clumped_sorted_POSfixed.vcf > core
awk 'NR <= 0 {print; next} {OFS="\t"; $3 = $1 ":" $2; print}' core > core_ID

cat header core_ID > herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed.vcf

#sort the association file
sort -k13,13V -k14,14g FDR_non_clumped_significant_sites_filtered.assoc.txt >> FDR_non_clumped_significant_sites_filtered_sorted.assoc.txt
#adding the header line for the association file :
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2
head -n 1 ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt > /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/header_line
cat header_line FDR_non_clumped_significant_sites_filtered_sorted.assoc.txt > FDR_non_clumped_significant_sites_filtered_sorted_header.assoc.txt


```

Create plink family files from the vcf :

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium
#make plink family files and create ID based on position information
module load plink
vcf=herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed.vcf

plink --vcf $vcf --out herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed --allow-extra-chr --recode --double-id 
```


Run plink to clump sites.



```

assoc=FDR_non_clumped_significant_sites_filtered_sorted_header.assoc.txt
plink --file herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed --clump $assoc --clump-p1 0.05 --clump-field FDR --clump-kb 100 --out herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed_100kb --allow-no-sex --allow-extra-chr --clump-snp-field ID

#worked!
Calculating allele frequencies... done.
Total genotyping rate is 0.926832.
92 variants and 108 people pass filters and QC.
Note: No phenotypes present.
--clump: 86 clumps formed from 92 top variants.
Results written to
herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed_100kb.clumped
```

### Step (2) : prep the input_file for ancestry_hmm

<ins>(1) Make bed file from clumped sites (86 sites) </ins>

```
#make bed file from the 86 clumped sites

awk 'NR <= 1 {next} {OFS="\t"; print $1,$4-1,$4}' herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed_100kb.clumped > herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed_100kb.bed

#remove the two last lines (crap)

head -n -2  herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed_100kb.bed > herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed_clean_100kb.bed
```

<ins>(2) Filter the herbarium vcf based on the bed file </ins>

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium
bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed_clean_100kb.bed

cd /cds3/kreiner/herbarium/
module load vcftools

vcf=merged_numericChr.vcf

vcftools --vcf $vcf --positions $bed --recode --stdout > herb_86clumped.vcf
```

<ins>(3) Get the read counts </ins>

Make sure the vcf file is sorted :
```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_prep_file
herb_86clumped.vcf
grep "^#" herb_86clumped.vcf > herb_86clumped_sorted.vcf
grep -v "^#" herb_86clumped.vcf | sort -k1,1V -k2,2g >> herb_86clumped_sorted.vcf

```
I ran the [vcf_to_read_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/vcf_to_read_counts.awk) script to convert the genotype information to read counts (will be the right side of the ancestry_hmm input file). 

The file is : /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_prep_file/herb_86clumped_allele_count.txt

<ins>(4) Calculate LD between those 86 sites </ins>

I use the script [calculate_ldhat_between_sites.Rmd](https://github.com/rozenn-pineau/Drought-paper/blob/main/calculate_ldhat_between_sites.Rmd) to calculate LD between each site in the bed file. 
(It is better to calculate LD before getting the read counts ready (step 3) because some sites that are outside of our LD map, such that they are removed from the analysis. In this case, everything was within the boundaries because the filtering had previously been done with the drought files).

<ins>(5) Filter ancestry panel for the herbarium sites </ins>

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_prep_file

bed=herb_86clump_for_ancestry.bed

vcf=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/vcfs/var_rudis_q75.vcf
vcftools --vcf $vcf --positions $bed --recode --stdout > var_rud_86clump.vcf

vcf=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/vcfs/var_tub_q75.vcf
vcftools --vcf $vcf --positions $bed --recode --stdout > var_tub_86clump.vcf
```

<ins>(6) Vcf to genotype counts for ancestry panel </ins>

I use the script [genotype_to_allele_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/genotype_to_allele_counts.awk) to calculate the number of genotypes in the ancestry panel :
```
/scratch/midway3/rozennpineau/drought/scripts/genotype_to_allele_counts.awk var_tub_86clump.vcf > allele_counts_86clump_var_tub.txt

/scratch/midway3/rozennpineau/drought/scripts/genotype_to_allele_counts.awk var_rud_86clump.vcf > allele_counts_86clump_var_rud.txt
```

<ins>(7) Assemble the input_file  </ins>

ancestry panel 0 :

/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_prep_file/allele_counts_86clump_var_tub.txt

ancestry panel 1  :

/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_prep_file/allele_counts_86clump_var_rud.txt

LD file :

/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_prep_file/herb_86clumped_ld.txt

sample counts :

/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_prep_file/herb_86clumped_allele_count.txt


```
#var rudis allele counts 
awk '{OFS="\t"; print $3,$4}' allele_counts_86clump_var_rud.txt > var_rud.allele_count

#rho column
tail -n 86  herb_86clumped_ld.txt | awk '{OFS="\t"; print $4}' > herb86_clumped.ld

#assemble chrom, pos, var tub allele counts, var rudis allele counts, rho, then sample read counts
paste allele_counts_86clump_var_tub.txt var_rud.allele_count herb86_clumped.ld herb_86clumped_allele_count.txt > herb86_clumped_input_file.txt
```


<ins>(7) Prep the sample file </ins>

The sample file is simply : samplename ploidy, in the same order as in the allele counts files.

```
bcftools query -l herb_893FDR_non_clumped_sorted_POSfixed_FORMATfixed_IDfixed.vcf > sample_file.tmp
awk '{OFS="\t" ; print $1, 2}' sample_file.tmp > sample_file.txt
```

### Step (3) : run ancestry_hmm

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/3_run/two_pulse

#!/bin/bash
#SBATCH --job-name=anc_hmm
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10GB

module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

ancestry_hmm -i herb86_clumped_input_file.txt -s sample_file.txt -a 2 0.33 0.67 -p 0 100000 0.33 -p 1 -9999 0.43 -p 1 -100 0.24 --tmax 10000
```
### Step (4) : process output from ancestry_hmm

At this step, we realized the probabilities were very low, even on an one-pulse model. This is most certainly because we are calling ancestry on 86 sites, which is not a lot. We thus decide to go back a few steps and to call ancestry on a bigger set of variants. 

## Calling ancestry on more sites


