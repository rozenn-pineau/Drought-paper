# Drought-paper
Suite of notes and scripts for the drought project.


# Outline
[Drought LD50 calculations](#Drought-LD50-calculations)

[Herbicide resistance : EPSPS depth estimations](#Herbicide-resistance-:-EPSPS-depth-estimations)

[Estimating global ancestry using ADMIXTURE](#Estimating-global-ancestry-using-ADMIXTURE)

[Drought dataset - ancestry mapping using ancestry hmm](#Drought-dataset---ancestry-mapping-using-ancestry-hmm)

[Processing results from ancestry_hmm output](#Processing-results-from-ancestry_hmm-output)

[GWAS on ancestry calls](#GWAS-on-ancestry-calls)

[GO enrichment analyses and sites of potential interest](#GO-enrichment-analyses-and-sites-of-potential-interest)
- snpEff

[Drought selection experiment trajectory analyses](#Drought-selection-experiment-trajectory-analyses)

[Herbarium dataset - ancestry mapping using ancestry hmm](#Herbarium-dataset---ancestry-mapping-using-ancestry-hmm)

[Calling ancestry on drought-informative sites only](#Calling-ancestry-on-drought-informative-sites-only)
    
[Calling ancestry on more sites](#Calling-ancestry-on-more-sites)

[Calculating Hazard Ratios on ancestry calls](#Calculating-Hazard-Ratios-on-ancestry-calls)

[Climate-Genotype correlation analyses](#Climate-Genotype-correlation-analyses)

# Drought LD50 calculations

# Herbicide resistance : EPSPS depth estimations
To quantify resistance in drought samples, we calculated the depth at EPSPS position in each sample. 

We first extracted the lines from the gff files that have "EPSPS" to identify the gene location in the genome: 

```
#extract lines with EPSPS
grep "EPSPS" Atub_193_hap2.all.sorted.gff > /scratch/midway3/rozennpineau/drought/intersect_gff.txt

```
There was one site within the genome 16 scaffolds : 
```
Scaffold_6      maker   mRNA    17472580        17484501
```
We further used mosdepth to extract the depth at this site for every sample (from the bam files)  : 

```
cd /cds3/kreiner/drought/bams
bed=/scratch/midway3/rozennpineau/drought/EPSPS/EPSPS.bed
for bam in *.bam
do
mosdepth -n -x --by $bed /scratch/midway3/rozennpineau/drought/EPSPS/mosdepth_output/$bam $bam
done

#From the mosdepth output, we extracted the depth into one file with the corresponding sample name.
for file in *regions.bed.gz
do
zcat $file | awk 'BEGIN {OFS="\t"}; {print $4}' >> depth.txt
done

#add the file name info in column 1 (extract name from file name?)
ls *regions.bed.gz > file_names.tmp #make a file of file names
cut -d "_" -f 1,2 file_names.tmp > samp_ID.tmp #extract important IF info
paste samp_ID.tmp depth.txt > ID_depth.txt #tab delim file
```
After extracting te depths for each sample at EPSPS gene, we need to extract the genome wide mean coverage for normalization : 

```
#make bed file from gff
awk 'BEGIN {OFS="\t"}; {print $1,$4,$5}' /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff >gff.bed
uniq gff.bed > gff_clean.bed

#extract depth for each coding site
bed=/scratch/midway3/rozennpineau/drought/EPSPS/gff_clean.bed

cd /cds3/kreiner/drought/bams

for bam in *.bam
do
mosdepth -n -x --by $bed /scratch/midway3/rozennpineau/drought/EPSPS/gff_depths/$bam $bam 
done

#gather in one file  : paste them together 

zcat AT27_827_193_2.dd.bam.regions.bed.gz | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}' > depth.txt
#410700 lines

for file in *regions.bed.gz
do
zcat $file | awk 'BEGIN {OFS="\t"}; {print $4}' > tmp
cp depth.txt tmp2
paste tmp2 tmp > depth.txt
done

#same file order as the EPSPS depth file

```


# Estimating global ancestry using ADMIXTURE
To calculate genome-wide ancestry, we used ADMIXTURE version 1.3.0 (Alexander et al., 2009). 
We tested different values for K and settled on K=2: 

```
#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=15:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100G   # memory per cpu-core

module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

my_bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/commongarden_allfiltsnps_193_hap2_numericChr_filt.bed

for K in 2;

do admixture --cv $my_bed $K ;

done
```



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
#Add header to bed file
awk '{print
ancestry_gwas_filtered_sites.bed


module load vcftools
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2
vcftools --vcf two_pulse_flexible_prop_2_values_ID.vcf --bed ancestry_gwas_filtered_sites_wheader.bed --out two_pulse_flexible_prop_2_values_ID_filtered --recode


```

output 

```
Parameters as interpreted:
	--vcf two_pulse_flexible_prop_2_values_ID.vcf
	--out two_pulse_flexible_prop_2_values_ID_filtered
	--recode
	--bed ancestry_gwas_filtered_sites_wheader.bed

After filtering, kept 282 out of 282 Individuals
Outputting VCF file...
	Read 37 BED file entries.
After filtering, kept 36 out of a possible 786261 Sites
Run Time = 3.00 seconds
````
36 loci in the output vcf file.


Getting the ancestry calls (in the form of genotypes) from the filtered vcf for more downstream analyses :
```
module load htslib
bgzip -f two_pulse_flexible_prop_2_values_ID_filtered.recode.vcf

tabix -f two_pulse_flexible_prop_2_values_ID_filtered.recode.vcf.gz

act-conda
bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' two_pulse_flexible_prop_2_values_ID_filtered.recode.vcf.gz > two_pulse_flexible_prop_2_values_ID_filtered_GT.txt #36 loci
```


## Checking : Does ancestry predict response to drought ?

Our expectation is that var. rudis ancestry is better adapted to drought than var. tuberculatus. 
Do we see this in our ancestry calls ? This is also simply a way to make sure that our pipeline was coded correctly. 

# GO enrichment analyses and sites of potential interest

We identified 36 sites significantly associated with drought adaptation. To test whether those sites are in functional regions of the genome with known effects, we extracted the regions from the annotated genome file.

The gff file has "Scaffold_" as chromosome names, so I need to update the bed chromosome names :

```
awk -F'\t' -vOFS='\t' '{ $1 = "Scaffold_" $1}1' ancestry_gwas_filtered_sites.bed > ancestry_gwas_filtered_sites_scaffold_names.bed 
```
Intersect the bed file with the gff file using bedtools :

```
bedtools intersect -b ancestry_gwas_filtered_sites_scaffold_names.bed -a /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff > intersect_FDR_gff_enriched_genes.txt
```
Extract the GO terms and Note field from the file : 

```
cut -d ";" -f 9 intersect_FDR_gff_enriched_genes.txt | grep GO > intersect_FDR_gff_GO_terms_bon.txt
cut -d ";" -f 10 intersect_FDR_gff_enriched_genes.txt | grep Note > intersect_FDR_gff_note.txt
```


### snpEff 

I followed the below steps to find the effect of the mutation on the protein function (thank you, Jake!) : 

```
#make a conda environment for snpeff
conda create -n snpeff

#activate that environment
conda activate snpeff

#install snpeff
conda install bioconda::snpeff

#move to the directory where conda installed snpeff
cd .conda/envs/snpeff/share/snpeff-5.2-1/

#make a new directory to store your database(s)
mkdir data

#change into there and make a directory for your database
cd data
mkdir Atub_193_hap2

#change into there and copy the genome sequence to a file called sequences.fa
cd Atub_193_hap2
cp /project/kreiner/data/genome/Atub_193_hap2.fasta sequences.fa

#either install gffread in this environment, or in my case, I use a different conda environment
#conda deactivate
#conda create -n gffread
#conda activate gffread
conda install bioconda::gffread

#   .. loaded 35100 genomic features from /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff

#use gffread to convert the gff file into gtf
gffread -E /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -T -o genes.gtf

#use gffread to create cds and protein sequence files
gffread -x cds.fa -g /project/kreiner/data/genome/Atub_193_hap2.fasta /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff
gffread -y protein.fa -g /project/kreiner/data/genome/Atub_193_hap2.fasta /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff

#move back over to your snpeff environment and add the database that you want to build to the snpeff config file
nano snpEff.config

# add the following lines to the config file right below the '# Non-standard Databases' lines, above 'Homo sapiens (hg19) (UCSC)'

# Atuberculatus genome, version 193_hap2
Atub_193_hap2.genome : Atub_193_hap2

#build the database
snpEff build -gtf22 -v Atub_193_hap2

#run snpeff on the vcf
#filter the vcf based on the bed file
bedtools intersect -header -b /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/FDR_significant_drought_sites.bed -a /scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/two_pulse_flexible_prop_2_values_ID.vcf > FDR_significant_drought_sites.vcf
#update scaffold names to include "Scaffold_"
awk -F'\t' -vOFS='\t' '{$1 = "Scaffold_" $1}1' FDR_significant_drought_sites.vcf > FDR_significant_drought_sites_scaffold_names.vcf
snpEff ann Atub_193_hap2 FDR_significant_drought_sites_scaffold_names.vcf > FDR_significant_drought_sites_ann.vcf

```

### CMH scan - GWAS output comparison
we want to compare the outputs from the CMH scans that identidied "agriculturally-adapted alleles" to the output of the GWAS on drought. I do this by comparing bed files both both (1) all loci, and (2) significant loci in both. 


On the cluster, I prepare the bed files for the CMH scans for comparison. 
```
#prepare the full CMH file bed (FDRdrought is the CMH scan output)
cat FDRdrought | \sed s/^\Scaffold_//g | awk -v OFS="\t" '{ print $1,$3,$3,$13}' > CMH_all.bed #63,979,747
#remove header
tail -n +2 CMH_all.bed > CMH_49338567.bed

#prepare the significant loci bed
cat FDRdrought | \sed s/^\Scaffold_//g | awk -v OFS="\t" '{ if ($13 <= 0.05) { print $1,$3,$3,$13} }' > CMH_383650.bed #383650

#I also subselect the sites that make through the Bonferonni threshold. 
cat FDRdrought | \sed s/^\Scaffold_//g | awk -v OFS="\t" '{ if ($14 <= 0.05) { print $1,$3,$3,$14} }' > CMH_BON.bed  #1419 loci

```

Then I use bedtools intersect to find intersections between beds. 

```
#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=5:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100GB

#conda environment
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

#load files
cd /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/
gwasbed=gwas_all.bed
cmhbed=CMH_all.bed
bedtools intersect -a $gwasbed -b $cmhbed -wa -wb -f 0.99 -r > intersect_all_drought_cmh_gwas.bed
#-r 1 -f 1 requires that there is at least 1 bp match, with 100% of the regions matching

gwasbed=gwas_893.bed
cmhbed=CMH_383650.bed
bedtools intersect -a $gwasbed -b $cmhbed -wa -wb -f 0.99 -r > intersect_significant_drought_cmh_gwas.bed

#for Bonferronni loci
gwasbed=gwas_893.bed
cmhbed=CMH_BON.bed 
bedtools intersect -a $gwasbed -b $cmhbed > intersect_BON_drought_cmh_gwas.bed #no intersection found

```


I further look for where in the genome the hits are : are there in genes ? I look into the gff file for this.
```
#prep bed file
awk '{OFS="\t";print $1, $2, $2}' common_significant_cmh_015.bed > cmh_015.bed
tail -n 13602 cmh_015.bed > cmh_015_noheader.bed

#add “Scaffold_” to overlap bed
awk -v OFS="\t" '{$1 = "Scaffold_" $1}1' cmh_015_noheader.bed > cmh_015_noheader_scaffold.bed
awk -v OFS="\t" '{$1 = "Scaffold_" $1}1' cmh_gwas_35.bed > cmh_gwas_35_scaffold.bed
awk -v OFS="\t" '{print $1, $2, $2}' common_significant_cmh_gwas_scaffold.bed > common_significant_cmh_gwas_scaffold_clean.bed
```

#bedtools intersect with gff
```
#remove header from bed file for bedtools intersect
bedtools intersect -a cmh_015_noheader_scaffold.bed -b /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -wo > overlap_cmh_015_gff.bed
bedtools intersect -a cmh_gwas_35_scaffold_noheader.bed -b /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -wo > overlap_gwas_cmh_35_gff.bed
bedtools intersect -a common_significant_cmh_gwas_scaffold_clean2.bed -b /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -wo > overlap_gwas_cmh_gff.bed
```

#extract gene names
```
cut -f 12 overlap_cmh_015_gff.bed | grep "Similar" | cut -d" " -f3 | cut -d":" -f1 | sort | uniq > gwas_cmh_genes_35.txt
cut -f 12 overlap_gwas_cmh_35_gff.bed | grep "Similar" | cut -d" " -f3 | cut -d":" -f1 | sort | uniq > gwas_cmh_genes_35.txt
cut -f 12 overlap_gwas_cmh_gff.bed | grep "Similar" | cut -d" " -f3 | cut -d":" -f1 | uniq > gwas_cmh_genes.txt
```

















# Drought selection experiment trajectory analyses

## upload scripts for trajectory analyses

### Compare to randomized distribution
Step (1) : generate list with allele frequency change from day 20 to day 1 for similar starting allele frequencies (0.001 precision) : [generate_null.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/generate_null.R).



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

### Step (1) : keep the variants common to the ancestry variant file and the herbarium variant file
To build the file that will be used to call ancestry at each of the herbarium site, we keep the intersection between the ancestry variant file and the herbairum variant file using bcftools isec.

```
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

#bgzip and tabix the vcf files

anc=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/fst/ancestry_fst_q75.vcf.gz

herb=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/merged_numericChr.vcf.gz

cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/

bcftools isec -p herb_ancestry $anc $herb

```

At this step, we have two variant files with the same variants for the two populations (976662 sites). I extract the position and chromosome information (bed) to calculate rho between each site. 

### Step (2) : calculate rho between each site

Extract chromosome and position information from the vcfs :

```
awk 'BEGIN {OFS="\t"} !/^#/ {print $1, $2-1, $2}' ancestry_common_976662.vcf  > ancestry_common_976662.be

```

Rscript to calculate rho between sites : [calculate_ldhat_between_sites.Rmd](https://github.com/rozenn-pineau/Drought-paper/blob/main/calculate_ldhat_between_sites.Rmd).

Because we exclude the sites that are outside of the region defined by the recombination map (the monotonic spline does not extrapolate outside of boundaries), we kept track of which sites to filter out of the variant files to further filter the vcf files.

```
#make bed file
awk '{print $1,$2,$3}' ancestry_herb_common_974134.ld | tail -n -974134 > ancestry_herb_common_974134.bed

#filter
module load vcftools
bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/ancestry_herb_common_974134.bed
vcf1=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/herb_ancestry/herb_common_976662.vcf
vcftools --vcf $vcf1 --bed $bed --out herb_common_974134.vcf --recode

vcf2=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/herb_ancestry/ancestry_common_976662.vcf
vcftools --vcf $vcf2 --bed $bed --out ancestry_common_974134.vcf --recode

```

### Step (3) : split ancestry vcf into tuberculatus versus rudis files

The ancestry file both has var rudis and var tuberculatus samples, that we will now split into two :

```
bcftools view -S var_rudis_samp.txt ancestry_common_974133.vcf > rudis_common_974133.vcf
bcftools view -S var_tub_samp.txt ancestry_common_974133.vcf > tub_common_974133.vcf

#check # of variants and samples
#44 samples for rudis and 21 samples for tuberculatus
#974133 variants for both
```


### Step (4) : extract allele counts from the ancestry variant files

[genotype_to_allele_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/genotype_to_allele_counts.awk)

```
/scratch/midway3/rozennpineau/drought/scripts/genotype_to_allele_counts.awk tub_common_974133.vcf > tub_common_974133.allele_counts

/scratch/midway3/rozennpineau/drought/scripts/genotype_to_allele_counts.awk rudis_common_974133.vcf > rudis_common_974133.allele_counts

```

### Step (5) : get genotypes for drought panel

2,0 for homs reference
1,1 for hets
0,2 for homs alternative

[vcf_to_read_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/vcf_to_read_counts.awk)

### Step (6) : put the file together

Paste columns together to make the full file, following instructions on ancestry_hmm github page. 

```


#var rudis allele counts 
awk '{OFS="\t"; print $3,$4}' tub_common_974133.allele_counts > var_rud.allele_counts

#rho column
awk '{OFS="\t"; print $4}' ancestry_herb_common_974133.ld > ancestry_herb_common_974133_rho.ld

#assemble chrom, pos, var tub allele counts, var rudis allele counts, rho, then sample read counts
paste tub_common_974133.allele_counts var_rud.allele_counts ancestry_herb_common_974133_rho.ld herb_common_974133.read_counts > herb_common_974133_input_file.txt


```

### Step (7) : make sample file


# Process output from ancestry_hmm

```
folder_name=$(basename "$PWD")
output_file="${folder_name}_values_09.txt"
cutoff=0.9

#initialize dataset with chrom and pos
awk 'BEGIN { OFS="\t" } {print $1,$2}' HB0900.posterior > $output_file

#loop through each individual and paste information to previous version of the file
for file in *.posterior; do

    header_name=$(basename "$file" .posterior)

    awk -v header="$header_name" 'NR == 1 { print header; }
    
    NR > 1 {
    
    result = "NA";
    if ($3 > $cutoff) {result = "0"}
    else if ($4 > $cutoff) {result = "1"}
    else if ($5 > $cutoff) {result = "1"}
    else if ($6 > $cutoff) {result = "2"}
    else if ($7 > $cutoff) {result = "2"}
    else if ($8 > $cutoff) {result = "2"}
    print result;
    }' $file > tmp

    paste $output_file tmp > tmp2

    mv tmp2 $output_file

done
```

This is for every site in common between the herbarium variant file and the ancestry panels. We are interest in the trajectories of drought-informative alleles. Let's filter the file for those sites only.  

### (1) filter ancestry call table for drought informative sites

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions

# Input files
bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/FDR_non_clumped_significant_sites.bed
assoc=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions/anc_calls_herbarium_all_sites_09.txt
out=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions/anc_calls_herbarium_893sites_09.txt

# Define column positions for chrom and pos in the TXT file (1-based index)
CHROM_COL=1
POS_COL=2

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

mv anc_calls_herbarium_893sites_09.txt anc_calls_herbarium_680sites_09.txt
#680 sites in common

#add header back to ancestry call file
```
We now have 680 sites with ancestry calls, that are drought informative sites. 

### (2) Filter association file for these 680 drought informative sites
```
#make bed from vcf file
vcf=anc_calls_herbarium_680sites_09.vcf 
bed=anc_calls_herbarium_680sites_09.bed

awk 'BEGIN {OFS="\t"} 
     !/^#/ {print $1, $2-1, $2}' "$vcf" > "$bed" #position in vcf is shifted by one

# Input files
bed=anc_calls_herbarium_680sites_09.bed
assoc=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt
out=FDR_non_clumped_680sites.assoc.txt

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
### (3) make vcf file for plink LD thinning and drop variants with more than 75% missing data

**make vcf file**

Plink takes in a vcf file, that we create based on the ancestry call table. 

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions/

# Input and output files
input_file="anc_calls_herbarium_680sites_09.txt"
output_file="anc_calls_herbarium_680sites_09.vcf"

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
awk 'NR <= 5 {print; next} {OFS="\t"; $3 = $1 ":" $2; print}' anc_calls_herbarium_680sites_09.vcf > anc_calls_herbarium_680sites_09_ID.vcf
```

**drop loci with more than 75% missing data**


### (4) Clumping based on LD

```
assoc=FDR_non_clumped_680sites_header.assoc.txt

#without missing loci threshold

plink --file anc_calls_herbarium_680sites_09_ID --clump $assoc --clump-p1 6.5830432740603e-05 --clump-field p_wald --clump-kb 100 --out two_pulse_flexible_prop_2_clumped_100kb --allow-no-sex --allow-extra-chr --clump-snp-field ID

#output
Logging to two_pulse_flexible_prop_2_clumped_100kb.log.
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --clump FDR_non_clumped_680sites_header.assoc.txt
  --clump-field p_wald
  --clump-kb 100
  --clump-p1 6.5830432740603e-05
  --clump-snp-field ID
  --file anc_calls_herbarium_680sites_09_ID
  --out two_pulse_flexible_prop_2_clumped_100kb

257091 MB RAM detected; reserving 128545 MB for main workspace.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (680 variants, 108 people).
--file: two_pulse_flexible_prop_2_clumped_100kb-temporary.bed +
two_pulse_flexible_prop_2_clumped_100kb-temporary.bim +
two_pulse_flexible_prop_2_clumped_100kb-temporary.fam written.
680 variants loaded from .bim file.
108 people (0 males, 0 females, 108 ambiguous) loaded from .fam.
Ambiguous sex IDs written to two_pulse_flexible_prop_2_clumped_100kb.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 108 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.677124.
680 variants and 108 people pass filters and QC.
Note: No phenotypes present.
--clump: 33 clumps formed from 680 top variants.
Results written to two_pulse_flexible_prop_2_clumped_100kb.clumped .

#--geno <threshold>: removes SNPs where the missing genotype proportion exceeds the provided threshold. --geno 0.25 removes loci with more than 25% missing genotype data, which corresponds to keeping loci with at least 75% non-missing genotype data.

```


### (5) extract clumped variants 
(1) make bed file to filter ancestry calls vcf file
```
awk -v OFS='\t' '{ print $1, $4, $4, $5}' two_pulse_flexible_prop_2_clumped_100kb.clumped > herbarium_filtered_sites.bed
```
(2) filter the vcf file based on the bed file
```
module load vcftools
vcftools --vcf anc_calls_herbarium_680sites_09_ID.vcf --bed herbarium_filtered_sites.bed --out anc_calls_herbarium_33sites --recode

Parameters as interpreted:
	--vcf anc_calls_herbarium_680sites_09_ID.vcf
	--out anc_calls_herbarium_33sites
	--recode
	--bed herbarium_filtered_sites.bed

After filtering, kept 108 out of 108 Individuals
Outputting VCF file...
	Read 34 BED file entries.
After filtering, kept 33 out of a possible 680 Sites
Run Time = 0.00 seconds
```
(3) getting ancestry calls for downstream analyses
```
module load htslib
bgzip -f anc_calls_herbarium_33sites.recode.vcf 
tabix -f anc_calls_herbarium_33sites.recode.vcf.gz

bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' anc_calls_herbarium_33sites.recode.vcf.gz > anc_calls_herbarium_33sites_GT.txt
```



## Estimating ancestry variation across the genome

To estimate the vaiation in ancestry across the genome, we called the mean ancestry per individual based on 1Mb windows across the genome using [estimate_ancestry_variation.sh](https://github.com/rozenn-pineau/Drought-paper/blob/main/estimate_ancestry_variation.sh). 



## Calculating Hazard Ratios on ancestry calls
To evaluate the advantage of bearing each drought-adapted site independently as well as the combination of them, we calculate the hazard ratios (Cox model) : 
- individual tests : [cox_ancestry_calls.Rmd](https://github.com/rozenn-pineau/Drought-paper/blob/main/cox_ancestry_calls.Rmd)
- test on polygenic scores : 


## Compare Drought sites from the Drought experiment and the sites identified as agriculturally adapted in Kreiner et al. 2022
We compared the sites that we identified in the drought experiment, to the sites that had been labeled as "agriculturally adapted" in the Kreiner 2022 paper.
There was 56 out of 893 sites in common with the drought dataset, and the drought only dataset where the CMH scan looked for agriculturally-enriched alleles (out of 383650 variants).

There was 54 sites in common between the drought dataset and the drought+common garden combined datasets where the CMH scan looked for agriculturally-enriched alleles (149883 variants). 

Example of how files were compared : 
```
#drought and common garden CMH scan comparison

cat FDRcmh_dcg | \sed s/^\"Scaffold_//g > tmp1

cat tmp1 | \sed s/\"//g > tmp2

awk '{ if ($13 <= 0.05) { print $1,$3-1,$3} }' tmp2 > FDR_significant_drought_commongarden_CMH_sites.bed #149883 sites

awk 'OFS="\t" {print $1,$2,$3}' FDR_significant_drought_commongarden_CMH_sites.bed > FDR_significant_drought_commongarden_CMH_sites_tab.bed # tab

bedtools intersect -a FDR_significant_drought_commongarden_CMH_sites_tab.bed  -b FDR_significant_drought_sites.bed > overlap_CHdrought_common_garden_drought.bed #54
```

In which genes are these variants situated ?


```
#update chromosome names
awk -F'\t' -vOFS='\t' '{ $1 = "Scaffold_" $1}1' FDR_significant_drought_commongarden_CMH_sites_tab.bed > FDR_significant_drought_commongarden_CMH_sites_scaffold_names.bed

#Intersect with bedtools :
bedtools intersect -b FDR_significant_drought_commongarden_CMH_sites_scaffold_names.bed -a /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff > intersect_CHMscan_drought_commongarden_gff.txt

#how many chromosomes are the genes in
grep Atub intersect_CHMscan_drought_commongarden_gff.txt | cut  -f1 | uniq # all 16 chromosomes

#grep unique gene names to count the number of genes
grep ID= intersect_CHMscan_drought_commongarden_gff.txt | cut -d'=' -f2 | grep hap2 | cut -d'-' -f1 | uniq > gene_names.tmp
cut -d';' -f1 gene_names.tmp | uniq > genes_commongarden_drought_gff.txt

```


# Herbarium and climate data analyses

### Extracting the climate data at each herbarium sample location
Step (1) - find the closest weather station with mean temperature and precipitation dating back to 1828, which is our earliest herbarium sample. 
Step (2) - extract the climate data for every year since 1828 for these weather stations. 


### Haplotype frequency change with climate change 
One we extracted the weather data for each sample at each date, we looked for consistent change through time in both climate and haplotype frequencies. 
Does an increase in temperature correlate with an increase in drought-adapted haplotype frequency ?

To do this, we calculated the physical distance between each pair of sample and kept the two samples that were the closest in space, with a maximum distance of 100 km. For each pwit, we calculated the allele frequency change (the time difference here can change according to the sample pair) and the climate change (total precipitation  and max temperature for the months of June, July, August combined). We ended up with 36 pairs for max temperature, 41 pairs for precipitation and 36 pairs for climate PC1). 

We identified a few different sites that track climate. 

# Climate-Genotype correlation analyses

[Figure 5 analyses - correlations between climate change and haplotype frequency changes]

Scrit for pairwise comparisons : order the pairwise distance between samples by ascending order; find the two closest samples, set aside; find the closest sample of the next one in the list, set aside; repeat until all samples have been used (23 pairs); calculate haplotype frequency change and climate change, test linear relationship, record if significant **and** postisive; go back to the list and start with the second sample this time, repeat (146 times for 146 possible pairs with distance in km >0 and < 100). 

For temperature max  : [run_pairwise_data.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/run_pairwise_data.R)

