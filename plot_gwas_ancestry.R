# Plot results from the gwas analysis done on ancestry mapped sites and day to full wilt with and without correction for population structure
rm(list= ls())

library(data.table)
library(ggplot2)
library(dplyr)
library(qqman)
library(tidyverse)
library(extrafont)


# Theoretical quantile function for -log10(p)
qlog10 <- function(p.values) {
  theoretical <- rank(p.values)/length(p.values)
  return(-log10(theoretical))
}
 
#this script is for the current set of filters that we are using - might be updated if we chose other combinations of filters. 
#It also does not include all previous analyses that led ud to the choice of these filters.

# ----------------------------------------------------------------------------- #
# read in regression results files
# ----------------------------------------------------------------------------- #

#new runs with -miss 0.2 and cutoff at 0.9
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/8.ancestry_hmm/gwas/gemma_cutoff09_miss02")
lmm_gemma <- read.table("result.assoc.txt", header= T, sep = "\t")
lmm_gemma_geno <- read.table("geno_corrected_gemma_gwas.assoc.txt", header= T, sep = "\t")
lmm_gemma_anc <- read.table("ancestry_corrected_gemma_gwas.assoc.txt", header= T, sep = "\t") 
#modify rs col to get chrom position information
lmm_gemma[c("chrom", "pos")] <- as.numeric(str_split_fixed(lmm_gemma$rs, ':', 2))
lmm_gemma_geno[c("chrom", "pos")] <- as.numeric(str_split_fixed(lmm_gemma_geno$rs, ':', 2))
lmm_gemma_anc[c("chrom", "pos")] <- as.numeric(str_split_fixed(lmm_gemma_anc$rs, ':', 2))

# ----------------------------------------------------------------------------- #
# correct p-values
# ----------------------------------------------------------------------------- #
df_vec <- c("lmm_gemma", "lmm_gemma_geno", "lmm_gemma_anc")

for (d in 1:length(df_vec)) {
  
  df <- get(df_vec[d])
  N <- dim(df)[1] 
  if (length(df$pval) >0 ){
  df$FDR <- p.adjust(df$pval, method = "fdr", n = N) 
  df$bon <- p.adjust(df$pval, method = "bonferroni", n = N)
  df$qlogp <- qlog10(df$pval)
  df <- df[order(-df$qlogp),]
  }
  if (length(df$p_wald) >0 ) {
    df$FDR <- p.adjust(df$p_wald, method = "fdr", n = N) 
    df$bon <- p.adjust(df$p_wald, method = "bonferroni", n = N)
    df$qlogp <- qlog10(df$p_wald)
    df <- df[order(-df$qlogp),]
  }
  
  #update dataset
  assign(df_vec[d], df)
  
}


# ----------------------------------------------------------------------------- #
# QQPlots
# ----------------------------------------------------------------------------- #

# QQplot with base 10 - lrt test for ancestry corrected gemma gwas
ggplot(lmm_gemma, aes(x = qlogp, y = -log10(p_wald))) + 
  geom_path() + ylim(0,10) + xlim(0,10) + 
  geom_abline(intercept = 0, slope = 1)

# QQplot with base 10 - lrt test for ancestry corrected gemma gwas
ggplot(lmm_gemma_anc, aes(x = qlogp, y = -log10(p_wald))) + 
  geom_path() + ylim(0,10) + xlim(0,10) + 
  geom_abline(intercept = 0, slope = 1)

# QQplot with base 10 - lrt test for genotype corrected gemma gwas
ggplot(lmm_gemma_geno, aes(x = qlogp, y = -log10(p_wald))) + 
  geom_path() + ylim(0,10) + xlim(0,10) + 
  geom_abline(intercept = 0, slope = 1)

# ----------------------------------------------------------------------------- #
# Genomic inflation factor
# ----------------------------------------------------------------------------- #
#inflation factor 

# Number of tests (GWAS SNPs)
sorted_pvals <- sort(lmm_gemma_anc$p_wald)
n <- length(sorted_pvals)

# Expected p-values under the null hypothesis
expected_pvals <- (1:n) / n

# Calculate the genomic inflation factor (lambda, λ)
observed_chi2 <- qchisq(1 - sorted_pvals, df=1)  # Chi-squared test statistic for each p-value
observed_lambda <- mean(observed_chi2) / qchisq(0.5, df=1)  # Mean observed chi-squared statistic divided by the expected (median)

# Output the inflation factor (lambda)
observed_lambda #2.356066
# Adjust p-values using genomic control
lmm_gemma_anc$inflation_p <- pmin(lmm_gemma_anc$p_wald * observed_lambda, 1)  # p-values can't exceed 1

#pmin = finds the element-wise minimum across multiple vectors; returns a vector of the same length as the input vectors, 
# where each element is the minimum of the corresponding elements across the input vectors

# Print adjusted p-values
hist(-log10(lmm_gemma_anc$inflation_p))

lmm_gemma_anc$FDR<-p.adjust(lmm_gemma_anc$inflation_p,method = "fdr")
lmm_gemma_anc$bon<-p.adjust(lmm_gemma_anc$inflation_p,method = "bonferroni")
lmm_gemma_anc$qlogp <- qlog10(lmm_gemma_anc$inflation_p)

# ----------------------------------------------------------------------------- #
# export association file with inflated corrected p-values
# ----------------------------------------------------------------------------- #
write.table(lmm_gemma_anc, "ancestry_corrected_inflated_gemma_gwas.assoc.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#export chromosome and positions that are significant FDR sites to filter the herbarium vcf
sites <- lmm_gemma_anc$FDR < 0.05 #893
export <- data.frame(chrom = lmm_gemma_anc$chrom[sites],
                     pos1 = lmm_gemma_anc$pos[sites]-1,
                     pos2 = lmm_gemma_anc$pos[sites])
write.table(export, "FDR_non_clumped_significant_sites.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# ----------------------------------------------------------------------------- #
# quick look at potentially significant SNPs
# ----------------------------------------------------------------------------- #

quick_stats <- matrix(NA,3,4)
quick_stats <- as.data.frame(quick_stats)
colnames(quick_stats) <- c("minFDR", "minBON", "numFDR", "numBON")
rownames(quick_stats) <- df_vec

for (d in 1:length(df_vec)) {
  
  df <- get(df_vec[d])

  quick_stats[d,1] <- min(df$FDR, na.rm = T)
  quick_stats[d,2] <- min(df$bon, na.rm = T) 
  quick_stats[d,3] <- length(which(df$FDR < 0.05))
  quick_stats[d,4] <- length(which(df$bon < 0.05))

}

quick_stats

# ----------------------------------------------------------------------------- #
# QQPlots
# ----------------------------------------------------------------------------- #

# QQplot with base 10 - lrt test for uncorrected manual gwas
ggplot(lmm_gemma_anc, aes(x = qlogp, y = -log10(inflation_p))) + 
  geom_path() + ylim(0,10) + xlim(0,10) + 
  geom_abline(intercept = 0, slope = 1)
# to test if there is enrichment in low p-values compared to the null, which are uniformly distributed under the null


# ----------------------------------------------------------------------------- #
# Manhattan plot
# ----------------------------------------------------------------------------- #

#GEMMA linear regression with correction
pdf(file="/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/0.Writing/1.Paper/1.Figures/Figure3/gwas_on_ancestry.pdf", 
    bg = "transparent", width=10, height=3, family = "Times New Roman")

idx <- which(lmm_gemma_anc$bon == max(lmm_gemma_anc$bon[which(lmm_gemma_anc$bon < 0.05)]))
idx2 <- which(lmm_gemma_anc$FDR == max(lmm_gemma_anc$FDR[which(lmm_gemma_anc$FDR < 0.05)]))

#colors
lmm_gemma_anc$cols <- "gray"
lmm_gemma_anc$cols[lmm_gemma_anc$FDR < 0.05] <- "red" #893
p1<- lmm_gemma_anc %>%
  ggplot(aes(pos/1000000,-log10(p_wald))) +
  geom_point(alpha=.5, col = lmm_gemma_anc$cols) +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  #facet_wrap(~ chrom, scales = "free_x")+
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,10) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = -log10(lmm_gemma_anc$p_wald[idx]),lty="dashed",lwd=.5, col = "#929292") +
  geom_hline(yintercept = -log10(lmm_gemma_anc$p_wald[idx2]),lty="dashed",lwd=.5, col = "#929292") 


p1 

dev.off()


#GEMMA linear regression with correction
pdf(file="/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/8.ancestry_hmm/gwas/gemma_cutoff09_miss02/gemma_geno_cutoff09_miss02.pdf", 
    bg = "transparent", width=8, height=4, family = "Times New Roman")

idx <- which(lmm_gemma_geno$bon == max(lmm_gemma_geno$bon[which(lmm_gemma_geno$bon < 0.05)]))
idx2 <- which(lmm_gemma_geno$FDR == max(lmm_gemma_geno$FDR[which(lmm_gemma_geno$FDR < 0.05)]))
p1<- lmm_gemma_geno %>%
  ggplot(aes(pos/1000000,-log10(p_wald))) +
  geom_point(alpha=.5, col = "steelblue") +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,10) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = -log10(lmm_gemma_geno$p_wald[idx]),lty="dashed",lwd=.5) +
  geom_hline(yintercept = -log10(lmm_gemma_geno$p_wald[idx2]),lty="dashed",lwd=.5, col = "gray") 


p1 

dev.off()

dim(lmm_gemma)
