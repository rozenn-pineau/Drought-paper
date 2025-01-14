library(tidyverse)
library(data.table)
rm(list= ls())

# values is a matrix with the ancestry information for each site (rows) and each individual (columns)
setwd("/scratch/midway3/rozennpineau/drought/ancestry_hmm/manual_gwas/two_pulse_flexible_prop_2")
values <- read.table("/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/two_pulse_flexible_prop_2_values.txt", sep = "\t", header = T)
#rm outliers
values <- values[,-c(which(colnames(values) == "P16_Nat_1_T" | colnames(values) == "P12_Nat_14_T"))] #282 samples
#check
dim(values)

samp_order <- read.table("../samp_order_noout.txt",  sep = "\t", header = F )

#keep chrom and pos info aside
pos_chrom <- values[,c(1,2)]
vals <- values[,-c(1,2)]

#turn genotypes into 0/0.5/1
vals <- vals/2

#upload survival info
phenos <- read.table("../drought_phenos.txt", sep = "\t", header = F)

#order phenos based on genotypes
phenos_ordered <- phenos[order(match(phenos$V1,samp_order$V1)),]

#upload covariate (PC1)
pc <- read.table("../PCA_output_noout.txt", sep = "\t", header = T)
#order pc based on samp_order
pc_ordered <- pc[order(match(pc$samp,samp_order$V1)),]


#rm file if already present
filename <- "logistic_covariate_regression_results.txt"
if (file.exists(filename)) {file.remove(filename)}
#make header
export <- c("chrom", "pos", "pval", "pvalcov", "zval", "zvalcov", "slope", "slopecov", "slope_err", "slope_err_cov")
write(export, filename,  sep = "\t", append = TRUE, ncolumns = 10)

#run logistic regression
for (i in 1:dim(vals)[1]) {

  #make dataset
  test <- data.frame(anc_value = t(vals[i,]),
                     survival = phenos_ordered$V2,
                     cov = pc_ordered$PC1)
  test <- test[complete.cases(test),]
  if (nrow(test) > 0 )
    {lm1 <- glm(data=test, test[,1] ~ test[,2] + test[,3])
  tmp <- as.data.frame(summary(lm1)$coefficients)

  pval <- tmp$`Pr(>|t|)`[2] #pvalue for phenotype
  pvalcov <- tmp$`Pr(>|t|)`[3] #pvalue for covariate
  zval <- tmp$`t value`[2] #test stat value for phenotype
  zvalcov <- tmp$`t value`[3] #test stat for covariate
  slope <- tmp$Estimate[2]
  slopecov <- tmp$Estimate
  slope_error <- tmp$`Std. Error`[2]
  slope_error_cov <- tmp$`Std. Error`[3]
  b0 <- tmp$Estimate[1]
  export <- c(pos_chrom[i,1], pos_chrom[i,2], pval, pvalcov, zval, zvalcov, slope, slopecov, slope_error, slope_error_cov)
  #save as I go
  write(export, filename,  sep = "\t", append = TRUE, ncolumns = 10)
  print(i)
    }
}
