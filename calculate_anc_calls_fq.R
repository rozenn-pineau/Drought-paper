#!/usr/bin/env Rscript 
  
library(data.table)
library(reshape2)

setwd("/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2")
anc <- read.table("two_pulse_flexible_prop_2_values.txt", sep = "\t", header = T)
anc_calls <- anc[,-c(1,2)]
samp_order <- colnames(anc)
idx <- c( which(samp_order == "P16_Nat_1_T") , which(samp_order == "P12_Nat_14_T") )
anc_calls <- anc_calls[,-idx] # 786261 x 280

#upload survival matrix
survival <- read.table("sample_survival_drought.txt", sep = "\t", header = F)

#calculate allele frequency
n_snp <- dim(anc_calls)[1]
n_samp <- dim(anc_calls)[2]
n_day <- 20
daily_af_all <- matrix(NA, n_day, n_snp)

for (day in 1:n_day) {

  for (mut in 1:n_snp) {

    af <- c()
    for (samp in 1:n_samp) {

      #if alive, count in AF
      if (survival[day,samp] == 1 & !is.na(anc_calls[mut,samp])) {

        if (anc_calls[mut,samp] == 0) {af <- c(af, 0)}
        if (anc_calls[mut,samp] == 1) {af <- c(af, 1)}
        if (anc_calls[mut,samp] == 2) {af <- c(af, 2)}

      }


    }

    daily_af_all[day, mut] <- sum(af)/(length(af)*2)

  }

}
dim(daily_af_all) #20 days times 35 SNPs

#store info
export <- cbind(anc[,c(1,2)], t(daily_af_all))
write.table(export, "anc_call_frequency_all.txt", sep = "\t", col.names = F, row.names = F )
