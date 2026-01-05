#!/usr/bin/env Rscript

#goal: to calculate the lengths of ancestry blocks along the genome

setwd("/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/randomization/")
ncol <- 9
nsamp <- 280 #nb of samples
#loop through files
vcfs <- Sys.glob("*.vcf")

for (f in 1:length(vcfs)) {

        #loop through samples
        anc <- read.table(vcfs[f], sep = "\t", header = F)      
        #for saving output file
        out <- sub('\\.vcf$', '', vcfs[f])
        
        for (samp in 1:nsamp) {
        
        #filter out missing calls
        
                anc_clean <- anc[anc[,ncol+samp] != "./.",]
                        
                if (dim(anc_clean)[1] > 0) { #if not empty 

                #store ancestry in two matrix to compare them
                chrom_len <- dim(anc_clean)[1]
                K1 <- anc_clean[1:(chrom_len-1),ncol+samp]
                K2 <- anc_clean[2:(chrom_len),ncol+samp]
                idx <- which(K1!=K2)

                #calculate length of blocks within chromosome for each sample
                if (length(idx) > 0) { #if not empty
                        ind_info <- matrix(NA,(length(idx)-1),4)
                        for (i in 2:length(idx)) {
                
                                start <- anc_clean[idx[i-1],2]
                                end <- anc_clean[idx[i],2]
                                L <- end - start

                                #record
                                ind_info[i-1,1] <- start
                                ind_info[i-1,2] <- end
                                ind_info[i-1,3] <- L
                                ind_info[i-1,4] <- samp
        
        
                                }
                        }
                }

        write.table(ind_info, paste("ancestry_block_info_append_chr", out, ".txt", sep = ""),, sep = "\t", col.names = F, row.names = F, quote =F, append = T)

        }
}
