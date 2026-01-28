#goal : to create bed files that have a similar distribution of ancestry block lengths as the observed drought bed file
  
#full <- read.table("../full_genome_anc_block_mean_len_tab_subset.bed", sep = "\t", header = F) #change to header = T with full file
full <- read.table("../ancestry_block_mean_786257.bed", sep = "\t", header = T)
obs <- read.table("observed_anc_block_len_distr.txt", sep = "\t", header = T) #plus or minus 100 bp
precision <- 500
precision <- 1000
precision <- 5000

for (perm in 1:100) {
        store_bed <- c()
        for (i in 1:length(obs$mid)) {
                idx_pool <- which(full[,4] > obs$mid[i]-precision & full[,4] < obs$mid[i]+precision)
                if (length(idx_pool)<1){break} #no matches
                idx_local <- sample(idx_pool,obs$count[i],replace=F)
                cur_bed <- full[idx_local,]
                store_bed <- rbind(store_bed, cur_bed)
                write.table(store_bed, paste("random_anc_block_", perm,".bed", sep = ""), sep = "\t", col.names = F, row.names = F, quote = F )
        }
}
