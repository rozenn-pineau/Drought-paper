#!/usr/bin/env Rscript

len <- read.table("ancestry_block_len_genome.txt", sep = " ", header = T)
#bed <- read.table("../gwas_893.bed", header = T)
bed <- read.table("two_pulse_flexible_prop_2_values_ID_header.bed", header = T)

nsamp <- 260
nlocus <- 786257
store <- matrix(NA, nlocus, nsamp)
for (samp in 246:nsamp) {
#per sample
  print(samp)

  for (locus in 1:nlocus) {
   #per locus
      #print(locus)

      chr <- bed$chrom[locus]
      pos <- bed$pos1[locus]

      subfile <- len[len$chr==chr & len$samp == samp,] #subset file by chromosome and sample
      #find location in position:

      #is subfile empty
      if (dim(subfile)[1] >0) {
              #is the block within boundaries
              if ( subfile$start[1] <= pos & subfile$end[dim(subfile)[1]] >= pos) {

                      #find block

                      block <- which(subfile$start >= pos)[1] - 1 #look at the block below the first block that is above the SNP start position

                      #if locus is at the very last spot

                      if (is.na(block) &  subfile$end[dim(subfile)[1]] >= pos) {store[locus,samp] <- subfile$len[dim(subfile)[1]]}

                      #if locus is below the first block
                      else if (block == 0){ store[locus, samp] <- NA }

                      #if locus is within the block
                      else if (subfile$start[block] <= pos & subfile$end[block] >= pos){ #if this is true, record block length

                        store[locus, samp] <- subfile$len[block]

                        }
                }
        }
   }

  write.table(t(store[,samp]), "ancestry_block_len_282x786257_append_246-260.txt", sep = "\t", row.names = F, col.names = F, quote = F, append = T)

}
