#!/usr/bin/env Rscript

  
find_site_up <- function(initial_anc, anc, ncol, samp,i) { # create a function 

        for (counter in 1:(i-1)) {
                #we changed chromosome
                if (anc$V1[i] != anc$V1[i-counter]) {return(counter);break}
                #we did not change chromosome
                if (anc[i-counter,ncol+samp] !=  initial_anc & anc$V1[i] == anc$V1[i-counter]) {return(counter);break}
        }
                #if end of chrom and still same ancestry
                if (anc[i-counter,ncol+samp] ==  initial_anc & anc$V1[i] == anc$V1[i-counter]){return(counter)}
}

find_site_down <- function(initial_anc, anc, ncol, samp,i) { # create a function

        #max position for this chromosome
        MAX <- which(anc$V2 == max(anc$V2[anc$V1 == anc$V1[i]]))- i
        for (counter in 1:MAX) {
                #we changed chromosome
                if (anc$V1[i] != anc$V1[i+counter]) {return(counter);break}
                #we did not change chromosome
                if (anc[i+counter,ncol+samp] !=  initial_anc & anc$V1[i] == anc$V1[i+counter]) {return(counter);break}
        }
        #if end of chrom and still same ancestry
        if (anc[i+MAX,ncol+samp] ==  initial_anc & anc$V1[i] == anc$V1[i+counter]){return(counter)}
}


#goal: to plot the distribution of ancestry block lengths around the drought adapted loci
setwd("/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought")
#anc <- read.table("anc_subset.vcf", sep = "\t", header = F)
anc <- read.table("anc_no_header.vcf", sep = "\t", header = F)
loci <- read.table("gwas_893.bed", sep = "\t", header = T)
#sort loci by chromosome
loci <- loci[order(loci$chrom),]

#which(anc$V1 %in% loci$chrom)

nind <- 280 #nb of samples
ncol <- 9 #where the ancestry calls start
locus <- 1
#find a locus that is in the ancestry file
locus_idx <- which(anc$V1 %in% loci$chrom & anc$V2 %in% loci$pos1)
nlocus <- length(locus_idx)
total_dist <- matrix(NA, nlocus, nind)

#this loop is for drought-adpated loci
for (i in locus_idx) {
                
                print(i)
                #this loop is for samples
                dist <- matrix(NA,1,nind)
                dist_above <- matrix(NA,1,nind)
                dist_below <- matrix(NA,1,nind)
                for (samp in 1:nind){
                        
                        print(samp)
                        initial_anc <- anc[i,ncol+samp] #this is the initial ancestry
                        
                        #go up
                        idx <- find_site_up(initial_anc, anc, ncol, samp, i)#idx of where change happens
                        L <- anc$V2[i] - anc$V2[i-idx]
                        dist_above[samp] <- L

                        #go down
                        idx <- find_site_down(initial_anc, anc, ncol, samp, i)#idx of where change happens
                        L <- anc$V2[i+idx] - anc$V2[i]
                        dist_below[samp] <- L

                        #total distance
                        dist <- dist_above + dist_below



                }
        
        total_dist[locus,] <- dist
        write.table(dist, "distance_distribution_append.txt", sep = "\t", col.names = F, row.names = F, quote =F, append = T)
        locus <- locus + 1
        
}

write.table(total_dist, "distance_distribution.txt", sep = "\t", col.names = F, row.names = F, quote =F)
