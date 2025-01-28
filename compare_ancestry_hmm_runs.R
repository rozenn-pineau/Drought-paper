#!/usr/bin/env Rscript 

#goal : to compare two independent runs from ancestry_hmm and see if they give us similar results

#load data
setwd("/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop")
file1 <- read.table("chrom9", sep = "\t")
summary(file1)
setwd("/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2")
file2 <- read.table("chrom9", sep = "\t")
summary(file2)

#rm chrom pos info
file1 <- file1[,-c(1,2)]
file2 <- file2[,-c(1,2)]


setwd("/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/cmp_runs/")
#I am afraid this is going to be a lot and we won't see much
#randomly sample 9 individuals and 100 sites and plot 3x3
inds <- sample(1:282,9)
sites <- sample(1:dim(file1)[1],100)

pdf("chrom9_test.pdf",
    bg = "white", width=10, height=10)

par(mfrow=c(3,3))
for (i in 1:9){
  #from table to vector
  vec1 <- as.vector(t(as.matrix(file1[sites,inds[i]])))
  vec2 <- as.vector(t(as.matrix(file2[sites,inds[i]])))
  plot(vec1, vec2, xlab = "run 1", ylab = "run 2", pch = 16)
  
}

dev.off()

setwd("/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/cmp_runs/")
#trying to get a better overview of the data by looking at the difference between each pair of samples
vec <- c()
sites <- sample(1:dim(file1)[1],1000)
for (x1 in 1:dim(file1)[2]) {
  
  vec <- c(vec, file1[sites,x1]-file2[sites,x1])
  write.table(vec, "chrom9_distance.txt", sep = "\t", row.names = F, col.names = F, append = T)   
  
}

#plot
pdf("chrom9_hist.pdf",
    bg = "white", width=10, height=10)

hist(vec, 6, xlab = "difference between runs", main = "chrom 9")

dev.off()
