#Goal : how close are global ancestry values as calculated by Structure to the mean ancestry as calculated by ancestry_hmm ?


library(extrafont)

# values is a matrix with the ancestry information for each site (rows) and each individual (columns)
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/8.ancestry_hmm/1.results/cutoff/")
values <- read.table("two_pulse_flexible_prop_2_values.txt", sep = "\t", header = T)
#rm outliers
values_fil <- values[,-c(which(colnames(values) == "P16_Nat_1_T" | colnames(values) == "P12_Nat_14_T"))] #280 samples
rm(values)

#testing on shorter dataset
values_fil <- values_fil[1:1000,] #when testing

#calculate mean ancestry per individual
values_fil <- values_fil/2
admx <- colMeans(values_fil[,-c(1,2)], na.rm=T)
samp_order_adm <- colnames(values_fil[-c(1,2)])
#make dataset
admx_dt <- data.frame(samp = samp_order_adm, admx = admx)

#upload global ancestry
k2 <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/2.admixture/k2_structure_res.csv", sep = ",", header = T)
#remove the two outliers
k2_fil <- k2[-c(which(k2$samp == "P16_Nat_1_T" | k2$samp == "P12_Nat_14_T")),] #280 samples


#order samples 
full <- merge(k2, admx_dt, by = "samp")

#compare

pdf(file="/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/2.admixture/k2_vs_admx.pdf", 
    bg = "transparent", width=5, height=5, family = "Times New Roman")

par(family = "Times New Roman", cex.lab = 1.5, cex.axis = 1.5)
plot(full$k2.V1, full$admx, pch =16, 
     xlab = "Global ancestry (ADMIXTURE)", ylab = "Mean ancestry (ancestry_hmm)")


dev.off()
