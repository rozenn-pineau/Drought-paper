#!/usr/bin/env Rscript 
  
#goal : to plot the ancestry per site per chromosome per individual
library(reshape2) #note ! use module load R/3.6.3 to be able to use reshape on midway
library(ggplot2)
library(plyr)
library(tidyverse)

#plot colors
tub_col <- "#76528BFF"
rud_col <- "#CBCE91FF"
het_col <- "#44A3BB"

# values is a matrix with the ancestry information for each site (rows) and each individual (columns)
setwd("/scratch/midway3/rozennpineau/drought/ancestry_hmm/manual_gwas/two_pulse_flexible_prop_2")
values <- read.table("/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/two_pulse_flexible_prop_2_values.txt", sep = "\t", header = T)
#rm outliers
values <- values[,-c(which(colnames(values) == "P16_Nat_1_T" | colnames(values) == "P12_Nat_14_T"))] #282 samples

#transform file from large to long format
df_long <- melt(values, id.vars = c( "chrom", "position"))

#upload df
#df_long <- melt(values[1:100,], id.vars = c( "chrom", "position")) #when testing
colnames(df_long) <- c("chrom", "pos", "variable", "GT")

#upload longitude info
long <- read.table("../k2_structure_res.csv", sep = ",", header = T)
long_info <- data.frame(variable = long$samp, long = long$long)
#remove the two outliers
long_info <- long_info[-c(which(long_info$variable == "P16_Nat_1_T" | long_info$variable == "P12_Nat_14_T")),] #280 samples

#subsample (more than 196 000 000 lines, the figure stops rendering at chromosome 3 if I do not subsample)
df_long <- df_long[sort(sample(dim(df_long)[1], 10000000, replace = F )),]

#merge 
dataset <- merge(df_long, long_info, by = "variable" )
#sort samples in values based on long
dataset <- as.data.frame(dataset[order(dataset$long),]) # 280 77011
dataset$GT <- as.character(as.numeric(dataset$GT))


#plot
pdf(file="ancestry_per_chrom_sub_reduced_new_cols.pdf", family = "Times New Roman",
    bg = "transparent", width=10, height=10)


dataset%>% #filter(chr=="1") %>%
  #ggplot(aes(pos, variable, color=GT, fill=GT)) + x = reorder(category, -value), y = value
  ggplot(aes(x = pos, y = reorder(variable, -long), color=GT, fill=GT)) +
  geom_tile() +
  facet_grid(~chrom,scales = "free_x",space = "free_x") +
  #scale_color_viridis_d(direction = -1) +
  #scale_fill_viridis_d(direction = -1) +
  scale_colour_manual(values = c(tub_col, het_col, rud_col)) +
  scale_fill_manual(values = c(tub_col, het_col, rud_col)) +
  theme(
    axis.text.y = element_blank(),  # Adjust size to your preference
    axis.text.x = element_text(size = 6,angle = 90),
    panel.spacing = unit(0, "lines"),           # No space between facets
    panel.border = element_rect(color = "black", fill = NA, size = .5)) + # Black outline
  labs(fill="Ancestry", y="Individual",x="Genomic Position") +
  guides(color="none")

dev.off()
