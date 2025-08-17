#!/usr/bin/env Rscript 

library(dplyr)
library(extrafont)
library(geosphere)
library(data.table)
library(car)
library(scales)
library(paletteer)
library(lsmeans)

tmax <- read.table("tmax_NOAA.txt", sep = "\t", header = T)
n_snp <- dim(tmax[,-c(1:6)])[2]

# Simple function to find closest pairs
find_closest_pairs <- function(data,itr) {

  # Remove samples with missing coordinates
  clean_data <- data %>%
    filter(!is.na(samp_lat), !is.na(samp_lon))

  # Calculate all pairwise distances
  pairs <- expand.grid(i = 1:nrow(clean_data), j = 1:nrow(clean_data)) %>%
    filter(i < j) %>%  # Only unique pairs
    mutate(
      sample1 = clean_data$sample[i],
      sample2 = clean_data$sample[j],
      year1 = clean_data$year[i],
      year2 = clean_data$year[j],
      year_diff = abs(clean_data$year[i] - clean_data$year[j]),
      distance_km = distVincentyEllipsoid(
        cbind(clean_data$samp_lon[i], clean_data$samp_lat[i]),
        cbind(clean_data$samp_lon[j], clean_data$samp_lat[j])
      ) / 1000
    ) %>%
    arrange(distance_km) #order by ascending distance

  # Select non-overlapping pairs (greedy) - get maximum possible
  # Exclude pairs from same year and beyond distance threshold
  pairs <- pairs %>% filter(year_diff > 0, distance_km <= 100) #146 options for pairs

  #re-order table (which pair gets chosen first)
  #itr <- itr + 1
  pairs <- pairs[c(itr:146,1:(itr-1)),]

  selected <- data.frame()
  used_samples <- c(0)

  for(i in 1:nrow(pairs)) {
    if(!pairs$i[i] %in% used_samples &
       !pairs$j[i] %in% used_samples) {
      selected <- rbind(selected, pairs[i,])
      used_samples <- c(used_samples, pairs$i[i], pairs$j[i])
    }
  }

  return(selected)
}
n_perm <- 1
n_ind <- dim(tmax)[1]

for (p in 1:n_perm) {

      # Loop through every possible first pair
    store <- matrix(NA,1,146)
    #tmax$tmax_rdn <- tmax$tmax[sample(1:n_ind)]
    count <- 1
    for (itr in 1:146){ #146 pairs, itr gets + 1 in the loop

        # Find pairs
        herb_minpairwisedistance <- find_closest_pairs(tmax, itr)
        n_pairs <- dim(herb_minpairwisedistance)[1]

        # Calculate AF change and climate change for this set of pairs
        n_pairs <- dim(herb_minpairwisedistance)[1]
        herb_minpairwisedistance$tmax_diff <- NA
        af_change <- matrix(NA, n_pairs, n_snp)
        for (s in 1:n_pairs) {

          #use the index columns (I checked and they correspond to the initial dataset)
          year_diff <- tmax$year[herb_minpairwisedistance$i[s]] - tmax$year[herb_minpairwisedistance$j[s]]

          #if the difference in years is negative, do j - i 
          if ( year_diff < 0 )
            {herb_minpairwisedistance$tmax_diff[s] <- tmax$tmax[herb_minpairwisedistance$j[s]] - tmax$tmax[herb_minpairwisedistance$i[s]]
            af_change[s,] <- unlist(tmax[herb_minpairwisedistance$j[s],c(7:39)] - tmax[herb_minpairwisedistance$i[s],c(7:39)])}

          #if the difference in years is positive, do i - j
          else {herb_minpairwisedistance$tmax_diff[s] <- tmax$tmax[herb_minpairwisedistance$i[s]] - tmax$tmax[herb_minpairwisedistance$j[s]]
          af_change[s,] <- unlist(tmax[herb_minpairwisedistance$i[s],c(7:39)] - tmax[herb_minpairwisedistance$j[s],c(7:39)])}

        }
        af_change_tmax <- af_change

        # Calculate linear regression between climate change and haplotype freq change
        sig <- matrix(NA, n_snp, 1)

        for (n in 1:n_snp) {

          #linear regression for each site separately - tmax
          if (length(which(is.na(af_change_tmax[,n]) == TRUE )) < (n_pairs-3)) {
          lm <- lm(af_change_tmax[,n] ~ herb_minpairwisedistance$tmax_diff)

          #store coefficient if relationship is positive
          if ( summary(lm)$coefficients[, 1][2] > 0 ) {
          sig[n] <- summary(lm)$coefficients[, 4][2]
          }
            }
        }

        store[count] <- length(which(sig <= 0.05))
        count <- count + 1
    }

    write.table(store, "herb_data_pos_corr.txt", sep = "\t", col.names = F, row.names = F, quote = F, append = T)
    print(p)
}


#output file is one line per 146 options for smallest distance, each column is the number of significant relationships obtained after
#randomizing the climate variable
