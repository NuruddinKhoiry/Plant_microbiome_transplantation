# This R script is written by Ahmad Nuruddin Khoiri 
# Email: nuruddinkhoiri34@gmail.com

# This R script is a part of the manuscript entitled 
# "Co-transplantation of phyllosphere and rhizosphere microbes promotes microbial colonization and enhances sugarcane growth"

setwd("~/DiskE/GitHub/Plant_microbiome_transplantation/")

# Library
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)


# R Function

# phy: phyloseq object
# group_column: The metadata column containing groups to compare.
# source_code: The group code for the microbiome source. It must be found in "group_column."
# control_code: The group code for the control in the experiment. It must be found in "group_column."
# sink_code: The group code for the treatment as the sink in the experiment. It must be found in "group_column."
# SampleID: The metadata column containing sample ID or code.
# pct_replicate: Filter taxa present in fewer than x% of replicates.
# alpha: The prior counts to be added to all of the sample counts. (Please refer to Alpha1 on https://github.com/biota/SourceTracker_rc/blob/master/ipynb/Sourcetracking%20using%20a%20Gibbs%20Sampler.ipynb)
# prob_cut_off: The minimum probability score to determine if the taxon likely originated from the source instead of the surrounding environment (control).
successive_taxa <- function(phy, group_column, source_code, control_code, sink_code, SampleID,
                            pct_replicate=50, alpha=0.001, prob_cut_off=0.9) {
  
  map <- data.frame(sample_data(phy))
  map[, group_column] <- as.character(map[, group_column])
  otu <- as.matrix(otu_table(phy))
  otu <- otu[, rownames(map)]
  
  # Filter taxa by their prevalence, i.e., present in at least x% of replicates
  mydata <- list()
  for (i in unique(as.character(map[, group_column]))){
    G1 <- otu[ ,map[, group_column] == i]
    G1[rowSums(G1 > 0) < (pct_replicate/100)*dim(G1)[2], ] <- 0
    mydata[[i]] <- G1
  }
  
  names(mydata) <- rep('', length(names(mydata)))
  df <- cbind.data.frame(mydata)
  df <- df[ ,rownames(map)]
  df_filt <- df[rowSums(df) > 0, ]
  
  # Compare control vs source
  map2 <- map[map[,group_column] %in% c(control_code, source_code), c(SampleID, group_column)]
  df_filt2 <- df_filt[ ,rownames(map2)]
  
  mydata2 <- list()
  for (i in c(control_code, source_code)){
    G1 <- apply(df_filt2[ ,map2[, group_column] == i], 1, function(x) mean(x)) %>% data.frame()
    colnames(G1) <- i
    mydata2[[i]] <- G1
  }
  
  names(mydata2) <- rep('', length(names(mydata2)))
  df2 <- cbind.data.frame(mydata2)
  df2 <- df2[rowSums(df2) > 0, ]
  
  # Add prior counts to all of the sample counts and then estimate probability score
  df2 <- df2 + alpha
  df3 <- df2 %>% apply(2, function(x) x/sum(x)) %>% apply(1, function(x) x/sum(x))
  df4 <- df3 %>% t()
  df5 <- df4[df4[,2] >= prob_cut_off, ]
  
  # Estimate sink
  mydata3 <- list()
  for (i in unique(as.character(map[, group_column]))){
    G1 <- df_filt[, rownames(map[map[,group_column] == i, c(SampleID, group_column)])]
    mydata3[[i]] <- data.frame(1*(rowSums(G1)>0))
    colnames(mydata3[[i]]) <- i
  }
  
  names(mydata3) <- rep('', length(names(mydata3)))
  venndata_df <- cbind.data.frame(mydata3)
  
  df6 <- venndata_df[,c(source_code, sink_code)]
  df7 <- df6[rowSums(df6) == 2, ]
  
  df8 <- df7[rownames(df7) %in% rownames(df5), ]
  
  df9 <- as.data.frame(df5) %>% rownames_to_column('Taxa') %>%
    select(Taxa, all_of(source_code))
  
  colnames(df9) <- c('Taxa', 'Source_prob')
  
  df10 <- df8 %>% rownames_to_column('Taxa') %>%
    left_join(df9, 'Taxa') %>%
    left_join(as.data.frame(tax_table(phy)@.Data) %>% rownames_to_column('Taxa'), 'Taxa')

  return(df10)
}

# Load samples
load('input_phyloseq_objects.RData')

# Identify successfully transplanted taxa
# 1.A Identify successfully transplanted taxa from P treatment on the phyllosphere niche
taxaP_P = successive_taxa(phy = phyllosphere, 
                          group_column = 'MT',
                          SampleID = 'SampleID', 
                          sink_code = 'P',
                          source_code = 'Farm', 
                          control_code = 'C',
                          pct_replicate=50,
                          alpha = 0.001,
                          prob_cut_off = .9
                          )

taxaP_P

# 1.B Identify successfully transplanted taxa from PR treatment on the phyllosphere niche
taxaP_PR = successive_taxa(phy = phyllosphere,
                           group_column = 'MT',
                           SampleID = 'SampleID',
                           sink_code = 'PR',
                           source_code = 'Farm',
                           control_code = 'C',
                           pct_replicate=50,
                           alpha = 0.001,
                           prob_cut_off = .9
                           )

taxaP_PR


# 2.A Identify successfully transplanted taxa from R treatment on the rhizosphere niche
taxaR_R = successive_taxa(phy = rhizosphere,
                          group_column = 'MT', 
                          SampleID = 'SampleID', 
                          sink_code = 'R',
                          source_code = 'Farm', 
                          control_code = 'C',
                          pct_replicate=50,
                          alpha = 0.001,
                          prob_cut_off = .9
                          )

taxaR_R

# 2.B Identify successfully transplanted taxa from PR treatment on the rhizosphere niche
taxaR_PR = successive_taxa(phy = rhizosphere,
                           group_column = 'MT', 
                           SampleID = 'SampleID',
                           sink_code = 'PR',
                           source_code = 'Farm', 
                           control_code = 'C',
                           pct_replicate=50,
                           alpha = 0.001,
                           prob_cut_off = .9
                           )

taxaR_PR
