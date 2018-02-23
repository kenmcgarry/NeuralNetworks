# neuralnet_proteins.R
# Predict protein type based on gene ontology mappings
# Neural Computing Applications
# started : 19/1/2018
# completed:

library(ROCR)
library(kernlab)
library("e1071")
library(caret)
library(sand)
library(igraph)
library(GO.db)
library(dplyr)
library(tidyr)
library(ggplot2)

memory.limit(1510241024*1024) # allocate RAM memory (15 GBs)
setwd("C:/R-files/neuralnet")  # now point to where the new code lives
load("complexnets_23rdFebruary2018.RData")
source("neuralnet_protein_functions.R")  # load in the functions required for this work. 

# restore mcrap data from original source, rather than reuse.
mcrap <- data.table::transpose(as.data.frame(mmt))
mcrap <- data.frame(mcrap)
colnames(mcrap) <- rownames(mmt)
rownames(mcrap) <- colnames(mmt)

positives <- mcrap[mcrap$targets == 1,]  # get all targets (1,443)
negatives <- mcrap[mcrap$targets == 0,]  # get all nontargets (11,554)
allnegatives <- data.frame(negatives)
negatives <- sample_n(negatives, nrow(positives)) # only use 1,443 of them to match positives
balanced_dat <- rbind(positives,negatives)

## Prepare a training and a test set ##
ntrain <- round(nrow(balanced_dat)*0.8) # number of training examples
tindex <- sample(nrow(balanced_dat),ntrain) # indices of training samples
xtrain <- data.frame(balanced_dat[tindex,])
xtest <-  data.frame(balanced_dat[-tindex,])

ytrain <- xtrain[,150]  # class labels for training data
ytest <- xtest[,150]   # class labels for test data
ytest <- as.factor(ytest)

# pretty plots for paper
bar_plot_gg2(drug_targets,1,"red")  # plot all target proteins
bar_plot_gg2(hubtargetlist,2,"blue")  # plot target



