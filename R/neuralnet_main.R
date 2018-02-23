# neuralnet_proteins.R
# Predict protein type based on gene ontology mappings
# Neural Computing Applications
# started : 19/1/2018
# completed:

library(ROCR)
library(kernlab)
library(randomForest)
library("e1071")
library(caret)
library(sand)
library(igraph)
library(GO.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(neuralnet)


memory.limit(1510241024*1024) # allocate RAM memory (15 GBs)
setwd("C:/R-files/NeuralNet")  # now point to where the new code lives
load("NCA-February23rd2018.RData")
source("neuralnet_proteins_functions.R")  # load in the functions required for this work. 

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


# Train SVM on data using CARET
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(3233)

svm_linear <- train(as.factor(ytrain) ~., data = (xtrain), method = "svmLinear",
                    trControl=trctrl,
                    #preProcess = c("center", "scale"),
                    tuneLength = 10)

# Train SVM on data using e1071
svm_model <- svm(as.factor(ytrain)~ ., data=xtrain)
#training set predictions
pred_train <-predict(svm_model,xtrain)
mean(pred_train==ytrain)

#test set predictions
pred_test <-predict(svm_model,xtest)
mean(pred_test==ytest)

# Train MLP on data, need to arrange outputs differently for targets and nontargets
nnet_train <- xtrain[,1:150]
nnet_train <- cbind(nnet_train, nnet_train$target == "1")
nnet_train <- cbind(nnet_train, nnet_train$target == "0")
names(nnet_train)[151] <- 'target'
names(nnet_train)[152] <- 'nontarget'

coln <- colnames(nnet_train[1:149]) # columns' name
a <- as.formula(paste('target + nontarget ~ ' ,paste(coln,collapse='+')))
mlp_model <- neuralnet(a, nnet_train,lifesign="full",hidden=c(30))
         
# Now predict MLP on test data
mlp_predict <- compute(mlp_model, xtest[-5])$net.result
# Put multiple binary output to categorical output
maxidx <- function(arr) {return(which(arr == max(arr))) }

idx <- apply(mlp_predict, c(1), maxidx)
prediction <- c('target', 'nontarget')[idx]
table(prediction, ytest)


# pretty plots for paper
bar_plot_gg2(drug_targets,1,"red")  # plot all target proteins
bar_plot_gg2(hubtargetlist,2,"blue")  # plot target



