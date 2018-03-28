# neuralnet_proteins.R
# Predict protein type based on gene ontology mappings
# for Neural Computing Applications or possibly BMC Bioinformatics
# started : 19/1/2018
# completed:


memory.limit(1510241024*1024) # allocate RAM memory (15 GBs)
setwd("C:/R-files/NeuralNet")  # now point to where the new code lives
load("NCA-27thMarch2018.RData")
#load("matrixdata.RData") # contains mmt and mcrap
source("neuralnet_functions.R")  # load in the functions required for this work. 
source("neuralnet_data.R")  # load in raw data, preprocess it. 


### Only uncomment this code if you intend building training matrix (mm & mt) from scratch - it takes a long time!
#pnames <- unique(c(unique(ppi$Gene_A),unique(ppi$Gene_B)))
#mm <- go_slim_annotation(pnames)
#mt <- give_classlabels_mm(mm) 

## Convert matrix to dataframe and balance out the data by undersampling
mnew <- data.table::transpose(as.data.frame(mt))
mnew <- data.frame(mnew)
colnames(mnew) <- rownames(mt)
rownames(mnew) <- colnames(mt)
positives <- mnew[mnew$targets == 1,]  # get all targets (1,449)
negatives <- mnew[mnew$targets == 0,]  # get all nontargets (11,567)
allnegatives <- data.frame(negatives)
negatives <- sample_n(negatives, nrow(positives)) # only use 1,449 of them to match positives
balanced_dat <- rbind(positives,negatives)


## Prepare a training and a test set 
ntrain <- round(nrow(bdata)*0.8) # number of training examples
tindex <- sample(nrow(bdata),ntrain) # indices of training samples
xtrain <- data.frame(bdata[tindex,])
xtest <-  data.frame(bdata[-tindex,])
ytrain <- xtrain[,150]  # class labels for training data (column 150 is class label)
ytest <- xtest[,150]   # class labels for test data
ytest <- as.factor(ytest)

## Train SVM on data using e1071 package
svm_model <- e1071::svm(as.factor(ytrain)~ ., data=xtrain[,1:149],scale=FALSE,cross=20)
#test set predictions
pred_test <-predict(svm_model,xtest[,1:149])

pred_svm <- ROCR::prediction(messyfactor2int(pred_test),ytest)
svm.roc <- ROCR::performance(pred_svm, "tpr", "fpr")
svm.pr <- ROCR::performance(pred_svm, "prec", "rec")

acc <- table(pred_test, ytest)
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\nSVM accuracy calculated by (TP+TN)/(TP+TN+FP+FN)= ",(TP + TN)/(TP + TN + FP + FN))
caret::confusionMatrix(data=pred_test,reference=ytest,positive="1")
auc.perf = ROCR::performance(pred_svm, measure = "auc")
auc.perf@y.values


# Train RBF on data using e1071 package
rbf_model <- e1071::svm(as.factor(ytrain)~., data=xtrain[,1:149],scale=FALSE,
                        kernel="radial", gamma=1, cost =1, cross=10)
pred_rbf <-predict(rbf_model,xtest[,1:149])
pred <- ROCR::prediction(messyfactor2int(pred_rbf),ytest)
rbf.roc <- ROCR::performance(pred, "tpr", "fpr")
rbf.pr <- ROCR::performance(pred, "prec", "rec")

acc <- table(pred_rbf, ytest)
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\nRBF accuracy calculated by (TP+TN)/(TP+TN+FP+FN)= ",(TP + TN)/(TP + TN + FP + FN))
caret::confusionMatrix(data=pred_rbf,reference=ytest,positive="1")
auc.perf = ROCR::performance(pred, measure = "auc")
auc.perf@y.values

# Train Random Forest on data 
rf_model <-randomForest(as.factor(ytrain) ~.,data=xtrain[,1:149],proximity=TRUE,keep.forest=TRUE,ntree=500,mtry=12)
predicted_rf <- predict(rf_model,newdata=xtest[,1:149],type = "prob")  
pred_rf <- ROCR::prediction((predicted_rf[,2]),ytest)
auc.perf = ROCR::performance(pred_rf, measure = "auc")
auc.perf@y.values

rf.roc <- ROCR::performance(pred_rf, "tpr", "fpr")
rf.pr <- ROCR::performance(pred_rf, "prec", "rec")
x <- rf.roc@x.values
y <- rf.roc@y.values
pred_rf <- stats::predict(rf_model,xtest)
caret::confusionMatrix(data=pred_rf,reference=ytest,positive="1")

order(importance(rf_model), 2)

acc<-table((predicted_rf[,2]), ytest)
print(rf_model)
#round(importance(rf_model), 2)
#varImpPlot(rf_model,main="",type=2,color="black",pch=16) 
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\nRandom Forest accuracy calculated by (TP+TN)/(TP+TN+FP+FN)= ",(TP + TN)/(TP + TN + FP + FN))


# Train MLP on data, need to arrange outputs differently for targets and nontargets
# using neuralnet package.
nnet_train <- xtrain[,1:150]
nnet_train <- cbind(nnet_train, nnet_train$target == "1")
nnet_train <- cbind(nnet_train, nnet_train$target == "0")
names(nnet_train)[151] <- 'target'
names(nnet_train)[152] <- 'nontarget'

coln <- colnames(nnet_train[1:149]) # columns' name
a <- as.formula(paste('target + nontarget ~ ' ,paste(coln,collapse='+')))
mlp_model <- neuralnet::neuralnet(a,nnet_train,lifesign="full",act.fct="logistic",stepmax = 20000,
                       algorithm = "rprop+",  threshold = 0.01,
                       hidden=c(80,20),linear.output=TRUE)  # train the MLP
         
# Now predict MLP on test data
mlp_predict <- neuralnet::compute(mlp_model, xtest[,1:149])$net.result
# Put multiple binary output to categorical output
maxidx <- function(arr) {return(which(arr == max(arr))) }
idx <- apply(mlp_predict, c(1), maxidx)
prediction <- c('target', 'nontarget')[idx]
acc <- table(prediction, ytest)
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\nMLP accuracy calculated by (TP+TN)/(TP+TN+FP+FN)= ",(TP + TN)/(TP + TN + FP + FN))
caret::confusionMatrix(data=mlp_predict,reference=ytest,positive="1")

prediction[prediction =="target"] <- 1
prediction[prediction =="nontarget"] <- 0
prediction <- factor2int(as.integer(prediction))
pred_mlp <- ROCR::prediction(prediction,ytest)
auc.perf = ROCR::performance(pred_mlp, measure = "auc")
auc.perf@y.values



# Now try nnet package in conjunction with caret - use 10-fold CV
train_control <- trainControl(method="repeatedcv", number = 10,repeats = 5, 
                           classProbs=TRUE,summaryFunction = twoClassSummary)
my_grid <- expand.grid(size = seq(from=80, to=120, by=10), decay=seq(from=0.1, to=0.5, by=0.1))
# The data needs to in yet another different format!
ytrain_nnet <- ytrain
ytrain_nnet[ytrain_nnet == 0]  <- "nontarget"
ytrain_nnet[ytrain_nnet == 1]  <- "target"
ytest_nnet <- as.numeric(ytest)
ytest_nnet[ytest_nnet == 2]  <- "nontarget"
ytest_nnet[ytest_nnet == 1]  <- "target"

# Use nnet package - 
nnet_model <- nnet(as.factor(ytrain)~., data=xtrain[,1:149],size=50,MaxNWts=100000,maxit=1000,decay=.001)
acc <- table(predict(nnet_model, xtest[,1:149], type = "class"),ytest)
nnet_pred <- predict(nnet_model, xtest[,1:149], type = "class")
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\nMLP accuracy calculated by (TP+TN)/(TP+TN+FP+FN)= ",(TP + TN)/(TP + TN + FP + FN))
caret::confusionMatrix(data=nnet_pred,reference=ytest,positive="1")

pred_mlp <- ROCR::prediction(factor2int(nnet_pred),ytest)
mlp.roc <- ROCR::performance(pred_mlp, "tpr", "fpr")
mlp.pr <- ROCR::performance(pred_mlp, "prec", "rec")

auc.perf = ROCR::performance(pred_mlp, measure = "auc")
auc.perf@y.values


# ROC plots of classifiers
attributes(mlp.roc)$roc_name <- "MLP"
attributes(rf.roc)$roc_name <- "RandomForest"
attributes(rbf.roc)$roc_name <- "RBF"
attributes(svm.roc)$roc_name <- "SVM"
roc_plot(rbf.roc,mlp.roc,svm.roc,rf.roc)

# PR plots of classifiers
attributes(mlp.pr)$pr_name <- "MLP"
attributes(rf.pr)$pr_name <- "RandomForest"
attributes(rbf.pr)$pr_name <- "RBF"
attributes(svm.pr)$pr_name <- "SVM"
pr_plot(rf.pr,rbf.pr,mlp.pr,svm.pr)

# pretty plots for paper
bar_plot_gg2(drug_targets,1,"red")  # plot all target proteins
bar_plot_gg2(hubtargetlist,2,"blue")  # plot target

########## select proteins that are nontargets but not in train or test set #########
# use the trained classifiers on new data "unknown" for potential targets

allnontargets <- mnew[mnew$targets == 0,]
unknown <- negatives[!rownames(allnontargets) %in% rownames(negatives),]
unknown <- data.frame(unknown)

# shape up new data
uindex <- base::sample(nrow(unknown),10) # indices of training samples
candidates <- unknown[sample(1:nrow(unknown), 1000,replace=FALSE),] 
candidates <- unknown[1:600,1:149]
candidates <- data.frame(candidates)

candidates <- unknown[1:1280,1:149]

# VENN: ensure models output in common format to allow comparisions with Venn diagram however
# a large amount of processing is required for the mlp, since neuralnet package is very different.
candidates_rf <-  stats::predict(rf_model,candidates)
candidates_svm <- stats::predict(svm_model,candidates)
candidates_rbf <- stats::predict(rbf_model,candidates)
candidates_mlp <- neuralnet::compute(mlp_model,candidates)$net.result
  maxidx <- function(arr) {return(which(arr == max(arr))) }
  idx <- apply(candidates_mlp, c(1), maxidx)
  idx <- as.integer(idx)
  prediction <- c('target', 'nontarget')[idx]
  prediction <- gsub('nontarget', 0, prediction)
  prediction <- gsub('target', 1, prediction)

# make a nice well structured data frame  
cnames <- rownames(candidates)
comp_models <- data.frame(candidates=cnames,
                          rf =as.vector(candidates_rf),
                          svm=as.vector(candidates_svm), 
                          rbf=as.vector(candidates_rbf),
                          mlp=as.vector(prediction))

# use data frame to fill classifier lists of what they think are target proteins
p_rf  <- comp_models$candidates[comp_models$rf  ==1]
p_svm <- comp_models$candidates[comp_models$svm ==1]
p_rbf <- comp_models$candidates[comp_models$rbf ==1]
p_mlp <- comp_models$candidates[comp_models$mlp ==1]

Reduce(intersect, list(p_rf,p_svm,p_rbf,p_mlp))  # the proteins identified by all four classifiers

# plot a Venn diagram to highlight commonly identifed protein targets between four classifiers
plot_venn(p_rf,p_svm,p_rbf,p_mlp)


#### LATENT NETWORK STAGE ####
ppi_net <- build_network(ppi)  # create an igraph object.

ppi_latent <- build_latent(ppi_net)  # make the latent model.








