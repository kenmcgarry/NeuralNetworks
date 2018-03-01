# neuralnet_proteins.R
# Predict protein type based on gene ontology mappings
# Neural Computing Applications
# started : 19/1/2018
# completed:


memory.limit(1510241024*1024) # allocate RAM memory (15 GBs)
setwd("C:/R-files/NeuralNet")  # now point to where the new code lives
load("NCA-1stMarch2018.RData")
load("matrixdata.RData") # contains mmt and mcrap
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



# Train SVM on data using e1071 package
svm_model <- e1071::svm(as.factor(ytrain)~ ., data=xtrain[,1:149],scale=FALSE,cross=20)
#test set predictions
pred_test <-predict(svm_model,xtest[,1:149])

pred_svm <- ROCR::prediction(messyfactor2int(pred_test),ytest)
svm.roc <- ROCR::performance(pred_svm, "tpr", "fpr")
svm.pr <- ROCR::performance(pred_svm, "prec", "rec")

acc <- table(pred_test, ytest)
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\nSVM accuracy calculated by (TP+TN)/(TP+TN+FP+FN)= ",(TP + TN)/(TP + TN + FP + FN))


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


# Train Random Forest on data 
rf_model <-randomForest(as.factor(ytrain) ~.,data=xtrain[,1:149],proximity=TRUE,keep.forest=TRUE)
predicted_rf <- predict(rf_model,newdata=xtest[,1:149],type = "prob")  # predict(m,newdata_matrix,type='prob')

#predicted_rf <- as.vector(predicted_rf)
pred_rf <- ROCR::prediction((predicted_rf[,2]),ytest)
rf.roc <- ROCR::performance(pred_rf, "tpr", "fpr")
rf.pr <- ROCR::performance(pred_rf, "prec", "rec")

x <- rf.roc@x.values
y <- rf.roc@y.values
y <- as.numeric(as.character(unlist(y[[i]])))
x <- as.numeric(as.character(unlist(x[[i]])))

rf.roc@x.values <- x
slot(rf.roc, "x.values") <- x 
plot(rf.roc)
plot(rf.pr)

acc<-table((predicted_rf[,2]), ytest)
print(rf_model)
round(importance(rf_model), 2)
varImpPlot(rf_model,main="",type=2,color="black",pch=16) 
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

pred_mlp <- ROCR::prediction(factor2int(nnet_pred),ytest)
mlp.roc <- ROCR::performance(pred_mlp, "tpr", "fpr")
mlp.pr <- ROCR::performance(pred_mlp, "prec", "rec")

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

allnontargets <- mcrap[mcrap$targets == 0,]
unknown <- negatives[!rownames(allnontargets) %in% rownames(negatives),]
unknown <- data.frame(unknown)

# shape up new data
uindex <- base::sample(nrow(unknown),10) # indices of training samples
candidates <- unknown[sample(1:nrow(unknown), 1000,replace=FALSE),] 
candidates <- unknown[1:600,1:149]
candidates <- data.frame(candidates)

# ensure models output in common format to allow comparisions with Venn diagram however
# a far amount of processing is required with the mlp, since neuralnet package is
# very different.
candidates_rf <-  predict(rf_model,candidates)
candidates_svm <- predict(svm_model,candidates)
candidates_rbf <- predict(rbf_model,candidates)
candidates_mlp <- neuralnet::compute(mlp_model,candidates)$net.result
  maxidx <- function(arr) {return(which(arr == max(arr))) }
  idx <- apply(candidates_mlp, c(1), maxidx)
  idx <- as.integer(idx)
  prediction <- c('target', 'nontarget')[idx]
  prediction <- gsub('nontarget', 0, prediction)
  prediction <- gsub('target', 1, prediction)

  # make a nice well structured data frame  
cnames <- names(candidates_svm)
comp_models <- data.frame(candidates=cnames,
                          rf=as.vector(candidates_rf),
                          svm=as.vector(candidates_svm), 
                          rbf=as.vector(candidates_rbf),
                          mlp=as.vector(prediction))

# use data frame to fill classifier lists of what they think are target proteins
p_rf  <- comp_models$candidates[comp_models$rf ==1]
p_svm <- comp_models$candidates[comp_models$svm ==1]
p_rbf <- comp_models$candidates[comp_models$rbf ==1]
p_mlp <- comp_models$candidates[comp_models$mlp ==1]


# Use Venn diagram to highlight commonly identifed protein targets between four classifiers
library(VennDiagram)
plot.new()
venn.plot <- venn.diagram(list(p_rf,p_svm,p_rbf,p_mlp), 
                          NULL, 
                          fill=c("red", "blue","green","pink"), 
                          alpha=c(0.5,0.5,0.5,0.5), 
                          cex = 2, 
                          cat.fontface=2, 
                          margins =c(10,10),
                          cat.cex=2,
                          #main = "Venn Diagram showing shared side effects for donepezil,galantamine,rivastigmine",
                          category.names=c("RandomForest", "SVM","RBF","MLP"))
grid.draw(venn.plot)








