# neuralnet_proteins.R
# Predict protein type based on gene ontology mappings
# To be submitted to Advances in data Analysis and Classification
# http://www.springer.com/statistics/journal/11634
# started : 19/1/2018
# completed:

memory.limit(1510241024*1024) # allocate RAM memory (15 GBs)
setwd("C:/common_laptop/R-files/NeuralNet")  # now point to where the new code lives
load("hybrid12thApril2018.RData")
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
bdata <- rbind(positives,negatives)


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
caret::confusionMatrix(data=pred_test,reference=ytest, mode = "prec_recall", positive="1")
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
caret::confusionMatrix(data=pred_rbf,reference=ytest, mode = "prec_recall", positive="1")
auc.perf = ROCR::performance(pred, measure = "auc")
auc.perf@y.values

# Train Random Forest on data 
rf_model <-randomForest(as.factor(ytrain) ~.,data=xtrain[,1:149],proximity=TRUE,
                        keep.forest=TRUE,ntree=500,mtry=12)
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
caret::confusionMatrix(data=pred_rf,reference=ytest, mode = "prec_recall", positive="1")

# examine most important variable by order using Gini Index
impvar <- data.frame(importance(rf_model))
impvar <- cbind(rownames(impvar),impvar)  # we need variable names
arrange(impvar, desc(MeanDecreaseGini))

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

temppredict <- as.factor(ones_twos(idx))

prediction <- c('target', 'nontarget')[idx]
acc <- table(prediction, ytest)
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\nMLP accuracy calculated by (TP+TN)/(TP+TN+FP+FN)= ",(TP + TN)/(TP + TN + FP + FN))

caret::confusionMatrix(data=temppredict,reference=ytest,positive="1")
caret::confusionMatrix(data=temppredict,reference=ytest, mode = "prec_recall", positive="1")

prediction[prediction =="target"] <- 1
prediction[prediction =="nontarget"] <- 0
prediction <- factor2int(as.integer(prediction))
pred_mlp <- ROCR::prediction(prediction,ytest)
auc.perf <- ROCR::performance(pred_mlp, measure = "auc")
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
caret::confusionMatrix(data=as.factor(nnet_pred),reference=ytest,positive="1")
caret::confusionMatrix(data=as.factor(nnet_pred),reference=ytest, mode = "everything", positive="1")

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
candidates <- na.omit(candidates)

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

# make a nice well structured data frame, for some reason svm and rbf only produce 1264 candidates from 1280 test data  
cnames <- rownames(candidates)
comp_model <- data.frame(candidates=cnames,
                          rf =as.vector(candidates_rf),
                          svm=as.vector(candidates_svm), 
                          rbf=as.vector(candidates_rbf),
                          mlp=as.vector(prediction))

# use data frame to fill classifier lists of what they think are target proteins
p_rf  <- na.omit(comp_model$candidates[comp_model$rf  ==1])
p_svm <- na.omit(comp_model$candidates[comp_model$svm ==1])
p_rbf <- na.omit(comp_model$candidates[comp_model$rbf ==1])
p_mlp <- na.omit(comp_model$candidates[comp_model$mlp ==1])

#commonproteins <- Reduce(intersect, list(p_rf,p_svm,p_rbf,p_mlp))  # the proteins identified by all four classifiers

# plot a Venn diagram to highlight commonly identifed protein targets between four classifiers
plot_venn(p_rf,p_svm,p_rbf,p_mlp)


#### LATENT NETWORK STAGE ####
# VERY MEMORY INTENSIVE AND SLOW, (16 hours for 500 proteins+crossvalidation) HENCE ONLY USE NECESSARY DATA
# START NEW R SESSION FROM THIS POINT
# http://ptrckprry.com/course/ssd/lecture/latentnet.html

library(eigenmodel)
library(latentnet)

setwd("C:/R-files/NeuralNet")  # now point to where the new code lives
source("neuralnet_functions.R")  # load in the functions required for this work. 
source("neuralnet_data.R")
load("latent.RData")
ppi_net <- build_network(ppi) # create an igraph object with attributes of hub, target, proteintype, k-core.
set.seed(42)

# create permutations of smaller (PSIZE) ppi network (otherwise gives memory allocation error : 1.8GB)
PSIZE <- 200
ppi_names <- V(ppi_net)$name
new_names <- sample(ppi_names,PSIZE)
new_ppi <- induced.subgraph(ppi_net,which(V(ppi_net)$name %in% new_names))
A <- igraph::get.adjacency(new_ppi,sparse = FALSE)

######################################
# extract a series of subgraphs
k_cores <- c("PSMA3","PPIG","KIR3DS1","SELE","FPR2","FXYD6","PTPRK","COL4A3BP","RPA2","SLC38A6",
             "SLC10A3","SLC16A14","UGDH","STEAP1","FZD4","SLC16A13","GRID1","GOT1L1","MBTPS2","EPHA10", 
             "AKAP8",  "GPR6",   "ADAMTS2", "CELA3A", "MAT1A","FXYD6","PPIG","FZD5","CRYZL1","ACP5",
             "PTPRK","SLC4A2","FPR2","CD47","BRCA2","SELE","KIR3DS1","PSMA3")
explore_subgraph <- induced.subgraph(graph=ppi_net,vids=unlist(neighborhood(graph=ppi_net,order=1,
                                                                            nodes=k_cores[34])))
length(V(explore_subgraph)) 
A <- igraph::get.adjacency(explore_subgraph,sparse = FALSE)

######################################
# fit ERGM to A
samp_fit <- ergmm(A ~ euclidean(d=2))
samp_fit <- ergmm(A ~ euclidean(d=2, G=2)) #cluster based
plot(samp_fit, pie=TRUE,labels=TRUE)
summary(samp_fit)
plot(samp_fit)
mcmc.diagnostics(samp_fit)
######################################

ppi.fit1 <- eigenmodel::eigenmodel_mcmc(A, R=2, S=11000,burn=10000)

#Creat hub covariate
hub.op <- vertex_attr(new_ppi)$hub %o%  vertex_attr(new_ppi)$hub
hub.effect <- matrix(as.numeric(hub.op %in% c(1, 4, 9)), PSIZE, PSIZE)
hub.effect <- array(hub.effect,dim=c(PSIZE, PSIZE, 1))
ppi.fit2 <- eigenmodel_mcmc(A, hub.effect, R=2,S=11000,burn=10000)

# Create type (protein type) covariate but convert from strings to factors
type.fac <- as.factor(vertex_attr(new_ppi)$type)
type.fac <- as.numeric(type.fac)
type.op <- type.fac %o% type.fac
type.effect <- matrix(as.numeric(type.op %in% c(1, 4, 9)),PSIZE,PSIZE)
type.effect <- array(type.effect,dim=c(PSIZE, PSIZE, 1))
ppi.fit3 <- eigenmodel_mcmc(A, type.effect,R=2, S=11000, burn=10000)

# Create core covariate
core.op <- vertex_attr(new_ppi)$coreness %o% vertex_attr(new_ppi)$coreness
core.effect <- matrix(as.numeric(core.op %in% c(1, 4, 9)),PSIZE,PSIZE)
core.effect <- array(core.effect,dim=c(PSIZE, PSIZE, 1))
ppi.fit4 <- eigenmodel_mcmc(A, core.effect,R=2, S=11000, burn=10000)

# Get eigenvectors
latent1 <- eigen(ppi.fit1$ULU_postmean)$vec[, 1:2]
latent2 <- eigen(ppi.fit2$ULU_postmean)$vec[, 1:2]
latent3 <- eigen(ppi.fit3$ULU_postmean)$vec[, 1:2]
latent4 <- eigen(ppi.fit4$ULU_postmean)$vec[, 1:2]

# Examine posterior means of latent variables
apply(ppi.fit1$L_postsamp, 2, mean)
apply(ppi.fit2$L_postsamp, 2, mean)
apply(ppi.fit3$L_postsamp, 2, mean)
apply(ppi.fit4$L_postsamp, 2, mean)

# Test model goodness of fit
perm.index <- sample(1:((PSIZE*(PSIZE-1))/2)) # permutations
nfolds <- 5
nmiss <- ((PSIZE*(PSIZE-1))/2)/nfolds
Avec <- A[lower.tri(A)]
Avec.pred1 <- numeric(length(Avec))


# RECREATE SAME MODELS BUT WITH LOOP TO IMPLEMENT n-fold CROSS-VALIDATION 
perf1 <- list() # create empty vector of perfs for ROC-Curves
auc1 <- list() # create empty vector of area under curves
pred1 <- list() # create empty vector of predictions for PR-Curves

perf2 <- list() # create empty vector of perfs for ROC-Curves
auc2 <- list() # create empty vector of area under curves
pred2 <- list() #

perf3 <- list() # create empty vector of perfs for ROC-Curves
auc3 <- list() # create empty vector of area under curves
pred3 <- list() #

perf4 <- list() # create empty vector of perfs for ROC-Curves
auc4 <- list() # create empty vector of area under curves
pred4 <- list() #

for (j in 1:4){
  if(j==1) modelversion <- NULL   # select only connectivity
  if(j==2) modelversion <- hub.effect  # select only hub variable
  if(j==3) modelversion <- core.effect  # select only coreness variable
  if(j==4) modelversion <- type.effect  # select only type variable
  # CHUNK 31 - cross validation.
  for(i in seq(1,nfolds)){
    # Index of missing values.
    miss.index <- seq(((i-1) * nmiss + 1),(i*nmiss), 1)
    A.miss.index <- perm.index[miss.index]
    
    # Fill a new Atemp appropriately with NA's.
    Avec.temp <- Avec
    Avec.temp[A.miss.index] <- rep("NA", length(A.miss.index))
    Avec.temp <- as.numeric(Avec.temp)
    Atemp <- matrix(0, PSIZE, PSIZE) # varying these two numbers leads to different accuracies
    Atemp[lower.tri(Atemp)] <- Avec.temp
    Atemp <- Atemp + t(Atemp)
    
    # Now fit model and predict.
    Y <- Atemp
    model1.fit <- eigenmodel_mcmc(Y, modelversion,R=2,S=11000, burn=10000) # change model here
    model1.pred <- model1.fit$Y_postmean
    model1.pred.vec <- model1.pred[lower.tri(model1.pred)]
    Avec.pred1[A.miss.index] <- model1.pred.vec[A.miss.index]
  }
  #pred1 <- floor(Avec.pred1)
  Avec[Avec==2] <- 0
  if(j==1){
    pred1 <- ROCR::prediction(Avec.pred1, Avec)
    perf1 <- ROCR::performance(pred1, "tpr", "fpr")
    auc1 <- ROCR::performance(pred1, "auc")}
  if(j==2){
    pred2 <- ROCR::prediction(Avec.pred1, Avec)
    perf2 <- ROCR::performance(pred2, "tpr", "fpr")
    auc2 <- ROCR::performance(pred2, "auc")}
  if(j==3){
    pred3 <- ROCR::prediction(Avec.pred1, Avec)
    perf3 <- ROCR::performance(pred3, "tpr", "fpr")
    auc3 <-  ROCR::performance(pred3, "auc")}
  if(j==4){
    pred4 <- ROCR::prediction(Avec.pred1, Avec)
    perf4 <- ROCR::performance(pred4, "tpr", "fpr")
    auc4 <-  ROCR::performance(pred4, "auc")}
}



# Plot ROC curves
# ROC PLOT COMPARING THREE MODELS
plot(perf1, col="red", lwd=2)
plot(perf2, add = TRUE, col="blue",lwd=2)
plot(perf3, add = TRUE, col="green",lwd=2)
plot(perf4, add = TRUE, col="orange",lwd=2) 

library(ROCR)

# CHUNK 33
slot(auc1, "y.values")
slot(auc2, "y.values")
slot(auc3, "y.values")
slot(auc4, "y.values")

mpred2 <- as.numeric(unlist(slot(pred2,"predictions")))
caret::confusionMatrix(data=round(mpred2),reference=Avec, mode = "everything", positive="1")

## ERGM CODE
#ppi_s <- network::as.network(as.matrix(A),directed=FALSE)
#network::set.vertex.attribute(ppi_s,"hub",v.attrs$hub)
#network::set.vertex.attribute(ppi_s,"type",v.attrs$type)
#network::set.vertex.attribute(ppi_s,"hub",v.attrs$)
#network::set.vertex.attribute(ppi_s,"hub",v.attrs$hub)
#ppi_ergm <- build_ergm(ppi_net)  # make the ergm model.


