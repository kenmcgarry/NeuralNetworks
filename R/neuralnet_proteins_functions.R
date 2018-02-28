# neuralnet_proteins_functions.R
# Predict protein type based on gene ontology mappings

library(ROCR)
library(plotROC)
library(kernlab)
library(randomForest)
library("e1071")
library(caret)
library(sand)
library(igraph)
library(GO.db)
library(plyr) 
library(dplyr)
library(reshape2)
library(tidyr)
library(boot)
library(ggplot2)
library(neuralnet)
library(nnet)
library(RSNNS)

# Calculate some statistics about the disease gene network
# returns a list: net and nodes
get_gstatistics <- function(gt) {
  net <- data.frame( 
    modu=igraph::modularity(gt, membership(cluster_walktrap(gt))),
    avepath=igraph::average.path.length(gt),
    nedges=igraph::ecount(gt),
    nverts=igraph::vcount(gt),
    transit=igraph::transitivity(gt),
    diam=igraph::diameter(gt,weights=NA),
    connect=igraph::is.connected(gt))
  
  nodes <- data.frame(   
    closeness=igraph::estimate_closeness(gt,mode="all",cutoff=3),
    degree=(igraph::degree(gt)),
    betweenness=igraph::estimate_betweenness(gt,directed=FALSE,cutoff=3),
    hubness=igraph::hub_score(gt)$vector,
    central=vector(mode="integer", length=net$nverts),
    comm=vector(mode="integer", length=net$nverts))
  
  tmp <- igraph::cluster_walktrap(gt)
  nodes$comm <- as.vector(membership(tmp))
  alpha <- igraph::alpha_centrality(ppi_net,alpha=0.1)  
  nodes$central <- as.vector(alpha)
  
  cat("\nOverall network statistics:")
  cat("\n   Modularity ",net$modu)
  cat("\n   Average path ",net$avepath)
  cat("\n   N edges ",net$nedges)
  cat("\n   N vertices ",net$nverts)
  cat("\n   Transitivity ",net$transit)
  cat("\n   Diameter ",net$diam)
  cat("\n   Is connected? ",net$connect)
  gstats = list(net=net, nodes=nodes)
  return(gstats)
}


# ggplot2 version of barplot - prettier?
bar_plot_gg2 <- function(dt,br,mycolor){
  tempnames <- sort(table(dt$TargetClass),decreasing = TRUE)
  Frequency <- as.vector(tempnames)  # counts for each protein type
  pnames <- names(tempnames)
  df <- data.frame(Frequency,pnames)
  df$Frequency <- as.numeric(df$Frequency)
  
  if(br==1){mybreaks = c(0,500,1000,1500,2000,2500,3000,3500,4000);limits<-max(mybreaks)}
  if(br==2){mybreaks = c(0,200,400,600,800,1000);limits<-max(mybreaks)}
  
  ggplot(df, aes(x = pnames,y=(Frequency))) + scale_x_discrete(limits = pnames) + 
    geom_bar(stat="identity",fill=mycolor)+
    theme(axis.text.x=element_text(face="bold",angle=40,hjust=1,size=12)) +
    theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
    ylab("Frequency count of proteins") + 
    xlab("")+
    scale_y_continuous(expand = c(0,0),breaks = mybreaks,limits = c(0,limits)) +
    theme(axis.title.y = element_text(color="black", size=14, face="bold"))
}

# factor2int will convert as suggested - the ROC::prediction function requires 
# continuous values
factor2int <- function(targettype){
  targettype <- as.numeric(targettype)
  #for (i in 1:length(targettype)){
  #  if(targettype[i] == 1) {  
  #    targettype[i] <- 0}}
  #for (i in 1:length(targettype)){
  #  if(targettype[i] == 2) {  
  #    targettype[i] <- 1}}
  
  for (i in 1:length(targettype)){
    if(targettype[i] == 1) {  
      targettype[i] <- runif(1, 7.0, 9.9)}}
  for (i in 1:length(targettype)){
    if(targettype[i] == 0) {  
      targettype[i] <- runif(1,0.01,0.02)}}
  return(targettype)
}

messyfactor2int <- function(targettype){
  targettype <- as.numeric(targettype)
  for (i in 1:length(targettype)){
    if(targettype[i] == 1) {  
      targettype[i] <- 0}}
  for (i in 1:length(targettype)){
    if(targettype[i] == 2) {  
      targettype[i] <- 1}}
  
  for (i in 1:length(targettype)){
    if(targettype[i] == 1) {  
      targettype[i] <- runif(1, 7.0, 9.9)}}
  for (i in 1:length(targettype)){
    if(targettype[i] == 0) {  
      targettype[i] <- runif(1,0.01,0.02)}}
  return(targettype)
}

# plot several ROC curves on one plot
roc_plot <- function(...){
  args = list(...)
  df <- data.frame(x=0,y=0,classifier="DECTREE",stringsAsFactors = FALSE)
  
  for (i in 1:length(args)){
    x <- attributes(args[[i]])$x.values
    y <- attributes(args[[i]])$y.values
    #if(attr(args[[i]], "roc_name") == "RandomForest"){
    #  y <- as.numeric((unlist(y[i])))
    #  x <- as.numeric((unlist(x[i])))}
    classifier <- rep(attributes(args[[i]])$roc_name,length(x))
    df_temp <- data.frame(x,y,classifier,stringsAsFactors = FALSE)
    names(df_temp) <- names(df) 
    df <- rbind(df, df_temp) 
    }
  
  df <- df[-1,]    # 1st entry is rubbish, so remove it
  df <- na.omit(df)
  ggplot(data=df, aes(x=x, y=y, group=classifier,colour=classifier)) + geom_line(size=1.5) +
    labs(x="False Positive Rate",y="True Positive Rate") +
    labs(color="") +
    theme(legend.position = "bottom", legend.direction = "horizontal")
 }   


# plot several PR curves on one plot
pr_plot <- function(...){
  args = list(...)
  df <- data.frame(x=0,y=0,classifier="DECTREE",stringsAsFactors = FALSE)
  
  for (i in 1:length(args)){
    x <- attributes(args[[i]])$x.values
    y <- attributes(args[[i]])$y.values
    #y <- as.numeric(as.character(unlist(y[[i]])))
    #x <- as.numeric(as.character(unlist(x[[i]])))
    classifier <- rep(attributes(args[[i]])$pr_name,length(x))
    df_temp <- data.frame(x,y,classifier,stringsAsFactors = FALSE)
    names(df_temp) <- names(df) 
    df <- rbind(df, df_temp) 
  }
  
  df <- df[-1,]    # 1st entry is rubbish, so remove it
  df <- na.omit(df)
  ggplot(data=df, aes(x=x, y=y, group=classifier,colour=classifier)) + geom_line(size=1.5) +
    labs(x="Recall",y="Precision") +
    labs(color="") +
    theme(legend.position = "bottom", legend.direction = "horizontal")
}




