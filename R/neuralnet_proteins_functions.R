# neuralnet_proteins_functions.R
# Predict protein type based on gene ontology mappings

library(ROCR)
library(kernlab)
library(randomForest)
library("e1071")
library(caret)
library(sand)
library(igraph)
library(GO.db)
library(plyr) 
library(dplyr)
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







