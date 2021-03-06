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
library(ontologySimilarity)
library(ontologyIndex)
library(GSEABase)

data(go)
data(gene_GO_terms)
data(GO_IC)

# remove isolated proteins if they are connected to the giant component
delete_isolates <- function(gt) {
  isol <- V(gt)[igraph::degree(gt)==0]
  gt <- igraph::delete.vertices(gt, isol)
  return(gt)
  
}

build_network <- function(ppt){
  # unsure if protein targets are part of the giant connected ppi network
  # assign 1=target; 0=non-target to each protein
  
  # remove multiple genes
  ppt$Gene_A <-gsub("\\|.*","",ppt$Gene_A)
  ppt$Gene_B <-gsub("\\|.*","",ppt$Gene_B)  # 
  
  un_targets <- (unique(drug_targets$Gene))        # 1,860 unique protein targets
  length(un_targets)
  un_ppi <- (unique(c(ppt$Gene_A,ppt$Gene_B)))      # 15,792 unique general proteins in ppi
  length(un_ppi)
  
  joint_ppi <- un_targets[un_targets %in% un_ppi]  # 1,293 targets are part of giant connected network (we lose 567 targets!)
  not_ppi <- un_targets[!un_targets %in% un_ppi]  # here are the 567 targets)
  
  # dataframe containing targets and non-target proteins. Annotate with:
  # 1. target; 2. hub; 
  # create ppi network (igraph object) and annotate with target or not target
  ppi_temp <- igraph::graph.data.frame(ppt)
  ppi_temp <- igraph::as.undirected(ppi_temp)
  ppi_temp <- igraph::simplify(ppi_temp)  # remove duplicates and self-loops
  ppi_temp <- delete_isolates(ppi_temp)
  igraph::delete.vertices(igraph::simplify(ppi_temp), igraph::degree(ppi_temp)==0)
  
  V(ppi_temp)[1:vcount(ppi_temp)]$target <- 0   # Intialise all to zeros
  V(ppi_temp)[1:vcount(ppi_temp)]$hub <- 0   # Intialise all to zeros
  V(ppi_temp)[1:vcount(ppi_temp)]$type <- "Unknown"   # Intialise protein "type" to unknown
  V(ppi_temp)[1:vcount(ppi_temp)]$coreness <- 0   # Intialise all to zeros
  
  # get main component only - ignore lessor weakly connected groups
  V(ppi_temp)$comp <- igraph::components(ppi_temp)$membership
  ppi_temp <- igraph::induced_subgraph(ppi_temp,V(ppi_temp)$comp==1)
  
  # remove from joint_ppi the lost nodes 
  survivors <- V(ppi_temp)$name
  joint_ppi <- un_targets[un_targets %in% survivors] 
  
  # determine if protein is a hub
  netstats <- get_gstatistics(ppi_temp)
  hubs <- find_hubs(netstats[[2]])
  
  # get the coreness for each protein # set_vertex_attr("label", value = LETTERS[1:10])
  coreness <- graph.coreness(as.undirected(ppi_temp))
  
  # Remove small quantity proteins; Adhesion; Nuclear Other; Antibody; CD Molecules; Ribosomal; Cytokine; Surface Antigen; Membrane other
  drug_targets <-  # Only keep protein target types with at least 50 occurences
    drug_targets %>%
    add_count(TargetClass,sort=TRUE) %>%
    filter(n > 50)
  
  # change names for drug_target and protein_class types
  drug_targets$TargetClass[drug_targets$TargetClass == "Ion channel"] <- "IC"
  drug_targets$TargetClass[drug_targets$TargetClass == "Nuc receptor"] <- "NR"
  drug_targets$TargetClass[drug_targets$TargetClass == "Transcription"] <- "TF"
  drug_targets$TargetClass[drug_targets$TargetClass == "Cytosolic other"] <- "Transporter"
  drug_targets$TargetClass[drug_targets$TargetClass == "Unclassified"] <- "Unknown"
  protein_class$TargetClass[protein_class$TargetClass == "oGPCR"] <- "GPCR"
  protein_class$TargetClass[protein_class$TargetClass == "TF; Epigenetic"] <- "TF"
  protein_class$TargetClass[protein_class$TargetClass == "Epigenetic"] <- "TF"
  
  # assumes "more_proteins" has beenloaded by neuralnet_data.R
  drug_targets2 <- drug_targets[,2:3]
  all_proteins <- rbind(more_proteins,protein_class,drug_targets2)
  all_proteins <- all_proteins %>% distinct(Gene,.keep_all = TRUE) # remove duplicates
    
  ppi_names <- V(ppi_temp)$name
  ppi_present <- all_proteins[all_proteins$Gene %in% ppi_names,]
    
  ppi_temp <- igraph::set_vertex_attr(ppi_temp,"type",ppi_present$Gene,ppi_present$TargetClass) # Now assign target types
  
  ppi_temp <- igraph::set_vertex_attr(ppi_temp,"coreness",names(coreness),coreness) # Now assign coreness
  ppi_temp <- igraph::set_vertex_attr(ppi_temp,"target",joint_ppi,1) # Now assign "1" if protein is a target (very neat coding!)
  ppi_temp <- igraph::set_vertex_attr(ppi_temp,"hub",hubs$genenames,1) # Now assign "1" if protein is a hub (very neat coding!)
  # vertex_attr(g)
  return(ppi_temp)
}


# annotate_go() will receive a list of proteins and annotate with GO terms it will return
# a matrix of terms and proteins. Uses the GO data by Daniel Greene.
# This fucntion is modified from complexnetworks.
annotate_with_go <- function(plist){
  category <- c("MF","BP","CC")
  nproteins <- length(plist)
 
  tempgo <- gene_GO_terms[plist]  # Annotate!!
  tempgo <- tempgo[!sapply(tempgo, is.null)]  # Not all proteins have GO annotations so remove them.
  go_proteins <- names(tempgo) # Unfortunately, we are left with only 13,417 proteins.
  # sort GO terms by the three categories, breakdown maybe useful at later date for summary statistics
  cc <- go$id[go$name == "cellular_component"]
  bp <- go$id[go$name == "biological_process"]
  mf <- go$id[go$name == "molecular_function"] 
  temp_cc <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=cc, x))
  temp_bp <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=bp, x))
  temp_mf <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=mf, x))

  return(tempgo)
}

# go_slim_annotation() reduces the complexities of numerous GO annoatations into a few key terms.
# http://www.geneontology.org/page/go-slim-and-subset-guide#On_the_web. This uses the GSEABase package.
# WARNING: this takes a looong time to compute..... approx 3 hours on my laptop
go_slim_annotation <- function(mylist){
  gostuff <- annotate_with_go(mylist)
  gostuff <-gostuff[names(which(lapply(gostuff, length) >1))]  # keep only genes with 1 or more GO annotations
  assign("go_mf", TRUE, env=globalenv())  # global variables to enable recovery from errors generated by goSlim()
  assign("go_cc", TRUE, env=globalenv())
  assign("go_bp", TRUE, env=globalenv())
  
  # create the matrix for classification algorithms, if a GO term is present mark it with by "1" in that column
  mm <- matrix(0, 149, length(names(gostuff)))  # mm=Number of GO terms in GoSlim (CC+BP+MF) x Number of genes
  colnames(mm) <- names(gostuff)
  rownames(mm) <- give_rownames_mm()
  
  for (i in 1:length(names(gostuff))){
    myCollection <- GOCollection(gostuff[[i]])
    genename <- names(gostuff[i])
    obo <- system.file("extdata","goslim_generic.obo", package="GSEABase") # generic terms by GO consortium
    #obo <- system.file("extdata","goslim_chembl.obo", package="GSEABase") # Chembl Drug Target developed by Mutowo and Lomax
    #obo <- system.file("extdata","goslim_pir.obo", package="GSEABase") #Protein Info Resource by Darren Natale
    slim <- getOBOCollection(obo)
    go_mf <- tryCatch(goSlim(myCollection, slim, "MF"),error=function(e) {go_mf <- error_go_mf()})
    go_cc <- tryCatch(goSlim(myCollection, slim, "CC"),error=function(e) {go_cc <- error_go_cc()})
    go_bp <- tryCatch(goSlim(myCollection, slim, "BP"),error=function(e) {go_bp <- error_go_bp()})
    
    # a lot of hard coded magic numbers, the goslim has 149 terms.
    if(length(go_mf) ==1) {mm[1:43,i]  <- as.vector(matrix(0,ncol=43))} else{  # fill with zeros if no annotations found
      go_mf[go_mf$Count != 0,]$Count <- 1; # convert non-zero numbers into 1's
      mm[1:43,i] <- go_mf$Count}            # found MF annotations, assign to matrix
    if(length(go_cc)==1) {mm[44:78,i] <- as.vector(matrix(0,ncol=35))} else{
      go_cc[go_cc$Count != 0,]$Count <- 1
      mm[44:78,i] <- go_cc$Count}
    if(length(go_bp) == 1) {mm[79:149,i]<-  as.vector(matrix(0,ncol=71))} else{
      go_bp[go_bp$Count != 0,]$Count <- 1
      mm[79:149,i] <- go_bp$Count}
  }
  return(mm)  # return matrix  
}

# Add the "drug target" class to the "mm" matrix, this will make it ready for classification algorithms.
# but first remove nodes that do have GO terms.
give_classlabels_mm <- function(mm){
  
  survivors <- colnames(mm)
  mmt <- matrix(0, 1, ncol(mm)) 
  rownames(mmt) <- "targets"
  # delete nodes without GO terms from ppi_net 
  #biglist <- V(ppi_net)$name
  #lostnodes <- survivors[!biglist %in% survivors]
  #lostnodes <- lostnodes[is.na(lostnodes)] <- 0
  #ppi_net <-igraph::delete_vertices(ppi_net,lostnodes[1:2020])
  
  targets <- drug_targets$Gene  # get the original drug target data
  
  for(i in 1:length(survivors)){
    temp <- grep(survivors[i],targets)
    if(length(temp) != 0){mmt[i] <- 1}
  }
  
  mmt <- rbind(mm,mmt)
  return(mmt)
}

give_rownames_mm <- function(){
  gonames <- c("GO:0000988:MF:transcription factor activity, prot",
               "GO:0003674:MF:molecular_function",
               "GO:0003677:MF:DNA binding",
               "GO:0003700:MF:transcription factor activity, sequ",
               "GO:0003723:MF:RNA binding",
               "GO:0003729:MF:mRNA binding",
               "GO:0003735:MF:structural constituent of ribosome",
               "GO:0003924:MF:GTPase activity",
               "GO:0004386:MF:helicase activity",
               "GO:0004518:MF:nuclease activity",
               "GO:0004871:MF:signal transducer activity",
               "GO:0005198:MF:structural molecule activity",
               "GO:0008092:MF:cytoskeletal protein binding",
               "GO:0008134:MF:transcription factor binding",
               "GO:0008135:MF:translation factor activity, RNA",
               "GO:0008168:MF:methyltransferase activity",
               "GO:0008233:MF:peptidase activity",
               "GO:0008289:MF:lipid binding",
               "GO:0008565:MF:protein transporter activity",
               "GO:0016301:MF:kinase activity",
               "GO:0016491:MF:oxidoreductase activity",
               "GO:0016746:MF:transferase activity, transferring",
               "GO:0016757:MF:transferase activity, transferring",
               "GO:0016765:MF:transferase activity, transferring",
               "GO:0016779:MF:nucleotidyltransferase activity",
               "GO:0016791:MF:phosphatase activity",
               "GO:0016798:MF:hydrolase activity, acting on glyco",
               "GO:0016810:MF:hydrolase activity, acting on carbo",
               "GO:0016829:MF:lyase activity",
               "GO:0016853:MF:isomerase activity",
               "GO:0016874:MF:ligase activity",
               "GO:0016887:MF:ATPase activity",
               "GO:0019843:MF:rRNA binding",
               "GO:0019899:MF:enzyme binding",
               "GO:0022857:MF:transmembrane transporter activity",
               "GO:0030234:MF:enzyme regulator activity",
               "GO:0030533:MF:triplet codon-amino acid adaptor",
               "GO:0030555:MF:RNA modification guide activity",
               "GO:0030674:MF:protein binding, bridging",
               "GO:0032182:MF:ubiquitin-like protein binding",
               "GO:0042393:MF:histone binding",
               "GO:0043167:MF:ion binding",
               "GO:0051082:MF:unfolded protein binding",
               "GO:0000228:CC:nuclear chromosome",
               "GO:0000229:CC:cytoplasmic chromosome",
               "GO:0005575:CC:cellular_component",
               "GO:0005576:CC:extracellular region",
               "GO:0005578:CC:proteinaceous extracellular matrix",
               "GO:0005615:CC:extracellular space",
               "GO:0005618:CC:cell wall",
               "GO:0005622:CC:intracellular",
               "GO:0005623:CC:cell",
               "GO:0005634:CC:nucleus",
               "GO:0005635:CC:nuclear envelope",
               "GO:0005654:CC:nucleoplasm",
               "GO:0005694:CC:chromosome",
               "GO:0005730:CC:nucleolus",
               "GO:0005737:CC:cytoplasm",
               "GO:0005739:CC:mitochondrion",
               "GO:0005764:CC:lysosome",
               "GO:0005768:CC:endosome",
               "GO:0005773:CC:vacuole",
               "GO:0005777:CC:peroxisome",
               "GO:0005783:CC:endoplasmic reticulum",
               "GO:0005794:CC:Golgi apparatus2",
               "GO:0005811:CC:lipid particle",
               "GO:0005815:CC:microtubule organizing center",
               "GO:0005829:CC:cytosol",
               "GO:0005840:CC:ribosome",
               "GO:0005856:CC:cytoskeleton",
               "GO:0005886:CC:plasma membrane",
               "GO:0005929:CC:cilium",
               "GO:0009536:CC:plastid",
               "GO:0009579:CC:thylakoid",
               "GO:0030312:CC:external encapsulating structure",
               "GO:0031410:CC:cytoplasmic vesicle",
               "GO:0043226:CC:organelle",
               "GO:0043234:CC:protein complex",
               "GO:0000003:BP:reproduction",
               "GO:0000278:BP:mitotic cell cycle",
               "GO:0000902:BP:cell morphogenesis",
               "GO:0002376:BP:immune system process",
               "GO:0003013:BP:circulatory system process",
               "GO:0005975:BP:carbohydrate metabolic process",
               "GO:0006091:BP:generation of precursor metabolites",
               "GO:0006259:BP:DNA metabolic process",
               "GO:0006397:BP:mRNA processing",
               "GO:0006399:BP:tRNA metabolic process",
               "GO:0006412:BP:translation",
               "GO:0006457:BP:protein folding",
               "GO:0006461:BP:protein complex assembly",
               "GO:0006464:BP:cellular protein modification process",
               "GO:0006520:BP:cellular amino acid metabolic process",
               "GO:0006605:BP:protein targeting",
               "GO:0006629:BP:lipid metabolic process",
               "GO:0006790:BP:sulfur compound metabolic process",
               "GO:0006810:BP:transport",
               "GO:0006913:BP:nucleocytoplasmic transport",
               "GO:0006914:BP:autophagy",
               "GO:0006950:BP:response to stress",
               "GO:0007005:BP:mitochondrion organization",
               "GO:0007009:BP:plasma membrane organization",
               "GO:0007010:BP:cytoskeleton organization",
               "GO:0007034:BP:vacuolar transport",
               "GO:0007049:BP:cell cycle",
               "GO:0007059:BP:chromosome segregation",
               "GO:0007155:BP:cell adhesion",
               "GO:0007165:BP:signal transduction",
               "GO:0007267:BP:cell-cell signaling",
               "GO:0007568:BP:aging",
               "GO:0008150:BP:biological_process",
               "GO:0008219:BP:cell death",
               "GO:0008283:BP:cell proliferation",
               "GO:0009056:BP:catabolic process",
               "GO:0009058:BP:biosynthetic process",
               "GO:0009790:BP:embryo development",
               "GO:0015979:BP:photosynthesis",
               "GO:0016192:BP:vesicle-mediated transport",
               "GO:0019748:BP:secondary metabolic process",
               "GO:0021700:BP:developmental maturation",
               "GO:0022607:BP:cellular component assembly",
               "GO:0022618:BP:ribonucleoprotein complex assembly",
               "GO:0030154:BP:cell differentiation",
               "GO:0030198:BP:extracellular matrix organization",
               "GO:0030705:BP:cytoskeleton-dependent intracellular",
               "GO:0032196:BP:transposition",
               "GO:0034330:BP:cell junction organization",
               "GO:0034641:BP:cellular nitrogen compound metabolism",
               "GO:0034655:BP:nucleobase-containing compound catalysis",
               "GO:0040007:BP:growth",
               "GO:0040011:BP:locomotion",
               "GO:0042254:BP:ribosome biogenesis",
               "GO:0042592:BP:homeostatic process",
               "GO:0043473:BP:pigmentation",
               "GO:0044281:BP:small molecule metabolic process",
               "GO:0044403:BP:symbiosis, encompassing mutualism",
               "GO:0048646:BP:anatomical structure formation involvement",
               "GO:0048856:BP:anatomical structure development",
               "GO:0048870:BP:cell motility",
               "GO:0050877:BP:neurological system process",
               "GO:0051186:BP:cofactor metabolic process",
               "GO:0051276:BP:chromosome organization",
               "GO:0051301:BP:cell division",
               "GO:0051604:BP:protein maturation",
               "GO:0055085:BP:transmembrane transport",
               "GO:0061024:BP:membrane organization",
               "GO:0065003:BP:macromolecular complex assembly",
               "GO:0071554:BP:cell wall organization or biogenesis",
               "GO:0071941:BP:nitrogen cycle metabolic process")
  
} 

# The following functions are to avoid the errors kicked out by goSlim (and crashing R) 
# when it cant find any terms
error_go_mf <- function(){
  return("MF: Error detected!")
  #go_default <- "error"
}
error_go_cc <- function(){
  return("CC: Error detected!")
  #go_default <- "error"
}
error_go_bp <- function(){
  return("BP: Error detected!")
  #go_default <- "error"
}

# Drug targets of the various drugs - usefully contains protein type (e.g. GPCR) as well.
# from http://drugcentral.org/download
# DrugCentral is a comprehensive drug information resource for FDA drugs and drugs approved outside USA. The 
# resources can be searched using: drug, target, disease, pharmacologic action, terms. 
load_drugtargets <- function(){
  drug_targets <- read.csv(file="C://common_laptop//R-files//disease//drug.target.interaction.tsv", header=TRUE, sep="\t",stringsAsFactors = FALSE)
  names(drug_targets)[names(drug_targets)=="DRUG_NAME"] <- "DrugName"
  names(drug_targets)[names(drug_targets)=="TARGET_CLASS"] <- "TargetClass"
  names(drug_targets)[names(drug_targets)=="GENE"] <- "Gene"
  drug_targets <- drug_targets[,c(1,4,6)]  # Only need three variables
  drug_targets$DrugName <- firstup(drug_targets$DrugName)   # convert first letter to uppercase to match existing data
  drug_targets <- na.omit(drug_targets)# remove NA's
  
  # now unlist special entries, I edited the original file and replaced "|" with "/"
  drug_targets<-
    drug_targets %>% 
    mutate(Gene=strsplit(as.character(Gene), "/")) %>%   # symbols=Gene
    unnest(Gene)
  drug_targets$Gene <- toupper(drug_targets$Gene)  # all to uppercase
  
  # shorten some names, for ease printing etc
  drug_targets$TargetClass <- gsub('Nuclear hormone receptor', 'Nuc receptor', drug_targets$TargetClass)
  drug_targets$TargetClass <- gsub('Transcription factor', 'Transcription', drug_targets$TargetClass)
  drug_targets$TargetClass <- gsub('Membrane receptor', 'Membrane', drug_targets$TargetClass)
  
  return(drug_targets)
}

# Makes first letter of string uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}


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

ones_twos <- function(targettype){
  targettype <- as.numeric(targettype)
  for (i in 1:length(targettype)){
    if(targettype[i] == 1) {  
      targettype[i] <- 0}}
  for (i in 1:length(targettype)){
    if(targettype[i] == 2) {  
      targettype[i] <- 1}}
  
  return(targettype)
}

# plot several ROC curves on one plot
roc_plot <- function(...){
  args = list(...)
  df <- data.frame(x=0,y=0,classifier="DECTREE",stringsAsFactors = FALSE)
  
  for (i in 1:length(args)){
    x <- attributes(args[[i]])$x.values
    y <- attributes(args[[i]])$y.values
    classifier <- rep(attributes(args[[i]])$roc_name,length(x))
    df_temp <- data.frame(x,y,classifier,stringsAsFactors = FALSE)
    names(df_temp) <- names(df) 
    df <- rbind(df, df_temp) 
    }
  
  df <- df[-1,]    # 1st entry is rubbish, so remove it
  df <- na.omit(df)
  ggplot(data=df, aes(x=x, y=y, group=classifier,colour=classifier)) + geom_line(size=1.0) +
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
    classifier <- rep(attributes(args[[i]])$pr_name,length(x))
    df_temp <- data.frame(x,y,classifier,stringsAsFactors = FALSE)
    names(df_temp) <- names(df) 
    df <- rbind(df, df_temp) 
  }
  
  df <- df[-1,]    # 1st entry is rubbish, so remove it
  df <- na.omit(df)
  ggplot(data=df, aes(x=x, y=y, group=classifier,colour=classifier)) + geom_line(size=1.0) +
    labs(x="Recall",y="Precision") +
    labs(color="") +
    theme(legend.position = "bottom", legend.direction = "horizontal")
}


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
  alpha <- igraph::alpha_centrality(gt,alpha=0.1)  
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

plot_venn <- function(p_rf,p_svm,p_rbf,p_mlp){
  # Use Venn diagram to highlight commonly identifed protein targets between four classifiers
  library(VennDiagram)
  plot.new()
  venn.plot <- venn.diagram(list(p_rf,p_svm,p_rbf,p_mlp), 
                            NULL, 
                            fill=c("red", "blue","green","yellow"), 
                            alpha=c(0.5,0.5,0.5,0.5), 
                            cex = 2, 
                            cat.fontface=2, 
                            margins =c(10,10),
                            cat.cex=2,
                            #main = "Venn Diagram showing shared side effects for donepezil,galantamine,rivastigmine",
                            category.names=c("RandomForest", "SVM","RBF","MLP"))
  grid.draw(venn.plot)  
  
  
}

# find_hubs() when presented with graph stats object will search for hubs and return list
# it will also add genenames to hublist. Need to create igraph object first then run get
# gs_statistics() to get the required data on "degree" for each protein.
find_hubs <- function(gstats){
  genenames <- as.character(rownames(gstats))
  hublist <- cbind(gstats,genenames)
  cutoff <- quantile(gstats$degree, probs = c(0.70, 0.75, 0.8, 0.85, 0.9, 0.99), na.rm = TRUE) 
  hublist <- filter(hublist,degree > cutoff[1])
  hublist <- data.frame(lapply(hublist, as.character), stringsAsFactors=FALSE)
  
  return(hublist)
}

# is_hub_target() receives a list of hubs, the drug_target data and the PPI network to see if
# these hub proteins are also targets.
is_hub_target <- function(hlist,dt,ppi){
  hub_targ_list <- dt[1,] # instantiate before use
  gnames <- hlist$genenames
  totalgenes <- c(ppi[,1],ppi[,2])
  cat("\nWe have ",length(unique(totalgenes))," unique genes in PPI network")
  cat("\nWe have ",length(unique(gnames))," unique HUB genes in PPI network")
  for (i in 1:length(gnames)){
    gene <- gnames[i] # get hub genes individually and see if they in lists of targets
    glist <- filter(dt, Gene == gene)  # This bit is OK
    if(nrow(glist) > 0){
      hub_targ_list <- rbind(hub_targ_list,glist) }
  }
  # Line below removes duplicates that appear in two variables
  hub_targ_list <-hub_targ_list[!(duplicated(hub_targ_list[c("DrugName","Gene")]) | duplicated(hub_targ_list[c("DrugName","Gene")], fromLast = TRUE)), ]
  hub_targ_list <- hub_targ_list[-1,]    # 1st entry is rubbish, so remove it
  cat("\nWe have ",length(unique(hub_targ_list$Gene))," unique genes that are hubs AND targets")
  cat("\nWe have ",length(unique(dt$Gene)) - length(unique(hub_targ_list$Gene)),   " unique genes that are targets but NOT hubs")
  cat("\nWe have ",length(unique(dt$Gene)),   " unique genes that are targets in total")
  
  return(hub_targ_list)
}





