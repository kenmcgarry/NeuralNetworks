# neuralnet_complex_latent.R

# build complex network subgraph from main ppi network containing the names of nontargets and targets 
# proteins from data used to build classifiers.

proteins <- rownames(balanced_dat)  # get the names of the proteins in our data
ids <- match(proteins, V(ppi_net)$name)  # convert these to numeric index
ids <- ids[!is.na(ids)]   # unfortunatley some do not exist in main component, so remove NA's
  
ppi_subgraph <- induced.subgraph(graph=ppi_net,vids=ids)  # create the subgraph
length(V(ppi_subgraph)) 

gs <- get_gstatistics(ppi_subgraph)
head(gs)

# how to access attributes
vertex_attr(ppi_subgraph, "target")
vertex_attr(ppi_subgraph, "hub")
vertex_attr(ppi_subgraph, "type")
vertex_attr(ppi_subgraph, "comp")

# Add "hub" and "type" to attributes of ppi_subgraph, we already have "target" attribute
# key structures : hublist, drug_targets, protein_class, 






