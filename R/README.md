Because of GitHub file size limitations, there are three RDATA files containing the data structures required for the work described.

`1. classifiers1stMarch2018.RData` This contains the SVM, RBF, and MLP models.
`2. RandomForest1stMarch2018.RData` This contains the RandomForest model - taking up about 55MB in the environment.
`3. NCA-February28th2018.RData` This contain everything else.


`neuralnet_main.R` is the starting point, it calls in the other R files, loading in functions and data files.

`neuralnet_proteins_functions.R` contains the workhorse functions and loads in R libraries.

`neuralnet_complex_latent.R` produces plots that appear in the paper.

