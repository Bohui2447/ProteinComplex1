# ProteinComplex1
This work was conducted under R version 3.5.1  
First download the data from "Link" and save them in the current working directory  

# Install functions used in this project  
devtools::install_github("paodan/funcTools")  
library(funcTools)  

#Install keras package # keras version 2.2.5.0  
install.packages("keras")  
or   
devtools::install_github("rstudio/keras")  
library(keras)   

# Download dataset for model training
https://uu.figshare.com/articles/dataset/Dataset_for_protein_complexes_project/22043474
This website contains three datasets (training dataset, testing dataset, and validation dataset) and a deep learning model (the best model used in the manuscirpt to predict the protein-protein interaction).

Then download and open the 30_00_model_training_evaluation.R in Rscript and run the code for reproducing  
