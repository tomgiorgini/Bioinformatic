# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
rm(list=ls())
library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)
################################################
setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic")
source("enrichR_UP_DOWN/getEnrichment.R")
source("enrichR_UP_DOWN/getEnrichmentPlot.R")
################################################
dataset <- "BRCA"

dirRes <- "Result/"

if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The directory",dirRes,"already exists"))
}

dirDataset <- paste0(dirRes,dataset,"/")

if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The directory",dirDataset,"already exists"))
}

dirEnrich <- paste0(dirDataset,"Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}


################################################
top_term <- 10
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirDataset,"DEG.txt")

dbs <- listEnrichrDbs() #lista di tutti i db

dbs <- c("DisGeNET","GO_Molecular_Function_2025", "GO_Biological_Process_2025", "KEGG_2021_Human", "TRANSFAC_and_JASPAR_PWMs")

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")
# input_list <- input_list$genes
list <- split(input_list$GeneSymbol,input_list$direction)

df <- lapply(list, function(x){
  enrichr(x, dbs)
})

getEnrichment(df$UP,"UP")
getEnrichment(df$DOWN,"DOWN")


