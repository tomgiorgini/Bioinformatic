# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
rm(list=ls())
library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)
################################################
setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic/Project")
source("Source/getEnrichment.R")
source("Source/getEnrichmentPlot.R")
################################################

dirRes <- "Results/"

dirC_AD <- paste0(dirRes, "C_AD", "/")

if(!dir.exists(dirC_AD)) {
  dir.create(dirC_AD)
} else {
  print(paste("The directory", dirC_AD, "already exists"))
}

dirEnrich <- paste0(dirC_AD,"Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}


################################################
top_term <- 10
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirC_AD,"DEG.txt")

dbs <- listEnrichrDbs() #lista di tutti i db

dbs <- c("DisGeNET","GO_Molecular_Function_2025", "GO_Biological_Process_2025", "KEGG_2021_Human", "TRANSFAC_and_JASPAR_PWMs")

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")
# input_list <- input_list$genes
list <- split(input_list$GeneSymbol,input_list$direction)

df <- lapply(list, function(x){
  enrichr(x, dbs)
})

groups <- names(df)

lapply(groups, function(g) {
  tag <- gsub("\\.", "_", g)
  getEnrichment(df[[g]], tag)
})
