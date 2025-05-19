# https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html
rm(list=ls())
library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)
################################################
setwd("C:/Users/tomma/OneDrive/Documenti/GitHub/Bioinformatic/Project")
source("getEnrichment.R")
source("getEnrichmentPlot.R")
################################################

dirRes <- "Results/"

dirC_MCI <- paste0(dirRes, "C_MCI", "/")

if(!dir.exists(dirC_MCI)) {
  dir.create(dirC_MCI)
} else {
  print(paste("The directory", dirC_MCI, "already exists"))
}

dirEnrich <- paste0(dirC_MCI,"Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}


################################################
top_term <- 10
thr_pval <- 0.05
################################################
file_input_list <- paste0(dirC_MCI,"DEG.txt")

dbs <- listEnrichrDbs() #lista di tutti i db

dbs <- c("DisGeNET","GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human", "TRANSFAC_and_JASPAR_PWMs")

input_list <- read.table(file_input_list, sep = "\t", header = T, check.names = F, quote = "")
# input_list <- input_list$genes
list <- split(input_list$GeneSymbol,input_list$direction)

df <- lapply(list, function(x){
  enrichr(x, dbs)
})

getEnrichment(df$UP,"UP")
getEnrichment(df$DOWN,"DOWN")
