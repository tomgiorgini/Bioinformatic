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

dirC_AD <- paste0(dirRes, "COMMON", "/")

if(!dir.exists(dirC_AD)) {
  dir.create(dirC_AD)
} else {
  print(paste("The directory", dirC_AD, "already exists"))
}

dirEnrich <- paste0("Functional_Enrichment/")

if (!dir.exists(dirEnrich)){
  dir.create(dirEnrich)
}else{
  print(paste("The directory",dirEnrich,"already exists"))
}

################################################
top_term <- 10
thr_pval <- 0.05
################################################

c_ad_up = read.table('Results/C_AD/genes_UP.txt', sep = "\t", header = F, check.names = F, quote = "")
c_ad_up = c_ad_up$V1
c_ad_down = read.table("Results/C_AD/genes_DOWN.txt", sep = "\t", header = F, check.names = F, quote = "")
c_ad_down = c_ad_down$V1
c_mci_up = read.table("Results/C_MCI/genes_UP.txt", sep = "\t", header = F, check.names = F, quote = "")
c_mci_up = c_mci_up$V1
c_mci_down = read.table("Results/C_MCI/genes_DOWN.txt", sep = "\t", header = F, check.names = F, quote = "")
c_mci_down = c_mci_down$V1

dbs <- listEnrichrDbs() #lista di tutti i db

dbs <- c("DisGeNET","GO_Molecular_Function_2025", "GO_Biological_Process_2025", "KEGG_2021_Human", "TRANSFAC_and_JASPAR_PWMs")

common_up = Reduce(intersect, list(c_ad_up,c_mci_up))
common_down =  Reduce(intersect, list(c_ad_down,c_mci_down))
c_ad_up = setdiff(c_ad_up,common_up)
c_ad_down = setdiff(c_ad_down,common_down)
c_mci_up = setdiff(c_mci_up,common_up)
c_mci_down = setdiff(c_mci_down,common_down)

list <- list(
  common_up   = common_up,
  common_down = common_down,
  c_ad_up     = c_ad_up,
  c_ad_down   = c_ad_down,
  c_mci_up    = c_mci_up,
  c_mci_down  = c_mci_down
)

df <- lapply(list, function(x){
  enrichr(x, dbs)
})

getEnrichment(df$common_up,"common_UP")
getEnrichment(df$common_down,"common_DOWN")
getEnrichment(df$c_ad_up,"c_ad_UP")
getEnrichment(df$c_ad_down,"c_ad_DOWN")
getEnrichment(df$c_mci_up,"c_mci_UP")
getEnrichment(df$c_mci_down,"c_mci_DOWN")

