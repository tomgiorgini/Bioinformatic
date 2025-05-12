rm(list=ls())


options(stringAsFactors = F)

setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic/GEO")
library("GEOquery")

dirRes <- "Results/"
if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The
directory",dirRes,"already exists"))
}
dataset <- "GEO_HCC"
dirDataset <- paste0(dirRes,dataset,"/")
if (!dir.exists(dirDataset)){
  dir.create(dirDataset)
}else{
  print(paste("The
directory",dirDataset,"already exists"))
}

###########################################
# STEP 1: Downloading data
###########################################

series = "GSE22058"

tmp <- getGEO(GEO = series) #show(tmp) mostra i metodi su tmp

set <- tmp$`GSE22058-GPL6793_series_matrix.txt.gz`

pData <- phenoData(set)

metadata <- pData@data

aData <- assayData(set)

matrix <- data.frame(aData$exprs)

rm(pData,aData,tmp)

###########################################
# STEP 2: Preparing data
###########################################
annotation <- fData(set)

geneSymbol <- annotation$GeneSymbol

matrix <- matrix[annotation$ID,]

matrix <-aggregate(matrix,list(geneSymbol),"mean")
ind <- which(matrix$Group.1 == "")
matrix <- matrix[-ind,]

rownames(matrix) <- matrix$Group.1
matrix <- matrix[,-1]

rm(set,annotation,ind,geneSymbol)


############################################
#STEP 3: Extracting tumor and normal samples
###########################################
metadata <-metadata[,c("geo_accession","individual:ch1","tissue:ch1")]
list <- split(metadata,metadata$`tissue:ch1`)

normal <- list$`adjacent liver non-tumor`
tumor <- list$`liver tumor`

normal <-normal[!duplicated(normal$`individual:ch1`),]
tumor <-tumor[!duplicated(tumor$`individual:ch1`),]
normal <-normal[order(normal$`individual:ch1`),"geo_accession"]
tumor <-tumor[order(tumor$`individual:ch1`),"geo_accession"]

data <- matrix[,c(normal,tumor)]
dataN <- matrix[,normal]
dataC <- matrix[,tumor]
metadata <- metadata[c(normal,tumor),]

rm(list,matrix)


############################################
#STEP 4: Export data
###########################################
write.table(data, paste0(dirDataset,"matrix.txt"), sep= "\t", col.names = NA, row.names = T, quote = F)

write.table(normal,paste0(dirDataset,"normal.txt"), sep= "\t", col.names = F,row.names = F, quote = F)

write.table(tumor,paste0(dirDataset,"tumor.txt"), sep= "\t", col.names = F, row.names = F, quote = F)

write.table(metadata,paste0(dirDataset,"metadata.txt"),sep= "\t", col.names = T, row.names =F, quote = F)

# Remove downloaded data
unlink(series, recursive = TRUE)

