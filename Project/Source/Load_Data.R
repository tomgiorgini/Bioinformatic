rm(list=ls())


options(stringAsFactors = F)

setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic/Project")
library("GEOquery")

dirRes <- "Results/"
if (!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste("The
directory",dirRes,"already exists"))
}
dataset <- "GEO_ALZHEIMER"
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

series = "GSE63060"

tmp <- getGEO(GEO = series) #show(tmp) mostra i metodi su tmp

set <- tmp$GSE63060_series_matrix.txt.gz

pData <- phenoData(set)

metadata <- pData@data

aData <- assayData(set)

matrix <- data.frame(aData$exprs)

matrix<- 2^matrix

rm(pData,aData,tmp)


###########################################
# STEP 2: Preparing data
###########################################

annotation <- fData(set)

geneSymbol <- annotation$ILMN_Gene

matrix <- matrix[annotation$ID,]

matrix <-aggregate(matrix,list(geneSymbol),"mean")

ind <- which(matrix$Group.1 == "")

if(length(ind) > 0) matrix <- matrix[-ind,]

rownames(matrix) <- matrix$Group.1

matrix <- matrix[,-1]

matrix[is.na(matrix)] <- 0

rm(set,annotation,ind,geneSymbol)

###########################################
#STEP 3: Extracting case and control samples
###########################################

list <- split(metadata$geo_accession,metadata$`status:ch1`)

# example: AD
control <- list$CTL
caseAD <- list$AD # for retrieving instead MCI
caseMCI <- list$MCI

dataN <- matrix[,control]
dataAD <- matrix[,caseAD]
dataMCI <- matrix[,caseMCI]

dataC_AD <- matrix[,c(control,caseAD)]
dataC_MCI <- matrix[,c(control,caseMCI)]
dataAD_MCI <- matrix[,c(caseAD,caseMCI)]
metadata = metadata[,c("geo_accession","status:ch1")]
rm(list,matrix)

############################################
#STEP 4: Export data
###########################################

write.table(dataC_AD, paste0(dirDataset,"matrixC_AD.txt"), sep= "\t", col.names = NA, row.names = T, quote = F)
write.table(dataC_MCI, paste0(dirDataset,"matrixC_MCI.txt"), sep= "\t", col.names = NA, row.names = T, quote = F)
write.table(dataAD_MCI, paste0(dirDataset,"matrixAD_MCI.txt"), sep= "\t", col.names = NA, row.names = T, quote = F)

write.table(control,paste0(dirDataset,"control.txt"), sep= "\t", col.names = F,row.names = F, quote = F)

write.table(caseAD,paste0(dirDataset,"caseAD.txt"), sep= "\t", col.names = F, row.names = F, quote = F)

write.table(caseMCI,paste0(dirDataset,"caseMCI.txt"),sep= "\t", col.names = T, row.names =F, quote = F)

write.table(metadata,paste0(dirDataset,"metadata.txt"),sep= "\t", col.names = T, row.names =F, quote = F)
