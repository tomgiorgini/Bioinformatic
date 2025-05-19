rm(list = ls())
options(stringAsFactors = F)

library(stringr)
library(pheatmap)

setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic")

dirRes <- "Result/"

if(!dir.exists(dirRes)) {
  dir.create(dirRes)
} else {
  print(paste("The directory", dirRes, "already exists"))
}

tumor <- "BRCA"

dirTumor <- paste0(dirRes, tumor, "/")

if(!dir.exists(dirTumor)) {
  dir.create(dirTumor)
} else {
  print(paste("The directory", dirTumor, "already exists"))
}

filename_in <- paste0("matrix_RNAseq_", tumor, ".txt")

filename_list_normal <- paste0(dirTumor, "normal.txt")
filename_list_tumor <- paste0(dirTumor, "tumor.txt")
filename_DEG <- paste0(dirTumor, "DEG.txt")
filename_matrix_DEG <- paste0(dirTumor, "matrix_DEG.txt")
filename_heatmap <- paste0(dirTumor, "heatmap.pdf")

##parameters

prc_IQR <- 0.1
thr_FC <- 2
thr_pval <- 0.05
paired <- TRUE


## STEP 0

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, 
                  row.names = 1, quote = "", nrows = 10)

classes <- sapply(tmp, class)

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, 
                  row.names = 1, quote = "", colClasses = classes)

genes <- rownames(tmp)
pz <- colnames(tmp)

pzN <- grep("TCGA-\\w+-\\w+-1\\d", pz, value = TRUE)
pzC <- grep("TCGA-\\w+-\\w+-0\\d", pz, value = TRUE)

unique_pzN <- pzN[!duplicated(str_extract(pzN, "TCGA-\\w+-\\w+"))]
unique_pzC <- pzC[!duplicated(str_extract(pzC, "TCGA-\\w+-\\w+"))]

nameN <- str_extract(unique_pzN, "TCGA-\\w+-\\w+")
nameC <- str_extract(unique_pzC, "TCGA-\\w+-\\w+")

common_pz <- intersect(nameN, nameC)

pzN_com <- sapply(common_pz, function(x){grep(x, unique_pzN, value = TRUE)})
pzC_com <- sapply(common_pz, function(x){grep(x, unique_pzC, value = TRUE)})

dataN <- tmp[, pzN_com]
dataC <- tmp[, pzC_com]

data <- cbind(dataN, dataC)
pz_com <- colnames(data)

#rm(tmp, pz, nameN, nameC, pzN, pzC, common_pz, classes)

## STEP 1

# 1.1 Mean of the data

ovr_mean <- rowMeans(data, na.rm = T)
ind <- which(ovr_mean == 0)

if (length((ind) > 0)) {
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
}

rm(ind) 

# 1.2 Log2 the data

data <- log2(data + 1)
dataN <- log2(dataN + 1)
dataC <- log2(dataC + 1)

# 1.2 IQR FILTERING

variation <- apply(data, 1, IQR)

thr_IQR <- quantile(variation, prc_IQR)
ind <- which(variation <= thr_IQR)

if (length((ind) > 0)) {
  dataN <- dataN[-ind,]
  dataC <- dataC[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
}

rm(ind)

#############################################################################
## STEP 2

# 2.1 LogFC computation

logFC <- rowMeans(dataC, na.rm = T) - rowMeans(dataN, na.rm = T)

hist(logFC, main = "FC (logarithmic) frequency distribution", breaks = 100, 
     xlab = "log2FC", ylab = "Frequency", col = "red")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = "grey")

# 2.2 P-Value computation

N <- ncol(dataN)
M <- ncol(dataC)

pval <- apply(data, 1, function(x) {
  res <- t.test(x[1:N], x[(N+1):(M+N)], paired = paired)
  pval <- res$p.value
}) 

# 2.3 P-value adjustment

adj_pval <- p.adjust(pval, method = "fdr")

# 2.4 Volcano plot before filtering

plot(logFC, -log10(adj_pval),
     main = "Volcano plot before filtering",
     xlab = "log2 Fold Change (FC)",
     ylab = "-log10 p-value"
)

abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = "blue")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = "red")

# 2.5 Filtering: removing the genes lower than the threshold

ind = which(abs(logFC) < log2(thr_FC))

if(length(ind)>0){
  data = data[-ind,]
  dataC = dataC[-ind,]
  dataN = dataN[-ind,]
  genes = genes[-ind]
  logFC = logFC[-ind]
  pval = pval[-ind]
  adj_pval = adj_pval[-ind]
}

rm(ind)

ind = which(adj_pval >= thr_pval)

if(length(ind)>0){
  data = data[-ind,]
  dataC = dataC[-ind,]
  dataN = dataN[-ind,]
  genes = genes[-ind]
  logFC = logFC[-ind]
  pval = pval[-ind]
  adj_pval = adj_pval[-ind]
}

rm(ind)

#plot after filter 

plot(logFC, -log10(adj_pval),
     main = "Volcano plot after filtering",
     xlab = "log2 Fold Change (FC)",
     ylab = "-log10 p-value"
)

abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = "blue")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = "red")


##############################################################################
#Exporting results

direction = ifelse(logFC>0, "UP", "DOWN")

DEG = data.frame(str_split_fixed(genes, "\\|",2))
colnames(DEG) = c("geneSymbol","ensembl_id")

result = data.frame(GeneSymbol = DEG$geneSymbol,
                    ensembl_id =DEG$ensembl_id,
                    pval = pval, adj_pval = adj_pval, 
                    logFC = logFC, direction = direction)
result = result[order(result$logFC, decreasing = T),]

write.table(result, file=filename_DEG, row.names = F,
            sep = "\t", quote = F)

write.table(data, file=filename_matrix_DEG, row.names = T, col.names = NA,
            sep = "\t", quote = F)

write.table(pzN_com, file=filename_list_normal, row.names = F, col.names = F,
            sep = "\t", quote = F)

write.table(pzC_com, file=filename_list_tumor, row.names = F, col.names = F,
            sep = "\t", quote = F)


##########################################################################
## PLOTTING

#BOX PLOT
index = which.max(logFC)

gene_id = genes[index]

par(mar=c(5,4,4,2)+0.1)
# 3. Fai il boxplot
boxplot(
  as.numeric(dataN[index, ]),
  as.numeric(dataC[index, ]),
  main  = paste0(
    DEG$geneSymbol[index], ", adjusted p-val = ",
    format(adj_pval[index], digits = 2)
  ),
  notch = TRUE,
  ylab  = "Gene expression value",
  xlab  = "Condition",
  names = c("normal", "cancer"),
  col   = c("green", "orange"),
  pars  = list(boxwex = 0.3, staplewex = 0.6),
  cex.lab  = 1.2,    # ingrandisce le label degli assi
  cex.axis = 1       # ingrandisce i tick/names
)

#PIE CHART  
count = table(result$direction)
pie(count, labels = paste(names(count)), col = c("blue","gold"))

#HEATMAP

test = grepl('TCGA-\\w+-\\w+-1\\d',colnames(data))
condition = ifelse(test, "normal", "cancer")
annotation = data.frame(condition=condition)
rownames(annotation) = colnames(data)
vect_col = c("green","orange")
names(vect_col)=unique(condition)
annotation_colors = list(condition =vect_col )

pheatmap(data,scale = "row", border_colors = NA, cluster_cols = T, cluster_rows = T,
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         clustering_method = "complete",
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         color = colorRampPalette(colors = c("blue", "blue3","black","yellow3","yellow"))(100),
         show_rownames = F,
         show_colnames = F,
         cutree_rows = 2,
         cutree_cols = 2,
         width= 10, height = 10,
         filename = filename_heatmap)
