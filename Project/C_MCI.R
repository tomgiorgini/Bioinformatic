rm(list = ls())
options(stringAsFactors = F)

library(stringr)
library(pheatmap)

setwd("C:/Users/tomma/OneDrive/Documenti/GitHub/Bioinformatic/Project")
dirRes <- "Results/"

dirC_MCI <- paste0(dirRes, "C_MCI", "/")

if(!dir.exists(dirC_MCI)) {
  dir.create(dirC_MCI)
} else {
  print(paste("The directory", dirC_MCI, "already exists"))
}

filename_in <- "Data/matrixC_MCI.txt"

filename_list_case <-"Data/caseMCI.txt"
filename_list_control <-"Data/control.txt"

filename_DEG <- paste0(dirC_MCI, "DEG.txt")
filename_matrix_DEG <- paste0(dirC_MCI, "matrix_DEG.txt")
filename_heatmap <- paste0(dirC_MCI, "heatmap.pdf")


###################################################################

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, 
                  row.names = 1, quote = "")

classes <- sapply(tmp, class)

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, 
                  row.names = 1, quote = "", colClasses = classes)

genes <- rownames(tmp)

control <- read.table(filename_list_control, header =F, sep = "\t", check.names = F, quote = "")
case <- read.table(filename_list_case, header = F, sep = "\t", check.names = F, quote = "")

col_names <- control[[1]]
dataN <- tmp[, col_names]

col_names <- case[[1]]
dataC <- tmp[, col_names]
data <- cbind(dataN, dataC)

rm(tmp)

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
prc_IQR <- 0.1

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

##parameters

thr_FC <- 1.15
thr_pval <- 0.05
paired <- FALSE

# 2.1 LogFC computation

logFC <- rowMeans(dataC, na.rm = T) - rowMeans(dataN, na.rm = T)

hist(logFC, main = "FC (logarithmic) frequency distribution", breaks = 1000, 
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

ind = which(adj_pval >= thr_pval | abs(logFC) < log2(thr_FC))

if(length(ind)>0){
  data_f = data[-ind,]
  dataC_f = dataC[-ind,]
  dataN_f = dataN[-ind,]
  genes_f = genes[-ind]
  logFC_f = logFC[-ind]
  pval_f = pval[-ind]
  adj_pval_f = adj_pval[-ind]
}

rm(ind)

#plot after filter 

plot(logFC_f, -log10(adj_pval_f),
     main = "Volcano plot after filtering",
     xlab = "log2 Fold Change (FC)",
     ylab = "-log10 p-value"
)

abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = "blue")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = "red")


data = data_f
dataC = dataC_f
dataN= dataN_f
genes= genes_f
logFC = logFC_f
pval = pval_f
adj_pval = adj_pval_f

rm(data_f,dataC_f,dataN_f,genes_f,logFC_f,pval_f,adj_pval_f)



#EXPORTING RESULTS
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


# Convertiamo in vettori di caratteri
control_samples <- as.character(control[[1]])
case_samples <- as.character(case[[1]])

# Creiamo i data frame con l’annotazione
df_control <- data.frame(sample = control_samples, condition = "control")
df_case    <- data.frame(sample = case_samples, condition = "case")

# Uniamo i due data frame
annotation <- rbind(df_control, df_case)

# Usiamo i nomi dei campioni come rownames
rownames(annotation) <- annotation$sample

samples = annotation$sample

annotation = data.frame(condition=annotation$condition,row.names = samples)
vect_col = c("green","orange")
names(vect_col)=unique(annotation$condition)
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
