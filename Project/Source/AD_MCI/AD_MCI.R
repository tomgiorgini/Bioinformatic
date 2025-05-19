rm(list = ls())
options(stringAsFactors = F)

library(stringr)
library(pheatmap)

setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic/Project")
dirRes <- "Results/"

dirAD_MCI <- paste0(dirRes, "AD_MCI", "/")

if(!dir.exists(dirAD_MCI)) {
  dir.create(dirAD_MCI)
} else {
  print(paste("The directory", dirAD_MCI, "already exists"))
}

filename_in <- "Data/matrixAD_MCI.txt"

filename_list_ad <-"Data/caseAD.txt"
filename_list_mci <-"Data/caseMCI.txt"

filename_DEG <- paste0(dirAD_MCI, "DEG.txt")
filename_matrix_DEG <- paste0(dirAD_MCI, "matrix_DEG.txt")
filename_heatmap <- paste0(dirAD_MCI, "heatmap.pdf")


###################################################################

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, 
                  row.names = 1, quote = "")

classes <- sapply(tmp, class)

tmp <- read.table(filename_in, header = T, sep = "\t", check.names = F, 
                  row.names = 1, quote = "", colClasses = classes)

genes <- rownames(tmp)

ad <- read.table(filename_list_ad, header =F, sep = "\t", check.names = F, quote = "")
mci <- read.table(filename_list_mci, header = F, sep = "\t", check.names = F, quote = "")

col_names <- ad[[1]]
dataAD <- tmp[, col_names]

col_names <- mci[[1]]
dataMCI <- tmp[, col_names]
data <- cbind(dataAD, dataMCI)

rm(tmp)

## STEP 1

# 1.1 Mean of the data

ovr_mean <- rowMeans(data, na.rm = T)
ind <- which(ovr_mean == 0)

if (length((ind) > 0)) {
  dataAD <- dataAD[-ind,]
  dataMCI <- dataMCI[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
}

rm(ind) 

# 1.2 Log2 the data

data <- log2(data + 1)
dataAD <- log2(dataAD + 1)
dataMCI <- log2(dataMCI + 1)

# 1.2 IQR FILTERING
prc_IQR <- 0.1

variation <- apply(data, 1, IQR)

thr_IQR <- quantile(variation, prc_IQR)
ind <- which(variation <= thr_IQR)

if (length((ind) > 0)) {
  dataAD <- dataAD[-ind,]
  dataMCI <- dataMCI[-ind,]
  data <- data[-ind,]
  genes <- genes[-ind]
}

rm(ind)



#############################################################################
## STEP 2

##parameters

thr_FC <- 1
thr_pval <- 0.05
paired <- FALSE

# 2.1 LogFC computation

logFC <- rowMeans(dataMCI, na.rm = T) - rowMeans(dataAD, na.rm = T)

hist(logFC, main = "FC (logarithmic) frequency distribution", breaks = 1000, 
     xlab = "log2FC", ylab = "Frequency", col = "red")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = "grey")

# 2.2 P-Value computation

N <- ncol(dataAD)
M <- ncol(dataMCI)

pval <- apply(data, 1, function(x) {
  res <- t.test(x[1:N], x[(N+1):(M+N)], paired = paired)
  pval <- res$p.value
}) 

# 2.3 P-value adjustment

adj_pval <- p.adjust(pval, method = "fdr")

# 2.4 Volcano plot before filtering

pdf("Results/AD_MCI/VOLCANO_PRE_FILTERING.pdf")

plot(logFC, -log10(adj_pval),
     main = "Volcano plot before filtering",
     xlab = "log2 Fold Change (FC)",
     ylab = "-log10 p-value"
)

abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = "blue")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = "red")

dev.off()
# 2.5 Filtering: removing the genes lower than the threshold

ind = which(adj_pval >= thr_pval | abs(logFC) < log2(thr_FC))

if(length(ind)>0){
  data = data[-ind,]
  dataMCI = dataMCI[-ind,]
  dataAD = dataAD[-ind,]
  genes = genes[-ind]
  logFC = logFC[-ind]
  pval = pval[-ind]
  adj_pval = adj_pval[-ind]
}

rm(ind)

#plot after filter 


pdf("Results/AD_MCI/VOLCANO_POST_FILTERING.pdf")

plot(logFC, -log10(adj_pval),
     main = "Volcano plot after filtering",
     xlab = "log2 Fold Change (FC)",
     ylab = "-log10 p-value"
)

abline(h = -log10(thr_pval), lty = 2, lwd = 2, col = "blue")
abline(v = c(-log2(thr_FC), log2(thr_FC)), lty = 2, lwd = 2, col = "red")

dev.off()

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
control_samples <- as.character(ad[[1]])
case_samples <- as.character(mci[[1]])

# Creiamo i data frame con l’annotazione
df_control <- data.frame(sample = control_samples, condition = "AD")
df_case    <- data.frame(sample = case_samples, condition = "MCI")

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

pdf("Results/AD_MCI/boxplot.pdf")

index = which.max(logFC)

gene_id = genes[index]

boxplot(
  as.numeric(dataMCI[index, ]),
  as.numeric(dataAD[index, ]),
  main  = paste0(
    DEG$geneSymbol[index], ", adjusted p-val = ",
    format(adj_pval[index], digits = 2)
  ),
  notch = TRUE,
  ylab  = "Gene expression value",
  xlab  = "Condition",
  names = c("Alzheimer Disease", "Mild Cognitive Impairment"),
  col   = c("green", "orange"),
  pars  = list(boxwex = 0.3, staplewex = 0.6),
  cex.lab  = 1.2,    # ingrandisce le label degli assi
  cex.axis = 1       # ingrandisce i tick/names
)

dev.off()

