rm(list=ls())


options(stringAsFactors = F)
library(stringr)
library(ggbiplot)
library(qcc)
library(ggpubr)
library(factoextra)
library(corrplot)
library(FactoMineR)
library(RColorBrewer)
library(survminer)



setwd("C:/Users/tomma/OneDrive/Documenti/GitHub/Bioinformatic/Project")
dirRes <- "Results/"
dirC_MCI <- paste0(dirRes, "C_MCI", "/")

filename_in <- "Results/C_MCI/matrix_DEG.txt"

filename_list_case <-"Data/caseMCI.txt"
filename_list_control <-"Data/control.txt"

file_score_plot = paste0(dirC_MCI, "score_plot.pdf")
file_pareto_scree_plot = paste0(dirC_MCI,"pareto_screen_plot.pdf")
file_contribution_plot = paste0(dirC_MCI,"PC_contibution_plot.pdf")


data = read.table(filename_in, header = T, sep = "\t", quote = "", check.names=F, row.names=1)

name_DEG = data.frame(str_split_fixed(rownames(data),"\\|" ,2))
rownames(data) = name_DEG$X1


list_control = read.table(filename_list_control, header = F, sep = "\t", quote = "", check.names=F)
list_control = list_control$V1

list_case = read.table(filename_list_case, header = F, sep = "\t", quote = "", check.names=F)
list_case = list_case$V1

data = t(data[,c(list_case,list_control)])

groups = c(rep("MCI", length(list_case)), rep("control", length(list_control)))


# 2. Apply PCA
# Rows of data correspond to observations (samples), columns to variables (geni)

pca <- prcomp(data, center = T, scale. = T, retx = T) 

################################################
# 3. Compute score and score plot
# (scores = the coordinates of old data (observations) in the new systems, that are the PCs)

# pca$x = t(data)*pca$rotation 
scores <- pca$x 

# pdf("pca.pdf",width=5, height=5) 
g <- ggbiplot(pca, 
              obs.scale = 1, 
              var.axes = F, 
              ellipse = T, 
              groups = groups)
print(g)

################################################
# 4. Compute eigenvalue
# eigenvalues of the covariance matrix ordered in decreasing order (from the largest to the smallest)
eigenvalue = pca$sdev^2 

# variance explained by each PC
varS <- round(eigenvalue/sum(eigenvalue)*100, 2)
names(varS) = paste0('PC', seq(1,length(varS)))

pdf(file_pareto_scree_plot)

# pareto chart
pareto.chart(varS[1:10])

# scree plot
fviz_eig(pca,addlabels = TRUE)

dev.off()
################################################
# 5. Compute loadings (coefficient of each PC)
# the matrix of variable loadings (a matrix whose columns contain the eigenvectors)
# matrice di trasformazione dalle vecchie alle nuove coordinate

loadings <- pca$rotation


################################################
# Contributions of variables
# The function fviz_contrib() creates a barplot of row/column contributions. 
# A reference dashed line corresponds to the expected value if the contribution were uniform. 
# For a given dimension, any row/column with a contribution above the reference line could be considered 
# as important in contributing to the dimension.

# Contributions of variables to the PCs
contrib_var <- get_pca_var(pca)$contrib 
colnames(contrib_var) <- paste0('PC', seq(1,ncol(contrib_var)))

contrib_var <- contrib_var[order(contrib_var[,"PC1"], decreasing = T),]

pdf(file_contribution_plot, width=10,height=7)

corrplot(contrib_var[1:10,1:10], is.corr=FALSE, 
         tl.col = "black", 
         method = "color",
         col = brewer.pal(n = 10, name = "BuPu"),
         addCoef.col = "black")

# to PC1
fviz_contrib(pca, choice = "var", axes = 1, top = 10)
# to PC2
fviz_contrib(pca, choice = "var", axes = 2, top = 10)
# to PC3
fviz_contrib(pca, choice = "var", axes = 3, top = 10)

# total (PC1 + PC2 + PC3)
fviz_contrib(pca, choice = "var", axes = 1:3, top = 15)

dev.off()

