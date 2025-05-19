library(ggbiplot)
library(qcc)
library(ggpubr)
library(factoextra)
library(corrplot)
library(FactoMineR)
library(RColorBrewer)
library(survminer)

options(stringsAsFactors = F)
################################################
# 1. Importing data
data("myeloma")
data <- myeloma
data <- data[,-c(2:5)]

condition1 <- "Hyperdiploid" # more similar to normal cell
condition2 <- "MMSET" # multimyeloid cells

data1 <- data[data$molecular_group == condition1,]
data2 <- data[data$molecular_group == condition2,]

data <- rbind(data1,data2)
groups <- as.character(data[,"molecular_group"])
data <- data[,-1]
################################################
# 2. Apply PCA
# Rows of data correspond to observations (samples), columns to variables (geni)

pca <- prcomp(data, center = T, scale. = T, retx = T) 
################################################
# 3. Compute score and score plot
# (scores = the coordinates of old data (observations) in the new systems, that are the PCs)

# pca$x = t(data)*pca$rotation 
scores <- pca$x 

# alternative for score computation
# scores <- get_pca_ind(pca)$coord
# colnames(scores) <- paste0('PC', seq(1,ncol(scores)))

# pdf("pca.pdf",width=5, height=5) 
g <- ggbiplot(pca, 
              obs.scale = 1, 
              var.axes = F, 
              ellipse = T, 
              groups = groups)
print(g)
# dev.off() 

# alternative score plot
fviz_pca_ind(pca,
             col.ind = groups, # "cos2"
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE,
             repel = TRUE, # Avoid text overlapping
)

################################################
# 4. Compute eigenvalue
# eigenvalues of the covariance matrix ordered in decreasing order (from the largest to the smallest)
eigenvalue = pca$sdev^2 

# variance explained by each PC
varS <- round(eigenvalue/sum(eigenvalue)*100, 2)
names(varS) = paste0('PC', seq(1,length(varS)))

# pareto chart
pareto.chart(varS)

# scree plot
fviz_eig(pca,addlabels = TRUE)

# alternative for computation of eigenvalues and explained and cumulative variance
df <- get_eigenvalue(pca)
rownames(df) <- paste0('PC', seq(1,nrow(df)))
################################################
# 5. Compute loadings (coefficient of each PC)
# the matrix of variable loadings (a matrix whose columns contain the eigenvectors)
# matrice di trasformazione dalle vecchie alle nuove coordinate

loadings <- pca$rotation

# loadings plot
# high cos2 = good representation of the variable on the principal component. 
# (i.e., the variable is positioned close to the circumference of the correlation circle)
# low cos2 = the variable is not perfectly represented by the PCs
# (i.e., the variable is close to the center of the circumference of the correlation circle)
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

fviz_pca_var(pca, col.var = "cos2", #"contrib" # Color by the quality of representation or by contribution of each avriable to the PC1-2
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
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

corrplot(contrib_var, is.corr=FALSE, 
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





