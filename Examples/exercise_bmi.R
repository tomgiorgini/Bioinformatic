

rm(list= ls())
options(StringsAsFactors = F)

###############################
source("computeBMI.R")
#############################
pz <- paste("PZ", 1:100, sep = "_")
height <- sample(150:190, 100, replace = T)
weight = sample(50:90, 100, replace = T)

BMI = computeBMI(height,weight)

df = data.frame(height, weight, BMI)
row.names(df) = pz

write.table(df,"df_patients.txt", sep = "\t", col.names = NA, row.names = T)