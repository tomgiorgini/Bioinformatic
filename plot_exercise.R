#plot with abline
plot(cars$speed,cars$dist, xlab = "Speed", ylab = "Distance", col = "blue")
title(main = "Scatter plot with best-fit line")
abline(lm(dist ~ speed, data= cars), col= "red") #fit-line

# different type of points
library(MASS)
plot(Cars93$Weight, Cars93$EngineSize,
     col = c("black","red","green","blue","cyan","magenta"),
     ylab = "EngineSize",
     xlab = "Weight",
     main = "Example of plot")

#barplot

barplot(table(Cars93$Type),col = c("black","red","green","blue","cyan","magenta"))

barplot(VADeaths, beside = T, legend = T)

#histogram
hist(Cars93$Weight, breaks = 5)

#boxplot

boxplot(Cars93$Horsepower)
table(Cars93)


