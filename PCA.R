library(ggplot2)
library(ggbiplot)

# read data
data <- read.table("TermProjectData.txt", sep = ",", header=TRUE)

# extract numerical data from the data set
electricityData <- data.frame(data$Global_active_power, 
                              data$Global_reactive_power, 
                              data$Voltage, 
                              data$Global_intensity, 
                              data$Sub_metering_1, 
                              data$Sub_metering_2, 
                              data$Sub_metering_3)

# remove observations with NA values in numerical data
electricityData = na.omit(electricityData)
electricityData = scale(electricityData)

# execute PCA on numerical data
pca <- prcomp(electricityData, scale=TRUE)

print(pca$rotation)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

# plot a scree plot for all 7 PCAs
ggscreeplot(pca)
# plot a biplot for PCA1 and PCA2 (remove alpha=0 to see the data points)
ggbiplot(pca, choices=c(1,2), alpha=0)

loading_scores <- pca$rotation[,1]
feature_scores <- abs(loading_scores)
feature_score_ranked <- sort(feature_scores, decreasing=TRUE)
