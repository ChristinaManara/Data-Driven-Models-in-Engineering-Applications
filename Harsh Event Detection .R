library(dplyr)
library(DMwR2)
library(DDoutlier)
library(Rcpp)
library(abodOutlier)
library(factoextra)
library(cluster)

df <- read.csv('/Users/christinamanara/Desktop/Data Science & Machine Learning/Data Driven /3rd Assigment/all_LAT_uns.csv')
#scaled_data<scale(data)
scaled_df <- as.data.frame(scale(df))

########### Elbow Method ##########

#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(scaled_df, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


###########  LOF ###############

# Find outliers by setting an optional k
outlier_score <- LOF(scaled_df, k=10)

# Sort and find index for most outlying observations
names(outlier_score) <- 1:nrow(scaled_df)
sort(outlier_score, decreasing = TRUE)

# Inspect the distribution of outlier scores
hist(outlier_score)


########## Kmeans ###############

# Make this example reproducible
set.seed(1)

# Perform k-means clustering with k = 4 clusters
km <- kmeans(df, centers = 10, nstart = 25)

# View results
km

########## KNN  ############
library(FNN)
knn <- get.knn(data = scaled_df, k = 5)
scaled_df$score <- rowMeans(knn$nn.dist)
head(scaled_df, 4)


