########################################################################
##################           QUESTION 1          #######################

library(scatterplot3d)

n <- 100
cntr.dist <- 4
x <- c(rnorm(n, cntr.dist), rnorm(n, -cntr.dist), rnorm(n, cntr.dist))
y <- c(rnorm(n, cntr.dist), rnorm(n, -cntr.dist), rnorm(n, 0))
z <- c(rnorm(n, cntr.dist), rnorm(n, -cntr.dist), rnorm(n, 0))

# Combine variables into a data frame
data <- data.frame(x, y, z)

# Visualize the original data with scatterplot3d
scatterplot3d(data, pch = 16, color = "blue", main = "Original 3D Data")

# Run PCA on the data
pca_result <- prcomp(data, center = TRUE, scale. = TRUE)

# Plot Principal Component Plot
plot(pca_result$x, pch = 16, col = "blue", main = "Principal Component Plot")

# Plot Histogram of Variances along Principal Components
barplot(summary(pca_result)$importance[2, ], main = "Variances along Principal Components", xlab = "Principal Components", ylab = "Variance")

# Display parameters of PCA transformation
summary(pca_result)

# Written description
cat("\nPrincipal Components Analysis Results:")
cat("\nProportion of Variance Explained by Each PC:")
print(summary(pca_result)$importance[2, ])
cat("\n\nRotation (Loadings) Matrix:")
print(pca_result$rotation)


########################################################################
##################           QUESTION 2          #######################

library(ALL)
data(ALL)

expression = exprs(ALL)

### Find highest variance
# Calculate variance for each gene
gene_variances <- apply(expression, 1, var)

# Find the indices of the top 100 genes with the highest variance
top100_indices <- rank(-gene_variances, ties.method = "min")<= 100

# Find the indices of the top 1000 genes with the highest variance
top1000_indices <- rank(-gene_variances, ties.method = "min")<= 1000


# Get the gene names for the top 100 and top 1000 genes
top100_genes <- rownames(expression)[top100_indices]
top1000_genes <- rownames(expression)[top1000_indices]

### Run PCA
# Run PCA on top 100 genes (transposed and scaled)
pca_top100 <- prcomp(scale(t(exprs(ALL)[top100_genes,])))

# Run PCA on top 1000 genes (transposed and scaled)
pca_top1000 <- prcomp(scale(t(exprs(ALL)[top1000_genes,])))

# Run PCA on all genes (transposed and scaled)
pca_all <- prcomp(scale(t(exprs(ALL))))

### Find number of principal components
# Function to find the number of principal components for a given percentage of variance
find_components <- function(pca_result, target_variance) {
  cumulative_variance <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
  num_components <- sum(cumulative_variance <= target_variance)
  return(num_components)
}

# Find the number of principal components for 50% and 75% variance for all three PCAs
num_components_50_all <- find_components(pca_all, 0.5)
num_components_75_all <- find_components(pca_all, 0.75)

num_components_50_top100 <- find_components(pca_top100, 0.5)
num_components_75_top100 <- find_components(pca_top100, 0.75)

num_components_50_top1000 <- find_components(pca_top1000, 0.5)
num_components_75_top1000 <- find_components(pca_top1000, 0.75)

# Create a data frame for easier plotting
df <- data.frame(
  Category = rep(c("Top 100", "Top 1000", "All Genes"), each = 2),
  Variance = rep(c("50%", "75%"), times = 3),
  Components = c(
    num_components_50_top100, num_components_75_top100,
    num_components_50_top1000, num_components_75_top1000,
    num_components_50_all, num_components_75_all
  )
)

# Set up a layout for side-by-side histograms
par(mfrow = c(1, 2))

# Plot histograms for 50% variance
barplot(df$Components[df$Variance == "50%"], names.arg = df$Category[df$Variance == "50%"],
        main = "50% Variance",
        xlab = "Gene Set",
        ylab = "Number of Principal Components",
        col = c("lightpink", "lightgreen", "lightblue"))

# Plot histograms for 75% variance
barplot(df$Components[df$Variance == "75%"], names.arg = df$Category[df$Variance == "75%"],
        main = "75% Variance",
        xlab = "Gene Set",
        ylab = "Number of Principal Components",
        col = c("lightpink", "lightgreen", "lightblue"))

# Reset the plotting layout
par(mfrow = c(1, 1))



#######
### Visualize according to B/T disease subtype
# Convert the "BT" attribute to a categorical variable with two levels ("B" and "T")
ALL$BT <- factor(substring(ALL$BT, 1, 1))

# Extract the expression data as a matrix
expression_data <- exprs(ALL)
bt_data <- ALL$BT  # Create a separate vector for BT


# Create a color palette for B and T disease types
colors <- c("B" = "blue", "T" = "red")

# Assign colors based on disease types
point_colors <- colors[bt_data]

# Plot PCA projection with colored points
plot(pca_all$x[, 1], pca_all$x[, 2], col = point_colors, pch = 16,
     main = "PCA Projection of ALL Dataset",
     xlab = "Principal Component 1",
     ylab = "Principal Component 2")

# Add legend
legend("topright", legend = levels(bt_data), fill = colors)

# Plot PCA projection with colored points for top 100 genes
plot(pca_top100$x[, 1], pca_top100$x[, 2], col = point_colors, pch = 16,
     main = "PCA Projection (Top 100 Genes)",
     xlab = "Principal Component 1",
     ylab = "Principal Component 2")

# Add legend
legend("topright", legend = levels(bt_data), fill = colors)

# Plot PCA projection with colored points for top 1000 genes
plot(pca_top1000$x[, 1], pca_top1000$x[, 2], col = point_colors, pch = 16,
     main = "PCA Projection (Top 1000 Genes)",
     xlab = "Principal Component 1",
     ylab = "Principal Component 2")

# Add legend
legend("topright", legend = levels(bt_data), fill = colors)



########################################################################
##################           QUESTION 3          #######################

# Function to perform hierarchical clustering and plot the dendrogram
perform_hierarchical_clustering <- function(gene_set, title) {
  # Subset expression data based on the selected genes
  subset_expression <- expression_data[gene_set, ]
  
  # Calculate Euclidean distances between samples
  distances <- dist(t(subset_expression))
  
  # Perform hierarchical clustering
  hclust_result <- hclust(distances)
  
  # Plot dendrogram
  plot(hclust_result, main = title, xlab = "Samples", ylab = "Distance")
}

# Run hierarchical clustering for top 100 genes
perform_hierarchical_clustering(top100_genes, "Hierarchical Clustering - Top 100 Genes")

# Run hierarchical clustering for top 1000 genes
perform_hierarchical_clustering(top1000_genes, "Hierarchical Clustering - Top 1000 Genes")

# Run hierarchical clustering for all genes
perform_hierarchical_clustering(rownames(expression_data), "Hierarchical Clustering - All Genes")


########################################################################
##################           QUESTION 4          #######################

# Set a random seed for reproducibility
set.seed(123)

# Function to calculate cluster strength for a given gene set using "ward" method
calculate_cluster_strength_ward <- function(gene_set) {
  # Subset expression data based on the selected genes
  subset_expression <- expression_data[gene_set, ]
  
  # Calculate Euclidean distances between samples
  distances <- dist(t(subset_expression))
  
  # Perform hierarchical clustering with the "ward" method
  hclust_result <- hclust(distances, method = "ward")
  
  # Cut the dendrogram to form clusters
  clusters <- cutree(hclust_result, k = length(unique(bt_data)))
  
  # Visualize the dendrogram
  plot(hclust_result, cex = 0.7)
  rect.hclust(hclust_result, 4)
  
  # Return the clusters
  return(clusters)
}

# Filter genes with standard deviation greater than one
high_var_genes <- rownames(expression_data)[apply(expression_data, 1, sd) > 1]

# Run hierarchical clustering and calculate cluster strength for the actual data
actual_clusters_ward <- calculate_cluster_strength_ward(high_var_genes)




# Number of permutations (change to 1000 later)
num_permutations <- 10

# Initialize a matrix to store permutation results
permutation_results_ward <- matrix(NA, nrow = num_permutations, ncol = length(actual_clusters_ward))

# Perform permutation exercise
for (i in 1:num_permutations) {
  # Shuffle disease type labels
  shuffled_labels <- sample(bt_data)
  
  # Run hierarchical clustering and store the clusters for shuffled data
  permutation_results_ward[i, ] <- calculate_cluster_strength_ward(high_var_genes)
}

# Plot the actual clusters
par(mfrow=c(1, 2))
plot(hclust(dist(t(expression_data[high_var_genes, ])), method = "ward"), cex = 0.7)

# Plot the distribution of permuted clusters
plot(1:length(actual_clusters_ward), permutation_results_ward[1, ], col = "gray", pch = 16, ylim = c(0, max(permutation_results_ward)), xlab = "Sample", ylab = "Cluster", main = "Permutation Exercise")
for (i in 2:num_permutations) {
  points(1:length(actual_clusters_ward), permutation_results_ward[i, ], col = "gray", pch = 16)
}

# Add actual clusters to the plot
points(1:length(actual_clusters_ward), actual_clusters_ward, col = "red", pch = 16)
legend("topright", legend = c("Actual", "Permutations"), col = c("red", "gray"), pch = 16)

# Reset plotting layout
par(mfrow=c(1, 1))


###########################################
# DEBUG 4

# Set a random seed for reproducibility
set.seed(123)

# Function to perform permutation exercise for dendrogram heights
permute_dendrogram_heights <- function(expression_data, bt_data, num_permutations = 100) {
  # Heights of the dendrogram from the original data
  ori_heights <- hclust(dist(t(expression_data)), method = "ward")$height
  
  # Initialize a vector to store heights from permutations
  rnd_heights <- numeric()
  
  # Perform permutation exercise
  for (i in 1:num_permutations) {
    # Generate a random expression matrix by permuting samples
    exprs_rnd <- t(apply(expression_data, 1, sample))
    
    # Calculate the dendrogram height for the permuted data
    hclust_rnd <- hclust(dist(t(exprs_rnd)), method = "ward")
    rnd_heights <- c(rnd_heights, hclust_rnd$height)
  }
  
  # Plot the distribution of dendrogram heights
  plot(sort(ori_heights), rank(sort(ori_heights))/length(ori_heights),
       col = "red", xlab = "Height", ylab = "F(Height)")
  
  points(sort(rnd_heights), rank(sort(rnd_heights))/length(rnd_heights), col = "blue")
  
  points(max(ori_heights), 0.9, col = "red")
  text(max(ori_heights), 0.9, "Original", col = "red", pos = 4)
  
  points(max(rnd_heights), 0.8, col = "blue")
  text(max(rnd_heights), 0.8, "Scrambled", col = "blue", pos = 4)
  
  abline(v = c(30, 42), lty = 2, col = "blue")
}

# Perform permutation exercise and plot dendrogram heights
permute_dendrogram_heights(expression_data[high_var_genes, ], bt_data, num_permutations = 100)


##################
# DEBUG Q4
exprs.sd.gt.1 <- exprs(ALL)[apply(exprs(ALL),1,sd)>1.0,]
hclust.sd.gt.1 <- hclust(dist(t(exprs.sd.gt.1)))
ori.heights <- hclust.sd.gt.1$height

rnd.heights <- numeric()
for ( i.sim in 1:100 ) {
  exprs.rnd <- t(apply(exprs.sd.gt.1,1,sample))
  hclust.rnd <- hclust(dist(t(exprs.rnd)), method = "ward.D2")
  
  rnd.heights <- c(rnd.heights,hclust.rnd$height)
}
plot(ori.heights,rank(ori.heights)/length(ori.heights),
     col= "red", xlab="height", ylab="F(height)")
points(rnd.heights,rank(rnd.heights)/length(rnd.heights),col="blue")
points(15,0.9,col="red")
text(15,0.9,"original",col="red",pos=4)
points(15,0.8,col="blue")
text(15,0.8,"scrambled",col="blue",pos=4)
abline(v=c(30,42),lty=2,col="blue")
