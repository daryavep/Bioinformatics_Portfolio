# I wanted to clean up the tables for easier data manipulation

# Table 3
table3 <- read.csv(file="table3.csv", header=TRUE)

# Join together sub-headers
colnames(table3)[4:15] <- paste("mRNA.expression.levels", colnames(table3)[4:15], sep = "_")
colnames(table3)[16:21] <- paste("Protein.expression.levels", colnames(table3)[16:21], sep = "_")
table3 <- table3[-c(1, 1), ]

# Change column names for easier analysis
colnames(table3) <- c("ID", "uniprot", "protein.name", "mRNAC1", "mRNAC2", "mRNAC3", "mRNAC4", "mRNAC5", "mRNAC6", "mRNAH1", "mRNAH2", "mRNAH3", "mRNAH4", "mRNAH5", "mRNAH6", "proteinc1", "proteinc2", "proteinc3", "proteinh1", "proteinh2", "proteinh3")

# Table 4
table4 <- read.csv(file="table4.csv", header=TRUE)

# Join together sub-headers
colnames(table4)[4:15] <- paste("mRNA.expression.levels", colnames(table3)[4:15], sep = "_")
colnames(table4)[16:21] <- paste("Protein.expression.levels", colnames(table3)[16:21], sep = "_")
table4 <- table4[-c(1, 1), ]

# Change column names for easier analysis
colnames(table4) <- c("ID", "uniprot", "protein.name", "mRNAC1", "mRNAC2", "mRNAC3", "mRNAC4", "mRNAC5", "mRNAC6", "mRNAH1", "mRNAH2", "mRNAH3", "mRNAH4", "mRNAH5", "mRNAH6", "proteinc1", "proteinc2", "proteinc3", "proteinh1", "proteinh2", "proteinh3")

# Convert columns to numeric
expression_columns <- 4:21
table3[, expression_columns] <- sapply(table3[, expression_columns], as.numeric)
table4[, expression_columns] <- sapply(table4[, expression_columns], as.numeric)

# I realized only after the fact that LATTE provided .txt files of the tables, but I suppose manipulating excel files was good practice.

######################################################################
#####################      QUESTION 1     ############################

# Define column ranges as chimp or human for both mRNA and protein expression levels
mRNA.chimp.columns <- c(4:9)
mRNA.human.columns <- c(10:15)
protein.chimp.columns <- c(16:18)
protein.human.columns <- c(19:21)

# Calculate mean values for chimp and human mRNA and protein expression levels for each row in "table3"
mean.mRNA.chimp.table3 <- rowMeans(table3[, mRNA.chimp.columns])
mean.mRNA.human.table3 <- rowMeans(table3[, mRNA.human.columns])
mean.protein.chimp.table3 <- rowMeans(table3[, protein.chimp.columns])
mean.protein.human.table3 <- rowMeans(table3[, protein.human.columns])

# Create scatter plots for "table3"
par(mfrow = c(1, 2))
plot(mean.mRNA.human.table3, mean.mRNA.chimp.table3, xlab = "Human mRNA", ylab = "Chimp mRNA", main = "Sample table3 - mRNA")
plot(mean.protein.human.table3, mean.protein.chimp.table3, xlab = "Human protein", ylab = "Chimp protein", main = "Sample table3 - Protein")


# Calculate mean values for chimp and human mRNA and protein expression levels for each row in "table4"
mean.mRNA.chimp.table4 <- rowMeans(table4[, mRNA.chimp.columns])
mean.mRNA.human.table4 <- rowMeans(table4[, mRNA.human.columns])
mean.protein.chimp.table4 <- rowMeans(table4[, protein.chimp.columns])
mean.protein.human.table4 <- rowMeans(table4[, protein.human.columns])

# Create scatter plots for "table4"
par(mfrow = c(1, 2))
plot(mean.mRNA.human.table4, mean.mRNA.chimp.table4, xlab = "Human mRNA", ylab = "Chimp mRNA", main = "Sample table4 - mRNA")
plot(mean.protein.human.table4, mean.protein.chimp.table4, xlab = "Human protein", ylab = "Chimp protein", main = "Sample table4 - Protein")


######################################################################
#####################      QUESTION 2     ############################

# Calculate correlations for "table3"
cor.mRNA.chimp.human.table3 <- cor(mean.mRNA.chimp.table3, mean.mRNA.human.table3)
cor.protein.chimp.human.table3 <- cor(mean.protein.chimp.table3, mean.protein.human.table3)

# Calculate correlations for "table4"
cor.mRNA.chimp.human.table4 <- cor(mean.mRNA.chimp.table4, mean.mRNA.human.table4)
cor.protein.chimp.human.table4 <- cor(mean.protein.chimp.table4, mean.protein.human.table4)

# Compare and describe the results
cat("Correlations for mean mRNA levels between chimps and humans in table3:", cor.mRNA.chimp.human.table3, "\n")
cat("Correlations for mean protein levels between chimps and humans in table3:", cor.protein.chimp.human.table3, "\n")

cat("Correlations for mean mRNA levels between chimps and humans in table4:", cor.mRNA.chimp.human.table4, "\n")
cat("Correlations for mean protein levels between chimps and humans in table4:", cor.protein.chimp.human.table4, "\n")

######################################################################
#####################      QUESTION 3     ############################

# Fit linear models for mean mRNA and mean protein relationships in table3
model.mRNA.table3 <- lm(mean.mRNA.human.table3 ~ mean.mRNA.chimp.table3)
model.protein.table3 <- lm(mean.protein.human.table3 ~ mean.protein.chimp.table3)

# Fit linear models for mean mRNA and mean protein relationships in table4
model.mRNA.table4 <- lm(mean.mRNA.human.table4 ~ mean.mRNA.chimp.table4)
model.protein.table4 <- lm(mean.protein.human.table4 ~ mean.protein.chimp.table4)

# Summary and model diagnostics for table3 mean mRNA relationship
summary(model.mRNA.table3)
par(mfrow = c(1, 2))
plot(model.mRNA.table3)
anova(model.mRNA.table3)

# Summary and model diagnostics for table3 mean protein relationship
summary(model.protein.table3)
plot(model.protein.table3)
anova(model.protein.table3)

# Summary and model diagnostics for table4 mean mRNA relationship
summary(model.mRNA.table4)
plot(model.mRNA.table4)
anova(model.mRNA.table4)

# Summary and model diagnostics for table4 mean protein relationship
summary(model.protein.table4)
plot(model.protein.table4)
anova(model.protein.table4)

# Reset the plotting parameters to default
par(mfrow = c(1, 1))

######################################################################
#####################      QUESTION 4     ############################

# Scatter plot for mean mRNA vs mean protein levels in table3 for chimp
par(mfrow = c(2, 2))
plot(mean.mRNA.chimp.table3, mean.protein.chimp.table3, 
     xlab = "Mean Chimp mRNA", ylab = "Mean Chimp Protein", 
     main = "Chimp - Sample table3", pch = 19, col = "blue")

# Scatter plot for mean mRNA vs mean protein levels in table3 for human
plot(mean.mRNA.human.table3, mean.protein.human.table3, 
     xlab = "Mean Human mRNA", ylab = "Mean Human Protein", 
     main = "Human - Sample table3", pch = 19, col = "red")

# Scatter plot for mean mRNA vs mean protein levels in table4 for chimp
plot(mean.mRNA.chimp.table4, mean.protein.chimp.table4, 
     xlab = "Mean Chimp mRNA", ylab = "Mean Chimp Protein", 
     main = "Chimp - Sample table4", pch = 19, col = "blue")

# Scatter plot for mean mRNA vs mean protein levels in table4 for human
plot(mean.mRNA.human.table4, mean.protein.human.table4, 
     xlab = "Mean Human mRNA", ylab = "Mean Human Protein", 
     main = "Human - Sample table4", pch = 19, col = "red")

# Reset the plotting parameters to default
par(mfrow = c(1, 1))


######################################################################
#####################      QUESTION 5     ############################

# Calculate differences in mean mRNA and protein levels for table3
mRNA.difference.table3 <- mean.mRNA.human.table3 - mean.mRNA.chimp.table3
protein.difference.table3 <- mean.protein.human.table3 - mean.protein.chimp.table3

# Calculate differences in mean mRNA and protein levels for table4
mRNA.difference.table4 <- mean.mRNA.human.table4 - mean.mRNA.chimp.table4
protein.difference.table4 <- mean.protein.human.table4 - mean.protein.chimp.table4

# Plot differences for table3
par(mfrow = c(1, 2))
plot(mRNA.difference.table3, protein.difference.table3, 
     xlab = "Difference in mRNA Levels", ylab = "Difference in Protein Levels", 
     main = "Sample table3", pch = 19, col = "red")

# Plot differences for table4
plot(mRNA.difference.table4, protein.difference.table4, 
     xlab = "Difference in mRNA Levels", ylab = "Difference in Protein Levels", 
     main = "Sample table4", pch = 19, col = "blue")

# Reset the plotting parameters to default
par(mfrow = c(1, 1))


######################################################################
#####################      QUESTION 6     ############################

# Calculate correlation coefficients for differences between human and chimp protein and mRNA levels in table3
cor.diff.table3 <- cor(mRNA.difference.table3, protein.difference.table3)

# Calculate correlation coefficients for differences between human and chimp protein and mRNA levels in table4
cor.diff.table4 <- cor(mRNA.difference.table4, protein.difference.table4)

# Display the correlation coefficients
cat("Correlation coefficient for differences in table3:", cor.diff.table3, "\n")
cat("Correlation coefficient for differences in table4:", cor.diff.table4, "\n")

######################################################################
#####################      QUESTION 7     ############################

# Fit a linear model for the relationship between differences in table3
model.table3 <- lm(mRNA.difference.table3 ~ protein.difference.table3)

# Fit a linear model for the relationship between differences in table4
model.table4 <- lm(mRNA.difference.table4 ~ protein.difference.table4)

# Summary and model diagnostics for table3
summary(model.table3)
par(mfrow = c(1, 2))
plot(model.table3)
anova(model.table3)

# Summary and model diagnostics for table4
summary(model.table4)
plot(model.table4)
anova(model.table4)

# Reset the plotting parameters to default
par(mfrow = c(1, 1))

######################################################################
#####################      QUESTION 8     ############################

# Function to perform t-tests for mRNA or protein levels
t.test.hc <- function(x, mRNA.human.columns, mRNA.chimp.columns, protein.human.columns, protein.chimp.columns, var.equal = TRUE) {
  p.value.mRNA <- t.test(x[mRNA.human.columns], x[mRNA.chimp.columns], var.equal = var.equal)$p.value
  p.value.protein <- t.test(x[protein.human.columns], x[protein.chimp.columns], var.equal = var.equal)$p.value
  return(c(p.value.mRNA, p.value.protein))
}

# Initialize vectors to store p-values
p.values.mRNA.table3 <- numeric(nrow(table3))
p.values.protein.table3 <- numeric(nrow(table3))
p.values.mRNA.table4 <- numeric(nrow(table4))
p.values.protein.table4 <- numeric(nrow(table4))

count.p.values.mRNA.table3 <- 0
count.p.values.protein.table3 <- 0
count.p.values.mRNA.table4 <- 0
count.p.values.protein.table4 <- 0

# Loop through the rows and calculate p-values for each gene in table3
for (i in 1:nrow(table3)) {
  p.values <- t.test.hc(table3[i, ], mRNA.human.columns, mRNA.chimp.columns, protein.human.columns, protein.chimp.columns, var.equal = TRUE)
  p.values.mRNA.table3[i] <- p.values[1]
  p.values.protein.table3[i] <- p.values[2]
  
  # Check if p-values are under 0.01 and increment the counters
  if (p.values[1] < 0.01) {
    count.p.values.mRNA.table3 <- count.p.values.mRNA.table3 + 1
  }
  if (p.values[2] < 0.01) {
    count.p.values.protein.table3 <- count.p.values.protein.table3 + 1
  }
}

# Loop through the rows and calculate p-values for each gene in table4
for (i in 1:nrow(table4)) {
  p.values <- t.test.hc(table4[i, ], mRNA.human.columns, mRNA.chimp.columns, protein.human.columns, protein.chimp.columns, var.equal = TRUE)
  p.values.mRNA.table4[i] <- p.values[1]
  p.values.protein.table4[i] <- p.values[2]
  
  # Check if p-values are under 0.01 and increment the counters
  if (p.values[1] < 0.01) {
    count.p.values.mRNA.table4 <- count.p.values.mRNA.table4 + 1
  }
  if (p.values[2] < 0.01) {
    count.p.values.protein.table4 <- count.p.values.protein.table4 + 1
  }
}

# Print the counts of p-values under 0.01
output <- paste("Number of p-values under 0.01 for mRNA in table3:", count.p.values.mRNA.table3, "\n",
                "Number of p-values under 0.01 for protein in table3:", count.p.values.protein.table3, "\n",
                "Number of p-values under 0.01 for mRNA in table4:", count.p.values.mRNA.table4, "\n",
                "Number of p-values under 0.01 for protein in table4:", count.p.values.protein.table4)

# Print the combined output
cat(output)

# Plot histograms of p-values for mRNA and protein levels in table3
par(mfrow = c(1, 2))
hist(p.values.mRNA.table3, main = "mRNA - table3", xlab = "p-value", col = "lightblue")
hist(p.values.protein.table3, main = "Protein - table3", xlab = "p-value", col = "lightgreen")

# Plot histograms of p-values for mRNA and protein levels in table4
hist(p.values.mRNA.table4, main = "mRNA - table4", xlab = "p-value", col = "lightblue")
hist(p.values.protein.table4, main = "Protein - table4", xlab = "p-value", col = "lightgreen")

# Create X-Y scatterplot of p-values (protein vs. mRNA) for table3
plot(p.values.mRNA.table3, p.values.protein.table3, xlab = "p-value mRNA", ylab = "p-value Protein", main = "P-Values - table3", pch = 19, col = "red")

# Create X-Y scatterplot of p-values (protein vs. mRNA) for table4
plot(p.values.mRNA.table4, p.values.protein.table4, xlab = "p-value mRNA", ylab = "p-value Protein", main = "P-Values - table4", pch = 19, col = "blue")

# Reset the plotting parameters to default
par(mfrow = c(1, 1))

######################################################################
#####################      QUESTION 9     ############################

# Function to perform t-tests for mRNA or protein levels
t.test.hc <- function(x, mRNA.human.columns, mRNA.chimp.columns, protein.human.columns, protein.chimp.columns, var.equal = FALSE) {
  p.value.mRNA <- t.test(x[mRNA.human.columns], x[mRNA.chimp.columns], var.equal = var.equal)$p.value
  p.value.protein <- t.test(x[protein.human.columns], x[protein.chimp.columns], var.equal = var.equal)$p.value
  return(c(p.value.mRNA, p.value.protein))
}

# Initialize vectors to store p-values
p.values.mRNA.table3 <- numeric(nrow(table3))
p.values.protein.table3 <- numeric(nrow(table3))
p.values.mRNA.table4 <- numeric(nrow(table4))
p.values.protein.table4 <- numeric(nrow(table4))

count.p.values.mRNA.table3 <- 0
count.p.values.protein.table3 <- 0
count.p.values.mRNA.table4 <- 0
count.p.values.protein.table4 <- 0

# Loop through the rows and calculate p-values for each gene in table3
for (i in 1:nrow(table3)) {
  p.values <- t.test.hc(table3[i, ], mRNA.human.columns, mRNA.chimp.columns, protein.human.columns, protein.chimp.columns, var.equal = FALSE)
  p.values.mRNA.table3[i] <- p.values[1]
  p.values.protein.table3[i] <- p.values[2]
  
  # Check if p-values are under 0.01 and increment the counters
  if (p.values[1] < 0.01) {
    count.p.values.mRNA.table3 <- count.p.values.mRNA.table3 + 1
  }
  if (p.values[2] < 0.01) {
    count.p.values.protein.table3 <- count.p.values.protein.table3 + 1
  }
}

# Loop through the rows and calculate p-values for each gene in table4
for (i in 1:nrow(table4)) {
  p.values <- t.test.hc(table4[i, ], mRNA.human.columns, mRNA.chimp.columns, protein.human.columns, protein.chimp.columns, var.equal = FALSE)
  p.values.mRNA.table4[i] <- p.values[1]
  p.values.protein.table4[i] <- p.values[2]
  
  # Check if p-values are under 0.01 and increment the counters
  if (p.values[1] < 0.01) {
    count.p.values.mRNA.table4 <- count.p.values.mRNA.table4 + 1
  }
  if (p.values[2] < 0.01) {
    count.p.values.protein.table4 <- count.p.values.protein.table4 + 1
  }
}

# Print the counts of p-values under 0.01
output <- paste("Number of p-values under 0.01 for mRNA in table3:", count.p.values.mRNA.table3, "\n",
                "Number of p-values under 0.01 for protein in table3:", count.p.values.protein.table3, "\n",
                "Number of p-values under 0.01 for mRNA in table4:", count.p.values.mRNA.table4, "\n",
                "Number of p-values under 0.01 for protein in table4:", count.p.values.protein.table4)

# Print the combined output
cat(output)

# Plot histograms of p-values for mRNA and protein levels in table3
par(mfrow = c(1, 2))
hist(p.values.mRNA.table3, main = "mRNA - table3", xlab = "p-value", col = "lightblue")
hist(p.values.protein.table3, main = "Protein - table3", xlab = "p-value", col = "lightgreen")

# Plot histograms of p-values for mRNA and protein levels in table4
hist(p.values.mRNA.table4, main = "mRNA - table4", xlab = "p-value", col = "lightblue")
hist(p.values.protein.table4, main = "Protein - table4", xlab = "p-value", col = "lightgreen")

# Create X-Y scatterplot of p-values (protein vs. mRNA) for table3
plot(p.values.mRNA.table3, p.values.protein.table3, xlab = "p-value mRNA", ylab = "p-value Protein", main = "P-Values - table3", pch = 19, col = "red")

# Create X-Y scatterplot of p-values (protein vs. mRNA) for table4
plot(p.values.mRNA.table4, p.values.protein.table4, xlab = "p-value mRNA", ylab = "p-value Protein", main = "P-Values - table4", pch = 19, col = "blue")

# Reset the plotting parameters to default
par(mfrow = c(1, 1))

#############################################################################
# ANOVA TESTING

# ANOVA function for mRNA
anova.hc.mRNA <- function(x) {
  mRNA.anova <- x[, c(mRNA.human.columns, mRNA.chimp.columns)]
  species <- factor(c(rep("Human",length(mRNA.human.columns)), rep("Chimp",length(mRNA.chimp.columns))))
  lm.model <- lm(c(unlist(mRNA.anova)) ~ species)
  anova(lm.model)[1, "Pr(>F)"]
}

# ANOVA function for protein
anova.hc.protein <- function(x) {
  protein.anova <- x[, c(protein.human.columns, protein.chimp.columns)]
  species <- factor(c(rep("Human",length(protein.human.columns)), rep("Chimp",length(protein.chimp.columns))))
  lm.model <- lm(c(unlist(protein.anova)) ~ species)
  anova(lm.model)[1, "Pr(>F)"]
}

# Initialize vectors to store p-values
p.values.mRNA.table3 <- numeric(nrow(table3))
p.values.protein.table3 <- numeric(nrow(table3))
p.values.mRNA.table4 <- numeric(nrow(table4))
p.values.protein.table4 <- numeric(nrow(table4))

count.p.values.mRNA.table3 <- 0
count.p.values.protein.table3 <- 0
count.p.values.mRNA.table4 <- 0
count.p.values.protein.table4 <- 0

# Loop through the rows and calculate p-values for each gene in table3
for (i in 1:nrow(table3)) {
  # Perform ANOVA on linear models
  p.values <- c(anova.hc.mRNA(table3[i, ]), anova.hc.protein(table3[i, ]))
  p.values.mRNA.table3[i] <- anova.hc.mRNA(table3[i, ])
  p.values.protein.table3[i] <- anova.hc.protein(table3[i, ])
  # Check if p-values are under 0.01 and increment the counters
  if (p.values[1] < 0.01) {
    count.p.values.mRNA.table3 <- count.p.values.mRNA.table3 + 1
  }
  if (p.values[2] < 0.01) {
    count.p.values.protein.table3 <- count.p.values.protein.table3 + 1
  }
}

# Loop through the rows and calculate p-values for each gene in table4
for (i in 1:nrow(table4)) {
  # Perform ANOVA on linear models
  p.values <- c(anova.hc.mRNA(table4[i, ]), anova.hc.protein(table4[i, ]))
  p.values.mRNA.table4[i] <- anova.hc.mRNA(table4[i, ])
  p.values.protein.table4[i] <- anova.hc.protein(table4[i, ])
  # Check if p-values are under 0.01 and increment the counters
  if (p.values[1] < 0.01) {
    count.p.values.mRNA.table4 <- count.p.values.mRNA.table4 + 1
  }
  if (p.values[2] < 0.01) {
    count.p.values.protein.table4 <- count.p.values.protein.table4 + 1
  }
}


# Print the counts of p-values under 0.01
output <- paste("Number of p-values under 0.01 for mRNA in table3:", count.p.values.mRNA.table3, "\n",
                "Number of p-values under 0.01 for protein in table3:", count.p.values.protein.table3, "\n",
                "Number of p-values under 0.01 for mRNA in table4:", count.p.values.mRNA.table4, "\n",
                "Number of p-values under 0.01 for protein in table4:", count.p.values.protein.table4)

# Print the combined output
cat(output)

# Plot histograms of p-values for mRNA and protein levels in table3
par(mfrow = c(1, 2))
hist(p.values.mRNA.table3, main = "mRNA - table3", xlab = "p-value", col = "lightblue")
hist(p.values.protein.table3, main = "Protein - table3", xlab = "p-value", col = "lightgreen")

# Plot histograms of p-values for mRNA and protein levels in table4
hist(p.values.mRNA.table4, main = "mRNA - table4", xlab = "p-value", col = "lightblue")
hist(p.values.protein.table4, main = "Protein - table4", xlab = "p-value", col = "lightgreen")

# Create X-Y scatterplot of p-values (protein vs. mRNA) for table3
plot(p.values.mRNA.table3, p.values.protein.table3, xlab = "p-value mRNA", ylab = "p-value Protein", main = "P-Values - table3", pch = 19, col = "red")

# Create X-Y scatterplot of p-values (protein vs. mRNA) for table4
plot(p.values.mRNA.table4, p.values.protein.table4, xlab = "p-value mRNA", ylab = "p-value Protein", main = "P-Values - table4", pch = 19, col = "blue")

# Reset the plotting parameters to default
par(mfrow = c(1, 1))


######################################################################
#####################      QUESTION 10     ###########################

# Create a list of common gene IDs between table3 and table4
common_genes <- intersect(table3$ID, table4$ID)

# Subset the data for common genes
common_genes_table3 <- table3[table3$ID %in% common_genes, ]
common_genes_table4 <- table4[table4$ID %in% common_genes, ]

# Calculate the mean mRNA and protein levels for human and chimp
human_mRNA_table3 <- rowMeans(common_genes_table3[, c(mRNA.human.columns)])
human_mRNA_table4 <- rowMeans(common_genes_table4[, c(mRNA.human.columns)])
chimp_mRNA_table3 <- rowMeans(common_genes_table3[, c(mRNA.chimp.columns)])
chimp_mRNA_table4 <- rowMeans(common_genes_table4[, c(mRNA.chimp.columns)])
human_protein_table3 <- rowMeans(common_genes_table3[, c(protein.human.columns)])
human_protein_table4 <- rowMeans(common_genes_table4[, c(protein.human.columns)])
chimp_protein_table3 <- rowMeans(common_genes_table3[, c(protein.chimp.columns)])
chimp_protein_table4 <- rowMeans(common_genes_table4[, c(protein.chimp.columns)])

# Create scatterplots
par(mfrow = c(2, 2))
plot(human_mRNA_table3, human_mRNA_table4, xlab = "Human Mean mRNA - Table 3", ylab = "Human Mean mRNA - Table 4", main = "Human Mean mRNA Levels")
plot(chimp_mRNA_table3, chimp_mRNA_table4, xlab = "Chimp Mean mRNA - Table 3", ylab = "Chimp Mean mRNA - Table 4", main = "Chimp Mean mRNA Levels")
plot(human_protein_table3, human_protein_table4, xlab = "Human Mean Protein - Table 3", ylab = "Human Mean Protein - Table 4", main = "Human Mean Protein Levels")
plot(chimp_protein_table3, chimp_protein_table4, xlab = "Chimp Mean Protein - Table 3", ylab = "Chimp Mean Protein - Table 4", main = "Chimp Mean Protein Levels")
par(mfrow = c(1, 1))

######################################################################
#####################      QUESTION 11     ###########################

# Create a list of common gene IDs between table3 and table4
common_genes <- intersect(table3$ID, table4$ID)

# Initialize vectors to store p-values
p_values_human_mRNA <- numeric(length(common_genes))
p_values_human_protein <- numeric(length(common_genes))
p_values_chimp_mRNA <- numeric(length(common_genes))
p_values_chimp_protein <- numeric(length(common_genes))

# Loop through the common genes and calculate p-values
for (i in 1:length(common_genes)) {
  # Extract the current common gene ID
  current_gene <- common_genes[i]
  
  # Find the row indices for the current gene in both tables
  index_table3 <- which(common_genes_table3$ID == current_gene)
  index_table4 <- which(common_genes_table4$ID == current_gene)
  
  p_values_human_mRNA[i] <- t.test(common_genes_table3[index_table3, mRNA.human.columns], common_genes_table4[index_table4, mRNA.human.columns])$p.value
  p_values_human_protein[i] <- t.test(common_genes_table3[index_table3, protein.human.columns], common_genes_table4[index_table4, protein.human.columns])$p.value
  p_values_chimp_mRNA[i] <- t.test(common_genes_table3[index_table3, mRNA.chimp.columns], common_genes_table4[index_table4, mRNA.chimp.columns])$p.value
  p_values_chimp_protein[i] <- t.test(common_genes_table3[index_table3, protein.chimp.columns], common_genes_table4[index_table4, protein.chimp.columns])$p.value
}

# Print human mRNA data for a specific row (replace 'index_table3' with the desired row index)
cat("Human mRNA data for row", index_table3, "in common_genes_table3:\n")
print(common_genes_table3[index_table3, mRNA.human.columns])


# Create 2x2 layout for the plots
par(mfrow = c(2, 2))

# Plot p-values for human mRNA expression in table 3 vs 4
hist(p_values_human_mRNA, main = "Human mRNA Expression - Table 3 vs 4", xlab = "p-value", col = "lightblue")

# Plot p-values for human protein expression in table 3 vs 4
hist(p_values_human_protein, main = "Human Protein Expression - Table 3 vs 4", xlab = "p-value", col = "lightgreen")

# Plot p-values for chimp mRNA expression in table 3 vs 4
hist(p_values_chimp_mRNA, main = "Chimp mRNA Expression - Table 3 vs 4", xlab = "p-value", col = "lightblue")

# Plot p-values for chimp protein expression in table 3 vs 4
hist(p_values_chimp_protein, main = "Chimp Protein Expression - Table 3 vs 4", xlab = "p-value", col = "lightgreen")

# Reset the plotting parameters to default
par(mfrow = c(1, 1))

######################################################################
#####################      QUESTION 12     ###########################

# PART A:

# Load ALL dataset
library(ALL)
data(ALL)

# Create a data frame with 'fusion protein' and '1970_s_at' columns (taken from notes part 1)
fp.df <- data.frame(
  fp = pData(ALL)[, "fusion protein"],
  g = exprs(ALL)["1970_s_at", ]
)

# Calculate the mean gene expression for each group
means.by.group <- tapply(fp.df$g, fp.df$fp, mean)

# Calculate the coefficients manually
coefficient.intercept <- means.by.group[1]
coefficient.p190 <- means.by.group[2] - coefficient.intercept
coefficient.p210 <- means.by.group[3] - coefficient.intercept

# Print the coefficients
cat("Manual Coefficients:\n")
cat("Intercept:", coefficient.intercept, "\n")
cat("Coefficient for p190:", coefficient.p190, "\n")
cat("Coefficient for p210:", coefficient.p210, "\n")



# Run ANOVA
anova_result <- anova(lm(g ~ fp, data = fp.df))
cat("ANOVA Result:\n")
print(anova_result)

# Extract p-value from ANOVA
p_value_anova <- anova_result[["Pr(>F)"]][1]
cat("\nP-Value from ANOVA:", p_value_anova, "\n")

# Pairwise t-tests comparing 'g' across levels of 'fp' with equal variances assumed
t_test_p190 <- t.test(fp.df$g[fp.df$fp == "p190"], fp.df$g[fp.df$fp != "p190"], var.equal = TRUE)
t_test_p210 <- t.test(fp.df$g[fp.df$fp == "p210"], fp.df$g[fp.df$fp != "p210"], var.equal = TRUE)

cat("\nT-Test Results for p190 (assuming equal variances):\n")
print(t_test_p190)

cat("\nT-Test Results for p210 (assuming equal variances):\n")
print(t_test_p210)
