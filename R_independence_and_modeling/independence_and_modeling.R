# Step 1: Create simulated dataset
set.seed(123) # For reproducibility

# Define x and y
x <- 1:20
y <- 20 + 10 * x + 5 * x^2 + rnorm(20, 0, 200)

# Combine x and y into a data frame
simulated.data <- data.frame(x, y)

# Step 2: Fit linear models
# Fit models
model1 <- lm(y ~ 1, data = simulated.data)
model2 <- lm(y ~ x, data = simulated.data)
model3 <- lm(y ~ x + I(x^2), data = simulated.data)

# Step 3: Obtain observed and fitted values using predict()
# Predict values from the models
simulated.data$y.pred1 <- predict(model1)
simulated.data$y.pred2 <- predict(model2)
simulated.data$y.pred3 <- predict(model3)

#Step 4: Obtain residual values
# Calculate residuals
simulated.data$residuals1 <- resid(model1)
simulated.data$residuals2 <- resid(model2)
simulated.data$residuals3 <- resid(model3)

# Step 5: Plot the observed and fitted data points as well as the residual values

# Plot the observed vs. fitted values for all models
plot1 <- ggplot(simulated.data, aes(x = x, y = y)) +
  geom_point(aes(color = "Observed"), size = 3) +
  geom_line(aes(x = x, y = y.pred1, color = "Model 1"), linewidth = 1) +
  geom_line(aes(x = x, y = y.pred2, color = "Model 2"), linewidth = 1) +
  geom_line(aes(x = x, y = y.pred3, color = "Model 3"), linewidth = 1) +
  labs(color = "Legend") +
  theme_minimal() +
  ggtitle("Observed vs. Fitted Values")
print(plot1)

# Plot the residual values for all models
plot2 <- ggplot(simulated.data, aes(x = x, y = residuals1)) +
  geom_point(aes(color = "Model 1 Residuals"), size = 3) +
  theme_minimal() +
  ggtitle("Model 1 Residuals")
print(plot2)

plot3 <- ggplot(simulated.data, aes(x = x, y = residuals2)) +
  geom_point(aes(color = "Model 2 Residuals"), size = 3) +
  theme_minimal() +
  ggtitle("Model 2 Residuals")
print(plot3)

plot4 <- ggplot(simulated.data, aes(x = x, y = residuals3)) +
  geom_point(aes(color = "Model 3 Residuals"), size = 3) +
  theme_minimal() +
  ggtitle("Model 3 Residuals")
print(plot4)

####################################################################
############## Question 2

library(ALL)
data(ALL)
library(ggplot2)

# Create an empty data frame to store p-values
p.values <- data.frame(Gene = character(0), PValue = numeric(0))
# Loop through genes and fit linear models, then extract p-values
for (gene in 1:length(rownames(exprs(ALL)))) {
  model <- lm(pData(ALL)$age ~ exprs(ALL)[gene,])
  anova.result <- anova(model)
  p.value <- anova.result$"Pr(>F)"[1]  # Extract p-value for gene effect on age
  p.values <- rbind(p.values, data.frame(Gene = gene, PValue = p.value))
}

# Plot the distribution of p-values
ggplot(p.values, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  labs(title = "Distribution of P-Values for Gene Effect on Age",
       x = "P-Value") +
  theme_minimal()

# Find the gene with the most and least significant p-values
most.significant.gene <- p.values[p.values$PValue == min(p.values$PValue), ]
least.significant.gene <- p.values[p.values$PValue == max(p.values$PValue), ]

# Fit linear models for the most and least significant genes
most.significant.model <- lm(age ~ exprs(ALL)[most.significant.gene$Gene, ], data = ALL)
least.significant.model <- lm(age ~ exprs(ALL)[least.significant.gene$Gene, ], data = ALL)

# Generate diagnostic plots for the most significant gene
par(mfrow = c(2, 2))
plot(most.significant.model)

# Generate diagnostic plots for the least significant gene
par(mfrow = c(2, 2))
plot(least.significant.model)

# Model summaries and ANOVA results for the most significant gene
summary(most.significant.model)
anova(most.significant.model)

# Model summaries and ANOVA results for the least significant gene
summary(least.significant.model)
anova(least.significant.model)

#####################################################################
############## Question 3

# Set a seed for reproducibility
set.seed(123)

# Number of samples to generate
n.samples <- 10000
sample.size <- 10

# Initialize an empty vector to store p-values
p.values3 <- numeric(n.samples)

# Perform t-test and store p-values
for (i in 1:n.samples) {
  sample1 <- rnorm(sample.size)
  sample2 <- rnorm(sample.size)
  t.test.result <- t.test(sample1, sample2)
  p.values3[i] <- t.test.result$p.value
}

# Plot the distribution of p-values
hist(p.values3, main = "Distribution of p-values from t-test", xlab = "p-value")

# Count how many p-values are below 0.05
significant.p.values3 <- sum(p.values3 < 0.05)
cat("Number of p-values below 0.05:", significant.p.values3, "\n")

#####################################################################
############ Question 4

# Load necessary items
library(ALL)
data(ALL)
library(ggplot2)

# Create empty data frames to store p-values
gene.p.values <- data.frame(Gene = character(0), ANOVA.PValue = numeric(0))
residual.p.values <- data.frame(Gene = character(0), Shapiro.Wilk.PValue = numeric(0))

# Calculate time to remission in days
ALL.pdat <- pData(ALL)
date.cr.chr <- as.character(ALL.pdat$date.cr)
diag.chr <- as.character(ALL.pdat$diagnosis)
date.cr.t <- strptime(date.cr.chr,"%m/%d/%Y")
diag.t <- strptime(diag.chr,"%m/%d/%Y")
days2remiss <- as.numeric(date.cr.t - diag.t)

# Loop through genes and fit linear models for gene effects on remission
for (gene in 1:length(rownames(exprs(ALL)))) {
  # Create a data frame for the current gene
  model_data <- data.frame(gene.expression = exprs(ALL)[gene, ], days2remiss = days2remiss)
  
  # Fit the linear model
  model <- lm(days2remiss ~ gene.expression, data = model_data)
  
  # Extract ANOVA p-value for gene effect
  anova_result <- anova(model)
  anova_p_value <- anova_result$`Pr(>F)`[1]
  gene.p.values <- rbind(gene.p.values, data.frame(Gene = gene, ANOVA.PValue = anova_p_value))
  
  # Shapiro-Wilk test for normality of residuals
  residuals <- residuals(model)
  shapiro_wilk_p_value <- shapiro.test(residuals)$p.value
  residual.p.values <- rbind(residual.p.values, data.frame(Gene = gene, Shapiro.Wilk.PValue = shapiro_wilk_p_value))
}

# Plot the distribution of ANOVA p-values for gene effects
ggplot(gene.p.values, aes(x = ANOVA.PValue)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  labs(title = "Distribution of ANOVA P-Values for Gene Effects on Remission",
       x = "ANOVA P-Value") +
  theme_minimal()

# Plot the distribution of Shapiro-Wilk test p-values for normality of residuals
ggplot(residual.p.values, aes(x = Shapiro.Wilk.PValue)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black") +
  labs(title = "Distribution of Shapiro-Wilk Test P-Values for Normality of Model Residuals",
       x = "Shapiro-Wilk P-Value") +
  theme_minimal()

#################################################################
########### DEBUG AREA ##########################################

# Load necessary items
library(ALL)
data(ALL)
library(ggplot2)

# Create empty data frames to store p-values
gene.p.values <- data.frame(Gene = character(0), ANOVA.PValue = numeric(0))
residual.p.values <- data.frame(Gene = character(0), Shapiro.Wilk.PValue = numeric(0))

# Loop through genes and fit linear models for gene effects on remission
for (gene in 1:length(rownames(exprs(ALL)))) {
  model_data <- ALL
  model_data$gene.expression <- exprs(ALL)[gene, ]
  
  # Remove rows with missing values in gene.expression and remission
  model_data <- na.omit(model_data)
  model_data$remission <- as.numeric(model_data$remission)
  
  # Check if there are enough data points for analysis
  if (nrow(model_data) >= 3) {
    model <- lm(remission ~ gene.expression, data = model_data)
    
    # Extract ANOVA p-value for gene effect
    anova_result <- anova(model)
    anova_p_value <- anova_result$`Pr(>F)`[1]
    gene.p.values <- rbind(gene.p.values, data.frame(Gene = gene, ANOVA.PValue = anova_p_value))
    
    # Shapiro-Wilk test for normality of residuals
    residuals <- residuals(model)
    shapiro_wilk_p_value <- shapiro.test(residuals)$p.value
    residual.p.values <- rbind(residual.p.values, data.frame(Gene = gene, Shapiro.Wilk.PValue = shapiro_wilk_p_value))
  }
}

# Plot the distribution of ANOVA p-values for gene effects
ggplot(gene.p.values, aes(x = ANOVA.PValue)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  labs(title = "Distribution of ANOVA P-Values for Gene Effects on Remission",
       x = "ANOVA P-Value") +
  theme_minimal()

# Plot the distribution of Shapiro-Wilk test p-values for normality of residuals
ggplot(residual.p.values, aes(x = Shapiro.Wilk.PValue)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black") +
  labs(title = "Distribution of Shapiro-Wilk Test P-Values for Normality of Model Residuals",
       x = "Shapiro-Wilk P-Value") +
  theme_minimal()

#################################################################
########### DEBUG AREA ##########################################

# Load necessary items
library(ALL)
data(ALL)
library(ggplot2)

# Calculate remission based on 
date.cr.chr <- as.character(ALL.pdat$date.cr)
diag.chr <- as.character(ALL.pdat$diagnosis)
date.cr.t <- as.Date(date.cr.chr, "%m/%d/%Y")
diag.t <- as.Date(diag.chr, "%m/%d/%Y")
days2remiss <- as.numeric(date.cr.t - diag.t)

# Create empty data frames to store p-values
gene.p.values <- data.frame(Gene = character(0), ANOVA.PValue = numeric(0))
residual.p.values <- data.frame(Gene = character(0), Shapiro.Wilk.PValue = numeric(0))

# Loop through genes and fit linear models for gene effects on remission
for (gene in 1:length(rownames(exprs(ALL))) ) {
  # Put current gene's expression levels and days-to-remission into a dataframe
  df.tmp <- data.frame(G = exprs(ALL)[gene, ], D2R = days2remiss)
  
  # Fit the model with those data
  lm.tmp <- lm(D2R ~ G, df.tmp)
  
  # Extract ANOVA p-value for gene effect
  p.anova <- anova(lm.tmp)["G", "Pr(>F)"]
  
  # Calculate and save Shapiro p-value for the distribution of the residuals
  p.shapiro <- shapiro.test(resid(lm.tmp))$p.value
  
  # Save the p-values in the data frames
  gene.p.values <- rbind(gene.p.values, data.frame(Gene = gene, ANOVA.PValue = p.anova))
  residual.p.values <- rbind(residual.p.values, data.frame(Gene = gene, Shapiro.Wilk.PValue = p.shapiro))
}

# Plot the distribution of ANOVA p-values for gene effects
ggplot(gene.p.values, aes(x = ANOVA.PValue)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  labs(title = "Distribution of ANOVA P-Values for Gene Effects on Remission",
       x = "ANOVA P-Value") +
  theme_minimal()

# Plot the distribution of Shapiro-Wilk test p-values for normality of residuals
ggplot(residual.p.values, aes(x = Shapiro.Wilk.PValue)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black") +
  labs(title = "Distribution of Shapiro-Wilk Test P-Values for Normality of Model Residuals",
       x = "Shapiro-Wilk P-Value") +
  theme_minimal()


#################################################################
########### DEBUG AREA ##########################################

ALL.pdat <- pData(ALL)
# the code below is explained in the notes; we are just calculating the values of
# time to remission in each patient here:
date.cr.chr <- as.character(ALL.pdat$date.cr)
diag.chr <- as.character(ALL.pdat$diagnosis)
date.cr.t <- strptime(date.cr.chr,"%m/%d/%Y")
diag.t <- strptime(diag.chr,"%m/%d/%Y")
days2remiss <- as.numeric(date.cr.t - diag.t) # done, we finally have days2remiss
# let’s now fit the model for each gene separately (note that we chose to use a loop;
# you could define an appropriate function and use apply() instead!):
p.anova <- numeric()
p.shapiro <- numeric()
for ( i.row in 1:dim(exprs(ALL))[[1]] ) { # for each row (i.e. gene)
  # put current gene’s expression levels and days-to-remission into a dataframe:
  df.tmp <- data.frame(G=exprs(ALL)[i.row,],D2R=days2remiss)
  lm.tmp <- lm(D2R~G,df.tmp) # fit the model with those data
  p.anova[i.row] <- anova(lm.tmp)["G","Pr(>F)"] # save anova p-value for that fit
  # calculate and save Shapiro p-value for the distribution of the residuals:
  p.shapiro[i.row] <- shapiro.test(resid(lm.tmp))$p.value
}
# draw plots:
old.par <- par(mfrow=c(1,3),ps=20)
plot(hist(p.anova,plot=F),main="ANOVA for gene effect",xlab="p-value")
plot(hist(p.shapiro,plot=F),main="Shapiro-Wilk on the residuals",xlab="p-value")
plot(p.shapiro,p.anova,log="xy",xlab="p(Shapiro-Wilk)",ylab="p(ANOVA)")



#################################################################
########## DEBUG AREA ###########################################
# Load necessary libraries
library(ALL)
library(ggplot2)

# Load the ALL dataset
data(ALL)

# Subset the data to the first 100 genes for debugging
subset.ALL <- ALL[, 1:5]

# Create histograms for the first 100 genes
for (gene in rownames(exprs(subset.ALL))) {
  data <- exprs(subset.ALL)[gene, ]
  
  # Create a histogram
  hist(data, breaks = "FD", col = "lightblue", main = paste("Histogram of Gene:", gene))
  
  # Print the plot to display it
  print(last_plot())
}

