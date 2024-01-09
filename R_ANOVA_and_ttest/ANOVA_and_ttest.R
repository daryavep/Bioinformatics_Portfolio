# Set up looping numbers
max.sample.size <- 20 
n.sim <- 1000  # number of re-samplings we are going to run for each sample size

# Set up matrices to save data
sd.ori <- matrix(0, nrow = n.sim, ncol = max.sample.size)
sd.n <- matrix(0, nrow = n.sim, ncol = max.sample.size)
sd.n.2 <- matrix(0, nrow = n.sim, ncol = max.sample.size)

# Perform computations in a loop
for (sample.size in 2:max.sample.size) { 
  for (i.sim in 1:n.sim) { 
    # Draw a new random sample of size sample.size from a standard normal distribution
    sample.data <- rnorm(sample.size, mean = 0, sd = 1)
    
    # Calculate "classical" standard deviation (s) for the sample
    sd.ori[i.sim, sample.size] <- sd(sample.data)
    
    # Calculate the first modified estimator (sqrt(n)) for the same sample
    sd.n[i.sim, sample.size] <- sqrt(sample.size) * sd(sample.data)
    
    # Calculate the second modified estimator (sqrt(n-2)) for the same sample
    sd.n.2[i.sim, sample.size] <- sqrt(sample.size - 2) * sd(sample.data)
  }
}

# Calculate the mean and standard deviation of each estimator for each sample size
mean.sd.ori <- colMeans(sd.ori)
sd.sd.ori <- apply(sd.ori, 2, sd)

mean.sd.n <- colMeans(sd.n)
sd.sd.n <- apply(sd.n, 2, sd)

mean.sd.n.2 <- colMeans(sd.n.2)
sd.sd.n.2 <- apply(sd.n.2, 2, sd)

# Plot the results
plot(1:max.sample.size, type = "n", xlab = "Sample Size (n)", ylab = "Estimator Value",
     main = "Performance of Standard Deviation Estimators", ylim = c(0, 8))
points(1:max.sample.size, mean.sd.ori, col = "red", pch = 1, lwd = 2, ylim = c(0, 1.5))
points(1:max.sample.size, mean.sd.n, col = "blue", pch = 2, lwd = 2)
points(1:max.sample.size, mean.sd.n.2, col = "green", pch = 3, lwd = 2)

legend("topright", legend = c("s (n-1)", "sqrt(n)", "sqrt(n-2)"),
       col = c("red", "blue", "green"), pch = 1:3, lwd = 2)

# Plot the standard deviations of the estimators
plot(1:max.sample.size, type = "n", xlab = "Sample Size (n)", ylab = "Standard Deviation",
     main = "Spread of Standard Deviation Estimators", ylim = c(0, 8))
points(1:max.sample.size, sd.sd.ori, col = "red", pch = 1, lwd = 2, ylim = c(0, 0.5))
points(1:max.sample.size, sd.sd.n, col = "blue", pch = 2, lwd = 2)
points(1:max.sample.size, sd.sd.n.2, col = "green", pch = 3, lwd = 2)

legend("topright", legend = c("s (n-1)", "sqrt(n)", "sqrt(n-2)"),
       col = c("red", "blue", "green"), pch = 1:3, lwd = 2)

########################################################################

# Set the random seed for reproducibility
set.seed(123)

# Define parameters for the normal distribution
mean1 <- 2
mean2 <- 2
variance1 <- 3
variance2 <- 3
n1 <- 5
n2 <- 7
N <- 104

# Create a matrix to store resampled data
resampling.matrix <- matrix(0, nrow = N, ncol = n1 + n2)

# Step 2: Define a function to calculate the mean difference
mean.diff <- function(x) mean(x[1:n1]) - mean(x[(n1 + 1):(n1 + n2)])

# Step 3: Simulate and calculate mean differences
for (i in 1:N) {
  sample1 <- rnorm(n1, mean = mean1, sd = sqrt(variance1))
  sample2 <- rnorm(n2, mean = mean2, sd = sqrt(variance2))
  resampling.matrix[i, ] <- c(sample1, sample2)
}

# Calculate the mean differences for all resamplings
sim.diff <- apply(resampling.matrix, 1, mean.diff)

# Calculate the brute force p-values using ranks
sim.p.rk <- rank(-abs(sim.diff)) / length(sim.diff)

# Step 4: Calculate the t-test p-values for the samples
t.test.p.values <- numeric(N)
for (i in 1:N) {
  t.test.p.values[i] <- t.test(resampling.matrix[i, 1:n1], resampling.matrix[i, (n1 + 1):(n1 + n2)])$p.value
}

# Plot the distributions and scatter plots
# Distribution of mean differences
hist(sim.diff, main = "Distribution of Mean Differences", xlab = "Mean Difference")

# Distribution of t-test p-values
hist(t.test.p.values, main = "Distribution of T-Test P-Values", xlab = "P-Value")

# Distribution of brute force p-values
hist(sim.p.rk, main = "Distribution of Brute Force P-Values", xlab = "P-Value")

# Scatter plot of t-test p-values vs. mean differences
plot(sim.diff, t.test.p.values, main = "T-Test P-Values vs. Mean Differences",
     xlab = "Mean Difference", ylab = "P-Value")

# Scatter plot of brute force p-values vs. mean differences
plot(sim.diff, sim.p.rk, main = "Brute Force P-Values vs. Mean Differences",
     xlab = "Mean Difference", ylab = "P-Value")

# Scatter plot of t-test p-values vs. brute force p-values
plot(t.test.p.values, sim.p.rk, main = "T-Test P-Values vs. Brute Force P-Values",
     xlab = "T-Test P-Value", ylab = "Brute Force P-Value")
