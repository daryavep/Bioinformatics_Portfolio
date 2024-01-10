######################################################################
####################        QUESTION 1        ########################

# Load the required data
library(ALL)

# Load expression data
data(ALL)

# Calculate remission time based on diagnosis and remission date
remission <- as.numeric(as.Date(ALL$date.cr, "%m/%d/%Y") - as.Date(ALL$diagnosis, "%m/%d/%Y"))

# Build data frame
d2r.34852 <- data.frame(
  expression = exprs(ALL)["34852_g_at", ],
  remission = remission
)

# Exclude rows with NA values in either variable
d2r.34852 <- na.omit(d2r.34852)

# Fit a linear model
model <- lm(remission ~ expression, data = d2r.34852)

# Perform ANOVA
anova.result <- anova(model)
print(anova.result)

# Get the observed values and predicted values
observed.values <- d2r.34852$remission
predicted.values <- predict(model)

# Calculate residuals
residuals.manual <- observed.values - predicted.values

# Calculate the sum of squares (model) manually
ss.model.manual <- sum((predicted.values - mean(predicted.values))^2, na.rm = TRUE)

# Calculate the sum of squares (residual) manually using residuals
ss.residual.manual <- sum(residuals.manual^2, na.rm = TRUE)

# Calculate the degrees of freedom manually
df.model.manual <- 1  # Model has 1 predictor
df.residual.manual <- length(remission) - 2

# Calculate mean sum of squares manually
mean.sq.model_manual <- ss.model.manual / df.model.manual
mean.sq.residual.manual <- ss.residual.manual / df.residual.manual

# Output the manual results
cat("Sum of Squares (Model):", ss.model.manual, "\n")
cat("Mean Sum of Squares (Model):", mean.sq.model.manual, "\n")
cat("Sum of Squares (Residual):", ss.residual.manual, "\n")
cat("Mean Sum of Squares (Residual):", mean.sq.residual.manual, "\n")

######################################################################
####################        QUESTION 2        ########################

# Build data frame
d2r.34852 <- data.frame(
  expression = exprs(ALL)["34852_g_at", ],
  remission = remission
)

# Exclude rows with NA
d2r.34852 <- na.omit(d2r.34852)

n.sim <- 100  # number of times we want to rerun the whole procedure (crossval or boot)
n.xval <- 5    # number of cross-validations (we want 5x crossval)
n.boots <- c(10, 100, 1000)   # number of bootstrap samples

mse.xval <- numeric(n.sim)  # Create a numeric vector to save MSEs for crossval
mse.boots <- matrix(NA, nrow = n.sim, ncol = length(n.boots))  # Create a matrix to save MSEs for different bootstrap sizes

# run a single loop here, n.sim times, in which we rerun both 5x crossval and boots
for (i.sim in 1:n.sim) { 
  # rerun the whole 5x-crossval and boot multiple times:
  
  s2.xval <- numeric(n.xval)  # Create an empty numeric vector
  xval.grps <- sample(1:n.xval, size = length(d2r.34852$remission), replace = TRUE)  # split data randomly into 5 groups for the current rerun of 5x crossval
  
  for (i.xval in 1:n.xval) { 
    # run crossval 5 times 
    test.df <- d2r.34852[xval.grps == i.xval, ]
    train.df <- d2r.34852[xval.grps != i.xval, ]
    
    lm.xval <- lm(remission ~ expression, data = train.df)  # fit the linear model on the training data only
    test.pred <- predict(lm.xval, newdata = test.df)  # Use the model to predict on test data
    s2.xval[i.xval] <- mean((test.pred - test.df$remission)^2)  # Perform squared error on the test and concatenate with s2.val
  }
  
  mse.xval[i.sim] <- mean(s2.xval)  # save mean error for the current rerun
  
  # Run bootstrap with different sample sizes
  for (i.boot_size in seq_along(n.boots)) {
    s2.boot <- numeric(n.boots[i.boot_size])  # Create an empty numeric vector
    n.obs <- dim(d2r.34852)[1]
    
    for (i.boot in 1:n.boots[i.boot_size]) {
      train.idx <- sample(1:n.obs, size = n.obs, replace = TRUE)  # sample indexes of the training set
      test.idx <- setdiff(1:n.obs, train.idx)  # get the indexes of the test data
      
      train.df <- d2r.34852[train.idx, ]  # extract the training data
      test.df <- d2r.34852[test.idx, ]    # Extract the test data
      
      lm.boot <- lm(remission ~ expression, data = train.df)  # fit the model on the training set
      test.pred <- predict(lm.boot, newdata = test.df)  # use the model to predict on test
      
      s2.boot[i.boot] <- mean((test.pred - test.df$remission)^2)  # calculate squared errors and concatenate with s2.boot
    }
    
    mse.boots[i.sim, i.boot_size] <- mean(s2.boot)  # save mean error for the current rerun
  }
}

# Draw boxplot to compare distributions of per-run MSE
boxplot(mse.xval, mse.boots[,1], mse.boots[,2], mse.boots[,3],
        names = c("Cross-validation", "Bootstrap n=10", "Bootstrap n=100", "Bootstrap n=1000"),
        main = "Comparison of Predictive Accuracy Assessment Methods",
        ylab = "Mean Squared Error")

######################################################################
####################        QUESTION 4        ########################

# Question 3 at the end of the code

# Build data frame
d2r.34852 <- data.frame(
  expression = exprs(ALL)["34852_g_at", ],
  remission = remission
)

# Exclude rows with NA
d2r.34852 <- na.omit(d2r.34852)

# Find the row with the smallest remission value
min.remission.row <- which.min(d2r.34852$remission)

# Remove the row with the smallest remission value
d2r.34852.filtered <- d2r.34852[-min.remission.row, ]

# Number of groups for cross-validation
n.xval <- 5
n.sim <- 100  # number of times we want to rerun the whole procedure

# Initialize vectors to store MSE values
mse.xval.wo.out <- numeric()

for (i.sim in 1:n.sim) { 
  xval.grps <- sample(1:n.xval, size = length(d2r.34852.filtered$remission), replace = TRUE)
  
  # Perform 5-fold cross-validation
  for (i.xval in 1:n.xval) {
    # Set group i aside as the 'test set'
    test.df <- d2r.34852.filtered[xval.grps == i.xval, ]
    train.df <- d2r.34852.filtered[xval.grps != i.xval, ]
  
    # Fit the model on the training set
    lm.xval <- lm(remission ~ expression, train.df)
  
    # Predict on the test set
    test.pred <- predict(lm.xval, test.df)
  
    # Calculate MSE with the record removed
    #mse.xval.wo.out <- c(mse.xval.wo.out, mean((test.pred - test.df$expression)^2, na.rm = TRUE))
    s2.xval[i.xval] <- mean((test.pred - test.df$remission)^2)
  }
  mse.xval.wo.out[i.sim] <-(s2.xval)
}

# Box plot with both MSE values
boxplot(list("With Record" = mse.xval, "Without Record" = mse.xval.wo.out), main = "MSE Comparison")

######################################################################
####################        QUESTION 3        ########################

library(limma);library(ALL); data(ALL)

ALL.pdat <- pData(ALL)

# Remission set
x.d2r <- as.numeric(as.Date(ALL$date.cr, "%m/%d/%Y") - as.Date(ALL$diagnosis, "%m/%d/%Y"))
d2r.tmp <- x.d2r[!is.na(x.d2r)]

# Expression set
#ALL.exprs <- exprs(ALL)[,!is.na(days2remiss)]
exprs.tmp <- exprs(ALL)[,!is.na(x.d2r)]

# Design matrix
design.matrix <- cbind(rep(1,length(d2r.tmp)),d2r.tmp)

n.sim <- 100 # number of times we want to rerun the whole procedure
n.xval <- 5 # number of cross validations
n.boot <-100
s2.xval <- numeric() # empty vector
s2.boot <- numeric () # empty vector
#xval.grps <- sample((1:dim(exprs.tmp)[2])%%n.xval+1) # Split groups

mse.xval <- numeric()  # Create a numeric vector to save MSEs
mse.boot.anew <- numeric()

# Vector to store the top genes selected by cross-validation and bootstrap
top.genes.xval <- character()
all.best.boot.genes <- character()


# Repeat process n.sim times
for (i.sim in 1:n.sim) {
  xval.grps <- sample((1:dim(exprs.tmp)[2])%%n.xval+1) # Split groups
  
  # Cross val loop
  for ( i.xval in 1:n.xval ) {
    exprs.tmp.train <- exprs.tmp[,xval.grps!=i.xval] # select training data
    design.matrix.train <- design.matrix[xval.grps!=i.xval,] #design matrix
    d2r.fit.train <- lmFit(exprs.tmp.train,design.matrix.train) # fit model to training data
    best.gene.xval <- rownames( topTable( eBayes(d2r.fit.train), # track best gene
                                          "d2r.tmp") )[1]
    top.genes.xval <- c(top.genes.xval, best.gene.xval)
    
    tmp.df <- data.frame(G=exprs.tmp[best.gene.xval,],D2R=d2r.tmp) # create data frame for best gene and D2R
    lm.xval <- lm(D2R~G,tmp.df[xval.grps!=i.xval,]) # fit model to D2R ~ best gene
    test.pred <- predict(lm.xval,tmp.df[xval.grps==i.xval,]) # test data
    s2.xval <- c(s2.xval,(test.pred- 
                            tmp.df[xval.grps==i.xval,"D2R"])^2) # MSE
  }
  
  mse.xval[i.sim] <- mean(s2.xval)
  
  # Bootstrap loop
  for ( i.boot in 1:100 ) {
    boot.idx <- sample(1:length(d2r.tmp),replace=T)
    exprs.boot <- exprs.tmp[,boot.idx]
    dm.boot <- design.matrix[boot.idx,]
    d2r.fit.boot <- lmFit(exprs.boot,dm.boot)
    
    best.gene <- rownames(topTable(eBayes(d2r.fit.boot),"d2r.tmp"))[1]
    all.best.boot.genes <- c(all.best.boot.genes,best.gene)
    df.tmp <- data.frame(G=exprs.tmp[best.gene,],D2R=d2r.tmp)
    boot.df <- df.tmp[boot.idx,]
    test.df <- df.tmp[!(1:dim(df.tmp)[1]%in%boot.idx),]
    best.lm.boot <- lm(D2R~G,boot.df)
    s2.boot <- c(s2.boot,(predict(best.lm.boot,test.df)-test.df$D2R)^2)
  }
  mse.boot.anew[i.sim] <- mean(s2.boot)
}  

# Create a table with counts of top genes selected by cross-validation
table.xval <- table(top.genes.xval)

# Create a table with counts of top genes selected by bootstrap
table.boot <- table(all.best.boot.genes)

# Print or manipulate the tables as needed
print("Counts of top genes selected by cross-validation:")
print(sort(table.xval, decreasing = TRUE)[1:10])

print("Counts of top genes selected by bootstrap:")
print(sort(table.boot, decreasing = TRUE)[1:10])

boxplot(list(mse.xval, mse.boot.anew), 
        names = c("Cross-Validation", "Bootstrap"),
        main = "Comparison of Predictive Accuracy Assessment Methods", 
        ylab = "Mean Squared Error")

######################################################################
####################        QUESTION 5        ########################


library(limma);library(ALL); data(ALL)

ALL.pdat <- pData(ALL)

# Remission set
x.d2r <- as.numeric(as.Date(ALL$date.cr, "%m/%d/%Y") - as.Date(ALL$diagnosis, "%m/%d/%Y"))
d2r.tmp <- x.d2r[!is.na(x.d2r)]

# Expression set
exprs.tmp <- exprs(ALL)[,!is.na(x.d2r)]

# Design matrix
design.matrix <- cbind(rep(1,length(d2r.tmp)),d2r.tmp)

n.sim <- 100 # number of times we want to rerun the whole procedure
n.xval <- 5 # number of cross validations
n.boot <-100
s2.xval <- numeric() # empty vector
s2.boot <- numeric () # empty vector

mse.xval <- numeric()  # Create a numeric vector to save MSEs
mse.boot.anew <- numeric()

# Vector to store the top genes selected by cross-validation and bootstrap
top.genes.xval <- character()
all.best.boot.genes <- character()


# Repeat process n.sim times
for (i.sim in 1:n.sim) {
  xval.grps <- sample((1:dim(exprs.tmp)[2])%%n.xval+1) # Split groups
  
  # Cross val loop
  for ( i.xval in 1:n.xval ) {
    exprs.tmp.train <- exprs.tmp[,xval.grps!=i.xval] # select training data
    design.matrix.train <- design.matrix[xval.grps!=i.xval,] #design matrix
    d2r.fit.train <- lmFit(exprs.tmp.train,design.matrix.train) # fit model to training data
    best.gene.xval <- rownames( topTable( eBayes(d2r.fit.train), # track best 2 genes
                                          "d2r.tmp") )[1:2]
    top.genes.xval <- c(top.genes.xval, best.gene.xval)
    
    # Separate the top two genes
    top_gene_1 <- best.gene.xval[1]
    top_gene_2 <- best.gene.xval[2]
    
    # Create a data frame
    tmp.df <- data.frame(D2R = d2r.tmp, G1 = exprs.tmp[top_gene_1, ], G2 = exprs.tmp[top_gene_2, ])
    lm.xval <- lm(D2R ~ G1 + G2, data = tmp.df[xval.grps != i.xval,])
    
    test.pred <- predict(lm.xval,tmp.df[xval.grps==i.xval,]) # test data
    s2.xval <- c(s2.xval,(test.pred- 
                            tmp.df[xval.grps==i.xval,"D2R"])^2) # MSE
  }
  
  mse.xval[i.sim] <- mean(s2.xval)
  
  # Bootstrap loop
  for ( i.boot in 1:100 ) {
    boot.idx <- sample(1:length(d2r.tmp),replace=T)
    exprs.boot <- exprs.tmp[,boot.idx]
    dm.boot <- design.matrix[boot.idx,]
    d2r.fit.boot <- lmFit(exprs.boot,dm.boot)
    
    best.gene <- rownames(topTable(eBayes(d2r.fit.boot),"d2r.tmp"))[1:2]
    all.best.boot.genes <- c(all.best.boot.genes,best.gene)
    
    # Separate the top two genes
    top_gene_1 <- all.best.boot.genes[1]
    top_gene_2 <- all.best.boot.genes[2]
    
    # Create a data frame 
    df.tmp <- data.frame(D2R = d2r.tmp, G1 = exprs.tmp[top_gene_1, ], G2 = exprs.tmp[top_gene_2, ])
    boot.df <- df.tmp[boot.idx,]
    test.df <- df.tmp[!(1:dim(df.tmp)[1]%in%boot.idx),]
    # Apply the lm function 
    best.lm.boot <- lm(D2R ~ G1 + G2, data = boot.df)
    
    s2.boot <- c(s2.boot,(predict(best.lm.boot,test.df)-test.df$D2R)^2)
  }
  mse.boot.anew[i.sim] <- mean(s2.boot)
}  

# Create a table with counts of top genes selected by cross-validation
table.xval <- table(top.genes.xval)

# Create a table with counts of top genes selected by bootstrap
table.boot <- table(all.best.boot.genes)

# Print or manipulate the tables as needed
print("Counts of top genes selected by cross-validation:")
print(sort(table.xval, decreasing = TRUE)[1:10])

print("Counts of top genes selected by bootstrap:")
print(sort(table.boot, decreasing = TRUE)[1:10])

boxplot(list(mse.xval, mse.boot.anew), 
        names = c("Cross-Validation", "Bootstrap"),
        main = "Comparison of Predictive Accuracy Assessment Methods", 
        ylab = "Mean Squared Error")
