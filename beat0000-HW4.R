#  Homework 4
#  October 30, 2016
#  Author:  Matthew J. Beattie
#  Filename:  beat0000-HW4

# LOAD IN COMMON LIBRARIES
library(VIM)
library(pls)
library(caret)
library(glmnet)
library(lars)
library(AppliedPredictiveModeling)

##### PROBLEM 6.1:  Meat Spectroscopy Data ############################

# 6.1.a:  Read in libraries and data
library(caret)
data(tecator)

# 6.1.b:  Perform PCA to reduce the number of factors in the model
# Create a data frame on the combination of absorp and endpoints
df.fullData <- data.frame(endpoints, absorp)
names(df.fullData)[1:3] <- c("water", "fat", "protein")


# Explore data first
# Check for missingness
aggr.fullData <- aggr(df.fullData)
summary(aggr.fullData)
cat("There are no missing data in the Absorp dataset")

# Examine absorp data - depict correlations graphically
cormatAbsorp <- cor(absorp)
heatmap(cormatAbsorp, main = "Heatmap of Correlations in Absorp")
cat("It appears that 20 or so variables that tightly correlated in Absorp")

# Conduct PCA using prcomp
pcaAbsorp <- prcomp(absorp, scale = TRUE)
summary(pcaAbsorp)                      
cat("One component explains 98.6% of the variation of the Absorp dataset")
head(pcaAbsorp$rotation[,1:3])
biplot(pcaAbsorp)
# Note:  The dimension of the transformed data is 1

# 6.1.c:  Build models on transformed data
# Combine transformed absorp data with endpoints data into a new dataframe
PC1 <- pcaAbsorp$x[,1]
pcaData <- data.frame(endpoints, PC1)
names(pcaData) <- c("water", "fat", "protein", "PC")

# Create test and training sets
# 75% of the sample size
smp_size <- floor(0.75 * nrow(pcaData))
set.seed(123)
train_ind <- sample(seq_len(nrow(pcaData)), size = smp_size)
lm.data <- df.fullData[, -c(1,3)]
lm.train <- lm.data[train_ind, ]
lm.test <- lm.data[-train_ind, ]
fullData.test <- df.fullData[-train_ind,]
pca.train <- pcaData[train_ind, ]
pca.test <- pcaData[-train_ind, ]

# Perform simple linear regression using only PCA component 1
lm1 <- lm(data=pca.train, fat ~ PC)
summary(lm1)                                
AIC(lm1)
BIC(lm1)
plot(lm1$fitted.values, lm1$residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residual Plot for the Absorp PCA Model")

# Examine PCA model predictions vs observed data in test dataset
pred.lm1 <- predict(lm1, pca.test)
cat("Test MSE for Absorp PCA model is: ", mean((pred.lm1 - fullData.test[,2])^2))
plot(pred.lm1, fullData.test[,2], xlab = "PCA Predicted Values", ylab = " Observed Values",
     main = "PCA Performance on Test Absorp Data - One Component", col = "red")
abline(0,1)


# Perform linear regression on full dataset
lm2 <- lm(data=lm.train, fat ~., y=TRUE)
summary(lm2)
AIC(lm2)
BIC(lm2)
plot(lm2$fitted.values, lm2$residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residual Plot for the Absorp Full OLM Model")

# Examine full model predictions vs observed data in test dataset
pred.lm2 <- predict(lm2, lm.test)
mean((pred.lm2 - fullData.test[,2])^2)
plot(pred.lm2, fullData.test[,2], xlab = "Full Predicted Values", ylab = "Observed Values",
     main = "Absorp Full Model Performance on Test Data", col = "red")
abline(0,1)

# Compare manual PCR against pcr()
set.seed(1)
lm1a=pcr(fat ~., data=lm.train, scale =TRUE, validation ="CV")
summary(lm1a)
validationplot(lm1a, val.type="MSEP", main = "Cross Validation Plot for Number of Components in Absorp PCR")

# Examine PCR predictions vs observed data in test dataset
pred.lm1a <- predict(lm1a, lm.test, ncomp=5)
cat("Test MSE for Absorp PCR model is" , mean((pred.lm1a - fullData.test[,2])^2))
plot(pred.lm1a, fullData.test[,2], xlab = "PCR Predicted Values", ylab = "Observed Values",
     main = "Absorp PCR Performance on Test Data - Five Components", col = "red")
abline(0,1)

# Perform linear regression using via PLS
set.seed(1)
lm3 <- plsr(data=lm.train, fat ~., ncomp = 10, validation="CV")
summary(lm3)
plot(lm3$fitted.values, lm3$residuals)
plot(lm3, ncomp=1:8, asp=1, line=TRUE)
validationplot(lm3, val.type="MSEP")

# Examine PLS predictions vs observed data in test dataset
pred.lm3 <- predict(lm3, lm.test, ncomp=5)
cat("Test MSE for Absorp PLS model is: ", mean((pred.lm3 - fullData.test[,2])^2))
plot(pred.lm3, fullData.test[,2], xlab = "PLS Predicted Values", ylab = "Observed Values",
     main = "PLS Performance on Test Data", col = "red")
abline(0,1)



##### PROBLEM 6.2:  PERMEABILITY  #######################
# 6.2.a:  Read in data
library(AppliedPredictiveModeling)
data(permeability)

# 6.2.b:  Remove near zero variance predictors
x <- nearZeroVar(fingerprints, freqCut=90/10, uniqueCut=10)
nz.fingerprints <- fingerprints[,-x]
cat("There are ", ncol(nz.fingerprints), " remaining predictors in the Permeability dataset")

# 6.2.c:  Split data into training and test, pre-process, and tune a PLS model
# Convert data into a data frame and create training and test sets at a 75/25 split
df.nz.fingerprints <- data.frame(nz.fingerprints)
df.permeability <- data.frame(permeability)
perm.data <- data.frame(df.permeability, df.nz.fingerprints)
smp_size <- floor(0.75 * nrow(nz.fingerprints))
set.seed(123)
train_ind <- sample(seq_len(nrow(nz.fingerprints)), size = smp_size)
perm.train <- perm.data[train_ind, ]
perm.test <- perm.data[-train_ind, ]

# Conduct a PLS model on the training dataset and calculate R-squared
set.seed(1)
perm.pls <- plsr(data=perm.train, permeability ~., ncomp = 50, validation="CV")
summary(perm.pls)
cat("Optimal number of components for the Permeability PLS model is 6 based upon RMSEP")

plot(perm.pls$fitted.values, perm.pls$residuals, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residual Plot for Permeability PLS Model")
plot(perm.pls, ncomp=1:8, asp=1, line=TRUE)
validationplot(perm.pls, val.type="MSEP", main = "Cross Validation Plot for Permeabilty PLS Model")

perm.pred.train <- predict(perm.pls, perm.train, ncomp=6)
mse.train <- mean((perm.pred.train - perm.train$permeability)^2)      # mean squared error
r2.train <- 1 - sum((perm.train$permeability - perm.pred.train)^2) / 
  sum((perm.train$permeability - mean(perm.train$permeability))^2)
cat("R squared for the Permeability PLS training model is: ", r2.train)


# 6.2.d:  Predict response from test set and compute R-squared
# Examine PLS predictions vs observed data in test dataset
perm.pred <- predict(perm.pls, perm.test, ncomp=6)
mse.test <- mean((perm.pred - perm.test$permeability)^2)      # mean squared error
r2.test <- 1 - sum((perm.test$permeability - perm.pred)^2) / 
             sum((perm.test$permeability - mean(perm.test$permeability))^2)
cat("R squared for the Permeability PLS test model is: ", r2.test)

plot(perm.pred, perm.test[,1], xlab = "PLS Predicted Values", ylab = "Observed Values",
     main = "PLS Performance on Permeability Test Data", col = "red")
abline(0,1)

# 6.2.e:  Try an alternative model to PLS - LASSO
x <- as.matrix(perm.train[,-1])
y <- perm.train[,1]
perm.lasso <- glmnet(x, y)
plot(perm.lasso, main = "Permeability LASSO Branching")                      # investigate the LASSO model with branching of coefficients
cv.perm.lasso <- cv.glmnet(x, y)      # conduct cross validation to determine optimal lambda
plot(cv.perm.lasso, main = "Permeability LASSO Cross Validation")           # plot cross validation
cat("The optimal lambda parameter for the Permeability LASSO model is: ", cv.perm.lasso$lambda.min) 

# count the number of non-zero coefficients in the model
coef.lasso <- coef(cv.perm.lasso, s="lambda.min")
index <- coef.lasso[,1] == 0
nonzero.coef.lasso <- data.frame(coef.lasso[!index,])
cat("The number of non-zero coefficients in the Permeability LASSO model is: ", nrow(nonzero.coef.lasso))

#  Use LASSO model to predict
perm.pred.lasso <- predict(cv.perm.lasso, as.matrix(perm.test[,-1]), s="lambda.min")
plot(perm.pred.lasso, perm.test[,1], xlab = "LASSO Predicted Values", ylab = "Observed Values",
     main = "LASSO Performance on Permeability Test Data", col = "red")
abline(0,1)

#  Manually calculate R2 from LASSO prediction
r2.test.lasso <- 1 - sum((perm.test$permeability - perm.pred.lasso)^2) / 
  sum((perm.test$permeability - mean(perm.test$permeability))^2)
cat("R squared for the Permeability LASSO test model is: ", r2.test.lasso)



##### PROBLEM 6.3:  Manufacturing Improvement  ####################
# 6.3.a:  Read in data
data(ChemicalManufacturingProcess)
rawdf <- ChemicalManufacturingProcess

# 6.3.b:  Examine missingness and check for skewness
summary(aggr(rawdf))                              #Note:  MF Process 3 had the most missings at 15
rawdf.kNN <- kNN(rawdf, k=5)
plot(rawdf.kNN$ManufacturingProcess03, rawdf.kNN$Yield, col = factor(rawdf.kNN$ManufacturingProcess03_imp), 
     main="k-Nearest Neighbors Imputation Sample", xlab="MF Process 3", ylab="Yield")
qplot(rawdf.kNN$ManufacturingProcess03, geom="density", xlab = "MF Process 3", main = "Density Plot of MF Process 3")

# Create imputed dataset by removing logical columns
rawdf.imp <- subset(rawdf.kNN, select = -c(Yield_imp:ManufacturingProcess45_imp))

# 6.3.c:  Split data into training and test and build a model
#  Examine correlations to determine any collinearity
cormatAbsorp <- cor(rawdf.imp)
heatmap(cormatAbsorp, main = "Heatmap of Correlations in Manufacturing Data")

#  Create training and test data
smp_size <- floor(0.75 * nrow(rawdf.imp))
set.seed(123)
train_ind <- sample(seq_len(nrow(rawdf.imp)), size = smp_size)
rawdf.train <- rawdf.imp[train_ind, ]
rawdf.test <- rawdf.imp[-train_ind, ]

# Use LASSO to zero out coefficients and preserve ability to tweak specific factors
rawdf.pred <- as.matrix(rawdf.train[,-1])
rawdf.resp <- rawdf.train[,1]
cv.rawdf.lasso <- cv.glmnet(rawdf.pred, rawdf.resp)      # conduct cross validation to determine optimal lambda
plot(cv.rawdf.lasso, main = "Cross Validation Plot for Manufacturing LASSO Model") 
cat("The optimal lambda parameter for the Manufacturing LASSO model is: ", cv.rawdf.lasso$lambda.min) 

#  Use LASSO model to predict
rawdf.pred.lasso <- predict(cv.rawdf.lasso, as.matrix(rawdf.test[,-1]), s="lambda.min")
plot(rawdf.pred.lasso, rawdf.test[,1], xlab = "LASSO Predicted Values", ylab = "Observed Values",
     main = "LASSO Performance on Manufacturing Test Data", col = "red")
abline(0,1)

#  6.3.d:  Calculate R2 from the test set and compare it to the training R2
#  Manually calculate R2 from LASSO prediction
r2.test.lasso <- 1 - sum((rawdf.test$Yield - rawdf.pred.lasso)^2) / 
  sum((rawdf.test$Yield - mean(rawdf.test$Yield))^2)                      
cat("R squared for the Manufacturing LASSO test model is: ", r2.test.lasso)


#  Manually calculate R2 from LASSO model training
rawdf.train.lasso <- predict(cv.rawdf.lasso, rawdf.pred, s="lambda.min")
r2.train.lasso <- 1 - sum((rawdf.train$Yield - rawdf.train.lasso)^2) / 
  sum((rawdf.train$Yield - mean(rawdf.train$Yield))^2)                  
cat("R squared for the Manufacturing LASSO training model is: ", r2.train.lasso)

#  6.3.e:  Identify which predictors are most important
rawdf.coefs.lasso <- predict(cv.rawdf.lasso, as.matrix(rawdf.test[,-1]), 
                             s="lambda.min", type="coefficients")
rawdf.coefs.lasso





