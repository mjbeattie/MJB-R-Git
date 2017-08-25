
# Checking correlation using the Spearman method
data("mtcars")
cor(mtcars$cyl, mtcars$mpg, method="spearman")


# Comparing histograms in hist() and qplot()
library(ggplot2)
data(iris)
hist(iris$Sepal.Length)
qplot(data=iris, Sepal.Length)

# Plotting using jitter to unhide masked variables
ggplot(data=iris, aes(x=Sepal.Length, y=Petal.Width)) +
    theme_bw() + 
    geom_point(aes(fill=Species),
               alpha=I(.75), position="jitter",
               colour="black",pch=21,size=5) +
    labs(y = "Petal Width (cm)", x = "Sepal Length (cm)") +
    theme(legend.key=element_blank(),
          axis.title = element_text(size = 14))


# Determining outliers on a graph
data(Animals)
ggplot(Animals, aes(x=brain, y=body, label=rownames(Animals))) + geom_text(size=3)

# Analysing missing values
install.packages("mice")   #install and load the MICE
library(mice)              #load the MICE package:  Multivariate Imputation by Chained Equations 
library(VIM)
data(mtcars)
mtcars <- mtcars[ -c(3, 5, 6, 8:11) ] # Drop all columns except mpg, cyl, hp, and qsec
mtcars[6:15,"mpg"] <- NA              # Create missing variables in mpg
mtcars[13:21,"hp"] <- NA              # Create missing variables in hp
mtcars[10:18,"qsec"] <- NA            # Create missing variables in qsec
head(mtcars)

fourcyl <- mtcars[mtcars$cyl==4,]
missing_fourcyl <- aggr(fourcyl)
summary(missing_fourcyl)


# Principal component analysis
data(mtcars)
apply(mtcars,2,sd)                    # Check standard deviations of variables
prcomp(mtcars,scale=TRUE)             # Run prcomp for PCA with z-score standardization
summary(prcomp(mtcars,scale=TRUE))    # Summary stats for PCA
biplot(prcomp(mtcars,scale=TRUE))     # Generate biplot of PCA vectors


