# Example code to demonstrate multivariate imputation by chained equations (mice)
# ISE 5103 Intelligent Data Analytics
# Charles Nicholson
# September 2015


# the package mice: multivariate imputation by chained equations
library(mice)


# create some random sample data 
#-------------------------------------------------
n=300   #n equals the number of observations


#four variables
x1<- 5*runif(n)-rexp(n)
x2<- rnorm(n) + runif(n) - 0.5*x1
x3<- sqrt(abs(x2)) + 2*rlnorm(meanlog=x1, sdlog=1, n) + rnorm(n)
x4<- 0.5*x3*runif(n) + rnorm(mean=x2,sd=1,n) + 2*runif(n) + rnorm(n)

# let y be some function of x1, x2, and x3
y<-5*x1*runif(n)-3*rnorm(mean=x2,sd=1,n)+2*sqrt(abs(x3))+rnorm(n)


# create a data frame from the vectors
df<-data.frame(y,x1,x2,x3,x4)

dfFull<-df  #save the full data for later use
#-------------------------------------------------


hist(df$y)
hist(df$x1)
hist(df$x2)
hist(df$x3)
hist(x3)
hist(df$x4)

# introduce some missingness in the data for multiple variables using different rules
#-------------------------------------------------
df[y<25 & runif(n)<0.65,"x1"]<-NA

hist(dfFull$x1)
hist(df$x1)

u<-runif(n)
df[u*y*y>80,"x2"]<-NA


hist(x3+x4)
df[x3+x4 > 250+runif(n)*200,"x3"]<-NA

u<-runif(n)
u1<-runif(n)
df[(y*u+x3+x2+u1) > 15,"x4"]<-NA

df[x3+x1<3,"y"]<-NA


u<-runif(n)
u1<-runif(n)
df[((x3+x1+x2)*u > 15 & (x3+x1+x4)*u1 < 50),"y"]<-NA

#check the percent missing per variable
myfun<-function(x) mean(is.na(x))
apply(df,2,myfun)
#-------------------------------------------------


# perform the first two steps of MI using the "mice" command 
# create m=6 data sets and impute missing values 
imp<-mice(df,m=7,maxit=50)

# the output object is quite complex!
str(imp)

#take a look at how the means and variances of the imputed values are (hopefully) converging 
imp$chainMean
imp$chainVar

#can plot those means and variances
plot(imp)

# perform the third step of MI using the "with" command
# to perform a standard analysis (in this case, a linear regression) on each data set 
fit<-with(imp, lm(y~x1+x2+x3+x4))

?with
#perfrom the fourth step of MI, recombination, using the "pool" command 
est<-pool(fit)


plot(dfFull)      #pairs plot of full data
plot(df)          #pairs plot of available cases
plot(na.omit(df)) #pairs plot for complete cases


#coefficient estimates based on full data (before creating missingness)
summary(fullfit<-lm(data=dfFull,y~x1+x2+x3+x4))

#coefficient estimates based on complete cases (no imputation)
summary(missfit<-lm(data=df,y~x1+x2+x3+x4))

ymissfit<-predict(missfit,newdata=dfFull)

#coefficient estimates based on MICE (recombined estimates)
summary(est)

est$qbar



yMIfit<-est$qbar[1]+dfFull[,2]*est$qbar[2]+dfFull[,3]*est$qbar[3]+dfFull[,4]*est$qbar[4]

library(ggplot2)

qplot(dfFull$y,ymissfit)+coord_cartesian(xlim=c(-50,150),ylim = c(-50, 150)) +geom_abline() 
qplot(dfFull$y,yMIfit)+coord_cartesian(xlim=c(-50,150),ylim = c(-50, 150)) +geom_abline() 
