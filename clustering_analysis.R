#Example of clustering for ISE 5103 Intelligent Data Analytics
#Charles Nicholson
#November 2015


library(ggplot2)    #for graphics and the dataset for this example session
library(cluster)    #provides more cluster algorithms than base R (e.g. PAM)
library(useful)     #provides a plot function for clusters and "FitKMeans" and "PlotHartigan" functions
library(NbClust)    #provides tons of tools for identifying the "right" number of clusters
library(rgl)        #for 3D rotating plots

#check out the input paramters for k-means and various options
?kmeans


#to start let's look at some really simple data: 2 dimensions 
data (cars,package="caret")
cars<-cars[,c("Price","Mileage")]   #price and mileage for several cars listed for resale in 2005 Kelly Blue Book
head(cars)
summary(cars)

#let's scale the data..
carScaled<-scale(cars)  #default is mean centering with scaling by standard deviation


#kmeans is a function in the standard R package "stats"
carsKM <- kmeans(carScaled,3, nstart=10)          

#note above: nstart=10   -- this performs the clustering 10 times with 10 different initial
#seeds; the clustering result with the minimum error is kept.


#we can take a look at the cluster results
carsKM$centers  # the centroids of the final clusers (remember, these are scaled)

carsKM$size #and the size of each cluster


# a ggplot of the original 2D data color-coded by the cluster number
p<-qplot(data=cars, x=Price, y=Mileage, color=factor(carsKM$cluster))  #plot the 2 variables and the cluster color
g <- guide_legend("Cluster")                  #retitle the legend...
p + guides(color = g, size = g, shape = g)    #retitle the legend...


#create our own plot function to look for "Within cluster Sum of square error 'elbow' plot"
#defaults to 15 as clusters max

wssplot <- function(data, nc=15){                    

  par(mfrow=c(1,2))
  
  wss <- NULL  
  pctExp <-NULL
  
  for (k in 1:nc)
  {
     kclus <- kmeans(data, centers=k)
     wss[k] <- kclus$tot.withinss      #store the total within SSE for given k
     pctExp[k] <- 1-wss[k]/kclus$totss
  }
  
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")

  plot(1:nc, pctExp, type="b", xlab="Number of Clusters",
       ylab="Pct Explained")
  
  par(mfrow=c(1,1))
}


#unfortunaely, there is no obivous "elbow"
wssplot(carScaled,nc=30)


# "Hartigan's rule is one statistical technique to "find" the elbow
clusFit<-FitKMeans(carScaled,max.clusters=30,nstart=20)   #evaluates k using the "Hartigan" rule
clusFit
PlotHartigan(clusFit)

#FYI: Hartigan rule is: ADD CLUSTER if (sum(k$withinss)/sum(kplus1$withinss)-1)*(nrow(x)-k-1) > 10

#unfortunately, with even with this technique, it seems to want to add too many clusters for this data




#demonstrate "NbClust"  
#since it takes a long time to run, demo on a small sample

# the NbClust package contains several techniques for evaluating the 'right' value for k
# the package identifies the value of k indicated by each statistical procedure
# some of the techniques agree, others don't

carScaledsamp<-carScaled[sample(nrow(carScaled),200),]

NbClust(carScaledsamp,method="kmeans")






############### irisdata ###################

#let's get a little bit more complex, and see how we do with "ground truth" available data

# we know that there really are 3 different species of in the Iris data
# let's see how many clusters are identified using these techniques

set.seed(42)
data(iris)
iris$Species <- factor(iris$Species,
                       levels = c("versicolor","virginica","setosa"))

iris[,1:4]<-scale(iris[,1:4])  # first scale the data
head(iris)

wssplot(iris[,1:4],nc=5)  #again, no obvious elbow...

#let's try NbClust...
NbClust(iris[,1:4],method="kmeans")

#we will use 3 clusters for k-means
cl <- kmeans(iris[,1:4],3)
iris$cluster <- as.factor(cl$cluster)
head(iris)


#to visualize this 4d data in 3d let's use PCA, and then color code base on clusters
pc <- princomp(iris[,1:4], cor=TRUE, scores=TRUE)
summary(pc)

plot3d(pc$scores[,1:3], col=iris$cluster, main="k-means clusters")

#and a quick check to see how our 'natural clusters' align with the species data
table(iris$Species, iris$cluster)



dev.off()

#let's try hiearchial clustering for the iris data
#hiearchial clustering requires a distance matrix

di <- dist(iris[,1:4], method="euclidean")   # with hiearchical clustering, only need distance matrix

hc <- hclust(di, method="ward")
plot(hc, labels=FALSE)

rect.hclust(hc, k=3, border="red")     #if we were to "cut" at k=3, what are the groups?
iris$hcluster <- as.factor(cutree(hc, k=3))   #cutting at k=3, here are the assignments

#and a quick check to see how our 'natural clusters' align with the species data
table( iris$Species, iris$hcluster)




############## mtcars data ##################


#Some easier data to play around with for hiearchial clustering
#use small data set for example...
#mtcars: 1974 Motor Trend US magazine data
# fuel consumption, 10 aspects of design and performance 

data(mtcars)   
head(mtcars)
MT<-mtcars[,c("mpg","cyl","disp","hp","wt","qsec")]
MT<-scale(MT)


#NOTE: the 'daisy" function in package "cluster" can handle mixed data for distances using the Gower's Distance method
#e.g.  d<-daisy(mtcars,stand=T)  

d<-daisy(mtcars,stand=T)

#different distance functions will produce different results



dev.off()


#look at the different shapes of dendrograms based on the linkage techniques

hclus<-hclust(d,method="single")   #notice the long chains (e.g., very unbalanced)
plot(hclus)

hclus<-hclust(d,method="complete")
plot(hclus)

hclus<-hclust(d,method="average")
plot(hclus)

hclus<-hclust(d,method="ward")  # notice how balanced the clusters are
plot(hclus)




#for fun and optional -- take a look at movie data from imdb.com

#################### movies data ##################
#we will use "movie rating data" for the example from the "ggplot2" library
library(ggplot2movies)
?movies
data(movies)
View(movies)


#but let's clean it up a bit first...
movies2<-movies[complete.cases(movies),c(1:6,17)]    #keep only complete cases; keep some variables
movies2<-movies2[movies2$votes>2500,]                #exclude movies with only a few ratings

mDD<-duplicated(movies2[,c("title")])
movies2<-movies2[!mDD,]

dat<-round(scale(movies2[,-c(1,7)]),3)               #scale the numeric data (and round the results)
row.names(dat)<-movies2$title                        #keep the movie titles as row.names
dat<-as.data.frame(dat)


#we will start by looking at partitions
wssplot(dat,nc=30)   #-- again, no clear elbow!  


set.seed(100)      #just so that we will all get the same results!

kclus<-kmeans(dat, 6, nstart=5)    #how about trying k=6

plot(kclus,data=dat) #this plot function comes from the "useful" libary uses PCA 

clusInfo<-data.frame(kclus$centers,kclus$size)
clusInfo

#by evaluating the centroids of each variable, we can discuss the 'meaning' of the cluster

#         year     length      budget     rating       votes kclus.size
#1 -2.22851913 -0.2312077 -0.91936612  0.9339563 -0.20458470        183
#2  0.20304230 -0.1427341 -0.45749396  0.3648474 -0.14853776        662
#3  0.01190123  0.9218765  0.27377778  1.3547160  3.47887654         81
#4  0.56625879  0.3389201  1.56101917 -0.1473259  0.17708626        313
#5  0.38614789 -0.5685563 -0.04270892 -1.2621502 -0.48717136        426
#6 -1.06658462  3.0500462 -0.33092308  0.9478769  0.09375385         65

# in this example, cluster #1 are the oldest movies, with smaller budgets, and above average ratings
# whereas, clsuter #4 are big budget movies that are newer but have lower ratings

#(your results might be different!)


movieClus <- data.frame(movies2, clust=kclus$cluster, dat)
head(movieClus[movieClus$clust==2,c("title","rating","budget","year","mpaa")])



#helper function to evaluate "representative data points within clusters"
#and possibly outliers
#-- calculate the distance between the given cluster centers
#-- and the data, and return the distance between the data and the closest cluster center

# NOTE to the progammatically inclined:
# this function is really poorly written in that it requires a "global" name for the cluster object (kclus)
# Sorry.

closest.cluster <- function(x) 
{
  cluster.dist <- apply(kclus$centers, 1, function(y) sqrt(sum((x-y)^2)))
  return(cluster.dist[which.min(cluster.dist)[1]])
}


movieClus$dist <- apply(dat, 1, closest.cluster)   #apply the "distance to nearest cluster function"
movieClus <-movieClus [order(movieClus$dist),]         #sort the data by this new distance

#now we can look at top 10 most representative titles in the cluster

head(movieClus[movieClus$clust==1,c("title","rating","budget","year","dist")],10) 
head(movieClus[movieClus$clust==2,c("title","rating","budget","year","dist")],10)  
head(movieClus[movieClus$clust==3,c("title","rating","budget","year","dist")],10) 
head(movieClus[movieClus$clust==4,c("title","rating","budget","year","dist")],10) 
head(movieClus[movieClus$clust==5,c("title","rating","budget","year","dist")],10) 
head(movieClus[movieClus$clust==6,c("title","rating","budget","year","dist")],10) 




