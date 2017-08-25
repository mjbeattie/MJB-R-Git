#  Function Name:  binaryEvaluator
#  Description:  Runs goodness tests on a logistic regression model and generates fitness plots
#  Arguments:  observedData - observed data from a classification model
#              classProb - probability of positive classification from the model
#              predValue - user defined threshold for positive classification
#  Returns:  List of key statistics
#  Version:  1.0
#  Version Date:  November 20, 2016
#  Author:  Matthew J. Beattie
binaryEvaluator=function(observedData, classProb, predValue)
{

  #  Generate confusion matrix
  workingData <- data.frame(trueVal = observedData, predProb = classProb)
  workingData$pred <- as.numeric(workingData$predProb > predValue)
 
  cm <- confusionMatrix(workingData$trueVal, workingData$pred)
  
  #  Generate ROC Curve and AUC
  pred <- prediction(workingData$predProb, workingData$trueVal)    #ROC curve for training data
  perf <- performance(pred,"tpr","fpr") 
  plot(perf,colorize=TRUE, main = "ROC Curve", print.cutoffs.at = c(0.15,0.20,0.25,0.5,0.75)); 
  abline(0, 1, col="red") 
  AUC <- performance(pred, "auc")@y.values
  
  #lift chart (added by Matt Beattie)
  perf2 <- performance(pred, "rpp", "lift")
  plot(perf2, colorize = TRUE, main = "Lift Curve", print.cutoffs.at = c(0.25,0.5,0.75))
  
  #  Generate concordant pairs analysis (COMMENT OUT FOR LARGE DATASETS)
#  Con_Dis_Data = cbind(workingData$trueVal, workingData$predProb) 
#  
#  ones = Con_Dis_Data[Con_Dis_Data[,1] == 1,]
#  zeros = Con_Dis_Data[Con_Dis_Data[,1] == 0,]
#  conc=matrix(0, dim(zeros)[1], dim(ones)[1])   #build a matrix of 0's 
#  disc=matrix(0, dim(zeros)[1], dim(ones)[1])
#  ties=matrix(0, dim(zeros)[1], dim(ones)[1])
#  
#  for (j in 1:dim(zeros)[1])
#  {
#    for (i in 1:dim(ones)[1])
#    {
#      if (ones[i,2]>zeros[j,2])
#      {conc[j,i]=1}
#      
#      else if (ones[i,2]<zeros[j,2])
#      {disc[j,i]=1}
#      
#      else if (ones[i,2]==zeros[j,2])
#      {ties[j,i]=1}
#    }
#  }
#  
#  Pairs = dim(zeros)[1]*dim(ones)[1]              #total number of pairs
#  PercentConcordance = (sum(conc)/Pairs)*100
#  PercentDiscordance = (sum(disc)/Pairs)*100
#  PercentTied = (sum(ties)/Pairs)*100

  #  Compute D statistic
  dStat.1<-workingData[workingData$trueVal==1,]
  dStat.0<-workingData[workingData$trueVal==0,]
  dStat <- mean(dStat.1$predProb) - mean(dStat.0$predProb)
  
  
  #K-S chart  (Kolmogorov-Smirnov chart) 
  # measures the degree of separation 
  # between the positive (y=1) and negative (y=0) distributions
  workingData$group<-cut(workingData$predProb,seq(1,0,-.1),include.lowest=T)
  xtab <- table(workingData$group,workingData$trueVal)
  
  #make empty dataframe
  KS <- data.frame(Group=numeric(10),
                 CumPct0=numeric(10),
                 CumPct1=numeric(10),
                 Dif=numeric(10))
  
  #fill data frame with information: Group ID, 
  #Cumulative % of 0's, of 1's and Difference
  for (i in 1:10) {
    KS$Group[i]<-i
    KS$CumPct0[i] <- sum(xtab[1:i,1]) / sum(xtab[,1])
    KS$CumPct1[i] <- sum(xtab[1:i,2]) / sum(xtab[,2])
    KS$Dif[i]<-abs(KS$CumPct0[i]-KS$CumPct1[i])
  }
  
  KSStat <- KS[KS$Dif==max(KS$Dif),]
  
  maxGroup<-KS[KS$Dif==max(KS$Dif),][1,1]
  
  #K-S chart using plot
  plot(KS$Group, y = KS$CumPct0, type = "l", col="red", ylim=c(0,1),
       xlab="Deciles", ylab="Cumulative Percent", main="K-S Chart")
  lines(KS$CumPct1, col="blue", ylim=c(0,1))
  segments(maxGroup, KS$CumPct0[maxGroup], maxGroup, KS$CumPct1[maxGroup])
  
  # Plot the distribution charts for the true positives and true negatives
  plot(density(workingData[workingData$trueVal==1,2]), xlab = "Probability",
       main = "Probability Density Plot for Positive Classifications")
  plot(density(workingData[workingData$trueVal==0,2]), xlab = "Probability",
       main = "Probability Density Plot for Negative Classifications")

  #  Return a list of key statistics (REMOVE CONC/DISC AND PAIRS FOR LARGE DATASETS)
#  return(list("Percent Concordance"=PercentConcordance,
#              "Percent Discordance"=PercentDiscordance,
#              "Percent Tied"=PercentTied,
#              "Pairs"=Pairs,
#              "AUC"=AUC,
#              "D Statistic"=dStat,
#              "KS Statistic"=KSStat,
#              "Confusion Matrix"=cm))

    return(list("AUC"=AUC,
              "D Statistic"=dStat,
              "KS Statistic"=KSStat,
              "Confusion Matrix"=cm))
  
}
#***END OF BINARY MODEL EVALUATION FUNCTION***#
