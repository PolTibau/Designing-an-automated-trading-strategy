## Set working directory to current script's location.
fileloc <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(fileloc)
rm(fileloc)

library(xts)
library(zoo)
library(quantmod) # data, plotting, quant modelling
library(PerformanceAnalytics) # performance and risk management
library(matrixStats)
library(doParallel)

# Load and set up data
##frequency of sampling
tau=1 #data is daily. Try tau=20 (month), tau=60 (quarterly)

#Used Google and Deutsche Bank markets to test the strategies
#x <- read.csv("C:/Users/polti/Downloads/GOOG.csv", sep = ",", header = TRUE)
x <- read.csv("C:/Users/polti/Downloads/DB.csv", sep = ",", header = TRUE)
#If is necessary to delete the first column of .csv data
x<-x[,-1]

#x <- read.csv("C:/Users/polti/Downloads/NNetDB.csv", sep = ",", header = TRUE)
#x <- read.csv("C:/Users/polti/Downloads/NNetGOOG.csv", sep = ",", header = TRUE)
#x <- read.csv("C:/Users/polti/Downloads/ArimaGOOG.csv", sep = ",", header = TRUE)
#x <- read.csv("C:/Users/polti/Downloads/ArimaDB.csv", sep = ",", header = TRUE)
#x <- read.csv("C:/Users/polti/Downloads/GP_GOOG.csv", sep = ",", header = TRUE)
#x <- read.csv("C:/Users/polti/Downloads/GP_DB.csv", sep = ",", header = TRUE)

##use Date as index
data <- as.xts(zoo(as.matrix(x[,-1]), as.Date(as.character(x[,1]))))

## Target :  Adj.Close Price 
target <- data$Adj.Close
#names(target) <- "DB"
#names(target) <- "GOOG"


##Performance Measures for Trading strategies
Performance <- function(x,ntrades=1,cost=0) {
  cumRetx = Return.cumulative(x,geometric = TRUE) -ntrades*cost
  #annRetx = Return.annualized(x, scale=252,geometric = T) -ntrades*cost
  annRetx = (cumRetx +1)^{252/length(x)} -1
  #sharpex = SharpeRatio.annualized(x, scale=252)
  sharpex = annRetx/sd.annualized(x,scale=252)
  winpctx = length(x[x > 0])/length(x[x != 0])
  annSDx = sd.annualized(x, scale=252)
  DDs <- findDrawdowns(x)
  maxDDx = min(DDs$return)
  maxLx = max(DDs$length)
  
  Perf = c(cumRetx, annRetx, sharpex, winpctx, annSDx, maxDDx, maxLx, ntrades)
  names(Perf) = c("Cumulative Return", "Annual Return","Annualized Sharpe Ratio","Win %", 
                  "Annualized Volatility", "Maximum Drawdown", "Max Length Drawdown","n.trades")
  return(Perf)
}


#Computation of pivot points:
PP = (data$High + data$Low + data$Close) / 3
S1 = (PP * 2)-data$High
S2 = PP-(data$High-data$Low)
S3 = data$Low-2*(data$High-PP)
R1 = (PP * 2)-data$Low
R2 = PP + (data$High-data$Low)
R3 = data$High + 2*(PP-data$Low)
names(PP) = "PP"; names(S1) = "S1"; names(S2) = "S2"
names(S3) = "S3"; names(R1) = "R1"; names(R2) = "R2"
names(R3) = "R3"
PivotPoints = cbind(PP,S1,S2,S3,R1,R2,R3)

#Computation of matrix T to see market tendencies
T <- vector(length = length(PivotPoints[,1]))
T <- lag(PivotPoints$PP,-1) - data$Open
T[length(PivotPoints[,1])] = 0
names(T) = "T"


##Average Performance of MA crossover over short (1 yr, 6 mon) fixed length periods
##using Rolling windows. Window size 252 (a year of daily data) or 252/2 for 6 mon
## step-size fixed to 1 day (ToDo: make stepsize variable; include transaction costs)
RollingTestMAStrategy1 <- function(myStock, ts, lim, longshort, Pl, Nl,w_size=252) {
  myPosition <- sig <- Lag(Sent_and_PP_1(ts,lim,longshort, Pl, Nl),1)
  bmkReturns <- dailyReturn(myStock, type = "arithmetic")
  myReturns <- bmkReturns*sig
  names(bmkReturns) <- 'BH'
  names(myReturns) <- 'PPandSent_1'  ## paste(names(ts),".on.",names(myStock),sep="")
  tt <- na.omit(merge(bmkReturns,myReturns))
  n_windows = nrow(tt) - w_size
  if(n_windows<1) stop("Window size too large")
  
  perform = foreach(i=1:n_windows, .combine = rbind) %do%{
    bhStra = tt$BH[i:(w_size+i-1),]
    myStra = tt$PPandSent_1[i:(w_size+i-1),]
    per=rbind(BH=Performance(bhStra),PPandSent_1=Performance(myStra))
    return(per)
  }
  
  bhindx = seq(1,2*n_windows,2); meindx = seq(2,2*n_windows,2)
  BHmeans = colMeans2(perform,rows = bhindx)
  MEmeans = colMeans2(perform,rows = meindx)
  MeanPerf=rbind(BHmeans,MEmeans)
  colnames(MeanPerf)=colnames(perform)
  rownames(MeanPerf)=c("BH","PPandSent_1")
  return(list("AvgPerf"=MeanPerf,"NumWindows"=n_windows))
}

RollingTestMAStrategy2 <- function(myStock, ts, lim, longshort, Pl, Nl,w_size=252) {
  myPosition <- sig <- Lag(Sent_and_PP_2(ts,lim,longshort, Pl, Nl),1)
  bmkReturns <- dailyReturn(myStock, type = "arithmetic")
  myReturns <- bmkReturns*sig
  names(bmkReturns) <- 'BH'
  names(myReturns) <- 'PPandSent_2'  ## paste(names(ts),".on.",names(myStock),sep="")
  tt <- na.omit(merge(bmkReturns,myReturns))
  n_windows = nrow(tt) - w_size
  if(n_windows<1) stop("Window size too large")
  
  perform = foreach(i=1:n_windows, .combine = rbind) %do%{
    bhStra = tt$BH[i:(w_size+i-1),]
    myStra = tt$PPandSent_2[i:(w_size+i-1),]
    per=rbind(BH=Performance(bhStra),PPandSent_2=Performance(myStra))
    return(per)
  }
  
  bhindx = seq(1,2*n_windows,2); meindx = seq(2,2*n_windows,2)
  BHmeans = colMeans2(perform,rows = bhindx)
  MEmeans = colMeans2(perform,rows = meindx)
  MeanPerf=rbind(BHmeans,MEmeans)
  colnames(MeanPerf)=colnames(perform)
  rownames(MeanPerf)=c("BH","PPandSent_2")
  return(list("AvgPerf"=MeanPerf,"NumWindows"=n_windows))
}
  
RollingTestMAStrategy3 <- function(myStock, ts, lim, longshort, Pl, Nl,w_size=252) {
  myPosition <- sig <- Lag(Sent_and_PP_3(ts,lim,longshort, Pl, Nl),1)
  bmkReturns <- dailyReturn(myStock, type = "arithmetic")
  myReturns <- bmkReturns*sig
  names(bmkReturns) <- 'BH'
  names(myReturns) <- 'PPandSent_3'  ## paste(names(ts),".on.",names(myStock),sep="")
  tt <- na.omit(merge(bmkReturns,myReturns))
  n_windows = nrow(tt) - w_size
  if(n_windows<1) stop("Window size too large")
  
  perform = foreach(i=1:n_windows, .combine = rbind) %do%{
    bhStra = tt$BH[i:(w_size+i-1),]
    myStra = tt$PPandSent_3[i:(w_size+i-1),]
    per=rbind(BH=Performance(bhStra),PPandSent_3=Performance(myStra))
    return(per)
  }
  
  bhindx = seq(1,2*n_windows,2); meindx = seq(2,2*n_windows,2)
  BHmeans = colMeans2(perform,rows = bhindx)
  MEmeans = colMeans2(perform,rows = meindx)
  MeanPerf=rbind(BHmeans,MEmeans)
  colnames(MeanPerf)=colnames(perform)
  rownames(MeanPerf)=c("BH","PPandSent_3")
  return(list("AvgPerf"=MeanPerf,"NumWindows"=n_windows))
}

RollingTestMAStrategy4 <- function(myStock, ts, lim, longshort, Pl, Nl,w_size=252) {
  myPosition <- sig <- Lag(Sent_and_PP_4(ts,lim,longshort, Pl, Nl),1)
  bmkReturns <- dailyReturn(myStock, type = "arithmetic")
  myReturns <- bmkReturns*sig
  names(bmkReturns) <- 'BH'
  names(myReturns) <- 'PPandSent_4'  ## paste(names(ts),".on.",names(myStock),sep="")
  tt <- na.omit(merge(bmkReturns,myReturns))
  n_windows = nrow(tt) - w_size
  if(n_windows<1) stop("Window size too large")
  
  perform = foreach(i=1:n_windows, .combine = rbind) %do%{
    bhStra = tt$BH[i:(w_size+i-1),]
    myStra = tt$PPandSent_4[i:(w_size+i-1),]
    per=rbind(BH=Performance(bhStra),PPandSent_4=Performance(myStra))
    return(per)
  }
  
  bhindx = seq(1,2*n_windows,2); meindx = seq(2,2*n_windows,2)
  BHmeans = colMeans2(perform,rows = bhindx)
  MEmeans = colMeans2(perform,rows = meindx)
  MeanPerf=rbind(BHmeans,MEmeans)
  colnames(MeanPerf)=colnames(perform)
  rownames(MeanPerf)=c("BH","PPandSent_4")
  return(list("AvgPerf"=MeanPerf,"NumWindows"=n_windows))
}  

##Indicators for use in trade-signals:
# Bearish sentiment indicators
negativeP= na.fill(data$negativePartscr,0)
names(negativeP)<- "negative"
uncertaintyP= na.fill(data$uncertaintyPartscr,0)
findownP= na.fill(data$findownPartscr,0)

## Bullish sentiment indicators
positiveP= na.fill(data$positivePartscr, 0)
names(positiveP)<- "positive"
certaintyP= na.fill(data$certaintyPartscr,0)
finupP= na.fill(data$finupPartscr,0)

##Combinations
BULL = .33*(positiveP+certaintyP+finupP); 
BEAR= .33*(negativeP+uncertaintyP+findownP); 
names(BULL)<-"BULL"; names(BEAR)<-"BEAR"
# Bull-Bear ratio: NAs are interpolated with  leftmost and rightmost non-NA value
##(NAs can be produced when BULL=BEAR=0)
#We define NAs as 50 because with our algorithm it will be the same as NULL information
BBr<- na.fill(100*BULL/(BULL+BEAR),50)
names(BBr)<-"BBr"
# BBlog<- 0.5*log((BULL+1)/(BEAR+1))

#PNr <-na.fill(100*positiveP/(positiveP+negativeP),"extend")
PNlog<- 0.5*log((positiveP+1)/(negativeP+1))
names(PNlog)<-"PNr"
RVT <- na.approx(data$RVT)
names(RVT)<-"RVT"

S <- cbind(BBr, RVT);
#S <- na.omit(S)
ts <- cbind(S,T)


testStrategy <-function(myStock, ts, lim, longshort, Pl, Nl, tcost){
  myPosition <- sig <- Lag(Sent_and_PP_4(ts,lim,longshort, Pl, Nl),1)
  ##compute runs of different signals (1,0 or -1). The number of these runs give number of trades
  runs<-rle(as.vector(na.omit(sig)))
  lruns <- length(runs$lengths)
  bmkReturns <- dailyReturn(myStock, type = "arithmetic")
  myReturns <- bmkReturns*sig
  names(bmkReturns) <- 'BH'
  names(myReturns) <- 'PP_And_Sent_3'
  tt <- na.omit(merge(bmkReturns,myReturns))
  ##Performance
  charts.PerformanceSummary(cbind(tt$PP_And_Sent_3, tt$BH))
  cbind(PPandSent=Performance(tt$PP_And_Sent_3,lruns,tcost),BH=Performance(tt$BH,2,tcost))
}

#Comparison of all methods planned
testStrategy_comp <-function(myStock, ts, lim, longshort, Pl, Nl, tcost){
  myPosition1 <- sig1 <- Lag(Sent_and_PP_1(ts,lim,longshort, Pl, Nl),1)
  myPosition2 <- sig2 <- Lag(Sent_and_PP_2(ts,lim,longshort, Pl, Nl),1)
  myPosition3 <- sig3 <- Lag(Sent_and_PP_3(ts,lim,longshort, Pl, Nl),1)
  myPosition4 <- sig4 <- Lag(Sent_and_PP_4(ts,lim,longshort, Pl, Nl),1)
  ##compute runs of different signals (1,0 or -1). The number of these runs give number of trades
  runs1<-rle(as.vector(na.omit(sig1)))
  runs2<-rle(as.vector(na.omit(sig2)))
  runs3<-rle(as.vector(na.omit(sig3)))
  runs4<-rle(as.vector(na.omit(sig4)))
  
  lruns1 <- length(runs1$lengths)
  lruns2 <- length(runs2$lengths)
  lruns3 <- length(runs3$lengths)
  lruns4 <- length(runs4$lengths)
  
  bmkReturns <- dailyReturn(myStock, type = "arithmetic")
  myReturns1 <- bmkReturns*sig1
  myReturns2 <- bmkReturns*sig2
  myReturns3 <- bmkReturns*sig3
  myReturns4 <- bmkReturns*sig4
  
  names(bmkReturns) <- 'BH'
  names(myReturns1) <- 'PPandSent_1'
  names(myReturns2) <- 'PPandSent_2'
  names(myReturns3) <- 'PPandSent_3'
  names(myReturns4) <- 'PPandSent_4'
  
  tt <- na.omit(merge(bmkReturns,myReturns1,myReturns2,myReturns3,myReturns4))
  ##Performance
  charts.PerformanceSummary(cbind(tt$PPandSent_1, tt$PPandSent_2, tt$PPandSent_3, tt$BH)) #tt$PPandSent_3
  cbind(PPandSent_1=Performance(tt$PPandSent_1,lruns1,tcost),PPandSent_2=Performance(tt$PPandSent_2,lruns2,tcost),PPandSent_3=Performance(tt$PPandSent_3,lruns3,tcost), PPandSent_4=Performance(tt$PPandSent_4,lruns4,tcost), BH=Performance(tt$BH,2,tcost))
}


#Using only Pivot points
Sent_and_PP_4 <- function(ts, lim, longshort, Pl, Nl){
  D <- ts$T
  c = 0
  p = 5
  for(i in 1:(length(D))){
    if(D[i] > 0.2 & c != 1){
      D[i]=p
      c=p
    }
    else if(D[i]< -0.2 & longshort==-1 & c!=-1){
      D[i]=-p
      c=-p
    }
    else if(D[i]< -0.2 & longshort==0 & c!=0){
      D[i]=0
      c=0
    }
    else{
      D[i]=p
    }
  }
  print(D)
  return(D)
}

#Similar to 2 but with less restrictions to buy
Sent_and_PP_3 <- function(ts, lim, longshort, Pl, Nl){
  D <- rep(0,length(ts[,1]))
  print(dim(ts))
  print(length(D))
  compte=0
  for(i in 1:(length(D))){
    if(ts$T[i]>0 & ts$BBr[i] > Pl & compte != 1){
      D[i] = 1
      compte=1
    }
    else if(ts$BBr[i] < Nl & ts$T[i]<0 & longshort==-1 & compte != -1){
      D[i] = -1
      compte=-1
    }
    else if(ts$BBr[i] < Nl & ts$T[i]<0 & longshort==0 & compte != 0){
      D[i] = 0
      compte=0
    }
    else D[i] = compte
  }
  print(D)
  return(D)
}

#Combination of PP and Sentimental
Sent_and_PP_2 <- function(ts, lim, longshort, Pl, Nl){
  D <- rep(0,length(ts[,1]))
  compte = 0
  for(i in 1:(length(D))){
    if(i > lim & ts$RVT[i] > EMA(ts$RVT, lim)[i]){
      if(ts$BBr[i] > Pl & ts$T[i] > 0 & compte != 1){
        D[i] = 1
        compte = 1
      }
      else if(ts$BBr[i] < Nl & ts$T[i] < 0 & longshort==-1 & compte != -1){
        D[i] = -1
        compte = -1
      }
      else if(ts$BBr[i] < Nl & ts$T[i] < 0 & longshort==0 & compte != 0){
        D[i] = 0
        compte = 0
      }
      else D[i] = compte
    }
    else D[i] = compte
  }
  print(D)
  return(D)
}

#Using only Sentimental
Sent_and_PP_1 <- function(ts, lim, longshort, Pl, Nl){
  D <- rep(0,length(ts[,1]))
  compte<-0
  for(i in 1:(length(D))){
    if(i > lim & ts$RVT[i] > 1.5*EMA(ts$RVT, lim)[i]){
      if(ts$BBr[i] > Pl & compte!=1){
        D[i] = 1
        compte=1
      }
      else if(longshort==-1 & ts$BBr[i] < Nl & compte!=-1){
        D[i] = -1
        compte=-1
      }
      else if(longshort==0 & ts$BBr[i] < Nl & compte !=0){
        D[i] = 0
        compte = 0
      }
      else D[i] = compte
    }
    else D[i] = compte
  }
  print(D)
  return(D)
}

##Full period test. Tune s, m, longshort=0,-1
res<-testStrategy(target, ts, lim=10, longshort=-1, Pl=55, Nl=45, tcost=0.01)
res
res<-testStrategy(target, ts, lim=10, longshort=-1, Pl=65, Nl=35, tcost=0.05)
res
res<-testStrategy_comp(target, ts, lim=10, longshort=-1, Pl=55, Nl=45, tcost=0.01)
res

##using Rolling windows. Window size 252 (a year of daily data) or 252/2 for 6 mon
##Number of windows = full_period - window_size
meanperf_1 = RollingTestMAStrategy1(target, ts, lim=10, longshort=-1, Pl=65, Nl=35,w_size=252/2)
meanperf_2 = RollingTestMAStrategy2(target, ts, lim=10, longshort=-1, Pl=65, Nl=35,w_size=252/2)
meanperf_3 = RollingTestMAStrategy3(target, ts, lim=10, longshort=-1, Pl=50, Nl=50,w_size=252)
meanperf_4 = RollingTestMAStrategy4(target, ts, lim=10, longshort=-1, Pl=50, Nl=50,w_size=252)

View(meanperf_1$AvgPerf)
View(meanperf_2$AvgPerf)
View(meanperf_3$AvgPerf)
View(meanperf_4$AvgPerf)

meanperf_1$NumWindows
meanperf_2$NumWindows
meanperf_3$NumWindows


