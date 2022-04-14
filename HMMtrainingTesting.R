library(lubridate)
library(TTR)
library(ggplot2)
library(zoo)
library("depmixS4")

data <- read.table("TermProjectData.txt", sep =",", header=TRUE)
data$Global_intensity <- scale(data$Global_intensity)
data$Global_active_power <- scale(data$Global_active_power)
data$Date <- as.POSIXct(data$Date, format="%d/%m/%Y")
data$Time <- as.POSIXct(data$Time, format="%H:%M:%S")
data$week <- wday(ymd(data$Date), week_start=1)
# Saturday will be our chosen day for the term project
saturdays <- data[data$week == 6,]  

#set window from 7am to 12pm (5hour)
saturdayMornings <- subset(saturdays, saturdays$Time > as.POSIXlt("7:00:00", format="%H:%M:%S") & saturdays$Time < as.POSIXlt("12:00:00", format="%H:%M:%S"))

# Partition data into training and testing
# find all rows that are from 2009
oneYearData <- c(which(format(as.Date(saturdayMornings$Date, format="%Y-%m-%d"),"%Y") == 2009))

# Training data is from 2006, 2007, 2008
trainingData <- saturdayMornings[-oneYearData,]
trainingData <- na.omit(trainingData)

# Testing data 2009
testData <- saturdayMornings[oneYearData,]
testData <- na.omit(testData)

# need to figure out the ntimes, the number of rows per window
ntimesVector_TrainingData <- aggregate(trainingData$Date, by=list(trainingData$Date), FUN=length)
ntimesVector_TestData <- aggregate(testData$Date, by=list(testData$Date), FUN=length)
colnames(ntimesVector_TrainingData) <- c("Date", "Observations")
colnames(ntimesVector_TestData) <- c("Date", "Observations")



################################################################################
#                                                                              #
#                         Training Multivariate HMM                            #
#                                                                              #
################################################################################
# training 11 models from 4 - 24 steps of 2
set.seed(1)

model1 <- depmix(list(Global_intensity~1, Global_active_power~1),
                         data =trainingData, nstates=4, 
                         family=list(gaussian(),gaussian()),
                         ntimes= ntimesVector_TrainingData$Observations)

fm1 <- fit(model1)
print(BIC(fm1))
logLik(fm1)


model2 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=6, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm2 <- fit(model2)
print(BIC(fm2))

model3 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=8, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm3 <- fit(model3)
print(BIC(fm3))

model4 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=10, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm4 <- fit(model4)
print(BIC(fm4))

model5 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=12, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm5 <- fit(model5)
print(BIC(fm5))

model6 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=14, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm6 <- fit(model6)
print(BIC(fm6))

model7 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=16, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm7 <- fit(model7)
logLik(fm7)
print(BIC(fm7))

model8 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=18, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm8 <- fit(model8)
print(BIC(fm8))

model9 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=20, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm9 <- fit(model9)
print(BIC(fm9))

model10 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=22, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm10 <- fit(model10)
print(BIC(fm10))

model11 <- depmix(list(Global_intensity~1, Global_active_power~1),
                 data =trainingData, nstates=24, 
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector_TrainingData$Observations)

fm11 <- fit(model11)
print(BIC(fm11))
logLik(fm11)


################################################################################
#                                                                              #
#                        Write data to csv                                     #
#                                                                              #
################################################################################


#write.csv(logLkDf, "loglikehood.csv")
#write.csv(bicDf, "bic.csv")

################################################################################
#                                                                              #
#                         Testing Multivariate HMM                             #
#                                                                              #
################################################################################

# Best two models
lengthTraining <- nrow(ntimesVector_TrainingData)
lengthTesting <- nrow(ntimesVector_TestData)

model6_test = depmix(list(Global_intensity~1, Global_active_power~1),
                     data =testData, nstates=14, 
                     family=list(gaussian(),gaussian()),
                     ntimes= ntimesVector_TestData$Observations)
model6_test = setpars(model6_test, getpars(fm6))
fwbw <- forwardbackward(model6_test)
training6 <- logLik(fm6)/lengthTraining
print(training6)
test6 <- logLik(model6_test)/lengthTesting
print(test6)

model7_test = depmix(list(Global_intensity~1, Global_active_power~1),
                     data =testData, nstates=16, 
                     family=list(gaussian(),gaussian()),
                     ntimes= ntimesVector_TestData$Observations)
model7_test = setpars(model7_test, getpars(fm7))
fwbw <- forwardbackward(model7_test)
training7 <- logLik(fm7)/lengthTraining
print(training7)
test7 <- logLik(model7_test)/lengthTesting
print(test7)


################################################################################
#                                                                              #
#                         Graph Log-likelihood and BIC                         #
#                                                                              #
################################################################################

loglikeliHoodVec = c(logLik(fm1),logLik(fm2),logLik(fm3),logLik(fm4),logLik(fm5),logLik(fm6),logLik(fm7), logLik(fm8),logLik(fm9),logLik(fm10),logLik(fm11))
bicVec =  c(BIC(fm1), BIC(fm2), BIC(fm3), BIC(fm4), BIC(fm5), BIC(fm6),BIC(fm7),BIC(fm8),BIC(fm9),BIC(fm10),BIC(fm11))

xAxis = seq(4, 24, by=2)
logLkDf <- data.frame(loglikeliHoodVec)
bicDf <- data.frame(bicVec)

ggplot() +
  geom_line(data = logLkDf, aes(x = xAxis, y = loglikeliHoodVec, color=" Training Log Likelihood")) +
  geom_line(data = bicDf, aes(x = xAxis, y = bicVec, color="BIC")) + 
  labs(y= "Values", x = "Number of states", title= "Values vs Number of States")

################################################################################
#                                                                              #
#                         Graph Log-Normalized Testing                         #
#                           and Training                                       #
#                                                                              #
################################################################################

loglikeliHoodVec = c(logLik(fm6)/104,logLik(fm7)/104)
testLikeliHoodVec = c(logLik(model6_test)/48,logLik(model7_test)/48)

xAxis = seq(14, 16, by=2)
logLkDf <- data.frame(loglikeliHoodVec)
logLikTestDf <-data.frame(testLikeliHoodVec)

ggplot() +
  geom_line(data = logLkDf, aes(x = xAxis, y = loglikeliHoodVec, color="Normalized Training Log Likelihood")) +
  geom_line(data = logLikTestDf, aes(x = xAxis, y = testLikeliHoodVec, color="Normalized Testing Log Likelihood")) + 
  labs(y= "Values", x = "Number of states", title="Normalized Testing vs Normalized Training Log likelihood")


################################################################################
#                                                                              #
#                         Anomaly Detection                                    #
#                                                                              #
################################################################################

testForAnomalies <- function(filename) {
  data <- read.table(filename, sep =",", header=TRUE)
  data$Global_intensity <- scale(data$Global_intensity)
  data$Global_active_power <- scale(data$Global_active_power)
  data$Date <- as.POSIXct(data$Date, format="%d/%m/%Y")
  data$Time <- as.POSIXct(data$Time, format="%H:%M:%S")
  data$week <- wday(ymd(data$Date), week_start=1)
  # Saturday will be our chosen day for the term project
  saturdays <- data[data$week == 6,]
  saturdayMornings <- subset(saturdays, saturdays$Time > as.POSIXlt("7:00:00", format="%H:%M:%S") & saturdays$Time < as.POSIXlt("12:00:00", format="%H:%M:%S"))
  saturdayMornings <- na.omit( saturdayMornings)
  ntimesVector <- aggregate(saturdayMornings$Date, by=list(saturdayMornings$Date), FUN=length)
  colnames(ntimesVector) <- c("Date", "Observations")
  set.seed(1)
  model = depmix(list(Global_intensity~1, Global_active_power~1),
                 data =saturdayMornings, nstates=16,
                 family=list(gaussian(),gaussian()),
                 ntimes= ntimesVector$Observations)
  
  model = setpars(model, getpars(fm7))
  fwbw <- forwardbackward(model)
  return(fwbw$logLike/nrow(ntimesVector))
}

logLk1 <- testForAnomalies("DataWithAnomalies1.txt")
print(logLk1)
print(training7)

logLk2 <- testForAnomalies("DataWithAnomalies2.txt")
print(logLk2)
print(training7)

logLk3 <- testForAnomalies("DataWithAnomalies3.txt")
print(logLk3)
print(training7)









