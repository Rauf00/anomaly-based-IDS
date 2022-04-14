library(lubridate)
library(TTR)
library(ggplot2)
library(zoo)
library("depmixS4")
library(hms)

################################################################################
#                                                                              #
#                         Graph Global Active                                  #
#                         and Global Intensity                                 #
#                                                                              #
################################################################################

data <- read.table("TermProjectData.txt", sep =",", header=TRUE)
data$Global_intensity <- scale(data$Global_intensity)
data$Global_active_power <- scale(data$Global_active_power)
data$Sub_metering_3 <- scale(data$Sub_metering_3)

data$Date <- as.POSIXct(data$Date, format="%d/%m/%Y")
data$Time <- as.POSIXct(data$Time, format="%H:%M:%S")
data$week <- wday(ymd(data$Date), week_start=1)
# Saturday will be our chosen day for the term project
saturdays <- data[data$week == 6,]
saturdays$Time  <- as_hms(ymd_hms(saturdays$Time))
saturdays <- na.omit(saturdays)

AggSaturday <- setNames(aggregate(saturdays$Global_active_power, by=list(saturdays$Time), mean), c('Time', 'Mean_Global_Active_Power'))
ggplot(AggSaturday, aes(Time)) + 
  geom_line(aes(y = Mean_Global_Active_Power, colour = "Saturday"), size = 0.5) +
  ggtitle("Mean Global Active Power Per Hour") +
  ylab("Mean Global Active Power") +
  xlab("24 Hour")

AggSaturday <- setNames(aggregate(saturdays$Global_intensity, by=list(saturdays$Time), mean), c('Time', 'Mean_Global_Intensity'))

ggplot(AggSaturday, aes(Time)) + 
  geom_line(aes(y = Mean_Global_Intensity, colour = "Saturday"), size = 0.5) +
  ggtitle("Mean Global Intensity Per Hour") +
  ylab("Mean Global Intensity") +
  xlab("24 Hour")


AggSaturday <- setNames(aggregate(saturdays$Sub_metering_3, by=list(saturdays$Time), mean), c('Time', 'Mean Submetering 3'))

ggplot(AggSaturday, aes(Time)) + 
  geom_line(aes(y =`Mean Submetering 3`, colour = "Saturday"), size = 0.5) +
  ggtitle("Mean Sub Metering Per Hour") +
  ylab("mean Submetering") +
  xlab("24 Hour")






