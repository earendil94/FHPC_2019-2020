####ASSIGNMENT 2#####

#installing
install.packages("ggpubr")
install.packages("gridExtra")


#library
library(gridExtra)
library(ggplot2)
library(dplyr)

#Set WD
setwd("C:/BashShared/FHPC_2019-2020/Assignements/Assignment02/R")

#Reading files
#df.01 <- read.csv("../01_array_sum.csv", sep = ";")
#df.06 <- read.csv("../06_touch_by_all.csv", sep = ";")
df.01 <- read.csv("../../ulyssesResults/01_array_sum.csv", sep = ";")
df.06 <- read.csv("../../ulyssesResults/06_touch_by_all.csv", sep = ";")
df.preSum <- read.csv("../../ulyssesResults/prefixSum.csv", sep = ";")
df.preSumNew <- read.csv("../../ulyssesResults/prefixSumNew.csv", sep = ";")

#Filtering my stuff

str(df.01)
str(df.06)

df.01


#StrongScaling
strongScaling <- function(df){
  
  df.times <- df %>%
    group_by(N, THREADS) %>%
    summarise(time = mean(TIME))
  
  K <- unique(df.times$N)
  TH <- unique(df.times$THREADS)
  
  #p <- list()
  par(mfrow=c(2,3))
  for (n in K){
    df.temp <- filter(df.times, N == n)
    T1 <- df.temp %>%
      filter(THREADS == 1) %>%
      pull(time)
    TP <- df.temp %>%
      filter(THREADS != 1) %>%
      pull(time)
    
    sp <- T1/TP
    #Let's add speedup for 1 thread manually
    sp <- c(1,sp)
    
    data_to_graph <- tibble(Thr= TH, Speedup = sp)
    plot(data_to_graph, type="l")
    #p[[n]] <- ggplot(data_to_graph, aes(Thr,Speedup)) + geom_point()

  }
  #do.call(grid.arrange, plots)
}

parallelOverhead <- function(df){
  
  df.times <- df %>%
    group_by(N, THREADS) %>%
    summarise(time = mean(TIME))
  
  K <- unique(df.times$N)
  TH <- unique(df.times$THREADS)
  
  par(mfrow=c(1,2))
  for (n in K){
    df.temp <- filter(df.times, N == n)
    T1 <- df.temp %>%
      filter(THREADS == 1) %>%
      pull(time)
    TP <- df.temp %>%
      filter(THREADS != 1) %>%
      pull(time)
    
    sp <- T1/TP
    #Let's add speedup for 1 thread manually
    sp <- c(1,sp)
    
    data_to_graph <- tibble(Thr= TH, Speedup = sp)
    plot(data_to_graph, type="l")
    #p[[n]] <- ggplot(data_to_graph, aes(Thr,Speedup)) + geom_point()
    
  }
  #do.call(grid.arrange, plots)
}

strongScaling(df.01) #Right now this is less than one
strongScaling(df.06) #Stops scaling (even though slowly) from 
strongScaling(df.preSum)
strongScaling(df.preSumNew)






?arrangeGrob

df.times <- df.01 %>%
  group_by(N, THREADS) %>%
  summarise(time = mean(TIME))

df.times.06 <- df.06 %>%
  group_by(N, THREADS) %>%
  summarise(time = mean(TIME))

df.times.preSum <- df.preSum %>%
  group_by(N, THREADS) %>%
  summarise(time = mean(TIME))


print(df.times.preSum)
print(df.times.06, n=120)


strongScaling(df.06)






