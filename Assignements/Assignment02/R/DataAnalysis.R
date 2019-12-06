####ASSIGNMENT 2#####

#library
library(ggplot2)
library(dplyr)

#Reading files
df.01 <- read.csv("../01_array_sum.csv", sep = ";")
df.06 <- read.csv("../06_touch_by_all.csv", sep = ";")

#Filtering my stuff

str(df.01)
str(df.06)
#StrongScaling

par(mfrow = c(1,2))

strongScaling <- function(df){
  
  df.times <- df %>%
    group_by(N, THREADS) %>%
    summarise(time = mean(TIME))
  
  K <- unique(df.times$N)
  TH <- unique(df.times$THREADS)
  

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
    print(ggplot(data_to_graph, aes(Thr,Speedup)) + geom_point())
  }
}

mf
strongScaling(df.01)
strongScaling(df.06)

df.01 %>%
  group_by(N, THREADS) %>%
  summarise(avg = mean(TIME))

df.06 %>%
  group_by(N, THREADS) %>%
  summarise(avg = mean(TIME))




