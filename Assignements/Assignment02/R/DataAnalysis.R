####ASSIGNMENT 2#####

#installing
install.packages("ggpubr")
install.packages("gridExtra")


#library
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)

#Set WD
setwd("C:/BashShared/FHPC_2019-2020/Assignements/Assignment02/R")

#Reading files
df.01 <- read.csv("../../ulyssesResults/01_array_sum.csv", sep = ";")
df.06 <- read.csv("../../ulyssesResults/06_touch_by_all.csv", sep = ";")
df.preSum <- read.csv("../../ulyssesResults/prefixSum.csv", sep = ";")
df.preSumNew <- read.csv("../../ulyssesResults/prefixSumNew.csv", sep = ";")
df.preSumNewNonLocal <- read.csv("../../ulyssesResults/prefixSumNonLocal.csv", sep = ";")
df.06New <- read.csv("../../ulyssesResults/06_touch_by_all_NEW.csv", sep = ";")
#Filtering my stuff

df.01_8 <- df.01 %>%
  filter(N==100000000)

df.06_8 <- df.06 %>%
  filter(N==100000000)

df.01_9 <- df.01 %>%
  filter(N==1000000000)

df.06_9 <- df.06 %>%
  filter(N==1000000000)

df.06_9_NEW <- df.06New %>%
  filter(N == 1000000000)

#Some useful functions


mean_std_df <- function(df){
  df.times <- df %>%
  group_by(N, THREADS) %>%
  summarise(time = mean(TIME))
  
  return(df.times)
}


#StrongScaling for single N

strongScaling <- function(df){
  
  df.times <- mean_std_df(df)
  
  TH <- unique(df.times$THREADS)
  
  T1 <- df.times %>%
    filter(THREADS == 1) %>%
    pull(time)

  TP <- df.times %>%
    pull(time)
    
  
  sp <- T1/TP
  #Let's add speedup for 1 thread manually
  #sp <- c(1,sp)
  p <- 1:20
  e <- (1/sp - 1/p)/(1-1/p)
  
  data_to_graph <- tibble(Threads = TH, Speedup = sp, e = e)
  #plot(data_to_graph, type="l")
  #p <- ggplot(data_to_graph, aes(Thr,Speedup)) + geom_point()
  return(data_to_graph)
}

df_01_8_sp <- strongScaling(df.01_8)
df_06_8_sp <- strongScaling(df.06_8)

df_01_8_sp <- df_01_8_sp %>%
  mutate(label = "01")
df_06_8_sp <- df_06_8_sp %>%
  mutate(label= "06")

df_8_sp <- bind_rows(df_01_8_sp,df_06_8_sp )
df_8_sp

#png(file="C:\\Users\\franc\\Documents\\TriesteUni\\HPC\\Speedup_0_8.png",
#    width = 400, height=350)
df_8_sp_plot <- ggplot(df_8_sp, aes(x=Threads,y=Speedup,color=label)) + 
  geom_line(size=2) +
  ggtitle("N = 10^8") +
  theme(plot.title = element_text(hjust = 0.5))
df_8_e_plot <- ggplot(df_8_sp, aes(x=Threads,y=e,color=label)) + 
  geom_line(size=2) +
  ggtitle("N = 10^9") +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

#And now for 9
df_01_9_sp <- strongScaling(df.01_9)
df_06_9_sp <- strongScaling(df.06_9)

df_01_9_sp <- df_01_9_sp %>%
  mutate(label = "01")
df_06_9_sp <- df_06_9_sp %>%
  mutate(label= "06")

df_9_sp <- bind_rows(df_01_9_sp,df_06_9_sp )
df_9_sp


#png(file="C:\\Users\\franc\\Documents\\TriesteUni\\HPC\\Speedup_0_9.png",
#    width = 400, height=350)
df_9_sp_plot <- ggplot(df_9_sp, aes(x=Threads,y=Speedup,color=label)) + 
  geom_line(size=2) +
  ggtitle("N = 10^9") +
  theme(plot.title = element_text(hjust = 0.5))
df_9_e_plot <- ggplot(df_9_sp, aes(x=Threads,y=e,color=label)) + 
  geom_line(size=2) +
  ggtitle("N = 10^9") +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

png(file="C:\\Users\\franc\\Documents\\TriesteUni\\HPC\\Speedup_0.png",
    width = 600, height=350)
figure <- ggarrange(df_8_sp_plot,df_9_sp_plot)
figure
dev.off()

png(file="C:\\Users\\franc\\Documents\\TriesteUni\\HPC\\ParallelOverhead_0.png",
    width = 600, height=350)
figure_e <- ggarrange(df_8_e_plot,df_9_e_plot)
figure_e
dev.off()

df.preSum.sp <- strongScaling(df.preSum)
df.preSumNew.sp <- strongScaling(df.preSumNew)
df.preSumNewNonLocal.sp <- strongScaling(df.preSumNewNonLocal)

df.preSum.sp <- df.preSum.sp %>%
  mutate(label = "Naive")
df.preSumNew.sp <- df.preSumNew.sp %>%
  mutate(label = "Optimized")
df.preSumNewNonLocal.sp <- df.preSumNewNonLocal.sp %>%
  mutate(label = "Non-local")

df.preSum_sp <- bind_rows(df.preSum.sp,df.preSumNew.sp, df.preSumNewNonLocal.sp )

png(file="C:\\Users\\franc\\Documents\\TriesteUni\\HPC\\PrefixSumSpeedup.png",
    width = 600, height=350)
ggplot(df.preSum_sp, aes(x=Threads,y=Speedup,color=label)) +
  geom_line(size=1.2) +
  ggtitle("Prefix Sum") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


####PERF#####

df.perf_01 <- read.csv("../perf_06_touch_by_all.csv ", sep= ";")
df.perf_06 <- read.csv("../perf_01_array_sum.csv  ", sep= ";")

df.perf_01
df.perf_06

df.perf_01 <- df.perf_01 %>%
  filter(event == "instructions-per-cycle" | event == "LLC-cache-miss-rate" | event == "L1-cache-miss-rate")

df.perf_06 <- df.perf_06 %>%
  filter(event == "instructions-per-cycle" | event == "LLC-cache-miss-rate" | event == "L1-cache-miss-rate")

df.perf_01 <- df.perf_01 %>%
  mutate(label="touch-first")

df.perf_06 <- df.perf_06 %>%
  mutate(label="touch-by-all")

df.perf <- bind_rows(df.perf_01, df.perf_06)
df.perf <- df.perf %>%
  mutate(Number=as.numeric(Number), err = as.numeric(err))
df.perf

ggplot(data = df.perf, aes(x=event, y = Number,
                           ymin = Number - err*Number, ymax = Number + err*Number)) +
  geom_pointrange() +
  facet_grid(cols = vars(label))


summary(df.perf)


###### NUMA OMp #######

strongScalingNUMA <- function(df){
  df.times <- mean_std_df(df)
  TH <- 10:20
    
  T1 <- df.times %>%
    filter(THREADS == 1) %>%
    pull(time)
    
  TP <- df.times %>%
    filter(THREADS > 9) %>%
    pull(time)
    
  sp <- T1/TP
  p <- 10:20
  e <- (1/sp - 1/p)/(1-1/p)
    
  data_to_graph <- tibble(Threads = TH, Speedup = sp, e = e)
}

df.06.9 <- strongScalingNUMA(df.06_9)
df.06.9.New <- strongScalingNUMA(df.06_9_NEW)

df.06.9
df.06.9.New

df.06.9 <- df.06.9 %>%
  mutate(label = "original")

df.06.9.New <- df.06.9.New %>%
  mutate(label = "modified")

df.06.NUMA.sp <- bind_rows(df.06.9, df.06.9.New)

df_6_NUMA_plot <- ggplot(df.06.NUMA.sp, aes(x=Threads,y=Speedup,color=label)) + 
  geom_line(size=1.3) +
  ggtitle("N = 10^9") +
  theme(plot.title = element_text(hjust = 0.5))

png(file="C:\\Users\\franc\\Documents\\TriesteUni\\HPC\\Speedup_NUMA.png",
    width = 600, height=350)
df_6_NUMA_plot
dev.off()

