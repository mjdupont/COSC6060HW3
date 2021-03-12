#Matthew Dupont and Louis Finney
#R Exercises Assignment
#2021-03-05

library(tidyverse)

queue.df <- read.csv("queue.csv")
                 
aggregated.df$Speedup_c = aggregated.df$MaxTime_c/aggregated.df$MaxTime_s
aggregated.df$Speedup_rr = aggregated.df$MaxTime_rr/aggregated.df$MaxTime_s
