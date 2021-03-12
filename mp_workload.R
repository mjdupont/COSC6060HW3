#Matthew Dupont
#R Exercises Assignment
#2021-03-12

library(tidyverse)

queue.df <- read.csv("results/queue.csv")

queue.df.50k <- queue.df %>% filter(Terms == 50000)
              
plot <- ggplot(queue.df.50k, aes(x=ChunkSize, y=MaxTime)) + 
  geom_point(aes(color=Procs)) + 
  labs(title = "Time per chunkSize, 50000 terms", x="Chunk Size", y="Max Time")

show(plot)
ggsave("images/TimePerChunkSize")
