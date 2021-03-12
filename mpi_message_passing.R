#Matthew Dupont
#R Exercises Assignment
#2021-03-12

library(tidyverse)

queue.df <- read.csv("results/queue.csv")
queue.df$Procs <- as.factor(queue.df$Procs)

queue.df <- queue.df %>% rename(MessageSize=ChunkSize)

queue.df.50k <- queue.df %>% filter(Terms == 50000)
queue.df.10k <- queue.df %>% filter(Terms == 10000)
queue.df.1k <- queue.df %>% filter(Terms == 1000)
queue.df.100 <- queue.df %>% filter(Terms == 100)
queue.df.10 <- queue.df %>% filter(Terms == 10)
              
plot50k <- ggplot(queue.df.50k, aes(x=MessageSize, y=MaxTime)) + 
  geom_point(aes(color=Procs)) + 
  labs(title = "Time per Message Size, 50000 terms", x="Message Size", y="Max Time")

show(plot50k)
ggsave("images/TimePerMessageSize50k.png")

plot10k <- ggplot(queue.df.10k, aes(x=MessageSize, y=MaxTime)) + 
  geom_point(aes(color=Procs)) + 
  labs(title = "Time per Message Size, 10000 terms", x="Message Size", y="Max Time")

show(plot10k)
ggsave("images/TimePerMessageSize10k.png")

plot1k <- ggplot(queue.df.1k, aes(x=MessageSize, y=MaxTime)) + 
  geom_point(aes(color=Procs)) + 
  labs(title = "Time per Message Size, 1000 terms", x="Message Size", y="Max Time")

show(plot1k)
ggsave("images/TimePerMessageSize1k.png")

plot100 <- ggplot(queue.df.100, aes(x=MessageSize, y=MaxTime)) + 
  geom_point(aes(color=Procs)) + 
  labs(title = "Time per Message Size, 100 terms", x="Message Size", y="Max Time")

show(plot100)
ggsave("images/TimePerMessageSize100.png")

plot10 <- ggplot(queue.df.10, aes(x=MessageSize, y=MaxTime)) + 
  geom_point(aes(color=Procs)) + 
  labs(title = "Time per Message Size, 10 terms", x="Message Size", y="Max Time")

show(plot10)
ggsave("images/TimePerMessageSize10.png")


queue.det.df <- read.csv("results/queue.csv")
queue.det.df$Procs <- as.factor(queue.det.df$Procs)

queue.det.df <- queue.det.df %>% rename(MessageSize=ChunkSize)

plotdet <- ggplot(queue.det.df, aes(x=MessageSize, y=MaxTime)) + 
  geom_point(aes(color=Procs)) + 
  labs(title = "Detailed Time per Message Size, 1000 terms", x="Message Size", y="Max Time")

show(plotdet)
ggsave("images/TimePerMessageSize1000Detailed.png")
