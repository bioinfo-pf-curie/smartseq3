#!/usr/bin/env Rscript

.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(plyr)

stats<-read.csv("table_mqc.stats", sep = "\t")

mean_stat<-stats[1,-1]
mean_stat[,1]<-c("All_samples")
mean_stat[,-1]<-apply(stats[,-c(1,2)], 2, mean )

# add sd ????

write.table(mean_stat, "mean_table_mqc.stats")