#
# plot mapbin code counts
#

library(ggplot2)

df <- read.table("code_counts",header=F,sep='')

ggplot(df,aes(x=V2,y=V3)) + geom_point()

ggplot(df,aes(x=V2)) + geom_histogram(bins=1000) + scale_x_log10()
