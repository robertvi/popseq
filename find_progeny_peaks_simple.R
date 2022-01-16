#
# plot all progeny kmer histograms
# find first 3 minima and first 2 peaks
# use simple method based only on position of first peak

library(ggplot2)

#progeny
samples = as.character(read.csv("sample_list",header=F)$V1)

all = NULL
stats = NULL

for(x in samples)
{
    df = read.csv(sprintf("%s_histo_k31.csv",x),header=F)
    colnames(df) = c("kmerct","freq")
    df$samp = x
    all = rbind(all,df)

    #find first minimum and first peak
    #min1 = df$kmerct[  df$freq == min(df$freq[df$kmerct<10])  ]
    #max1 = df$kmerct[  df$freq == max(df$freq[df$kmerct>10&df$kmerct<30])  ]

    prev = df$freq[-9999:-10000]
    val = df$freq[2:9999]
    nxt = df$freq[-1:-2]
    minima = which(prev>val&nxt>val)[1] + 1 # first minimum
    maxima = which(prev<val&nxt<val)[1] + 1 # first maximum

    #extrapolate from first peak to remaining positions
    maxima[2] = round(2.0* maxima[1])
    minima[2] = round((maxima[1]+maxima[2])/2.0)
    minima[3] = round(maxima[2] + (maxima[2] - minima[2])*1.5)

    nxt = data.frame(samp=x,min1=minima[1],max1=maxima[1],min2=minima[2],max2=maxima[2],min3=minima[3])

    stats = rbind(stats,nxt)

    ggplot(df,aes(x=kmerct,y=freq)) + geom_line() + xlim(0,80) + ylim(0,2.5e7) + geom_vline(xintercept=c(minima,maxima))
    ggsave(sprintf("%s_kmer_histograms.png",x))

    #readline(prompt="press enter")
}



ggplot(all,aes(x=kmerct,y=freq)) + geom_line(aes(color=samp)) + xlim(0,80) + ylim(0,2.5e7)
ggsave("all_progeny_kmer_histograms.png")

write.table(stats,file="peak_stats_simple",quote=F,sep=" ",row.names=F)
