#
# plot all redgauntlet, hapil and progeny kmer histograms
# find first minimum and first peak
#

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
    minima = which(prev>val&nxt>val)[1:10] + 1 # first ten minima
    maxima = which(prev<val&nxt<val)[1:10] + 1 # first ten maxima

    #limit to kmer count less than 60
    minima = minima[minima<60]
    maxima = maxima[maxima<60]

    #RGXH056a has no true second peak, approximate to most similar other sample (RGXH052a)
    if(x == "RGXH056a")
    {
        #manually set approx positions of second "minimum" and "peak"
        minima[2] = 26
        maxima[2] = 32
    }

    #estimate minimum 3 as halfway to expected positino of peak 3 from peak 2
    minima[3] = round(maxima[2] + 0.5 * (maxima[2] - maxima[1]))

    if(length(minima) == 3 && length(maxima) == 2)
    {
        nxt = data.frame(samp=x,min1=minima[1],max1=maxima[1],min2=minima[2],max2=maxima[2],min3=minima[3])
    }
    else
    {
        print("error")
    }

    stats = rbind(stats,nxt)

    g=ggplot(df,aes(x=kmerct,y=freq)) + geom_line() + xlim(0,100) + scale_y_log10() + geom_vline(xintercept=c(minima,maxima))
    plot(g)

    #readline(prompt="press enter")
}



#ggplot(all,aes(x=kmerct,y=freq)) + geom_line(aes(color=samp)) + scale_x_log10() + scale_y_log10()

write.table(stats,file="peak_stats_auto",quote=F,sep=" ",row.names=F)
