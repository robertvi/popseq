#
# plot all redgauntlet, hapil and progeny kmer histograms
# find first minimum and first peak
#

library(ggplot2)

#hapil
df <- read.csv("ha_histo_k31.csv",header=F)
colnames(df) = c("kmercount","freq")
#ggplot(df,aes(x=kmercount,y=freq)) + geom_point() + scale_x_log10() + scale_y_log10()
df[df$freq==min(df$freq[df$kmercount<25]),]
df[df$freq==max(df$freq[df$kmercount>10&df$kmercount<40]),]
df[df$freq==min(df$freq[df$kmercount>25&df$kmercount<50]),]
df[df$freq==max(df$freq[df$kmercount>40&df$kmercount<75]),]
df[df$freq==min(df$freq[df$kmercount>70&df$kmercount<100]),]
g=ggplot(df,aes(x=kmercount,y=freq)) + geom_point() + xlim(0,100) + ylim(0,2e7) + geom_vline(xintercept=c(10,25,40,52,91))
ggsave("hapil_kmer_histo.png",plot=g)

#plot with classification intervals
frc = 0.6 #fraction of best guess intervals to use
max0 = round((10-1)*frc+1) #max kmercount to be assigned confidently as 0-copy
min1 = round((10-25)*frc+25) #min kmercount to be assigned confidently as 1-copy
max1 = round((40-25)*frc+25) #max kmercount to be assigned confidently as 1-copy
min2 = round((40-52)*frc+52) #min kmercount to be assigned confidently as 2-copy
g=ggplot(df,aes(x=kmercount,y=freq)) + geom_point() + xlim(0,100) + ylim(0,2e7) + geom_vline(xintercept=c(max0,min1,max1,min2))
ggsave("hapil_kmer_histo_intervals.png",plot=g)
#max0 min1 max1 min2
#6 16 34 45

#redgauntlet
df <- read.csv("rg_histo_k31.csv",header=F)
colnames(df) = c("kmercount","freq")
#ggplot(df,aes(x=kmercount,y=freq)) + geom_point() + scale_x_log10() + scale_y_log10()
df[df$freq==min(df$freq[df$kmercount<100]),]
df[df$freq==max(df$freq[df$kmercount>50&df$kmercount<150]),]
df[df$freq==min(df$freq[df$kmercount>100&df$kmercount<200]),]
df[df$freq==max(df$freq[df$kmercount>150&df$kmercount<250]),]
df[df$freq==min(df$freq[df$kmercount>250&df$kmercount<300]),]
g=ggplot(df,aes(x=kmercount,y=freq)) + geom_point() + xlim(0,350) + ylim(0,6e6) + geom_vline(xintercept=c(41,100,144,207,283))
ggsave("redgauntlet_kmer_histo.png",plot=g)
#sum(df$freq[df$kmercount>41&df$kmercount<144])
#==> 279675725 candidate 1-copy kmers
#10*279675725 ==> 2,796,757,250 bytes 2.8GiB (kmer + count for all of those)

frc = 0.6 #fraction of best guess intervals to use
max0 = round((41-1)*frc+1) #max kmercount to be assigned confidently as 0-copy
min1 = round((41-100)*frc+100) #min kmercount to be assigned confidently as 1-copy
max1 = round((144-100)*frc+100) #max kmercount to be assigned confidently as 1-copy
min2 = round((144-207)*frc+207) #min kmercount to be assigned confidently as 2-copy
g=ggplot(df,aes(x=kmercount,y=freq)) + geom_point() + xlim(0,300) + ylim(0,6e6) + geom_vline(xintercept=c(max0,min1,max1,min2))
ggsave("redgauntlet_kmer_histo_intervals.png",plot=g)
#max0 min1 max1 min2
#25 65 126 169
#sum(df$freq[df$kmercount>65&df$kmercount<126])
#==> 231,842,633 high confidence 1-copy kmer anchors

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
    minima[3] = maxima[2] + 0.5 * (maxima[2] - maxima[1])

    if(length(minima) == 3 && length(maxima) == 2)
    {
        nxt = data.frame(samp=x,min1=minima[1],max1=maxima[1],min2=minima[2],max2=maxima[2],min3=minima[3])
    }
    else
    {
        print("error")
    }

    stats = rbind(stats,nxt)
}

#ggplot(all,aes(x=kmerct,y=freq)) + geom_point(aes(color=samp)) + scale_x_log10() + scale_y_log10()

#intercept = stats$max1

ggplot(all,aes(x=kmerct,y=freq)) + geom_line(aes(color=samp)) + xlim(0,75) + ylim(0,3e7) + geom_vline(xintercept=stats$max1)
ggsave("progeny_kmer_histo.png")

#load base count info
bases = read.table("base_counts",sep='',header=T)

stats = merge(stats,bases)
ggplot(stats,aes(x=bases,y=max1)) + geom_point()

write.table(stats,file="peak_stats_tmp",quote=F,sep=" ",row.names=F)
