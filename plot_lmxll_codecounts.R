library(ggplot2)

#get chromosome lengths
basedir = "/home/vicker/rjv_mnt/cluster/octoploid_mapping/consensus_map4/popn_RGxHA/map/haplotypes"
chrlen = NULL
for(chrno in c("1","2","3","4","5","6","7"))
{
    for(sub in c("A","B","C","D"))
    {
        chrm = paste0(chrno,sub)
        fname = sprintf("%s/%s_haplotypes.csv",basedir,chrm)
        tmp = read.table(fname,header=T,sep=',')
        maxcm = max(tmp$cM)
        chrlen = rbind(chrlen,data.frame(chrm=chrm,chrno=chrno,sub=sub,len=maxcm))
    }
}

chrlen$chrno = as.factor(chrlen$chrno)
chrlen$sub = as.factor(chrlen$sub)

df = read.table("best_lmxll_code_counts",header=F,sep='',colClasses = c("character","character","integer","integer","integer","double","double","double"))
colnames(df) = c("code","chrm","count0","count1","phase0","mincm","maxcm","meancm")

#count total markers
df$sum = df$count0+df$count1

#classify as high confidence if > 2 markers
df$conf = 0
df$conf[df$sum>2] = 1
df$conf = as.factor(df$conf)

#split chromosome into number and subgenome letter
df$chrm = as.factor(df$chrm)
df$chrno = substr(df$chrm,1,1)
df$sub = substr(df$chrm,2,2)

#order bins and number them
df = df[order(df$chrm,df$meancm),] 
df$bin = 1:nrow(df)
df$subbin = 0

#number bins within each chromosome
for(chrno in c("1","2","3","4","5","6","7"))
{
    for(sub in c("A","B","C","D"))
    {
        count = nrow(df[df$chrno==chrno&df$sub==sub,])
        df$subbin[df$chrno==chrno&df$sub==sub] = 1:count
    }
}

#show distribution of counts
ggplot(df,aes(x=count0,y=count1)) + geom_point()

#plot position of bin centres by chromosome
ggplot(df,aes(x=chrm,y=meancm)) + geom_point()

#box plot to show ranges of bins
ggplot(df,aes(x=chrm,y=meancm)) + geom_point() + geom_segment(aes(x=chrm,xend=chrm,y=mincm,yend=maxcm))

ggplot(df,aes(y=subbin,x=meancm)) + geom_point(aes(color=conf)) +
    geom_segment(aes(y=subbin,yend=subbin,x=mincm,xend=maxcm,color=conf)) +
    geom_segment(data=chrlen,aes(y=0,yend=0,x=0.0,xend=len)) +
    facet_grid(chrno~sub,scales="free") +
    labs(x="centimorgan",y="bin")+
    theme(legend.position="none")

