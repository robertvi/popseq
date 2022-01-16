library(ggplot2)

#inpname = "best_matchup_counts"
#outname = "kmers2bins.png"

inpname = "best_fuzzy_counts"
outname = "kmers2bins_fuzzy.png"

#========================chromosome lengths========================
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


#========================bin positions========================
pos = read.table("best_lmxll_code_counts",header=F,sep='',colClasses = c("character","character","integer","integer","integer","double","double","double"))
colnames(pos) = c("code","chrm","count0","count1","phase0","mincm","maxcm","meancm")

#count total markers
pos$sum = pos$count0+pos$count1

#classify as high confidence if > 2 markers
pos$conf = 0
pos$conf[pos$sum>1] = 1
pos$conf = as.factor(pos$conf)

#split chromosome into number and subgenome letter
pos$chrm = as.factor(pos$chrm)
pos$chrno = substr(pos$chrm,1,1)
pos$sub = substr(pos$chrm,2,2)

#order bins and number them
pos = pos[order(pos$chrm,pos$meancm),]
pos$bin = 1:nrow(pos)
pos$subbin = 0

#number bins within each chromosome
for(chrno in c("1","2","3","4","5","6","7"))
{
    for(sub in c("A","B","C","D"))
    {
        count = nrow(pos[pos$chrno==chrno&pos$sub==sub,])
        pos$subbin[pos$chrno==chrno&pos$sub==sub] = 1:count
    }
}

pos$subbin = -(pos$subbin %% 4 + 1) * 0.3

#====================kmer matches to bins=====================
#chrm phase meancm kmermatches
df = read.table(inpname,header=F,sep='',colClasses = c("character","integer","double","integer"),stringsAsFactors=F)
colnames(df) = c("chrm","phase","meancm","matches")
df$chrm = as.factor(df$chrm)

#encode phase as the sign of the match count
#add 1 to the match count to ensure that +1 and -1 do not both map to zero
#df$phaselogmatch = (2*df$phase - 1) * log10(df$matches+1)
df$logmatch = log10(df$matches+1)
df$phase = as.factor(df$phase)

#split chromosome into number and subgenome letter
df$chrm = as.factor(df$chrm)
df$chrno = substr(df$chrm,1,1)
df$sub = substr(df$chrm,2,2)

g = ggplot(df) +
    geom_segment(data=chrlen,alpha=0.5,aes(y=0,yend=0,x=0.0,xend=len)) +
    geom_segment(data=pos,alpha=0.5,aes(y=subbin,yend=subbin,x=mincm,xend=maxcm,colour=conf)) +
    geom_point(data=pos,alpha=0.5,aes(y=subbin,x=meancm,color=conf)) +
    geom_col(width=2,position="identity",alpha=0.4,aes(x=meancm,y=logmatch,color=phase,fill=phase)) +
    facet_grid(chrno~sub,scales="free_x") +
    labs(title="popseq: matches of RG unique kmers to genetic map bins\ncircle=mean marker position in bin, pink circle=bin with only 1 map marker, bar colour indicates phase (haplotype) of matching kmer",x="centimorgan",y="bin positions (bottom) + log10 kmer matches (transparent overlapping bars)") +
    theme(legend.position="none")

plot(g)

ggsave(outname,plot=g,width=17.5,height=10) #sizes are inches
