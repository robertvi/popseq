#load each linkage group file and find max cm value
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
        chrlen = rbind(chrlen,data.frame(chrm=chrm,len=maxcm))
    }
}

