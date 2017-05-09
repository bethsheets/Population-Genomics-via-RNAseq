library(hierfstat)

snps<-read.delim('<file>.012',header=F)
snps<-snps[,-1]
pos<-read.delim('<file>.012.pos',header=F)
indv<-read.delim('<file>.012.indv',header=F)
meta<-read.delim("meta.txt",header=T)

colnames(snps)<-paste(pos[,1],pos[,2],sep='-')
rownames(snps)<-meta$Sample
snps<-as.matrix(snps)

hf<-snps
hf[hf==0]<-11
hf[hf==1]<-12
hf[hf==2]<-22
hf[hf==-1]<-NA
pop<-as.numeric(meta$Pop)
hf2<-as.data.frame(cbind(pop,hf))

write.bayescan(dat=hf2, diploid=TRUE, fn='<outfile>.bsc')


