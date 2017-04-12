library(WGCNA)
library(data.table)
library(gdsfmt)
library(SNPRelate)

setwd('<directory of your file')

#make a matrix of your 012 SNP data
#first column holds index, not genotype info, also there is no header

data<-data.frame(fread('<file>.012',na.strings='-1',header=F))
snps<-read.delim('<file>.012.pos',header=F)
indv<-read.delim('<file>.012.indv',header=F)
rownames(data)<-indv[,1]
data<-data[,-1]

#add meta data
meta<-read.delim("meta.txt",header=T)
#this file contains the fields: Sample and Pop; you will need to use these terms to match the script below
#your meta file must be ordered the same way as your 012 matrix, you can generate a meta file in the correct order easily from your 012.indv file

snpgdsCreateGeno("<input>.gds",genmat=as.matrix(data),sample.id=meta$Sample,snp.chromosome=snps[,1],snp.position=snps[,2],snpfirstdim=F)

gfile<-snpgdsOpen("<input>.gds")

#create a PCA
pca<-snpgdsPCA(gfile,num.thread=2,missing.rate=1)  #a missing rate of 1 means you're allowing all missing data
pc.perc<-pca$varprop*100

plot(pca$eigenvect[,2],pca$eigenvect[,1],col=meta$Pop,cex=1.3,cex.lab=1.2,xlab="PC2",ylab="PC1",lwd=2,pch=20,axes=T)
legend("topleft",bty="n",cex=1.2,legend=unique(meta$Pop),fill=unique(meta$Pop))

#Calculate global FST
snpgdsFst(gfile,population=as.factor(meta$Pop),method="W&C84",autosome.only=F)
