library(DESeq2)

#read in your counts file
counts<-read.delim("merged_counts.txt")

#add meta file with population information (with header Pop)
meta<-read.delim('meta.txt', header=TRUE)

#make contig labels into row names and remove that column
rownames(counts)<-counts[,1]
counts<-as.matrix(counts[,-1])
head(counts)

#normalize read counts and filter for high counts
normCounts<-t(counts)/estimateSizeFactorsForMatrix(counts)

#PCA of normalized gene expression
pc.out<-prcomp(normCounts)
plot(pc.out$x[,1],pc.out$x[,2],col=meta$Pop,pch=19,ylab="PC1", xlab="PC2", cex.lab=0.8,cex=1.2)
legend("topleft",bty="n",cex=1.2,legend=unique(meta$Pop),fill=unique(meta$Pop))