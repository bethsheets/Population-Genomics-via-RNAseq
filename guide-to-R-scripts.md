General guide to R scripts for analyzing your transcriptome data

## Gene expression analysis

### Compare gene expression using DESeq2
```
library(DESeq2) 
counts<-read.delim(‘counts.txt,row.names=1) 
cds<-DESeqDataSetFromMatrix(counts.txt,meta,~cluster) 
cds<-DESeq(cds) 
res<-results(cds) 
sig<-res[which(res$padj<0.05),] 
write.table(sig,file=‘DEcontigs.txt’,quote=F,sep=‘\t') 
```

### WGCNA 
- in program R
- [website](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)
- script TBD

### Identify expression SNPs
- see Noah Rose's R package: https://github.com/noahrose/vcf2eqtl


## SNP ANALYSIS 

### A simple way to format your 012 SNP matrix
```
snps<-read.delim('file.012', header=F)
pos<-read.delim('file.012.pos',header=F)
indv<-read.delim<-('file.012.indv',header=F)

colnames(snps)<-paste(pos[,1],pos[,2],sep='-')
rownames(snps)<-indv[,1]
snps<-as.matrix(snps)

#PCA of SNPs
pc.out<-prcomp(snps)
summary(pc.out)
plot(pc.out$x[,1],pc.out$x[,2])	#PC1 v PC2

```

### Add meta data to your SNP matrix
- 	make a meta data file with info about individuals (location, date, etc.)
- 	make sure your meta file is ordered the same as your vcfs! (i.e. ls your samples in the terminal to see their order)
- script TBD


### Allow for missing SNP data with SNPrelate
- 

### Test for loci under selection using BayeScan
- [download program](http://cmpg.unibe.ch/software/BayeScan/download.html)
- 	identifies putative loci under selection




