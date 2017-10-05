General guide to R scripts for analyzing your transcriptome data

## Gene expression analysis

### Normalize gene expression (DeSeq)
```
library(DESeq2) 
#get merged_counts.txt from get-bam-counts.sh script
counts<-read.delim("merged_counts.txt")

#make a meta file with samples in the same order as how they are listed in your directory
meta<-read.delim("meta.txt")

#make contig labels into row names and remove that column
rownames(counts)<-counts[,1]
colnames(counts)<-meta$sample
counts<-as.matrix(counts[,-1])

#normalize read counts and filter for high counts (i)e. more than 10 reads/site)
normCounts<-t(counts)/estimateSizeFactorsForMatrix(counts)
normCounts_10<-normCounts[,colMeans(normCounts)>10]

transposed<-t(normCounts_10)
write.table(transposed,file="normCounts_10.txt", quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
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
#create 012 files from a your vcf with vcftools-012genotype-matrix.sh
snps<-read.delim('file.012', header=F)
pos<-read.delim('file.012.pos',header=F)
indv<-read.delim<-('file.012.indv',header=F)

#the order of samples in your meta file should match how they are listed in your computer's directory. this example has a meta file with: sample, pop
meta<-read.delim('meta.txt', header=TRUE)
rownames(snps)<-meta$sample
colnames(snps)<-paste(pos[,1],pos[,2],sep='-')

#create a PCA of SNPs
pc.out<-prcomp(snps)
summary(pc.out)

#plot with samples colored by population
plot(pc.out$x[,1],pc.out$x[,2], col=meta$pop,main="Title of your PCA", xlab="PC1", ylab="PC2", pch=16, cex=1.5)

#add sample IDs to points
text(pc.out$x[,1],pc.out$x[,2],labels=meta$id,pos=4,cex=0.7, offset=0.1)

#add a legend box
legend(x="topright",legend=unique(meta$pop),fill=unique(meta$pop))
```


###Detect outliers with Bayescan
- this is an FST based outlier detection method
- warning: we do not use the FST outputs from this program


### Detect outliers with PCAdapt
- this detects outliers off the first and second principal component
- it may be usefule to use both Bayescan and PCAdapt to remove potential loci under selection to create a neutral SNP dataset



### Look at ancestry with admixture
- Make a plink file from your VCF
	`bash plink-from-vcf.sh infile.vcf outfile`
- Run admixture on the cluster






