source('~/Scripts/Bayescan_plot_R.r', chdir = TRUE)

results<-plot_bayescan("~/Desktop/bayescan/Bg_SNPs_10persite_fst.txt", FDR=0.05, add_text=FALSE)
outliers<-results$outliers
results$nb_outliers


setwd('~/Desktop/bgSNPs/')

#make a matrix of your 012 SNP data
#first column holds index, not genotype info, also there is no header
meta<-read.delim("~/Desktop/bgSNPs/meta.txt",header=T)
snps<-read.delim('~/Desktop/bgSNPs/Bg_HMSORLB_filtered_10persite.012',na.strings=-1,header=F)
snps<-snps[,-1]
pos<-read.delim('~/Desktop/bgSNPs/Bg_HMSORLB_filtered_10persite.012.pos',header=F)
indv<-read.delim('~/Desktop/bgSNPs/Bg_HMSORLB_filtered_10persite.012.indv',header=F)
colnames(snps)<-paste(pos[,1],pos[,2],sep='-')
rownames(snps)<-meta$Sample
snps<-as.matrix(snps)

###########Format SNP Table with SNPs in Rows and Samples in Columns#########
snps2<-t(snps) #transposes the SNP table to Samples as Columns and Contigs as rows

##########Subset Outliers################### 
nosnps-c(1:16839)
snps3<-as.data.frame(cbind(nosnps,snps2))
outSNPs<-snps3[match(results$outliers, snps3$nosnps),]
outSNPs2<-outSNPs[,-1]


#write Outlier SNPs to Excel
write.table(outSNPs, file = "Bg-Bayescan-SNPs.txt")


############PCA of outliers############
outPCA<-t(outSNPs2)
#meta<-read.delim('conch4_all_snps_meta.txt', header = T)
#outPCA<-as.data.frame(cbind(meta$Number,outPCA))

pc.out<-prcomp(outPCA)
summary(pc.out)
plot(pc.out$x[,1],pc.out$x[,2],col=meta$Pop,xlab='PC1',ylab='PC2',pch=19)
legend('topright',legend=unique(meta$Pop),fill=unique(meta$Pop))