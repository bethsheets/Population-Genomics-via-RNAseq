library(RColorBrewer)
colors <- brewer.pal(12,"Set3")[c(1:8,10:12)]


setwd("<path>")

###Read in Admixture files
meta <- read.delim("meta.txt",header=T)
inds <- read.delim("inds.txt",header=F)

K2 <- read.delim("<file>.2.Q",header=F,sep=" ")

K3 <- read.delim("<file>.3.Q",header=F,sep=" ")

K4 <- read.delim("<file>.4.Q",header=F,sep=" ")

K5 <- read.delim("<file>.5.Q",header=F,sep=" ")

#Make bar plots for K2-5 results

#set the margins
par(mfrow=c(4,1),mar=c(2,3,1,0.2),mgp=c(1,0.2,0),tck=0,cex.lab=1.5,cex.axis=0.8,fg="black",col.lab="black",col.axis="black",oma=c(7,1,0,0))

#K2 bar plot
barplot(t(K2),col=colors[c(1:5)],axes=F,ylab="Ancestry",ylim=c(0,1),axisnames=F,border="gray30")
axis(1,at=###,labels=inds,las=2)
axis(2,line=-1,cex.axis=1.1,las=2)

#K3 bar plot
barplot(t(K3),col=colors[c(1:5)],axes=F,ylab="Ancestry",ylim=c(0,1),axisnames=F,border="gray30")
axis(2,line=-1,cex.axis=1.1,las=2)
abline(v=splits,lwd=3)

#K4 bar plot
barplot(t(K4),col=colors[c(1:5)],axes=F,ylab="Ancestry",ylim=c(0,1),axisnames=F,border="gray30")
axis(2,line=-1,cex.axis=1.1,las=2)
abline(v=splits,lwd=3)

#K5 bar plot
barplot(t(K5),col=colors[c(1:5)],axes=F,ylab="Ancestry",ylim=c(0,1),axisnames=F,border="gray30")
axis(1,at=mids,labels=locs,las=2,cex.axis=1.3)
axis(2,line=-1,cex.axis=1.1,las=2)



