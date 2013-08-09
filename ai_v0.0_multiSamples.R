##################################################################
## ASE model for use on multiple samples 
##
## Created: 08/07/2013
##
## Author: CTH
##
## Version0.0: to be run interactively
##################################################################
## Allele specific funtions. 
source('/nfs/hpnfs/groups/piquelab/charvey/ASE/src/fitAseModels.v2.R')
## Some ploting functions
source('/nfs/hpnfs/groups/piquelab/gmb/AI/results/qqman.r')
##################################################################
setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/cleanBedfiles")
##################################################################
pileupFile='HUVEC_0_combined.gz'
unionFile='HUVEC_union.bed.gz'

command=paste("gunzip -cd ",pileupFile)
unionCommand=paste("gunzip -cd ",unionFile)

pileupDat <- read.table(file=pipe(command),header=F,quote="",comment.char="",as.is=T,sep="\t") 
names(pileupDat) <- c("Ref", "Alt", "Err")

pileupChrom <- read.table(file=pipe(unionCommand),header=F,quote="",comment.char="",as.is=T,sep="\t") 
names(pileupChrom) <- c("chr", "pos-1", "pos","af")

pileup<-cbind(pileupChrom, pileupDat)

chrList <- unique(pileup$chr)
chrList <- paste("chr",sort(as.numeric(gsub("chr","",chrList))),sep="")
factor(chrList,levels=chrList)
pileup$chr <- factor(pileup$chr,levels=chrList)
pileup <- pileup[order(pileup$chr,pileup$pos),];

#####################################################################
#system(paste("mkdir -p",sName))
## Write the original data just in case we want to view later
#write.table(pileup,file=gzfile(paste(fileName,'.loci.gz',sep='')),row.names=F,col.names=T,sep='\t',quote=F)
##write.table(file=paste(fileName,'.stats.txt',sep=''),saveData,row.names=T,col.names=T,quote=F,sep='\t') # Must be data frame
##saveData$mean.coverage <- mean(pileup$num.reads)

####################################
oraf <- pileup$Ref/(pileup$Ref + pileup$Alt)
omac <- pmin(pileup$Alt,pileup$Ref)
otac <- (pileup$Ref + pileup$Alt)

tr.l <- 0		##15
tr.u <- 4e+13	##400
ind <- otac>tr.l & otac<tr.u;

names(ind) <- pileup$rsID

ref=pileup$Ref[ind];
alt=pileup$Alt[ind];
phi <- pileup$af[ind]

####################################

##bed <- pileup[otac>65,1:3]
##bed$chr <- sub("chr","",bed$chr)
##tmpFile <- paste(Sys.getenv("TMPDIR","/tmp"),"/",sName,'.bed.gz',sep='')
##write.table(bed,file=gzfile(tmpFile),row.names=F,col.names=F,sep='\t',quote=F)

##source('/wsu/home/groups/piquelab/testCentipede/AI/results/fitAseModels.v2.R')
##aux <- fitAse(ref,alt);

asenull <- fitAseNull(ref,alt);
##asenull <- fitAseNull(ref,alt,max.it=1,eps=0.1);

## Use allel frequencies for genotype prior
gmat <- cbind(g0=(1-phi)^2,g1=2*phi*(1-phi),g2=phi^2)
ase2 <- fitAseNull(ref,alt,log.gmat=log(gmat))

eps <- ase2$eps;
##ase3 <- fitAseNull(ref,alt,log.gmat=log(gmat),fixGprior=F)
##ase2 <- fitAse(ref,alt,log.gmat=asenull$log.gt,eps=0.1,fixEps=TRUE,fixT1=TRUE,max.it=300,t1=1E-6);
##ase2 <- fitAse(ref,alt,eps=0.1,fixEps=TRUE,fixT1=TRUE,max.it=300,t1=1E-6);

##eps <- 0.2

## Calculate rho using eps
rho <- comp.rho(ref,alt,eps)

## Simple rho estimate
rho0 <- plogis(log(ref)-log(alt));

##logliktop = ref * log((1 - eps) * rho + eps * (1 - rho)) + alt * log((1 - eps) * (1 - rho) + eps * rho);
##logRatio <- logliktop-asenull$loglik
##pval <- (1-pchisq(2*logRatio,df=1))

lrt <- lrtEpsRhoBinom(ref,alt,eps,rho)
pval <- (1-pchisq(2*lrt,df=1))

qv <- qvalue(pval,pi0.method="bootstrap")
sum(qv$qv<0.1)

log.het <- pmin(ase2$log.gt[,2],asenull$log.gt[,2])
het <- exp(log.het)

hetInd <- het > 0.99
sum(hetInd)

## ## Find MLE for Beta-Binomial dispersion parameter
## aux <- optim(1,fn=fnLogLikBetaBinomialNull,gr=grLogLikBetaBinomialNull,R=ref[hetInd],A=alt[hetInd],method="L-BFGS-B",lower=1,upper=10,hessian=T)
## se <- sqrt(1/aux$hessian)
## Dmax <- exp(aux$par-0*se)
## Dmax
## str(aux)

############################################

D <- 1:100;
D <- exp((0:500)/50)
##D <- seq(1,40,0.01)
aux <- sapply(D,function(D){
  sum(logLikBetaBinomialRhoEps(0.5,eps,D,ref[hetInd],alt[hetInd]))
  ##fnLogLikBetaBinomialNull(log(D),ref[hetInd],alt[hetInd])
  ##grLogLikBetaBinomialNull(log(D),ref[hetInd],alt[hetInd])
})
D[which.max(aux)]
##plot(D,aux,log="x",pch='.',cex=3)

Dmax <- min(D[which(aux>=max(aux)-2)])
Dmax <- D[which.max(aux)]
Dmax



## Find MLE for rho from the Beta-Binomial. 
Dmax2 <- Dmax
aux <- optim(rep(0,sum(hetInd)),fn=logLikBetaBinomial2,gr=gLogLikBetaBinomial,D=Dmax2,R=ref[hetInd],A=alt[hetInd],method="L-BFGS-B")
rho3 <- plogis(aux$par)
## Could I use rho2 instead?


## Recaluclate Het LRT using the beta-bionomial 
lrt3 <- logLikBetaBinomialRhoEps(rho3,eps,Dmax2,ref[hetInd],alt[hetInd]) - logLikBetaBinomialRhoEps(0.5,eps,Dmax,ref[hetInd],alt[hetInd])
pval3 <- (1-pchisq(2*lrt3,df=1))

rho2 <- rho
rho2[hetInd] <- rho3

aux <- pmax(logLikBetaBinomialRhoEps(0.0,eps,Dmax,ref,alt),
            logLikBetaBinomialRhoEps(1.0,eps,Dmax,ref,alt),
            logLikBetaBinomialRhoEps(0.5,eps,Dmax,ref,alt))
lrt2 <- logLikBetaBinomialRhoEps(rho2,eps,Dmax2,ref,alt) - aux;
pval2 <- (1-pchisq(2*lrt2,df=1))



qv2 <- qvalue(pval2[hetInd],pi0.method="bootstrap")
sum(qv2$qv<0.1)
100-qv2$pi0*100

####################################

aux <- (pbeta(1-eps,1+alt,1+ref)-pbeta(eps,1+alt,1+ref))/(1-2*eps)
##bf <- aux/exp(asenull$loglik)
bf <- aux/exp(ase2$loglik)
##bf <- aux/exp(pmax(log.gt[,"g0"],log.gt[,"g1t0"],log.gt[,"g2"]))

##hist(log(bf[bf>4]),breaks=100)
##plot(-log10(pval2),log(bf[hetInd]))

aux <- cumsum(sort(bf))
aux <- aux/(1:length(aux))
p0 <- mean(aux<1.1)

####################################
####################################
####################################
####################################

qval <- het*0+1
qval[hetInd] <- qv2$qval

pileup2 <- cbind(pileup[ind,c("chr","pos-1","pos","ref","alt","rsID","af","errors","ref.matches","alt.matches")],otac=otac[ind],het=het,ase2$gt,pval=pval2,bf=bf,qval)
pileup2[pileup2$qval<0.1,]

write.table(file=gzfile(paste(fileName,'.all.bed.gz',sep='')),pileup2,row.names=F,col.names=F,quote=F,sep='\t')

write.table(file=gzfile(paste(fileName,'.sig.bed.gz',sep='')),pileup2[pileup2$qval<0.2,],row.names=F,col.names=F,quote=F,sep='\t')


##hist(lrt,breaks=200) 

##X11(type="Xlib",display="localhost:10.0")
setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/out")

## Plot quality values to see if they make any sense
png(file=paste(fileName,'.plots.quality.png',sep=''))

plot(as.numeric(names(qual.table)),as.numeric(qual.table),xlab="Base calling qualities",ylab="Freq.",pch='.',cex=6,log='y')
title(sName)
abline(h=0,lty=3)

dev.off()

png(file=paste(fileName,'.plots.ecdf.png',sep=''))
plot(ecdf(log10(otac[otac>0])),main=sName,xlab=expression(log[10](x)),ylab=expression(F[x>0](x)))
abline(v=log10(c(tr.l,tr.u)),lty=3)
dev.off()

##chrCol <- rainbow(length(chrList))
chrCol <- rep(c("orange","darkblue"),length(chrList))[1:length(chrList)]
names(chrCol) <- chrList

png(file=paste(fileName,'.plots.oraf.png',sep=''),width=800,height=400)
layout(t(c(1,1,1,2)))
par(cex=1.0)
oldmar <- par("mar")
mar <- oldmar
mar[4] <- 0.0
par(mar=mar)
ind2 <- ind & (abs(pileup$af-0.5)<0.4)
##x <- 1:sum(ind2)
plot(oraf[ind2],xlab="Chromosome order",ylab="Obs. reference allele Freq.",pch='.',cex=3,col=chrCol[pileup$chr[ind2]],axes=F)
axis(2)
x.at <- c(which(!duplicated(pileup$chr[ind2])),sum(ind2))
axis(1,at=x.at,labels=FALSE,cex=0.3)
abline(v=x.at,lty=3)
text((x.at[-1]+x.at[-length(x.at)])*0.5, -0.15, labels = chrList, srt = 45, pos =3, xpd = TRUE,cex=0.7)
title(sName)
abline(h=0.5,lty=3)
mar <- oldmar
mar[2] <- 0.0
par(mar=mar)
aux <- hist(jitter(oraf[ind2]),breaks=101,plot=F)
barplot(aux$counts+1,horiz=TRUE,log="x")
title(xlab="Freq.")
par(mar=oldmar)
dev.off()

####################################
png(file=paste(fileName,'.plots.rho.png',sep=''))

plot(rho0,rho,pch='.',cex=2,xlab=expression(rho[0]),ylab=expression(hat(rho)),main=sName)
points(rho0[hetInd],rho3,pch='.',cex=2,col='blue')
abline(0,1,col='red',lty=2)

dev.off()


## ## Plot to see how different are the call on heterozygozity
## layout(1)
## smoothScatter(jitter(pmin(phi,1-phi)),jitter(ase2$log.gt[,2]-asenull$log.gt[,2]),pch='.',cex=2)
## ##smoothScatter(jitter(otac[ind]),jitter(ase2$log.gt[,2]-asenull$log.gt[,2]),pch='.',cex=2)
## smoothScatter(jitter(otac[ind]),jitter(ase2$log.gt[,2]-asenull$log.gt[,2]),pch='.',cex=2)
## ##smoothScatter(jitter(otac[ind]),jitter(het),pch='.',cex=2)

##abline(0,1,col='red')
#
##hist(ase2$gt[ase2$gt>0.001 & ase2$gt<0.999],breaks=1000)

png(file=paste(fileName,'.plots.qqplot.png',sep=''))

qq(pval[hetInd])
qqp <- qqplot(-log10(ppoints(sum(hetInd))),-log10(pval3),plot.it=F)
points(qqp,pch='.',cex=4,col='blue')

qqp <- qqplot(-log10(ppoints(sum(hetInd))),-log10(pval2[hetInd]),plot.it=F)
points(qqp,pch='.',cex=7,col='darkgreen')
legend("topleft",c(paste("D = Inf, eps =",round(eps*100,digits=2),"%"),
                   paste("Het 99%, D = ",round(Dmax,digits=2)),
                   paste("Het unc., pi0 = ",round(qv2$pi0*100,digits=2),"%")),fill=c("black","blue","darkgreen"))

dev.off()

sName
## THE-END
