##################################################################
## Filters a pilup file from a single sample
## 
## Created: 08/06/2013
## 
## Author: RPR, CH
##
## Version0.0: hacked from 'ai_v3.6_noImpute' to be run interactively
##################################################################
qual <- function(x) { strtoi(charToRaw(x),16L)-33 }

filt <- function(r,q,t) {
	ret=''
	for (i in 1:nchar(r)) {
		if (substr(r,i,i) != '$' & substr(r,i,i) != '^') {
			if (qual(substr(q,i,i)) >= t) { ret <- paste(ret,substr(r,i,i),sep='') }
		}
	}
	ret
}
##################################################################
##setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/bedfiles")
##################################################################
#pileupFile="HUVEC_EtOH_6hr.pileup.bed.gz"
#pileupFile="HUVEC_EtOH_24hr.pileup.bed.gz"
#pileupFile="HUVEC_Dex_6hr.pileup.bed.gz"
#pileupFile="HUVEC_Dex_24hr.pileup.bed.gz"
##################################################################
thresh <- 13  # Base call quality threshold  NOT used
mincov <- 4   # Minimum coverage
maxcov <- 200000 # Max coverage
##outDir <- 'res/test'
##outDir <- gsub("pileup/","res/",gsub(".pileup.bed.gz","",pileupFile))

##################################################################
cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
  pileupFile <-cargs[1];
if(length(cargs)>=2)
  outDir <-cargs[2];
#if(length(cargs)>=3)
#  outLociFile <- cargs[3];
##################################################################

## Isolate the file name for the output files
##fileName = unlist(strsplit(pileupFile,'/'))
##fileName = fileName[length(fileName)] # remove the path to the pileup
##fileName = paste(outDir,fileName,sep='/')
##################################################################
## Output to make sure everything is as expected
#cat("#Pileup File:",pileupFile,"\n")
#cat("#Threshold:",thresh,"\n")
###cat("#fileName:",fileName,"\n")
sName <- gsub("pileup/","",gsub(".pileup.bed.gz","",pileupFile));
cat("#sName:",sName,"\n")
cat("outDir:",outDir,"\n")
##################################################################
## Create an object to save data in
saveData <- data.frame(pileup.file=c(pileupFile)) #,impute.file=c(imputeFile))
##################################################################
## QC STEPS
# Remove all loci whose read coverage is below 4 (our picked threshold)
# cov = read coverage per allele to define heterozygote
# RPR: We could pipe the data in and apply these filters, to make all this go faster
## pileup <- pileup[pileup$num.reads >= mincov,]
## pileup <- pileup[pileup$num.reads <= maxcov,]

command=paste("less ",pileupFile,
		"| awk ' $5 >=",mincov," && $5 <=",maxcov,"'")
#cat("#pipe: ",command,"\n")

## RPR: pipe speeds up, as.is=T speeds up and avouids having to use character conversion later...
pileup <- read.table(file=pipe(command),header=F,quote="",comment.char="",as.is=T,sep="\t") 
##names(pileup) <- c("chr","pos-1","pos","ref","num.reads","read.alleles","read.quality","xchr","xpos-1","xpos","rsID","TKG.Ref","alt","xcov","xoverlap")
names(pileup) <- c("chr","pos-1","pos","ref","num.reads","read.alleles","read.quality","xchr","xpos-1","xpos","rsID","TKG.Ref","alt","d1","af","d2")


# See if the ref allels match, then discard uncesessary columns
##saveData$ref.match <- length(which(toupper(as.character(pileup$ref)) == as.character(pileup$TKG.Ref))) == dim(pileup)[1]
indMatch <- (toupper(pileup$ref) == pileup$TKG.Ref)
saveData$ref.match <- sum(!indMatch)==0;
saveData$ref.missmatches <- sum(!indMatch);
pileup <- pileup[indMatch,c("chr","pos-1","pos","ref","alt","rsID","num.reads","read.alleles","read.quality","af")]
## I'm also discarding mismatching rows.
stopifnot(mean(indMatch)>0.8) ##RPR stop if too many errors?
rm(indMatch)

## ## Plot coverage histogram
## hist(log2(pileup$num.reads),breaks=200)

# Duplicates arise from the impute data for 3 (so far) determined reasons: tri+ alleleic SNPs,
# indels (should already be filtered), and incongruencies between illumina and 1KG data
d1 <- duplicated(paste(pileup$chr,pileup$pos,sep=":"))
d2 <- duplicated(paste(pileup$chr,pileup$pos,sep=":"),fromLast=T)
saveData$num.duplicates <- length(which(d1)) + length(which(d2))
pileup <- pileup[!(d1 | d2),]
rm(d1,d2)
## RPR: Position did not inlude chr

## Clean up the read data
## TODO - after proper filtering, this first line shoudlnt matter anymore. Confirm.

## Filter the alleles to remove those at the beginning and end of a mapped read
## Unnecessary if data loaded with: as.is=T
## pileup$read.alleles <- as.character(pileup$read.alleles)
## pileup$read.quality <- as.character(pileup$read.quality)
pileup$read.alleles.filt <- mapply(gsub,'[a-zA-Z.,]\\$','$',pileup$read.alleles) ## Remove start base 
pileup$read.alleles.filt <- mapply(gsub,'\\^[[:punct:][:alnum:]][a-zA-Z.,]','^',pileup$read.alleles.filt)  ## Remove end base

##head((strsplit(pileup$read.quality,'')))
## RPR alternative way to extract base quality numbers. 
qual <- sapply(1:nrow(pileup),function(ii){qual(pileup$read.quality[ii])})
qual <- unlist(qual)
qtr <- quantile(qual,seq(0,1,0.1))
qual.table <- table(qual);
qtr
qual.table
##RPR: We could decide the the threshold based on the quantile...?
##RPR: For the moment I remove the function...



##pileup$read.alleles.filt <- mapply(filt,pileup$read.alleles.filt,pileup$read.quality,MoreArgs=list(t=thresh))
pileup <- pileup[nchar(pileup$read.alleles.filt)>0,]
## RPR: we could do some of the filters as a pipe. 

## Should be no more odd chars - no Ns, no +s, -s, ^s, $s, etc (still... assert that there are only ACGTs left?)
##pileup$read.alleles.clean <- gsub('[^AGCTagct.,]',"",pileup$read.alleles)
pileup$read.alleles.clean <- mapply(gsub,'[\\.\\,]',pileup$ref,pileup$read.alleles.filt)
pileup$read.alleles.clean <- toupper(pileup$read.alleles.clean)
pileup$ref <- toupper(pileup$ref) # for downstream analysis

# For each read at each location, calculate the number of ref and alt allele matches (alt defined by 1KG)
pileup$ref.matches <- as.integer(nchar(mapply(gsub,paste('[^',as.character(pileup$ref),']',sep=""),'',pileup$read.alleles.clean)))
pileup$alt.matches <- as.integer(nchar(mapply(gsub,paste('[^',as.character(pileup$alt),']',sep=""),'',pileup$read.alleles.clean)))

## RPR: I'm trying to see if we can calculate errors by counting the alleles not matching the ref,
## but it is not possible if the mapper does not allow errors. We could hash all alternate bases at the SNP location, this would allow
## to estimate the error rate, by counting the number of bases not matching REF/ALT alleles.
## This is relatively easy to do and could help to estimate accuracy in determining Hets.
pileup$errors <- as.integer(nchar(mapply(gsub,paste('[^ACGT]',sep=""),'',pileup$read.alleles.clean))) - (pileup$alt.matches + pileup$ref.matches)
#sum(pileup$errors)
eps0 <- sum(pileup$errors)/sum(pileup$ref.matches+pileup$alt.matches)
#eps0

## Reorder by position so chr names appear right

chrList <- unique(pileup$chr)
chrList <- paste("chr",sort(as.numeric(gsub("chr","",chrList))),sep="")
factor(chrList,levels=chrList)
pileup$chr <- factor(pileup$chr,levels=chrList)
pileup <- pileup[order(pileup$chr,pileup$pos),];

#####################################################################
##
##outFile=gsub("pileup.bed.gz", "pileup.clean.bed.gz", pileupFile)
system(paste("mkdir -p",outDir))
fileName <- paste(outDir,"/",sName,"",sep="");
outFile <- paste(fileName,".pileup.clean.bed.gz",sep="");
## Write the original data just in case we want to view later
write.table(pileup[,c(1:6,10,13:15)],file=gzfile(paste(outFile,sep='')),row.names=F,col.names=F,sep='\t',quote=F)
##write.table(file=paste(fileName,'.stats.txt',sep=''),saveData,row.names=T,col.names=T,quote=F,sep='\t') # Must be data frame

#####################################################################
#####################################################################
#####################################################################
#####################################################################

##source("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/src/fitAseModels.v2.R")
source("/nfs/hpnfs/groups/piquelab/charvey/ASE/src/fitAseModels.v2.R")
source('/nfs/hpnfs/groups/piquelab/gmb/AI/results/qqman.r')
library(qvalue)

####################################
oraf <- pileup$ref.matches/(pileup$ref.matches + pileup$alt.matches)
omac <- pmin(pileup$alt.matches,pileup$ref.matches)
otac <- (pileup$ref.matches + pileup$alt.matches)

tr.l <- 15
tr.u <- 400
ind <- otac>tr.l & otac<tr.u;

names(ind) <- pileup$rsID

ref=pileup$ref.matches[ind];
alt=pileup$alt.matches[ind];
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
cat("#Dmax:\t",Dmax,"\n");


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










cat("End:",sName,"\n")
## THE-END
