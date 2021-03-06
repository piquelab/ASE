## Examination of allele-specific effects using ENCODE data and Illumina/impute2 genotyping data
## Gregory Moyerbrailean, Wayne State Unviersity 2013
# To look for SNPs causing allele-specific effects, we will compare DNase- and/or CHiP-seq data 
# (read data from samtools mpileup) to SNPs, both genotyped and predicted by impute2.
# In this way, we can assess, for a given heterozygous SNP position whether or not we see 
# an allele-specific effect, possibly indicated by an imbalance of reads for the two alleles

## Version 3.4
## Roger is making some changes
## 1) as.is=T when opening file
## 2) added some plotting functions
## 3) Removed some steps and plots. 
## 4) Added new method for genotyping and to make ASE calls


## New in version 3.3
# Elminated reading of SNP allels (handled elsewhere in pipeline)
# Added steps for filtering out the first and last bases from reads (WORKING)

## New in version 3.2
# Added minimum coverage constraint of 20x to heterozygous SNPs
# Added maximum coverage constraint of 100x to filter possible unannotated repeats
# Added ability to read in alternate alleles
# Added density plots for visualization of results

## New in version 3.1
# q-value calculation added
# Significance now defined as q-value of less than 0.1
# Graph is now colored by q-value significance
# Graph now plots log10((ref+1)/(alt+1)) to allow plotting of "0" values
# Het data output is now true bed format, i.e. no header and '.bed' suffix

library(qvalue)

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

## Allele specific funtions. 
source('/nfs/hpnfs/groups/piquelab/gmb/AI/results/fitAseModels.v2.R')

## Some ploting functions
source('/nfs/hpnfs/groups/piquelab/gmb/AI/results/qqman.r')

# pileupFile="../pileups/huvec.wg.cmycPeaks.impute.pileup"
pileupFile="pileup/UwDgfAoaf.pileup.bed.gz"
pileupFile="pileup/UwDgfHuvec.pileup.bed.gz"
pileupFile="pileup/OpenChromDnaseLncapAndro.pileup.bed.gz"
##pileupFile="HUVEC_Dex_6hr.pileup.bed.gz"
# pileupFile="pileup/OpenChromDnaseMixGm10266Gm12878.pileup.bed.gz"
# setwd("/wsu/tmp/pipeline2/")
# pileupFile="pileup/UwDgfHuvec.pileup.bed.gz"
# OpenChromDnaseLncapAndro

## -- 
## setwd("~/piquelab/130702_Data/analyses/pileup")

##RPR: Moving some options up so it is easier to change them later
thresh <- 13  # Base call quality threshold  NOT used
mincov <- 4   # Minimum coverage
maxcov <- 20000 # Max coverage

outDir <- 'res/test'
outDir <- gsub("pileup/","res/",gsub(".pileup.bed.gz","",pileupFile))

cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
  pileupFile<-cargs[1];
if(length(cargs)>=2)
  outDir<-cargs[2];
if(length(cargs)>=3)
  thresh<-as.integer(cargs[3]);

## Isolate the file name for the output files
fileName = unlist(strsplit(pileupFile,'/'))
fileName = fileName[length(fileName)] # remove the path to the pileup
fileName = paste(outDir,fileName,sep='/')

## Output to make sure everything is as expected
cat("#Pileup File:",pileupFile,"\n")
cat("#Threshold:",thresh,"\n")
cat("#fileName:",fileName,"\n")
sName <- gsub("pileup/","",gsub(".pileup.bed.gz","",pileupFile));
cat("#sName:",sName,"\n")

## Create an object to save data in
saveData <- data.frame(pileup.file=c(pileupFile)) #,impute.file=c(imputeFile))

## QC STEPS
# Remove all loci whose read coverage is below 4 (our picked threshold)
# cov = read coverage per allele to define heterozygote
# RPR: We could pipe the data in and apply these filters, to make all this go faster
## pileup <- pileup[pileup$num.reads >= mincov,]
## pileup <- pileup[pileup$num.reads <= maxcov,]

command=paste("less ",pileupFile,
  "| awk ' $5 >=",mincov," && $5 <=",maxcov,"'")
cat("#pipe: ",command,"\n")

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
sum(pileup$errors)
eps0 <- sum(pileup$errors)/sum(pileup$ref.matches+pileup$alt.matches)
eps0

## Reorder by position so chr names appear right

chrList <- unique(pileup$chr)
chrList <- paste("chr",sort(as.numeric(gsub("chr","",chrList))),sep="")
factor(chrList,levels=chrList)
pileup$chr <- factor(pileup$chr,levels=chrList)
pileup <- pileup[order(pileup$chr,pileup$pos),];


#####################################################################
system(paste("mkdir -p",sName))
## Write the original data just in case we want to view later
write.table(pileup,file=gzfile(paste(fileName,'.loci.gz',sep='')),row.names=F,col.names=T,sep='\t',quote=F)
##write.table(file=paste(fileName,'.stats.txt',sep=''),saveData,row.names=T,col.names=T,quote=F,sep='\t') # Must be data frame
##saveData$mean.coverage <- mean(pileup$num.reads)

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

##hist(pval2[hetInd & pval2<0.9],breaks=100)
##hist(pval3,breaks=1000)

##qplot(qv2)


####################################


##hist(pval2[pval2<0.99],breaks=200)

## names(ind) <- pileup$rsID
## names(pval2) <- names(ind[ind])
## aux <- qv2$qv<0.1;
## aux <- names(pval2)[aux]
## qq(pval2,annotate=aux)

##qq(pval[hetInd])


## ## Plot the starting density of ref alleles / total
## pdf(file=paste(fileName,'.plots.startDens.pdf',sep=''))
## plot(density(pileup[!(pileup$ref.matches==0 & pileup$alt.matches==0),]$ref.matches/(pileup[!(pileup$ref.matches==0 & pileup$alt.matches==0),]$ref.matches+pileup[!(pileup$ref.matches==0 & pileup$alt.matches==0),]$alt.matches)),main="Ref/Total - All Loci")
## abline(v=0.5);
## title(sName)
## dev.off()



## tr.otac <- 20
## ##b <- c(0:21-0.5)/20;
## b <- c(0:tr.otac)/tr.otac;
## oraf <- pileup$ref.matches/(pileup$ref.matches + pileup$alt.matches)
## omac <- pmin(pileup$alt.matches,pileup$ref.matches)
## otac <- (pileup$ref.matches + pileup$alt.matches)

## bHist <- hist(oraf[omac>=2 & otac>tr.otac],breaks=b,plot=F);
## tr.omac <- c(2,3,4,5,6,8,10,12,15,20,25)
## aux <- sapply(tr.omac,function(tr1){
##   bHist <- hist(oraf[omac>=tr1 & otac>=tr.otac],breaks=b,plot=F);
##   bHist$counts;#/sum(bHist$counts);
## })

## cpal <- colorRampPalette(c("purple","darkblue","blue","green","orange","red","darkred"))
## mcol <- cpal(length(tr.omac))

## pdf(file=paste(fileName,'.plots.allHist.pdf',sep=''))

## plot(NA,xlim=c(0,1),ylim=c(0,max(aux)),col='lightblue',main=sName,xlab="Observed reference allele freq.","Frequency");
## sapply(1:ncol(aux),function(ii){lines(bHist$mids,aux[,ii],lwd=3,col=mcol[ii]);ii})
## abline(v=0.5,lwd=4,col='darkblue',lty=2)
## abline(v=c(0.25,0.75),lty=2)
## legend("topleft",paste("",tr.omac),lwd=3,col=mcol,bty="n",title="Het. tresh >=")

## dev.off()


##cat("#bHist:",sName,bHist$counts,"\n",sep="\t")

## ## -- Filter for the predicted hets using the coverage threshold --
## hetcov <- 2
## pileup.hets <- pileup[pileup$ref.matches >= hetcov & pileup$alt.matches >= hetcov,]
## pileup.hets <- pileup.hets[pileup.hets$num.reads >= 20,]
## pileup.hets <- pileup.hets[order(pileup.hets$chr,pileup.hets$pos),]
## pileup.hets$prop <- pileup.hets$ref.matches/(pileup.hets$ref.matches + pileup.hets$alt.matches)

## write.table(file=gzfile(paste(fileName,'.hets.bed.gz',sep='')),pileup.hets,row.names=F,col.names=F,quote=F,sep='\t')

## if (dim(pileup.hets)[1] > 2) {
## 	# Plot the density & histogram of ref/total after filtering for (predicted) hets
## 	pdf(file=paste(fileName,'.plots.hetDens.pdf',sep=''))
## 	plot(density(pileup.hets$prop),main=sName)
## 	abline(v=0.5)
## 	dev.off()
	
## 	pdf(file=paste(fileName,'.plots.hetHist.pdf',sep=''))

## 	b <- c(0:51)/51
## 	hetHist <- hist(pileup.hets$prop,breaks=b,plot=F);
##         plot(hetHist,xlim=c(0,1),col='lightblue',main=sName,xlab="Ref/Total");
## 	abline(v=0.5)
##         abline(v=c(0.25,0.5,0.75))
        
## 	dev.off()
##         cat("#hetHist:",sName,hetHist$counts,"\n",sep="\t")
	

## 	# Run a binomial test on the heterozygous sites
## 	sig <- 0.1	# 0.001 for Huvec DgfDNase data
## 	pileup.hets$binom <- apply(data.frame(a=pileup.hets$ref.matches,b=pileup.hets$alt.matches), 1, function(x) binom.test(c(x[1],x[2]),alternative="two.sided")$p.value)
## 	pileup.hets$qvalue <- qvalue(pileup.hets$binom)$qvalues
## 	pileup.hets$sig <- unlist(lapply(pileup.hets$qvalue, function(x,sig) ifelse(x < sig, 1, 0), sig))

##         x <- 1:length(pileup.hets$binom)
##         x <- x/max(x);
##         #qqplot(-log10(x),-log10(pileup.hets$binom))
        
## 	write.table(file=paste(fileName,'.hets.sig.bed',sep=''),pileup.hets[pileup.hets$sig == 1,],row.names=F,col.names=F,quote=F,sep='\t')
## 	out <- c(fileName,dim(pileup)[1],dim(pileup.hets)[1],dim(pileup.hets[pileup.hets$sig==1 & pileup.hets$prop > .5,])[1],dim(pileup.hets[pileup.hets$sig==1 & pileup.hets$prop < .5,])[1])
## 	write.table(file=paste(fileName,'.counts.txt',sep=''),t(out),row.names=F,col.names=F,quote=F,sep='\t')
	
## 	# Plot the density of ref/total for significantly imbalanced hets
## 	if (sum(pileup.hets$sig) > 1) {
## 		pdf(file=paste(fileName,'.plots.sigDens.pdf',sep=''))
                
## 		plot(density(pileup.hets[pileup.hets$sig==1,]$prop),main="Ref/Total - Imbalanced Hets")
## 		abline(v=0.5)
                
## 		dev.off()
		
## 		pdf(file=paste(fileName,'.plots.sigHist.pdf',sep=''))
                
## 		plot(hist(pileup.hets[pileup.hets$sig==1,]$prop,breaks=b,xlim=c(0,1),col='lightblue',main="Predicted Hets",xlab="Ref/Total"))
## 		abline(v=0.5)
                
## 		dev.off()
## 	}

## 	## Save a plot showing the distribution of read alleles
## 	pdf(file=paste(fileName,'.plots.mirMan.pdf',sep=''))
        
## 	plot(1:length(pileup.hets$ref.matches), (pileup.hets$ref.matches) / (pileup.hets$ref.matches + pileup.hets$alt.matches), 
## 		 pch=19, col=pileup.hets$sig+1, xlab="Position", ylab="Reference Reads / Total Reads", 
## 		 main="Read Alleles For Predicted Heterozygous Loci", ylim=c(0,1))
## 	#points(1:length(pileup.hets[pileup.hets$sig==1,]$sig), pileup.hets[pileup.hets$sig==1,]$prop, pch=19,col='red')
## 	abline(h=0.5,lwd=3)
## 	dev.off()
## }
## write.table(file=paste(fileName,'.stats.txt',sep=''),saveData,row.names=T,col.names=T,quote=F,sep='\t') # Must be data frame

sName
## THE-END
