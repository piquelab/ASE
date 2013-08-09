##################################################################
## Filters a pilup file from a single sample
## 
## Created: 08/06/2013
##
## Author: RPR, CH
##
## Version0.1: set up for a make file 
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
setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/bedfiles")
##################################################################
#pileupFile="HUVEC_EtOH_6hr.pileup.bed.gz"
#pileupFile="HUVEC_EtOH_24hr.pileup.bed.gz"
#pileupFile="HUVEC_Dex_6hr.pileup.bed.gz"
pileupFile="HUVEC_Dex_24hr.pileup.bed.gz"
##################################################################
thresh <- 13  # Base call quality threshold  NOT used
mincov <- 1   # Minimum coverage
maxcov <- 200000 # Max coverage
##outDir <- 'res/test'
##outDir <- gsub("pileup/","res/",gsub(".pileup.bed.gz","",pileupFile))

##################################################################
#cargs<-commandArgs(trail=TRUE);
#if(length(cargs)>=1)
#  pileupFile <-cargs[1];
#if(length(cargs)>=2)
#  outBedFile <-cargs[2];
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
#cat("#sName:",sName,"\n")
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
quality<- sapply(1:nrow(pileup),function(ii){qual(pileup$read.quality[ii])})
quality <- unlist(quality)
qtr <- quantile(quality,seq(0,1,0.1))
quality.table <- table(quality);
#qtr
#quality.table
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
outFile=gsub("pileup.bed.gz", "pileup.clean.bed.gz", pileupFile)
setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/cleanBedfiles")
#system(paste("mkdir -p",sName))
## Write the original data just in case we want to view later
write.table(pileup[,c(1:6,10,13:15)],file=gzfile(paste(outFile,sep='')),row.names=F,col.names=F,sep='\t',quote=F)
##write.table(file=paste(fileName,'.stats.txt',sep=''),saveData,row.names=T,col.names=T,quote=F,sep='\t') # Must be data frame

##cat("End:",sName,"\n")
## THE-END
