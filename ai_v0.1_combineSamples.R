##################################################################
## Combines pileups accross samples y first creating the union 
## of all loci, then combining reads.
##
## Created: 08/06/2013
##
## Author: RPR
##
## Version0.1: set up for a make file 
##################################################################
setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/cleanBedfiles")
##################################################################  
fileList <- dir("./"," *pileup.clean.bed.gz$")
#fileList <- dir("./",".*Dex.*pileup.clean.bed.gz$")
## grep -v grep non matching lines -w matching whole word 
##system("less *Dex*.pileup.clean.bed.gz | grep -v -w '^chr' | cut -f 1-3 | bedSort stdin stdout | uniq | gzip > union.bed.gz")
system("less *.pileup.clean.bed.gz | grep -v -w '^chr' | cut -f 1-3 | bedSort stdin stdout | uniq | gzip > HUVEC_union.bed.gz")

sNames = gsub(".pileup.clean.bed.gz","",fileList)
sNames
##################################################################  
##################################################################  
ExtractFields <- function(fn){			
	##fn <- fileList[1]
	cat("Processing:",fn,"\n")
	command=paste("intersectBed -a union.bed.gz -b ",fn," -wao | cut -f 12-14 ",sep="")
	aa <- read.table(pipe(command),sep="\t",as.is=T,na.strings=".")
	aa[is.na(aa)] <- 0
	aa			
}

aux <- sapply(fileList,ExtractFields)

Ref <- as.matrix(do.call(cbind,aux[1,]))
colnames(Ref) <- sNames
Alt <- as.matrix(do.call(cbind,aux[2,]))
colnames(Alt) <- sNames
Err <- as.matrix(do.call(cbind,aux[3,]))
colnames(Err) <- sNames

allRef<-apply(Ref, MARGIN=1, sum)
allAlt<-apply(Alt, MARGIN=1, sum)
allErr<-apply(Err, MARGIN=1, sum)

outData<-as.matrix(cbind(allRef, allAlt, allErr))
outFile<-"HUVEC_0_combined.gz"

write.table(outData,file=gzfile(paste(outFile,sep='')),row.names=F,col.names=F,sep='\t',quote=F)
