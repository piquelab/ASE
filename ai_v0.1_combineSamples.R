##################################################################
## Combines pileups accross samples by first creating the union 
## of all loci, then combining reads.
##
## Created: 08/06/2013
##
## Author: RPR
##
## Version0.1: use to determine if samples are from the same individual
##################################################################
##################################################################  
UnionExtractFields <- function(fileList, combine=FALSE){			
	#browser()
	tmpFile <- scan(pipe("mktemp -t"),character(0))
	system(paste("less ",fileList," | grep -v -w '^chr' | cut -f 1-3,7 | bedSort stdin stdout | uniq | gzip > ",tmpFile))
	sNames <- strsplit(gsub("/nfs/rprscratch/pipeline3/res/","",gsub(".pileup.clean.bed.gz","",fileList)),"/")
	sNames<-sapply(sNames, `[`, 1)	
	
	anno <- read.table(gzfile(tmpFile),sep="\t",as.is=T)	
	
	aux <- sapply(fileList,function(fn){			
				#fn <- fileList[1]
				cat("Processing:",fn,"\n")
				command=paste("intersectBed -a ",tmpFile," -b ",fn," -wao | cut -f 1-3,12-14 ",sep="")
				aa <- read.table(pipe(command),sep="\t",as.is=T,na.strings=".")
				aa[is.na(aa)] <- 0
				stopifnot(identical(aa[,1:3],anno[,1:3]))
				aa[,-(1:3)]
			})
	colnames(anno) = c("chr","pos0","pos","af")
	
	Ref <- as.matrix(do.call(cbind,aux[1,]))
	colnames(Ref) <- sNames
	Alt <- as.matrix(do.call(cbind,aux[2,]))
	colnames(Alt) <- sNames
	Err <- as.matrix(do.call(cbind,aux[3,]))
	colnames(Err) <- sNames
	
	return.list<-list(ref=Ref,alt=Alt,err=Alt,anno=anno);
	
	if(combine==TRUE){
		allRef<-apply(Ref, MARGIN=1, sum)
		allAlt<-apply(Alt, MARGIN=1, sum)
		allErr<-apply(Err, MARGIN=1, sum)
		
		return.list$all<-as.matrix(cbind(allRef, allAlt, allErr))
		
	}
	return(return.list)
}
##################################################################  
#ase.dat<-UnionExtractFields(pileups,combine=TRUE)
##################################################################  