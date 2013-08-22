##################################################################
## ASE model for use on multiple samples 
##
## Created: 08/18/2013
##
## Author: CTH
##
## Version0.2: Use the bayes factor
##################################################################
## Data aggregation functions
source('/nfs/hpnfs/groups/piquelab/charvey/ASE/src/ai_v0.1_combineSamples.R')
## ASE model fitting functions funtions 
source('/nfs/hpnfs/groups/piquelab/charvey/ASE/src/ai_v0.0_fitAseModels.R')
## Some ploting functions
source('/nfs/hpnfs/groups/piquelab/gmb/AI/results/qqman.r')
#setwd("/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/cleanBedfiles")
#pileups <- dir("./"," *pileup.clean.bed.gz$")
##################################################################
setwd("/nfs/rprscratch/pipeline3/res")
##################################################################
pileups <- as.character(scan("selectByDmax.txt", what = character(), strip.white=TRUE))
pileups<-paste("/nfs/rprscratch/pipeline3/res/",pileups,"/",pileups,".pileup.clean.bed.gz",sep="")
ase.dat<-UnionExtractFields(pileups)
##################################################################  
ref<-ase.dat[[1]]
alt<-ase.dat[[2]]
err<-ase.dat[[3]]
sample.ids<-colnames(ase.dat[[1]])
##################################################################
GetId<-function(sample.list,sample.number){sample.list[sample.number]}
GetSampleBF<-function(bf.matrix,sample.string){bf.matrix[sample.string,]}
##################################################################
n.samples<-dim(ref)[2]

tr<-12
indMat <- (ref+alt)>tr
selRows<- rowSums(indMat)>=70;
ref <- ref[selRows,]
alt <- alt[selRows,]
err <- err[selRows,]
indMat <- indMat[selRows,]

rhoall <- ref/(alt+ref)
simmat <- cor(rhoall)

heatmap(1-abs(simmat))

aux <- t(sapply(1:(n.samples-1),function(ii){
					cat("Processing ii:",ii,"\n")
					aux <- rep(0,n.samples)
					aux[(ii+1):n.samples] <- sapply((ii+1):n.samples,function(jj){
								#ii<-22; jj<-124
	
								indVec <- indMat[,ii] & indMat[,jj]
								##sum(indVec)
								
								#row.1<-ref[indVec,ii]*(ref[indVec,ii]+alt[indVec,ii])^-1
								#row.2<-ref[indVec,jj]*(ref[indVec,jj]+alt[indVec,jj])^-1
								
								#cor(row.1,row.2)
								all.ref <- cbind(ref[indVec,ii],ref[indVec,jj])
								all.alt <- cbind(alt[indVec,ii],alt[indVec,jj])
								
																
								llkMulti <- fitAseNullMulti(all.ref,all.alt,verbose=FALSE,fixGprior=FALSE)$logliksum
								llk.1 <- fitAseNull(ref[indVec,ii],alt[indVec,ii],verbose=FALSE,fixGprior=FALSE)$logliksum
								llk.2 <- fitAseNull(ref[indVec,jj],alt[indVec,jj],verbose=FALSE,fixGprior=FALSE)$logliksum		
								
								bf <- (llkMulti-llk.1-llk.2)
								#bf <- exp(llk.1+llk.2-llkMulti)
								
							})
					aux
				}))


##################################################################  
setwd("/nfs/hpnfs/groups/piquelab/charvey/ASE/out")
image.directory <- '/wsu/home/fl/fl97/fl9788/piquelab/charvey/ASE/out/'
load("aux_08202013")
load("ase.dat_08212013")
##################################################################  
#aux <- rbind(aux,rep(0,n.samples))
#dim(aux)
#aux <- aux + t(aux)
#diag(aux)<-0

png(file=paste(image.directory,'bayesFactor.complete.png',sep=''))

avgcov <- round(colMeans(ref+alt),digits=2)

colnames(aux) <- paste("a",avgcov,colnames(ref),sep=":")
rownames(aux) <- colnames(aux)
image(aux)

dev.off()
##################################################################  
# range of aux -14762,92 to 13358



vmax <- apply(aux,1,max)
vmin <- apply(aux,1,min)

plot(density(vmax));
aux2 <- aux[vmax>600,vmax>600];

aux2[aux2<1] <- 1

dmat <- 1/aux2

dmat[dmat>70000] <- 70000
diag(dmat) <- 0

png(file=paste(image.directory,'bayesFactor.bfs.png',sep=''),width=2000,height=2000)
heatmap(dmat,margins=c(40,40),cexCol=2,cexRow=2)
#heatmap(dmat,margins=c(40,40),cexCol=2,cexRow=2)
dev.off()

plot(GetSampleBF(dmat,"EpiUwRmapDNaseFetal.Intestine.Small.SRX201835"))
##diag(aux2) <-1
##heatmap(1-aux2)

## need to identify what is causing come stronf, and some possiby spurious associations,
## could be differences in informedness of sites
## 1: for starters, lets compare 4363 to 1835
comp.strong<-cbind(ref[,"EpiUwRmapDNaseFetal.Intestine.Large.SRX204363"],ref[,"EpiUwRmapDNaseFetal.Intestine.Small.SRX201835"])
indVec1<-indMat[,"EpiUwRmapDNaseFetal.Intestine.Large.SRX204363"] & indMat[,"EpiUwRmapDNaseFetal.Intestine.Small.SRX201835"]
comp.strong<-comp.strong[indVec1,]

comp.weak<-cbind(ref[,"EpiUwRmapDNaseFetal.Muscle.Trunk.SRX214045"],ref[,"EpiUwRmapDNaseFetal.Lung.Right.SRX089270"])
indVec2<-indMat[,"EpiUwRmapDNaseFetal.Muscle.Trunk.SRX214045"]&indMat[,"EpiUwRmapDNaseFetal.Lung.Right.SRX089270"]
comp.weak<-comp.weak[indVec2,]

##aux <- do.call(rbind,aux)
aux[1:10,1:10]
str(aux)

##################################################################  
