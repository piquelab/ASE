##################################################################
## Added model fitting for multiple samples with individual epsilon 
##
## Created: 08/09/2013
##
## Author: RPR, CH
##
## Version0.0: to be run interactively
##################################################################
## 
## Updated from /nfs/hpnfs/groups/piquelab/gmb/AI/results/fitAseModels.R
##

logLikRhoEps <- function(rho,eps,ref,alt){
  ref * log((1 - eps) * rho + eps * (1 - rho)) + alt * log((1 - eps) * (1 - rho) + eps * rho)
}

## Beta function maximization.
logLikBetaBinomial <- function(R,A,p,D){
  aux <- (lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D) + lgamma(D)) ##+ lgamma(R+A+1) - lgamma(A+1) - lgamma(R+1)
  aux
}

logLikBetaBinomialRhoEps <- function(rho,eps,D,R,A){
  p <- (rho*(1-eps)+(1-rho)*eps)
  aux <- (lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D) + lgamma(D)) ##+ lgamma(R+A+1) - lgamma(A+1) - lgamma(R+1)
  aux
}

fnLogLikBetaBinomialNull <- function(log.D,R,A){
  D <- exp(log.D);
  pD <- D*0.5;
  aux <- (lgamma(R+pD) + lgamma(A + pD) - lgamma(R+A+D) - 2*lgamma(pD) + lgamma(D));
  -sum(aux)
}
grLogLikBetaBinomialNull <- function(log.D,R,A){
  D <- exp(log.D);
  pD <- 0.5*D;
  aux <- (0.5* digamma(R+pD) + 0.5 * digamma(A + pD) - digamma(R+A+D) - digamma(pD) + digamma(D));
  -sum(aux)
}

logLikBetaBinomial2 <- function(logit.p,D,R,A){
  p <- plogis(logit.p)
  aux <- (lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D) + lgamma(D));
  -sum(aux)
}

gLogLikBetaBinomial <- function(logit.p,D,R,A){
  p <- plogis(logit.p)
  aux <- D * (digamma(R+p*D) - digamma(A + (1-p)*D) - digamma(p*D) + digamma((1-p)*D)) ;
  as.matrix(-aux)
}

comp.rho <- function(ref,alt,eps){
  rho <- (ref * (1 - eps) - eps * alt) / ((ref + alt) * (1 - 2 * eps))
  rho[rho<0] <- 0
  rho[rho>1] <- 1
  rho
}

lrtEpsRhoBinom <- function(ref,alt,eps,rho){
  log.gt <- cbind(g0   = ref * log(1 - eps) + alt * log(eps),
                  g1t0 = (ref + alt) * log(0.5),
                  g1t1 = ref * log((1 - eps) * rho + eps * (1 - rho)) + alt * log((1 - eps) * (1 - rho) + eps * rho),
                  g2   = ref * log(eps) + alt * log(1 - eps)
                  )
  lrt <- log.gt[,"g1t1"]-pmax(log.gt[,"g0"],log.gt[,"g1t0"],log.gt[,"g2"])
}


## aux <- optim(0,fn=fnLogLikBetaBinomialNull,gr=grLogLikBetaBinomialNull,R=ref,A=alt,method="BFGS")
## D <- exp(aux$par)

## aux <- optim(0.2,fn=logLikBetaBinomial2,gr=gLogLikBetaBinomial,D=0.1,R=100,A=50,method="BFGS",control=list(reltol=1E-10))
## plogis(aux$par)

## aux <- optim(0.2,fn=logLikBetaBinomial2,gr=gLogLikBetaBinomial,D=1,R=73,A=1,method="BFGS")

## L <- length(ref);

## aux <- optim(rep(0,L),fn=logLikBetaBinomial2,gr=gLogLikBetaBinomial,D=0.5,R=ref,A=alt,method="BFGS")
## rho3 <- plogis(aux$par)

## pini <- rep(0,L)
## logLikBetaBinomial2(pini,D=1,R=ref,A=alt)

fitAse <- function(ref,alt,eps=0.1,rho,t1=1E-100,log.gmat,max.it=100,tol=1E-16,fixEps=FALSE,fixRho=FALSE,fixT1=FALSE){
  L <- length(ref);
  stopifnot(L == length(alt));
  ## Parameter initialization
  ## eps=0.1;
  if(missing(rho))
    rho=rep(0.5,L);
  logit.rho <- qlogis(rho);
  ##log.gmat <- matrix(log(1/3),L,3)
  if(missing(log.gmat)){
    log.gmat <- matrix(log(c(0.25,0.5,0.25)),1,3)
    colnames(log.gmat) <- c("g0","g1","g2");
  }
  logit.t1 <- qlogis(t1);
  it.num <- 0
  logit.eps <- qlogis(0.1);
  converged <- FALSE;
  while(it.num <= max.it){
############ E-step #############
    log.gt <- cbind(g0   = ref * log(1 - eps) + alt * log(eps) + log.gmat[,"g0"],
                    g1t0 = (ref + alt) * log(0.5) + log.gmat[,"g1"] + log(1-t1),
                    g1t1 = ref * log((1 - eps) * rho + eps * (1 - rho)) + alt * log((1 - eps) * (1 - rho) + eps * rho) + log.gmat[,"g1"] + log(t1),
                    g2   = ref * log(eps) + alt * log(1 - eps) + log.gmat[,"g2"]
                    )
    ## This normalizes marginal posterior probabilities to add 1
    ##aux <- log.gt %*% rep(1/4,4)
    ##log.gt <- apply(log.gt,2,function(col){col-aux})
    ## the previous step not really necessary but may help stabilize large numbers, IT DOES NOT WORK
    aux <- log(exp(log.gt) %*% rep(1,4));
    log.gt <- apply(log.gt,2,function(col){col-aux})
    gt <- exp(log.gt);
    if(it.num == max.it)
      break;
    ## rowSums(head(exp(aux2))) ## should be ones.
    ## Constructing the W matrix
    ## 
    log.W <- cbind(w1  = log.gt[,"g0"]   ,
                   w2  = log.gt[,"g0"]   ,
                   w4  = log.gt[,"g1t0"] + log(1-eps),
                   w6  = log.gt[,"g1t0"] + log(eps),
                   w8  = log.gt[,"g1t1"] + log(1-eps) + log(1-rho) - log((1 - eps) * (1 - rho) + eps * rho),
                   w9  = log.gt[,"g1t1"] + log(1-eps) + log(rho)   - log((1 - eps) * rho + eps * (1 - rho)),
                   w10 = log.gt[,"g1t1"] + log(eps)   + log(rho)   - log((1 - eps) * (1 - rho) + eps * rho),
                   w11 = log.gt[,"g1t1"] + log(eps)   + log(1-rho) - log((1 - eps) * rho + eps * (1 - rho)),               
                   w12 = log.gt[,"g2"]   ,
                   w15 = log.gt[,"g2"]
                   )
    W <- (exp(log.W))
    if((it.num == max.it) | (converged==TRUE))
      break;
###### M-STEP  #####
    converged <- TRUE;
    ## w6=w7 and w5=w4 , w2=w1, 
    ## epsilon
    if(!fixEps){
      aux.num <- (sum(rowSums(W[,c("w6","w11","w15")]) * ref) + sum(rowSums(W[,c("w2","w6","w10")]) * alt));
      aux.den <- (sum(rowSums(W[,c("w1","w4","w9")]) * ref) + sum(rowSums(W[,c("w4","w8","w12")]) * alt));      
      new.logit.eps <- log(aux.num)-log(aux.den);
      if(abs(logit.eps-new.logit.eps)>tol)
        converged <- FALSE;
      logit.eps <- new.logit.eps;
      eps <- plogis(logit.eps)
    }
    ## rho vector
    if(!fixRho){
      new.logit.rho <- log(W[,"w9"] * ref + W[,"w10"] * alt + 0) - log(W[,"w11"] * ref + W[,"w8"] * alt + 0);
      new.logit.rho[new.logit.rho>100] <- 100
      new.logit.rho[new.logit.rho< -100] <- -100      
      rho.change <-  max(abs(new.logit.rho-logit.rho));
      cat("rho.change",rho.change,rho.change > tol,"\n");
      if(rho.change > tol){
        converged <- FALSE
      }
      logit.rho <- new.logit.rho;
      rho <- plogis(logit.rho)
    }
    ## gmat
    ## aux <- colSums(gt);
    ## aux <- aux/sum(aux);
    ## log.gmat[,"g0"] <- log(aux["g0"])
    ## log.gmat[,"g1"] <- log(aux["g1t0"]+aux["g1t1"])
    ## log.gmat[,"g2"] <- log(aux["g2"])
    ## ## opt 2
    ## log.gmat[,"g0"] <- log.gt[,"g0"]
    ## log.gmat[,"g1"] <- log(gt[,"g1t0"]+gt[,"g1t1"]) 
    ## log.gmat[,"g2"] <- log.gt[,"g2"]
    ## # t1
    if(!fixT1){
      new.logit.t1 <- log(sum(gt[,"g1t1"])) - log(sum(gt[,"g1t0"]))
      if(abs(logit.t1-new.logit.t1)>tol)
        converged <- FALSE;
      logit.t1 <- new.logit.t1;
      t1 <- plogis(logit.t1)
    }
    ##
    it.num <- it.num +1;
    cat("#it:",it.num,"eps=",eps,"t1=",t1,"Post:",colMeans(gt),"rhoS:",summary(rho),"\n");
    stopifnot(!is.na(eps))
  }
  invisible(list(gt=gt,rho=rho,eps=eps))
}


fitAseNull <- function(ref,alt,eps=0.1,log.gmat,max.it=100,tol=1E-16,fixGprior=TRUE){
  L <- length(ref);
  stopifnot(L == length(alt));
  ## Parameter initialization
  ##log.gmat <- matrix(log(1/3),L,3)
  if(missing(log.gmat)){
    log.gmat <- matrix(log(c(0.25,0.5,0.25)),1,3)
    colnames(log.gmat) <- c("g0","g1","g2");
  }
  it.num <- 0
  logit.eps <- qlogis(0.1);
  converged <- FALSE;
  logliksum <- 0;
  while(it.num <= max.it){
  ############ E-step #############
    log.gt <- cbind(g0   = ref * log(1 - eps) + alt * log(eps) + log.gmat[,"g0"],
                    g1 = (ref + alt) * log(0.5) + log.gmat[,"g1"],
                    g2   = ref * log(eps) + alt * log(1 - eps) + log.gmat[,"g2"]
                    )
    log.gt[log.gt < (-200)] <- (-200)
    ##browser()
    ## This normalizes marginal posterior probabilities to add 1
    loglik <- log(exp(log.gt) %*% rep(1,3));
    new.logliksum <- sum(loglik);
    if(abs(new.logliksum-logliksum)<tol)
      converged <- TRUE;
    log.gt <- apply(log.gt,2,function(col){col-loglik})
    gt <- exp(log.gt);
    if((it.num == max.it) | (converged==TRUE))
      break;
    ###### M-STEP  #####
    ## w6=w7 and w5=w4 , w2=w1, 
    ## epsilon
    converged <- TRUE;
    new.logit.eps <- log(sum( gt[,"g2"] * ref + gt[,"g0"] * alt)) - log(sum( gt[,"g0"] * ref + gt[,"g2"] * alt))
   ## browser()
    if((logit.eps-new.logit.eps)>tol)
      converged <- FALSE;
    logit.eps <- new.logit.eps;
    eps <- plogis(logit.eps)
    ## ## ## opt 2
    if(!fixGprior){
      log.gmat[,"g0"] <- log.gt[,"g0"]
      log.gmat[,"g1"] <- log.gt[,"g1"]
      log.gmat[,"g2"] <- log.gt[,"g2"]
    }
    ## ##
    it.num <- it.num +1;
    cat("#it:",it.num,"eps=",eps,"Post:",colMeans(gt),"loglik",logliksum,"DeltaLogLik",abs(new.logliksum-logliksum),"\n");
    stopifnot(!is.na(eps))
    logliksum <- new.logliksum;
  }
  invisible(list(gt=gt,log.gt=log.gt,eps=eps,loglik=loglik,logliksum=logliksum))
}

fitAseNullMulti <- function(ref,alt,eps=rep(0.1,ncol(ref)),log.gmat,max.it=100,tol=1E-16,fixGprior=TRUE){
	L <- nrow(ref);
	S <- ncol(ref);
	##browser()
	stopifnot(L == nrow(alt));
	## Parameter initialization
	##log.gmat <- matrix(log(1/3),L,3)
	if(missing(log.gmat)){
		log.gmat <- matrix(log(c(0.25,0.5,0.25)),1,3)
		colnames(log.gmat) <- c("g0","g1","g2");
	}
	it.num <- 0
	logit.eps <- qlogis(eps);
	converged <- FALSE;
	logliksum <- 0;
	rs <- rowSums((ref + alt)) * log(0.5);
	while(it.num <= max.it){
		############ E-step #############
		log.gt <- cbind(ref %*% log(1 - eps) + alt %*% log(eps) + log.gmat[,"g0"],
				        rs  + log.gmat[,"g1"],
				        ref %*% log(eps) + alt %*% log(1 - eps) + log.gmat[,"g2"]
		)
		colnames(log.gt) <- c("g0","g1","g2");
		log.gt[log.gt < (-200)] <- (-200)
		##browser()
		## This normalizes marginal posterior probabilities to add 1
		loglik <- log(exp(log.gt) %*% rep(1,3));
		new.logliksum <- sum(loglik);
		if(abs(new.logliksum-logliksum)<tol)
			converged <- TRUE;
		## Normalize such that rows of exp(log.gt) add to 1 
		log.gt <- apply(log.gt,2,function(col){col-loglik})
		gt <- exp(log.gt);
		if((it.num == max.it) | (converged==TRUE))
			break;
		###### M-STEP  #####
		## w6=w7 and w5=w4 , w2=w1, 
		## epsilon
		converged <- TRUE;
		num <- t(gt[,"g2"]) %*% ref + t(gt[,"g0"]) %*% alt;
		den <- t(gt[,"g0"]) %*% ref + t(gt[,"g2"]) %*% alt;
		new.logit.eps <- as.vector(log(num) - log(den));
		## browser()
		if(max(abs(logit.eps-new.logit.eps))>tol)
			converged <- FALSE;
		logit.eps <- new.logit.eps;
		eps <- plogis(logit.eps)
		## ## ## opt 2
		if(!fixGprior){
			log.gmat[,"g0"] <- log.gt[,"g0"]
			log.gmat[,"g1"] <- log.gt[,"g1"]
			log.gmat[,"g2"] <- log.gt[,"g2"]
		}
		## ##
		it.num <- it.num +1;
		cat("#it:",it.num,"eps=",eps,"Post:",colMeans(gt),"loglik",logliksum,"DeltaLogLik",abs(new.logliksum-logliksum),"\n");
		stopifnot(!is.na(eps))
		logliksum <- new.logliksum;
	}
	invisible(list(gt=gt,log.gt=log.gt,eps=eps,loglik=loglik,logliksum=logliksum))
}



#########################

testPval <- function(eps){
  log.gt <- cbind(g0   = ref * log(1 - eps) + alt * log(eps),
                  g1 = (ref + alt) * log(0.5),
                  g2   = ref * log(eps) + alt * log(1 - eps)
                  )
  aux <- apply(log.gt,1,which.max)
  hetInd <- (aux==2);
  rho <- comp.rho(ref,alt,eps)
  lrt <- lrtEpsRhoBinom(ref,alt,eps,rho)
  pval <- (1-pchisq(2*lrt,df=1))
  pval[hetInd]
}

## qqplot(-log10(runif(length(pval[hetInd]),0.0,1)),-log10(pval[hetInd]))
## abline(0,1)

## sapply(c(0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3),function(eps){
##   aux <- testPval(eps)
##   qq <- qqplot(-log10(runif(length(aux),0.0,1)),-log10(aux),plot.it=F)
##   points(qq,pch='.',cex=3,col='red')
##   eps
## })



## p <- seq(0,1,0.01)
## plot(p,logLikBetaBinomial(60,30,p,100),type="l")
## lines(p,logLikBetaBinomial(60,30,p,10000))
## lines(p,logLikBetaBinomial(60,30,p,0.5))
## lines(p,logLikRhoEps(p,0.1,60,30),col='red',lty=2)
## lines(p,logLikRhoEps(p,0.01,60,30),col='red')
## lines(p,logLikBetaBinomialRhoEps(p,0.1,100,60,30),col='blue')

## plot(p,logLikBetaBinomial(50,4,p,100),type="l",ylim=c(-40,-15))
## lines(p,logLikBetaBinomial(50,4,p,10000))
## lines(p,logLikBetaBinomial(50,4,p,0.5))
## lines(p,logLikRhoEps(p,0.3,50,4),col='red',lty=3)
## lines(p,logLikRhoEps(p,0.1,50,4),col='red',lty=2)
## lines(p,logLikRhoEps(p,0.01,50,4),col='red')
## lines(p,logLikBetaBinomialRhoEps(p,0.01,10,50,4),col='blue')
## lines(p,logLikBetaBinomialRhoEps(p,0.01,30,50,4),col='blue',lty=3)

## plot(p,logLikBetaBinomial(56,11,p,10000),type="l",ylim=c(-100,-30))
## lines(p,logLikBetaBinomial(56,11,p,10000))
## lines(p,logLikRhoEps(p,eps,56,11),col='red')
## lines(p,logLikRhoEps(p,0.1,56,11),col='red',lty=3)
## lines(p,logLikBetaBinomialRhoEps(p,eps,1000,56,11),col='blue',lty=3)
## lines(p,logLikBetaBinomialRhoEps(p,eps,30,56,11),col='blue')

## lines(p,logLikBetaBinomial(10,10,p,1))

## plot(p,logLikBetaBinomial(200,100,p,1000),type="l",lwd=4)
## lines(p,logLikBetaBinomial(200,100,p,100),lwd=3)
## lines(p,logLikBetaBinomial(200,100,p,10),lwd=2)
## lines(p,logLikBetaBinomial(200,100,p,1),lwd=1)

## plot(p,logLikBetaBinomial(4,2,p,10),type="l",lwd=4)
## lines(p,logLikBetaBinomial(20,10,p,10),lwd=1)
## lines(p,logLikBetaBinomial(200,100,p,10),lwd=1)
## lines(p,logLikBetaBinomial(2000,1000,p,10),lwd=1)

#########################
#########################
#########################
#########################

## qqp <- qqplot(-log10(runif(length(pval2[pval2<0.999]),0.0,0.999)),-log10(pval2[pval2<0.999]),plot.it=F)
## plot(qqp,pch='.',cex=2)
## ##points(qqp,pch='.',cex=2,col='blue')
## abline(0,1)

#########################

## size.h0 <- ref+alt
## aux <- apply(ase3$gt,1,which.max)
## eps.h0 <- 0.05
## myprob <- c(1-eps.h0,0.5,eps.h0)
## ref.h0 <- rbinom(length(ref),size=size.h0,prob=myprob[aux])
## ase.h0 <- fitAseNull(ref.h0,size.h0-ref.h0)
## rho.h0 <- comp.rho(ref.h0,size.h0-ref.h0,ase.h0$eps)
## lrt.h0 <- lrtEpsRhoBinom(ref.h0,size.h0-ref.h0,ase.h0$eps,rho.h0)
## pval.h0 <- (1-pchisq(2*lrt.h0,df=1))
## hetInd.h0 <- ase.h0$gt[,2]>0.9
## ##hist(pval.h0[hetInd.h0],breaks=200)
## qq2 <- qqplot(-log10(runif(length(pval.h0[pval.h0<0.999 & hetInd.h0]),0.0,0.999)),-log10(pval.h0[pval.h0<0.999 & hetInd.h0]),plot.it=F)
## ##plot(qq,pch='.',cex=3)
## points(qq,pch='.',cex=3)
## abline(0,1)

#########################

## D <- 1:100;
## D <- exp((0:500)/50)
## ##D <- seq(1,40,0.01)
## aux <- sapply(D,function(D){
##   sum(logLikBetaBinomialRhoEps(0.5,eps,D,ref[hetInd],alt[hetInd]))
##   ##fnLogLikBetaBinomialNull(log(D),ref[hetInd],alt[hetInd])
##   ##grLogLikBetaBinomialNull(log(D),ref[hetInd],alt[hetInd])
## })
## D[which.max(aux)]
## plot(D,aux,log="x",pch='.',cex=3)
##Dmax <- min(D[which(aux>=max(aux)-2)])
##Dmax <- D[which.max(aux)]
##Dmax

#########################

## TRY:
## Use allele frequencies, pop-gen. 

## t1 <- 0.1
## log.gt <- cbind(g0   = ref * log(1 - eps) + alt * log(eps) + log.gmat[,"g0"],
##                 g1t0 = (ref + alt) * log(0.5) + log.gmat[,"g1"] + log(1-t1),
##                 g1t1 = ref * log((1 - eps) * rho + eps * (1 - rho)) + alt * log((1 - eps) * (1 - rho) + eps * rho) + log.gmat[,"g1"] + log(t1),
##                 g2   = ref * log(eps) + alt * log(1 - eps) + log.gmat[,"g2"]
##                 )
## aux <- log(exp(log.gt) %*% rep(1,4));
## log.gt <- apply(log.gt,2,function(col){col-aux})
## gt <- exp(log.gt);

## lrt <- log.gt[,"g1t1"]-pmax(log.gt[,"g0"],log.gt[,"g1t0"],log.gt[,"g2"])

## pval <- chisq.test(2*lrt)

#########################

## a <- 1
## b <- 1
## N <- 10
## compS <- function(N,a,b){
##   p <- rbeta(N,a,b)
##   X <- (runif(N)<=p)
##   sum(X)
## }
## p <- rbeta(N,a,b)
## mean(p)
## var(p)

## a/(a+b)*N
## a*b/(a+b)^2*N
## S <- replicate(100000,compS(N,a,b))
## mean(S)
## var(S)

## ##qqplot(rbinom(10000,N,a/(a+b)),S)
## ##abline(0,1)

## aux <- table(S)
## x <- as.numeric(names(aux))
## y <- aux/sum(aux)
## exp.y <- dbinom(x,N,a/(a+b),log=TRUE);

## plot(-exp.y,-log(y))
## abline(0,1)



#####################################

## colSums(gt)

## tr.asb <- 0.5
## idx <- (gt[,"g1t1"]>tr.asb)
## aux <- (sum(rho[idx]>0.5)/sum(rho[(gt[,"g1t1"]>=0.0)]>0.5))/
## (sum(rho[idx]<0.5)/sum(rho[(gt[,"g1t1"]>=0.0)]<0.5))
## aux
## hist(rho[idx],breaks=200)
## sum(idx)
## min(ref[idx])
## min(alt[idx])
## max(rho[idx])
## min(rho[idx])

## head(gt[ which(ref==alt),])

## round(colSums(gt))

## head(rho[ which(ref==alt)])
## head(rho)

## cat("\n") 
## cat("\n") 


