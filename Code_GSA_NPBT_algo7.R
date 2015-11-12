#########################################################################
#####            Algo7
#########################################################################
require(MCMCpack)
require(matlab)
require(mnormt)
### Likehood function
Loglikhodfunc_knSig<- function(val,i,X,Y,Sig,delt0){ ## 
  n1 <- nrow(X)
  n2 <- nrow(Y)
  delt <- delt0
  delt[i] <- val 
  n0 <- (1/n1 +1/n2)^{-1}
  tp <- solve(Sig)%*%(n0*tcrossprod(matrix(colMeans(Y)-colMeans(X)-delt,ncol=1)))
  as.numeric(-.5*sum(diag(tp)))
}

Loglikhodfunc_knSig <- Vectorize(Loglikhodfunc_knSig,c("val")) ## vectorize the function with respect to its arguments

### Simulate a new Phic
# Mean $\tilde{\mu}_0
Post_delt_KnSig <- function(delti,X,Y,parm0,Sig0,delt){ ## delti is a single value of delta, 
  n1 <- nrow(X)
  n2 <- nrow(Y)
  n0 <- (1/n1 +1/n2)^{-1}
  mu0 <- parm0[1]
  tau0 <- parm0[2]
  p0 <- parm0[3]
  Sig <- Sig*(1/n0)
  id0 <- (1:length(delt))[delt==delti]  ## find the Cii with same Cii values
  id1 <- (1:length(delt))[delt!=delti]  ## find the Cii with same Cii values
  D <- colMeans(Y)-colMeans(X)
  invSig <- chol2inv(Sig[id1,id1])
  mu0 <- sum(matrix(D[id0],ncol=1) + Sig[id0,id1]%*%(invSig%*%matrix(delt[id1] - D[id1],ncol=1)))
  on0 <- matrix(rep(1,length(id0)),ncol=1)
  sig0 <- t(on0)%*%(Sig[id0,id0] - Sig[id0,id1]%*%(invSig%*%Sig[id1,id0]))%*%on0 ## 
  sigs <- 1/(1/tau0 + 1/sig0) 
  mus <- mu0*sigs/sig0
  
  pis <- (1 + ((1-p0)/p0)*dnorm(0,sd=sqrt(tau0))/dnorm(0,mean=mus,sd=sqrt(sigs)))^{-1}
  
  roll <- rbinom(n=1,size=1,prob=pis)
  (1-roll)*rnorm(n=1,mean=mus,sd=sqrt(sigs))
}

Post_delt_KnSig <- function(delti,X,Y,parm0,Sig0,delt){ ## delti is a single value of delta, 
  n1 <- nrow(X)
  n2 <- nrow(Y)
  n0 <- (1/n1 +1/n2)^{-1}
  mu0 <- parm0[1]
  tau0 <- parm0[2]
  p0 <- parm0[3]
  Sig <- Sig*(1/n0)
  id0 <- (1:length(delt))[delt==delti]  ## find the Cii with same Cii values
  if(length(id0) < length(delt)){
    id1 <- (1:length(delt))[delt!=delti]  ## find the Cii with same Cii values
    D <- colMeans(Y)-colMeans(X)
    invSig <- solve(Sig[id1,id1])
    mu0 <- sum(matrix(D[id0],ncol=1) + Sig[id0,id1]%*%(invSig%*%matrix(delt[id1] - D[id1],ncol=1)))
    on0 <- matrix(rep(1,length(id0)),ncol=1)
    sig0 <- t(on0)%*%(Sig[id0,id0] - Sig[id0,id1]%*%(invSig%*%Sig[id1,id0]))%*%on0 ## 
  } else{
    D <- colMeans(Y)-colMeans(X)
    mu0 <- sum(D)
    on0 <- matrix(rep(1,length(id0)),ncol=1)
    sig0 <- t(on0)%*%(Sig[id0,id0])%*%on0
  }
  sigs <- 1/(1/tau0 + 1/sig0) 
  mus <- mu0*sigs/sig0
  
  pis <- (1 + ((1-p0)/p0)*dnorm(0,sd=sqrt(tau0))/dnorm(0,mean=mus,sd=sqrt(sigs)))^{-1}
  
  roll <- rbinom(n=1,size=1,prob=pis)
  (1-roll)*rnorm(n=1,mean=mus,sd=sqrt(sigs))
}

Post_delt_KnSig <- function(delti,X,Y,parm0,Sig0,delt){ ## delti is a single value of delta, 
  n1 <- nrow(X)
  n2 <- nrow(Y)
  n0 <- (1/n1 +1/n2)^{-1}
  mu0 <- parm0[1]
  tau0 <- parm0[2]
  p0 <- parm0[3]
  Sig <- Sig*(1/n0)
  id0 <- (1:length(delt))[delt==delti]  ## find the Cii with same Cii values
  if(length(id0) < length(delt)){
    id1 <- (1:length(delt))[delt!=delti]  ## find the Cii with same Cii values
    D <- colMeans(Y)-colMeans(X)
    invSig <- solve(Sig[id1,id1])
    mu0 <- mean(matrix(D[id0],ncol=1) + Sig[id0,id1]%*%(invSig%*%matrix(delt[id1] - D[id1],ncol=1)))
    on0 <- matrix(rep(1,length(id0)),ncol=1)/length(id0)
    sig0 <- t(on0)%*%(Sig[id0,id0] - Sig[id0,id1]%*%(invSig%*%Sig[id1,id0]))%*%on0 ## 
  } else{
    D <- colMeans(Y)-colMeans(X)
    mu0 <- mean(D)
    on0 <- matrix(rep(1,length(id0)),ncol=1)/length(id0)
    sig0 <- t(on0)%*%(Sig[id0,id0])%*%on0
  }
  sigs <- 1/(1/tau0 + 1/sig0) 
  mus <- mu0*sigs/sig0
  
  pis <- (1 + ((1-p0)/p0)*dnorm(0,sd=sqrt(tau0))/dnorm(0,mean=mus,sd=sqrt(sigs)))^{-1}
  
  roll <- rbinom(n=1,size=1,prob=pis)
  (1-roll)*rnorm(n=1,mean=mus,sd=sqrt(sigs))
}
######################################################################### 
## Prog1
UpdatCi_alg7_1 <- function(Cii_nw,X,Y,phiC_nw,alpha0,parm0,Sig){ 
  n <- length(Cii_nw)
  Cii <- Cii_nw
  phiC <- rep(NA,2*length(phiC_nw))
  phiC[1:n] <- phiC_nw
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  p0 <- parm0[3]
  
  cp <- 0
  for(i in 1:n){
    if(sum(Cii[i]==Cii[-i]) > 0){## Ci is not a singleton
      ival <- unique(Cii[-i])
      ci_nw <- n+i ## newly created state
      #phiC_tp <- phiC_nw ## create a copy to delt
      phiC_tp <- (1-rbinom(n=M,size=1,prob=p0))*rnorm(n=M,mean=mu0,sd=sqrt(sig0)) # simulate a new phiC 
      r <- log(alpha0)+ Loglikhodfunc_knSig(phiC_tp,i,X,Y,Sig,phiC_nw) - log(n-1) - Loglikhodfunc_knSig(phiC_nw[i],i,X,Y,Sig,phiC_nw)
      acep_prop <- min(c(1, ifelse(r >= 0,1,exp(r))))
      roll <- rbinom(n=1,size=1,prob=acep_prop)
      if(roll){
        Cii[i] <- ci_nw 
        phiC[n+i]<- phiC_tp[i]}
    }else{
      ival <- unique(Cii[-i]) ##distinct states except i
      ci_nw <- ival[rmultinom(n=1,size=1,prob=as.vector(table(factor(Cii[-i],levels = unique(Cii[-i])))/(n-1)))==1] ## chose a value of Cj
      r <- log(n-1)+ Loglikhodfunc_knSig(phiC_nw[ci_nw],i,X,Y,Sig,phiC_nw) - log(alpha0) - Loglikhodfunc_knSig(phiC_nw[i],i,X,Y,Sig,phiC_nw)
      acep_prop <- min(c(1, ifelse(r >= 0,1,exp(r))))
      roll <- rbinom(n=1,size=1,prob=acep_prop)
      Cii[i] <- roll*ci_nw +(1-roll)*Cii[i]
    }
  }
  ival <- unique(Cii) ## look at the unique Cii's
  phiC <- phiC[Cii]   ## Store the delta_is
  for(i in 1:length(ival)){ Cii[Cii == ival[i]] <- i}
  return(list("ci"=Cii,"Phic"=phiC))
} 

## Prog1 a
UpdatCi_alg7_10 <- function(X,Y,phiC,alpha0,parm0,Sig){ 
  n <- length(phiC)
  phiC_nw <- phiC
  #phiC <- rep(NA,2*length(phiC_nw))
  #phiC[1:n] <- phiC_nw
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  p0 <- parm0[3]
  
  cp <- 0
  for(i in 1:n){
    if(sum(phiC_nw[i]==phiC_nw[-i]) > 0){## Ci is not a singleton
      #phiC_tp <- phiC_nw ## create a copy to delt
      phiC_tp <- (1-rbinom(n=1,size=1,prob=p0))*rnorm(n=1,mean=mu0,sd=sqrt(sig0)) # simulate a new phiC 
      r <- log(alpha0)+ Loglikhodfunc_knSig(phiC_tp,i,X,Y,Sig,phiC_nw) - log(n-1) - Loglikhodfunc_knSig(phiC_nw[i],i,X,Y,Sig,phiC_nw)
      acep_prop <- min(c(1, ifelse(r >= 0,1,exp(r))))
      roll <- rbinom(n=1,size=1,prob=acep_prop)
      if(roll){
         phiC_nw[i] <- phiC_tp}
    }else{
      
      phiC_tp <- (unique(phiC_nw[-i]))[rmultinom(n=1,size=1,prob=as.vector(table(factor(phiC_nw[-i],levels = unique(phiC_nw[-i]))))/(n-1))==1] ## chose a value of Cj
      r <- log(n-1)+ Loglikhodfunc_knSig(phiC_tp,i,X,Y,Sig,phiC_nw) - log(alpha0) - Loglikhodfunc_knSig(phiC_nw[i],i,X,Y,Sig,phiC_nw)
      acep_prop <- min(c(1, ifelse(r >= 0,1,exp(r))))
      roll <- rbinom(n=1,size=1,prob=acep_prop)
      phiC_nw[i] <- roll*phiC_tp +(1-roll)*phiC_nw[i]
    }
  }
  
  return(phiC_nw)
} 

###### Algo7 prog 2
UpdatCi_alg7_2 <- function(Cii_nw,X,Y,phiC,alpha0,parm0,Sig){
  n <- length(Cii_nw)
  Cii <- Cii_nw
   phiC_nw <-  phiC
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  p0 <- parm0[3]
  cp <- 0
  for(i in 1:n){
    if(sum(Cii[i]==Cii[-i])> 0){## Ci is not a singleton
      ival <- unique(Cii[-i])
      prob0 = log(as.vector(table(factor(Cii[-i],levels = unique(Cii[-i])))/(n-1))) + Loglikhodfunc_knSig(unique(phiC_nw),i,X,Y,Sig,phiC_nw) 
      # cat("***prob0=",prob0,"\n")
      Cii[i] <- ival[rmultinom(n=1,size=1,prob=prob0)==1]
    }
  }
  return(Cii)
}

## prog 2a
UpdatCi_alg7_20 <- function(X,Y,phiC,alpha0,parm0,Sig){
  n <- length(phiC)
  phiC_nw <-  phiC
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  p0 <- parm0[3]
  cp <- 0
  for(i in 1:n){
    if(sum(phiC_nw[i]==phiC_nw[-i]) > 0){## Ci is not a singleton
      prob0 = log(as.vector(table(factor(phiC_nw[-i],levels = unique(phiC_nw[-i])))/(n-1))) + Loglikhodfunc_knSig(unique(phiC_nw[-i]),i,X,Y,Sig,phiC_nw) 
      prob0 <- prob0 - max(prob0)
      # cat("***prob0=",prob0,"\n")
      phiC_tp <- phiC_nw[-i]
      phiC_nw[i] <- (unique(phiC_nw[-i]))[rmultinom(n=1,size=1,prob=exp(prob0))==1]
    }
  }
  return(phiC_nw)
}

###%%%%% prog 1
McMUp_alg7 <- function(Cii,Y,phiC,Nsam,alpha0,parm0){
  out_phic <- NULL
  out_cii <- NULL
  Cii_new <- Cii
  phiC_new <- phiC 
  
  for(m in 1:Nsam){
    #cat("m=",m,"--")
    out<- UpdatCi_alg7_1(Cii_new,Y,phiC_new,alpha0,parm0)
    Cii_new <- out[["ci"]]
    #cat("Ci=",Cii_new,"\n")
    phiC_new <- out[["Phic"]]
    #cat("PhiC=",phiC_new,"\n")
    out1 <- UpdatCi_alg7_2(Cii_new,Y,phiC_new,alpha0,parm0)
    Cii_new <- out1
    
    #c_val <- unique(Cii_new) ## values of phiC non na!
    #for(l0 in 1:length(c_val)){Cii_new[Cii_new==c_val[l0]] <- l0}
    c_val <- unique(Cii_new) ## values of phiC non na!
    for(lx in 1:length(c_val)){
      phiC_new[Cii_new==c_val[lx]] <- SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0)}
    
    out_phic <- rbind(out_phic, phiC_new)
    out_cii <- rbind(out_cii,Cii_new)
  }
  return(list("Phic"=out_phic,"ci"=out_cii))
}

###%%%%% prog 2
McMUp_alg70 <- function(X,Y,phiC,Nsam,alpha0,parm0,Sig){
  out_phic <- NULL
  phiC_new <- phiC 
  
  for(m in 1:Nsam){
    #cat("m=",m,"--")
    phiC_new1  <- UpdatCi_alg7_10(X,Y,phiC_new,alpha0,parm0,Sig)
    #cat("PhiC=",phiC_new,"\n")
    phiC_new2  <- UpdatCi_alg7_20(X,Y,phiC_new1,alpha0,parm0,Sig)
    
    phiC_new20 <- phiC_new2
    ival <- unique(phiC_new2)
    for(lx in 1:length(ival)){
      id00 <- (1:length(phiC_new2))[phiC_new2 == ival[lx]] 
      phiC_new20[id00] <- Post_delt_KnSig(ival[lx],X,Y,parm0,Sig,phiC_new2) #SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0)}
      #phiC_new2[id00] <- phiC_new20[id00]
    }
    phiC_new2 <- phiC_new20
    out_phic <- rbind(out_phic, phiC_new2)
    #out_cii <- rbind(out_cii,Cii_new)
  }
  return(out_phic)
}


set.seed(1)
Nsam <- 10000
M=3
alpha0 <- 1
p = 30
n1 = 20 
n2 = 20
parm0 <-c(0,5,.95)
phiC_nw <- rnorm(p)
phiC_nw <- rep(0,p)
phiC_new <- rnorm(p)
Sig<- diag(rep(1,p))
delt0 <- rep(0,p)

X <- rmnorm(n=n1,mean=rep(0,p),varcov = Sig)
meany <- c(10,10,rep(0,p-2))
Y <- rmnorm(n=n2,mean=meany,varcov = Sig)


UpdatCi_alg7_10(X,Y,phiC_nw,alpha0,parm0,Sig)
UpdatCi_alg7_20(X,Y,phiC_nw,alpha0,parm0,Sig)

out <- McMUp_alg70(X,Y,phiC_nw,Nsam,alpha0,parm0,Sig)

mean(rowSums(out==0)==ncol(out)) 
mean(rowSums(cbind(out[,1]!=0,out[,2]!=0,out[,-c(1:2)]==0))==ncol(out))

colMeans(out==0)
colMeans(out==0)
#out <- McMUp_alg7(Cii,Y,phiC_new,Nsam,alpha0,parm0)
