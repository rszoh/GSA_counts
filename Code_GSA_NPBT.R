### Test an example of gene set testing  using 
### DP mixture with atom at 0
### We use the algorigthm 
require(MCMCpack)
require(matlab)
require(mnormt)


Loglikhodfunc<- function(val,i,X,Y,Sig,delt0){ ## 
  n1 <- nrow(X)
  n2 <- nrow(Y)
  delt <- delt0
  delt[i] <- val 
  tp <- solve(Sig)%*%(n2*tcrossprod(matrix(colMeans(Y)-colMeans(X)-delt,ncol=1)) + crossprod(scale(X,scale=F)) + crossprod(scale(Y,scale=F)))
  as.numeric(-.5*(n1+n2-1)*(unlist(determinant(Sig,logarithm=T))[1]) - .5*sum(diag(tp)))
}

Loglikhodfunc <- Vectorize(Loglikhodfunc,c("val")) ## vectorize the function with respect to its arguments


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

#***************************************************************************************************
#***************************************************************************************************
#***************************************************************************************************
## Unknown Sigma
UpdatCi_alg8 <- function(Cii_nw,X,Y,phiC_nw,alpha0,parm0,M,Sig){ 
  p <- ncol(Y)
  Cii <- Cii_nw
  #phiC <- numeric(n)
  phiC <- phiC_nw ## vector of length p
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  p0 <- parm0[3]
  cp <- 0
  for(i in 1:p){
    ival <- unique(Cii[-i]) ## unique delta
    K_m <- length(ival)
    h <- K_m + M
    phiC_tp <- numeric(h)
    phiC_tp[1:K_m] <- phiC[ival] ## new phiC_tp with respect to new Ciis
    Cil <- Cii[-i] ## label for the ci
    for(l0 in 1:K_m){ Cil[Cil == ival[l0]] <- l0} ## relabel the Cis
    if(sum(Cii[i]== Cil) > 0){## Ci is not a singleton
      roll <- rbinom(n=M,size=1,prob=p0)
      phiC_tp[(K_m+1):h] <- (1-roll)*rnorm(n=M,mean=mu0,sd=sqrt(sig0))} ##  simulate new values from the prior
    else{
      roll <- rbinom(n=M,size=1,prob=p0)
      Cii[i] <- K_m+1
      phiC_tp[(K_m+1):h] <- (1-roll)*rnorm(n=M,mean=mu0,sd=sqrt(sig0))}
    Cii[-i]<- Cil
    prob0 <- numeric(h)
#     cat("delta:",phiC_tp[Cii],"\n")
#     cat("Cii:",Cii,"\n")
    cat("Diff:",length(Loglikhodfunc(unique(Cii[-i]),i,X,Y,Sig,phiC_tp[Cii])) - length(log(as.vector(table(Cil))*(1/(p-1+alpha0)))),"---")
    cat("Diff1-2:",length(Loglikhodfunc(unique(Cii[-i]),i,X,Y,Sig,phiC_tp[Cii]))- K_m,"\n")
    #cat("Loglike:",length(Loglikhodfunc(unique(Cii[-i]),i,X,Y,Sig,phiC_tp[Cii])),"\n")
  #cat("Other",length(log(as.vector(table(Cil))*(1/(p-1+alpha0)))),"\n")
#     cat("CSt:",log(as.vector(table(Cii[-i]))*(1/(p-1+alpha0))),"\n")
    prob0[1:K_m] <- log(as.vector(table(Cil))*(1/(p-1+alpha0))) + Loglikhodfunc(phiC_tp[1:K_m],i,X,Y,Sig,phiC_tp[Cii])#   pnorm(Y[i]-phiC_tp[1:K_m])*
    prob0[(K_m+1):h] <- log((alpha0/M/(p-1+alpha0))) + Loglikhodfunc(phiC_tp[(K_m+1):h],i,X,Y,Sig,phiC_tp[Cii])   #pnorm(Y[i]-phiC_tp[(K_m+1):h]) 
    #cat("probo=",prob0,"\n")
    prob0 <- exp(prob0-max(prob0)) ## adjust the sacel a little bit
#     cat("probo=",prob0,"\n")
#     cat("length(prob0)=",length(prob0),"\n")
#     cat("h:",h,"\n")
    Cii[i] <- (1:h)[rmultinom(n=1,size=1,prob=prob0)==1] 
    phiC <- phiC_tp[Cii]
    ival <- unique(Cii)
    for(j in 1:length(ival)){ Cii[Cii == ival[j]] <- j}
  }
  
  #phiC <- phiC[Cii]
  #for(i in 1:length(ival)){ Cii[Cii == ival[i]] <- i}
  return(list("ci"=Cii,"Phic"=phiC))
} 

## known Sigma
UpdatCi_alg8_KnSig <- function(Cii_nw,X,Y,phiC_nw,alpha0,parm0,M,Sig){ 
  p <- ncol(Y)
  Cii <- Cii_nw
  phiC <- phiC_nw  ##vector of varying length unique \tau_s
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  p0 <- parm0[3]
  cp <- 0
  for(i in 1:p){
    ival <- unique(Cii[-i]) ## unique delta except i
    K_m <- length(ival)
    h <- K_m + M
    phiC_tp <- numeric(h)
    phiC_tp[1:K_m] <- unique(phiC[Cii[-i]]) # new phiC_tp with respect to new Ciis
    Cil <- Cii[-i] ## label for the ci
    for(l0 in 1:K_m){ Cil[Cii[-i] == ival[l0]] <- l0} ## relabel the Cis
    
    if(sum(Cii[i]== Cil) > 0){## Ci is not a singleton
      roll <- rbinom(n=M,size=1,prob=p0)
      phiC_tp[(K_m+1):h] <- (1-roll)*rnorm(n=M,mean=mu0,sd=sqrt(sig0)) ##  simulate new values from the prior
        }else{
      roll <- rbinom(n=M-1,size=1,prob=p0)
      phiC_tp[(K_m+1)] <- phiC[Cii[i]]
      Cii[i] <- K_m+1
      phiC_tp[(K_m+2):h] <- (1-roll)*rnorm(n=M-1,mean=mu0,sd=sqrt(sig0))}
    Cii[-i] <- Cil
    prob0 <- numeric(h)
    #     cat("delta:",phiC_tp[Cii],"\n")
    #     cat("Cii:",Cii,"\n")
    cat("Diff:",length(Loglikhodfunc_knSig(phiC_tp[1:K_m],i,X,Y,Sig,phiC_tp[Cii])) - length(as.vector(table(factor(Cil,levels=unique(Cil))))*(1/(p-1+alpha0))),"---")
    cat("Diff1-2:",length(Loglikhodfunc_knSig(phiC_tp[1:K_m],i,X,Y,Sig,phiC_tp[Cii]))- K_m," --")
    cat("unique CI",length(unique(Cil)),"--K_m:",K_m,"--i=",i,"\n")
    if(K_m!=length(unique(Cii[-i]))) {stop("violation \n")}
    #cat("Loglike:",length(Loglikhodfunc(unique(Cii[-i]),i,X,Y,Sig,phiC_tp[Cii])),"\n")
    #cat("Other",length(log(as.vector(table(Cil))*(1/(p-1+alpha0)))),"\n")
    #     cat("CSt:",log(as.vector(table(Cii[-i]))*(1/(p-1+alpha0))),"\n")
    prob0[1:K_m] <- log(as.vector(table(factor(Cil,levels=unique(Cil))))*(1/(p-1+alpha0))) + Loglikhodfunc_knSig(phiC_tp[1:K_m],i,X,Y,Sig,phiC_tp[Cii])#   pnorm(Y[i]-phiC_tp[1:K_m])*
    prob0[(K_m+1):h] <- log((alpha0/M/(p-1+alpha0))) + Loglikhodfunc_knSig(phiC_tp[(K_m+1):h],i,X,Y,Sig,phiC_tp[Cii])   #pnorm(Y[i]-phiC_tp[(K_m+1):h]) 
    #cat("probo=",prob0,"\n")
    prob0 <- exp(prob0-max(prob0)) ## adjust the sacel a little bit
    Cii[i] <- (1:h)[rmultinom(n=1,size=1,prob=prob0)==1] 
    #phiC <- (phiC_tp[unique(Cii)])
    if(sum(phiC_tp[Cii] == 0.0)> 1)
    {id01 <- (1:length(Cii))[phiC_tp[Cii] == 0.0] ## index of all the Ciis pointing to zero
     if(length(unique(Cii[id01])) > 1){ Cii[id01] <- min(Cii[id01])}
    }
    phiC <- unique(phiC_tp[Cii]) ## unique values of delta
    ival <- unique(Cii)
    Ci0 <- Cii
    for(j in 1:length(ival)){ Ci0[Cii == ival[j]] <- j}
    Cii <- Ci0
    
    
    phiC_nw00 <- phiC_tp
    Cii00 <- Cii
    
  }
  
  #phiC <- phiC[Cii]
  #for(i in 1:length(ival)){ Cii[Cii == ival[i]] <- i}
  return(list("ci"=Cii,"Phic"=phiC,"Cii00"= Cii00,"phiC00"=phiC_nw00))
} 


## known Sigma new
UpdatCi_alg8_KnSig2 <- function(X,Y,phiC_nw,alpha0,parm0,M,Sig){ 
  p <- ncol(Y)
  phiC <- phiC_nw  ##vector of varying length unique \tau_s
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  p0 <- parm0[3]
  cp <- 0
  for(i in 1:p){
    ival <- unique(phiC[-i]) ## unique delta except i
    K_m <- length(ival)
    h <- K_m + M
    phiC_tp <- numeric(h)
    phiC_tp[1:K_m] <- unique(phiC[-i]) # new phiC_tp with respect to new Ciis

    if(sum(phiC[i]== phiC) > 0){## Ci is not a singleton
      roll <- rbinom(n=M,size=1,prob=p0)
      phiC_tp[(K_m+1):h] <- (1-roll)*rnorm(n=M,mean=mu0,sd=sqrt(sig0)) ##  simulate new values from the prior
    }else{
      roll <- rbinom(n=M-1,size=1,prob=p0)
      phiC_tp[(K_m+1)] <- phiC[i]
      phiC_tp[(K_m+2):h] <- (1-roll)*rnorm(n=M-1,mean=mu0,sd=sqrt(sig0))}
      
    prob0 <- numeric(h)
    #     cat("delta:",phiC_tp[Cii],"\n")
    #     cat("Cii:",Cii,"\n")
    #cat("Diff:",length(Loglikhodfunc_knSig(phiC_tp[1:K_m],i,X,Y,Sig,phiC_tp[Cii])) - length(as.vector(table(factor(Cil,levels=unique(Cil))))*(1/(p-1+alpha0))),"---")
    #cat("Diff1-2:",length(Loglikhodfunc_knSig(phiC_tp[1:K_m],i,X,Y,Sig,phiC_tp[Cii]))- K_m," --")
    #cat("unique CI",length(unique(Cil)),"--K_m:",K_m,"--i=",i,"\n")
    if(K_m!=length(unique(phiC[-i]))) {stop("violation \n")}
    #cat("Loglike:",length(Loglikhodfunc(unique(Cii[-i]),i,X,Y,Sig,phiC_tp[Cii])),"\n")
    #cat("Other",length(log(as.vector(table(Cil))*(1/(p-1+alpha0)))),"\n")
    #     cat("CSt:",log(as.vector(table(Cii[-i]))*(1/(p-1+alpha0))),"\n")
    prob0[1:K_m] <- log(as.vector(table(factor(phiC[-i],levels=unique(phiC[-i]))))*(1/(p-1+alpha0))) + Loglikhodfunc_knSig(phiC_tp[1:K_m],i,X,Y,Sig,phiC)#   pnorm(Y[i]-phiC_tp[1:K_m])*
    prob0[(K_m+1):h] <- log((alpha0/M/(p-1+alpha0))) + Loglikhodfunc_knSig(phiC_tp[(K_m+1):h],i,X,Y,Sig,phiC)   #pnorm(Y[i]-phiC_tp[(K_m+1):h]) 
    #cat("probo=",prob0,"\n")
    prob0 <- exp(prob0-max(prob0)) ## adjust the sacel a little bit
    phiC[i] <- phiC_tp[rmultinom(n=1,size=1,prob=prob0)==1] 
  }
  return(phiC)
} 
#**************************************************************************************************
# Mean $\tilde{\mu}_0
Post_delt <- function(X,Y,delti,Cii,parm0,Sig,delt){ ## delti is a single value of delta, 
  mu0 <-parm0[1]
  tau0 <- parm0[2]
  p0 <- parm0[3]
  
  id0 <- (1:length(Cii))[Cii==delti]  ## find the Cii with same Cii values
  id1 <- (1:length(Cii))[Cii!=delti]  ## find the Cii with same Cii values
  D <- colMeans(Y)-colMeans(X)
mu0 <- sum(D[id0] + Sig[id0,id1]%*%(solve(Sig[id1,id1])%*%matrix(delt[id1] - D[id1],ncol=1)))
on0 <- matrix(rep(1,length(id0)),ncol=1)
 sig0 <- t(on0)%*%(Sig[id0,id0] - Sig[id0,id1]%*%(solve(Sig[id1,id1])%*%Sig[id1,id0]))%*%on0 ## 
 sigs <- 1/(1/tau0 + 1/sig0) 
 mus <- mu0*sigs/sig0
 
 pis <- (1 + ((1-p0)/p0)*dnorm(0,sd=sqrt(tau0))/dnorm(0,mean=mus,sd=sqrt(sigs)))^{-1}
 
 roll <- rbinom(n=1,size=1,prob=pis)
 (1-roll)*rnorm(n=1,mean=mus,sd=sqrt(sigs))
}

# Mean $\tilde{\mu}_0
Post_delt_KnSig_old <- function(delti,X,Y,parm0,Sig0,delt){ ## delti is a single value of delta, 
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
##%%%%% MCMC steps
McMUp_alg8 <- function(Cii,X,Y,phiC,Nsam,alpha0,parm0,M,Sig0){
  out_phic <- NULL
  out_cii <- NULL
  out_sig <- NULL
  Cii_new <- Cii
  phiC_new <- phiC 
  Sig <- Sig0
  n2 <- nrow(Y)
  n1 <- nrow(X)
  p <- ncol(X)
  for(m in 1:Nsam){
    #cat("m=",m,"--")
    out<- UpdatCi_alg8(Cii_new,X,Y,phiC_new,alpha0,parm0,M,Sig)
    Cii_new <- out[["ci"]]
    phiC_new <- out[["Phic"]]
    
    c_val <- unique(Cii_new) ## values of phiC non na!
    for(lx in 1:length(c_val)){
      phiC_new[Cii_new==c_val[lx]] <- Post_delt(X,Y,c_val[lx],Cii,parm0,Sig/n2,phiC_new) # SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0)}
    }
    #SS0 <- n2*tcrossprod(matrix(colMeans(Y)-colMeans(X)- phiC_new,ncol=1)) + crossprod(scale(X,scale=F)) + crossprod(scale(Y,scale=F))
      #n2*tcrossprod(colMeans(Y)-colMeans(X)- phiC_new) + tcrossprod(X) +tcrossprod(Y)
    #Sig<- riwish(v=p+n1+n2,S = SS0)   
    Sig <- diag(rep(1,p))
    out_phic <- rbind(out_phic, phiC_new)
    out_cii <- rbind(out_cii,Cii_new)
    out_sig <- rbind(out_sig,lower.tri(Sig,diag=T))
  }
  return(list("Phic"=out_phic,"ci"=out_cii))
}

##%%%%% MCMC steps (known covariance function)
McMUp_alg8KnSig <- function(Cii,X,Y,phiC,Nsam,alpha0,parm0,M,Sig0){
  out_phic <- NULL
  out_cii <- NULL
  out_sig <- NULL
  Cii_new <- Cii
  phiC_new <- phiC 
  Sig <- Sig0
  n2 <- nrow(Y)
  n1 <- nrow(X)
  p <- ncol(X)
  for(m in 1:Nsam){
    #cat("m=",m,"--")
    
    cat("\n Starting \n  Cii",unique(Cii_new),"-- delta",phiC_new,"\n \n")
    
    if(sum(phiC_new[Cii_new] == 0.0)> 1) ## multiple 0
    {id01 <- (1:length(Cii_new))[phiC_new[Cii_new] == 0.0] ## index of all the Ciis pointing to zero
    if(length(unique(Cii_new[id01])) > 1){ Cii_new[id01] <- min(Cii_new[id01])}
    }
    phiC_new <- unique(phiC_new[Cii_new]) ## unique values of delta
    #ival <- unique(Cii)
    
    cat("\n Starting \n  Cii",length(unique(Cii_new)),"-- delta",length(phiC_new),"\n \n")
    #ival <- unique(phiC_new[Cii_new])
    #for(i in 1:length(ival)){ Cii_new[phiC_new[Cii_new] == ival[i]] <- i}
    out<- UpdatCi_alg8_KnSig(Cii_new,X,Y,phiC_new,alpha0,parm0,M,Sig)
    Cii_new <- out$ci
    phiC_new <- out$Phic
  cat("\n \n Cii",unique(Cii_new),"-- delta",phiC_new,"\n \n")
     for(lx in 1:length(phiC_new)){
       phiC_new[lx] <- Post_delt_KnSig(phiC_new[lx],X,Y,Cii_new,parm0,Sig,phiC_new[Cii_new]) # SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0)}
     }
  
    Sig <- diag(rep(1,p))
    out_phic <- rbind(out_phic, phiC_new[Cii_new])
    out_cii <- rbind(out_cii,Cii_new)
    out_sig <- rbind(out_sig,lower.tri(Sig,diag=T))
  }
  return(list("Phic"=out_phic,"ci"=out_cii))
}


##%%%%% MCMC steps (known covariance function) ## new methods
McMUp_alg8KnSig2 <- function(X,Y,phiC,Nsam,alpha0,parm0,M,Sig0){
  out_phic <- NULL
  out_cii <- NULL
  out_sig <- NULL
  phiC_new0 <- phiC 
  Sig <- Sig0
  n2 <- nrow(Y)
  n1 <- nrow(X)
  p <- ncol(X)
  for(m in 1:Nsam){
    phiC_new <- UpdatCi_alg8_KnSig2(X,Y,phiC_new0,alpha0,parm0,M,Sig)
    phiC_new00 <- phiC_new
    phiC_new01 <- phiC_new
    ival <- unique(phiC_new)
    for(lx in 1:length(ival)){
      id00 <- (1:length(phiC_new01))[phiC_new01 == ival[lx]]
      phiC_new00[id00] <- Post_delt_KnSig(ival[lx],X,Y,parm0,Sig,phiC_new) # SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0)}
      #phiC_new <- phiC_new0
      }
    phiC_new <- phiC_new00
    Sig <- diag(rep(1,p))
    out_phic <- rbind(out_phic, phiC_new)
    #out_cii <- rbind(out_cii,Cii_new)
    #out_sig <- rbind(out_sig,lower.tri(Sig,diag=T))
  }
  return(out_phic)
}

#####################################################
#Examples
set.seed(1)
  Nsam <- 5000
  M=5
  alpha0 <- 1
  p = 30
  n1 = 10 
  n2 = 10
  parm0 <-c(0,3,.977)
  phiC_nw <- rep(0,p)
  phiC_new <- rnorm(p)
  Cii <- 1:p
  Sig<- diag(rep(1,p))
  delt0 <- rep(0,p)
  X <- rmnorm(n=n1,mean=rep(0,p),varcov = Sig)
  q <- 2
    mu2 <-2 
  meany <- c(rep(mu2,q),rep(0,p-q))
  Y <- rmnorm(n=n2,mean=meany,varcov = Sig)
  
  out <- Loglikhodfunc_knSig(1,1,X,Y,Sig,delt0) 
  out0 <- UpdatCi_alg8_KnSig2(X,Y,phiC_new,alpha0,parm0,M=3,Sig)
  out0
  Post_delt_KnSig(out0[1],X,Y,parm0,Sig,out0)
  out1 <- McMUp_alg8KnSig2(X,Y,out0,Nsam,alpha0,parm0,M,Sig) 
 mean(rowSums(out1==0)==ncol(out1)) 
 mean(rowSums(cbind(out1[,1]!=0,out1[,2]!=0,out1[,-c(1:2)]==0))==ncol(out1))
  colMeans(out1)
  plot(out1[1:100,1],type="s")
  
out <- Loglikhodfunc_knSig(1,1,X,Y,Sig,delt0) 
out0 <- UpdatCi_alg8_KnSig(Cii,X,Y,phiC_new,alpha0,parm0,M=3,Sig)
out0$ci
out0$Phic
UpdatCi_alg8_KnSig(out0$ci,X,Y,out0$Phic,alpha0,parm0,M=3,Sig)

Post_delt_KnSig(X,Y,out0$Phic[1],out0$ci,parm0,Sig,out0$Phic[out0$ci])
out1 <- McMUp_alg8KnSig(Cii,X,Y,phiC_new[Cii],Nsam,alpha0,parm0,M,Sig) 

out <- Loglikhodfunc(1,1,X,Y,Sig,delt0) 
out0 <- UpdatCi_alg8(Cii,X,Y,phiC_nw,alpha0,parm0,M=3,Sig)
out1 <- McMUp_alg8(Cii,X,Y,phiC_nw,Nsam,alpha0,parm0,M,Sig) 

