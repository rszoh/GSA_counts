
#######################################################################
### Algortihm 5 and 7 Neal 2000
#######################################################################
###Simulate a new state for the latent variables
UpdatCi_alg5 <- function(i,Cii_nw,Y,phiC_nw,alpha0,parm0,R){ 
  Cii <- Cii_nw
  phiC <- phiC_nw
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  n <- length(Y)
  cp <- 0

  replicate(n=R,{
    ival <- unique(Cii)
    ival <- unique(Cii) ## unique phiC
    prob0 <- numeric(length(ival))    ## all currents states
    for(k in 1:length(ival)){prob0[k] <- sum(Cii[-i]==ival[k])/(n-1+alpha0) ## eq 5.4 (a) of Neal(2000)
    }
    prob0 <- c(prob0, alpha0/(n-1+alpha0)) ## eq 5.4 (b) of Neal(2000)
    Cii_temp <- rmultinom(n=1,size=1,prob=prob0) ## standardize the probs
    #Ci_pro <- (1:length(prob0))[Cii_temp==1] ## return what i of c_i is simulated (proposed ci)
    Ci_pro <- (c(ival,max(ival)+1))[Cii_temp==1] ## return what i of c_i is simulated (proposed ci)
    if(sum(Ci_pro == Cii)==0){ ## new state
      phiC_tp <- rnorm(n=1,mean=mu0,sd=sig0) # simulate a new phiC 
      acep_prop <- min(c(1, pnorm(Y[i]- phiC_tp)/pnorm(Y[i]-phiC[Cii[i]])))
      #cat("\n**",acep_prop,"***\n")
      roll <- rbinom(n=1,size=1,prob=acep_prop)
      #cat("\n**",roll,"***\n")
        Cii[i] <- roll*Ci_pro +(1-roll)*Cii[i] 
    }
    })
  return(list("ci"=Cii,"Phic"=phiC))
} 


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 McMUp_alg5 <- function(Cii,Y,phiC,Nsam,alpha0,parm0,R){
  out_phic <- NULL
  out_cii <- NULL
  Cii_new <- Cii
  phiC_new <- phiC 
  
  for(m in 1:Nsam){
    
    for(i in 1:length(Y)){
      out<- UpdatCi_alg5(i,Cii_new,Y,phiC_new,alpha0,parm0,R)
      Cii_new <- out[["ci"]]
      phiC_new <- out[["phic"]]
    }
    
    c_val <- unique(Cii_new) ## values of phiC non na!
    for(l0 in 1:length(c_val)){Cii_new[Cii_new==c_val[l0]] <- l0}
    c_val <- unique(Cii_new) ## values of phiC non na!
    for(lx in c_val){
      phiC_new[Cii_new==c_val[lx]] <- SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0)}
    
    out_phic[[m]] <- phiC_new
    out_cii <- rbind(out_cii,Cii_new)
  }
  return(list("Phic"=out_phic,"ci"=out_cii))
}

#########################################################################
#####            Algo7
#########################################################################
 ### Simulate a new Phic
 
 SimphiC <- function(c,Cii,Y,phiC,alpha0,parm0){ ## 
   mu0 <- parm0[1]
   sig0 <- parm0[2] 
   Ytp <- Y[Cii==c]
   ni <- length(Ytp)
   rnorm(n=1,mean = ((ni*(sig0^2))*mean(Ytp) + mu0)/(1+ni*(sig0^2)), sd=sqrt(sig0^2/((1+ni*(sig0^2))))) 
 }
 
 ######################################################################### 
  
UpdatCi_alg7_1 <- function(Cii_nw,Y,phiC_nw,alpha0,parm0){ 
  Cii <- Cii_nw
  phiC <- numeric(2*length(phiC_nw))
  phiC[1:n] <- phiC_nw
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  n <- length(Y)
  cp <- 0
  for(i in 1:n){
  if(sum(Cii[i]==Cii[-i]) > 0){## Ci is not a singleton
    ival <- unique(Cii[-i])
    ci_nw <- n+i ## newly created state
    phiC_tp <- rnorm(n=1,mean=mu0,sd=sig0) # simulate a new phiC 
    acep_prop <- min(c(1, alpha0*pnorm(Y[i]- phiC_tp)/((n-1)*pnorm(Y[i]-phiC[Cii[i]]))))
    roll <- rbinom(n=1,size=1,prob=acep_prop)
    if(roll){
    Cii[i] <- ci_nw 
    phiC[n+i]<- phiC_tp}
    
  }
  else{
    ival <- unique(Cii[-i]) ##distinct states except i
    ci_nw <- ival[rmultinom(n=1,size=1,prob=as.vector(table(Cii[-i]))/(n-1))==1] ## chose a value of Cj
    acep_prop <- min(c(1, (n-1)*pnorm(Y[i]- phiC[ci_nw])/(alpha0*pnorm(Y[i]-phiC[Cii[i]]))))
    roll <- rbinom(n=1,size=1,prob=acep_prop)
    Cii[i] <- roll*ci_nw +(1-roll)*Cii[i]
  }
  }
  ival <- unique(Cii)
  phiC <- phiC[Cii]
  for(i in 1:length(ival)){ Cii[Cii == ival[i]] <- i}
  return(list("ci"=Cii,"Phic"=phiC))
} 

###### Algo7
UpdatCi_alg7_2 <- function(Cii_nw,Y,phiC_nw,alpha0,parm0){
  Cii <- Cii_nw
  phiC <- phiC_nw
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  n <- length(Y)
  cp <- 0
  for(i in 1:n){
    if(sum(Cii[i]==Cii[-i])> 0){## Ci is not a singleton
      ival <- unique(Cii[-i])
      prob0 = as.vector(table(Cii[-i]))*pnorm(Y[i]- phiC[ival])/(n-1)
     # cat("***prob0=",prob0,"\n")
     Cii[i] <- ival[rmultinom(n=1,size=1,prob=prob0)==1]
  }
  }
  return(Cii)
}

###%%%%%
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

alpha0 <- 5
Nsam=10000
parm0 <-c(0,2)
R=1
out <- McMUp_alg7(Cii,Y,phiC,Nsam,alpha0,parm0)


##################################################################
###%%%%   Algo. 8
#################################################################
UpdatCi_alg8 <- function(Cii_nw,Y,phiC_nw,alpha0,parm0,M){ 
  n <- length(Y)
  Cii <- Cii_nw
  #phiC <- numeric(n)
  phiC <- phiC_nw
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  
  cp <- 0
  
  for(i in 1:n){
    ival <- unique(Cii[-i])
    K_m <- length(ival)
    h <- K_m + M
    phiC_tp <- numeric(h)
    phiC_tp[1:K_m] <- phiC[ival]
    Cil <- Cii[-i] ## label for the ci
    for(l0 in 1:K_m){ Cil[Cil == ival[l0]] <- l0} ## relabel the Cis
    
    if(sum(Cii[i]== Cil) > 0){## Ci is not a singleton
      phiC_tp[(K_m+1):h] <- rnorm(n=M,mean=mu0,sd=sig0)}
      else{
        Cii[i] <- K_m+1
        phiC_tp[(K_m+2):h] <- rnorm(n=M-1,mean=mu0,sd=sig0)}
    prob0 <- numeric(h)
    prob0[1:K_m] <- as.vector(table(Cii[-i]))*pnorm(Y[i]-phiC_tp[1:K_m])*(1/(n-1+alpha0))
    prob0[(K_m+1):h] <- (alpha0/M/(n-1+alpha0))*pnorm(Y[i]-phiC_tp[(K_m+1):h]) 
    Cii[i] <- (1:h)[rmultinom(n=1,size=1,prob=prob0)==1] 
    Cii[-i]<- Cil
    phiC <- phiC_tp[Cii]
  }
    
  ival <- unique(Cii)
  phiC <- phiC[Cii]
  for(i in 1:length(ival)){ Cii[Cii == ival[i]] <- i}
  return(list("ci"=Cii,"Phic"=phiC))
} 

##%%%%%
McMUp_alg8 <- function(Cii,Y,phiC,Nsam,alpha0,parm0,M){
  out_phic <- NULL
  out_cii <- NULL
  Cii_new <- Cii
  phiC_new <- phiC 
  
  for(m in 1:Nsam){
    #cat("m=",m,"--")
    out<- UpdatCi_alg8(Cii_new,Y,phiC_new,alpha0,parm0,M)
    Cii_new <- out[["ci"]]
    phiC_new <- out[["Phic"]]
    
    c_val <- unique(Cii_new) ## values of phiC non na!
    for(lx in 1:length(c_val)){
      phiC_new[Cii_new==c_val[lx]] <- SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0)}
    
    out_phic <- rbind(out_phic, phiC_new)
    out_cii <- rbind(out_cii,Cii_new)
  }
  return(list("Phic"=out_phic,"ci"=out_cii))
}



#####################################################
n <- 20  # sample size
set.seed(2003)
#Clus <- c(-2,0,2) 
#mu <- rep(Clus,each=floor(n/length(Clus)))
p0 = rmultinom(n,size=1,prob=rep(1,3))
Y <- p0[1,]*rnorm(n,mean=-2) + p0[2,]*rnorm(n,mean=0) + p0[3,]*rnorm(n,mean=2)
hist(Y)
### Algorithm #2
Cii <- 1:length(Y) ## initial states
phiC <- rnorm(length(Y))

alpha0 <- 5
parm0 <- c(0,5)

alpha0 <- 1
Nsam=5000
parm0 <-c(0,2)
M=5
out <- McMUp_alg8(Cii,Y,phiC,Nsam,alpha0,parm0,M)
tail(out$Phic)
tail(out$ci)

