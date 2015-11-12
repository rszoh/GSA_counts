#######################################################################
### Algortihm 2 Neal 2000
#######################################################################
###Simulate a new state for the latent variables

SimCi <- function(i,Cii,Y,phiC,alpha0,parm0){
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  n <- length(Y)
  prob0 <- numeric(length(phiC))    ## all possible states
  ival <- (1:length(phiC))[!is.na(phiC)] ## number of states for phiC not na (available) 
  for(k in ival)
  {
    prob0[k] <- (sum(Cii[-i]==k)/(n-1+alpha0))*pnorm((Y[i]-phiC[k])) ## eq 3.6 (a) of Neal(2000)
  }
  ival0 <- (1:length(phiC))[is.na(phiC)] ## get the unoccupied states  
  #prob0[ival0[1]] <- c(prob0, alpha0*pnorm((Y[i]-mu0)/sqrt(1+sig0^2))/(n-1+alpha0)) ## eq 3.6 (a) of Neal(2000)
  prob0[ival0[1]] <-  alpha0*pnorm((Y[i]-mu0)/sqrt(1+sig0^2))/(n-1+alpha0) ## eq 3.6 (a) of Neal(2000)
  prob00 <- prob0[!is.na(prob0)] ## get ride of the unaccopied states
  prob00 <- prob00/sum(prob00) ## standardize the probs
  Cii_temp <- rmultinom(n=1,size=1,prob=prob00) ## standardize the probs
  ((1:length(prob0))[!is.na(prob0)])[Cii_temp==1] ## return what i of c_i is simulated
}

#SimCi(1,Cii,Y,phiC,alpha0,parm0)

### Simulate a new Phic

SimphiC <- function(c,Cii,Y,phiC,alpha0,parm0){ ## 
  mu0 <- parm0[1]
  sig0 <- parm0[2] 
  Ytp <- Y[Cii==c]
  ni <- length(Ytp)
  rnorm(n=1,mean = ((ni*sig0^2)*mean(Ytp) + mu0)/(1+ni*sig0^2), sd=sqrt(sig0^2/((1+ni*sig0^2)))) 
}

#SimphiC(10,Cii,Y,phiC,alpha0,parm0)
####%%%%%%%%% update function

McMUp <- function(Cii,Y,phiC,Nsamp,alpha0,parm0){
  out_phic <- NULL
  out_cii <- NULL
  for(m in 1:Nsamp){
    for(i in 1:length(Y)){
      if(sum(Cii[i]==Cii[-i])==0){
        phiC[Cii[i]] <- NA 
        Cii[i] <- SimCi(i,Cii,Y,phiC,alpha0,parm0) ## Simulate a new Ci
        if(sum(Cii[i]==Cii[-i])==0){ ## if simulated Ci is a new state
          phiC[Cii[i]] <- rnorm(n=1,mean = ((sig0^2)*mean(Y[i]) + mu0)/(1+sig0^2), sd=sqrt(sig0^2/(1+sig0^2))) }  
      }
    }
    c_val <- (1:length(phiC))[!is.na(phiC)] ## values of phiC non na!
    for(lx in c_val){
      phiC[lx] <- SimphiC(lx,Cii,Y,phiC,alpha0,parm0)}
    out_phic <-rbind(out_phic, phiC)
    out_cii <- rbind(out_cii,Cii)
  }
  return(list("Phic"=out_phic,"ci"=out_cii))
}

######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#######################################################################
### Algortihm 5 Neal 2000
#######################################################################
UpdatCi_alg5 <- function(i,Cii_nw,Y,phiC_nw,alpha0,parm0,R){ 
  Cii <- Cii_nw
  phiC <- phiC_nw
  mu0 <- parm0[1]
  sig0 <- parm0[2]
  n <- length(Y)
  cp <- 0
 # ival <- unique(Cii) ### unique state of Ci
     ### reoder phiC 
  #phiC <- unique(phiC)
   ## change the ci to 1, 2, .. the Cii and phiC
  #ival <- (1:length(phiC))[!is.na(phiC)] ## number of states for phiC not na (available) 
  for(lx in 1:R){
    #phiC <- phiC[Cii]
    ival <- unique(Cii)
    for(l0 in 1:length(ival)){Cii[Cii == ival[l0]] <- l0}
    ival <- unique(Cii) ## unique phiC
    prob0 <- numeric(length(ival))    ## all currents states
    for(k in 1:length(ival)){prob0[k] <- sum(Cii[-i]==ival[k])/(n-1+alpha0) ## eq 5.4 (a) of Neal(2000)
    }
    
    cat("\n**",length(Cii),"***\n")
    cat("\n**",sum(!is.na(Cii)),"***")
    prob0 <- c(prob0, alpha0/(n-1+alpha0)) ## eq 5.4 (b) of Neal(2000)
    cat("\n",prob0,"***")
    cat("\n",length(prob0),"***")
    #prob0 <- prob0/sum(prob0) ## standardize the probs
    Cii_temp <- rmultinom(n=1,size=1,prob=prob0) ## standardize the probs
    Ci_pro <- (1:length(prob0))[Cii_temp==1] ## return what i of c_i is simulated (proposed ci)
    cat("\n",length(phiC),"***")
    cat("\n Cii:",Cii,"***")
    cat("\n phiC:",phiC,"***")
    if(sum(Ci_pro == Cii)==0){ ## new state
      phiC_tp <- rnorm(n=1,mean=mu0,sd=sig0) # simulate a new phiC 
      cat("\n phiC_tp:",phiC_tp,"***")
      cat("\n phiC[Cii[i]]:",phiC[Cii[i]],"***")
      cat("\n length(phiC):",length(phiC),"***")
      cat("\n [Cii[i]:",Cii[i],"***")
      acep_prop <- min(c(1, pnorm(Y[i]- phiC_tp)/pnorm(Y[i]-phiC[Cii[i]])))
      cat("\n**",acep_prop,"***\n")
      roll <- rbinom(n=1,size=1,prob=acep_prop)
      cat("\n**",roll,"***\n")
      if(roll){
        Cii[i] <- Ci_pro 
        #phiC <- c(phiC,phiC_tp)}
       ## save the new state
    }
    #else{
    #acep_prop <- min(c(1, pnorm(Y[i]-phiC[ival==Ci_pro])/pnorm(Y[i]-phiC[Cii[i]])))
    #roll <- rbinom(n=1,size=1,prob=acep_prop)
    #Cii[i] <- Ci_pro
    }
  }
  
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
  for(lx in c_val){
    phiC_new <- c(phiC_new,SimphiC(c_val[lx],Cii_new,Y,phiC_new,alpha0,parm0))}
  
    out_phic[[m]] <- phiC_new
     out_cii <- rbind(out_cii,Cii_new)
  }
  return(list("Phic"=out_phic,"ci"=out_cii))
}

#McMUp_alg5(Cii,Y,phiC,Nsam,alpha0,parm0,R)
#######################################################################
### Algortihm 8 Neal 2000
#######################################################################




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(2003)
M <- 1000000 
yi<- rnorm(n=1)
mu0 <- 2
sig0 <- 1
val1 <- rnorm(M,mean=mu0,sd=sig0)
mean(pnorm(yi-val1))
pnorm((yi-mu0)/sqrt(1+sig0^2))

n <- 30  # sample size
set.seed(2003)
Clus <- c(-1,0,1) 
mu <- rep(Clus,each=floor(n/length(Clus)))
Y <- rnorm(n,mean = mu)

### Algorithm #2
Cii <- 1:length(mu) ## initial states
phiC <- rnorm(length(mu))

alpha0 <- 5
parm0 <- c(0,5)

#Nsam=50000
#out <- McMUp(Cii,Y,phiC,Nsamp=Nsam,alpha0,parm0)

dim(out[["ci"]])
sort(unique(Cii))
(out[["Phic"]])[Nsam,]
table(out[["ci"]][Nsam,])

############ 
alpha0 <- 5
Nsam=100
R=1
out5<- McMUp_alg5(Cii,Y,phiC,Nsam,alpha0,parm0,R=1)

out5[["ci"]][100,]
out5[["phiC"]][[100]]
