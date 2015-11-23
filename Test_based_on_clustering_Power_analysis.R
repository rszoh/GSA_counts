#!/share/apps/R-latest/bin/Rscript
.libPaths("/data/rszoh/R_LIBS_311")

cmdArgs = commandArgs(trailingOnly = TRUE)
idf = as.numeric(cmdArgs)

set.seed(idf)
#% Simulation Code
#% For the Group 
#%
########################################
source("Code_GSA_NPBT.R")

## Proportion of elements of true mean that are zero

p0 <- c(0,.01,.05,.5)
p <- 200
n1 <- 50
n2 <- 50

### Sigma1 
B <- .85*eye(25,25) + .15*ones(25,25)
Sig1 <- as.matrix(bdiag(B,B,B,B))
Sig1 <- kronecker(eye(8),B)

N = 100 ## number of samples
## We run the simulation for N=100 samples
#setup parallel backend to use 8 processors
cl<-makeCluster(8)
registerDoParallel(cl)

#start time
strt<-Sys.time()

res <- foreach(icount(N),.combine=cbind)%dopar%{
  Nsam <- 20000
  M = 5
  alpha0 <- 1
  parm0 <-c(0,10,.977) ## sig0  = parm0[2], where sig0 is the variance
  phiC_nw <- rep(0,p)
  phiC_new <- rnorm(p)
  Cii <- 1:p
  Sig<- diag(rep(1,p))
  delt0 <- rep(0,p)
  X <- rmnorm(n=n1,mean=rep(0,p),varcov = Sig)
  q <- floor(p*(p0[idf]))
  mu2 <- 5 
  if(q>0){meany <- c(rep(mu2,q),rep(0,p-q))} else{
    meany <- rep(0,p)
  }
  Y <- rmnorm(n=n2,mean=meany,varcov = Sig)
  out1 <- UpdatCi_alg8_KnSig2(X,Y,phiC_new,alpha0,parm0,M=3,Sig)
  if(q > 0){c(p0,mean(rowSums(cbind(out1[,1:q]!=0,out1[,-c(1:q)]==0))==ncol(out1)),colMeans(out1))} else {
    c(p0,mean(rowSums(out1==0)==ncol(out1)),colMeans(out1))
  } 
}


cat("\n Simulation Donce \n")
colnames(res) <- c("P0","Prob_true", paste("parm",1:p,sep=""))

path_maj <-"/data/rszoh/"

dir.create(paste(path_maj,"Test_Based_on_Clustering",sep=""),showWarnings=FALSE)
loc <- paste(path_maj,"Test_Based_on_Clustering",sep="") 
write.csv(res)

cat("\n ** Done saving output ** \n")

