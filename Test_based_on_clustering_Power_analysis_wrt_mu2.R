#!/share/apps/R-latest/bin/Rscript
.libPaths("/data/rszoh/R_LIBS_311")

cmdArgs = commandArgs(trailingOnly = TRUE)
idf = as.numeric(cmdArgs)

#set.seed(idf)
#% Simulation Code
#% For the Group 
#%
#idf = 1
########################################
path <- "/share_home/rszoh/test_R_script/GSA_counts"
setwd(path)
source("Code_GSA_NPBT.R")

## Proportion of elements of true mean that are zero

p0 <- c(.05)
mu2vl <- sqrt(10)*c(.1,.5,1.5,3,6)
p <- 200
n1 <- 50
n2 <- 50

### Sigma1 
B <- .85*eye(25,25) + .15*ones(25,25)
#Sig <- as.matrix(bdiag(B,B,B,B))
Sig <- kronecker(eye(8),B)

N = 50 ## number of samples
## We run the simulation for N=100 samples
#setup parallel backend to use 8 processors
nc <- min(c(detectCores(),8))
cl <- makeCluster(nc)
registerDoParallel(cl)

#start time
strt<-Sys.time()

res <- foreach(icount(N),.packages=c("matlab","mnormt","Matrix","MCMCpack"),.combine=rbind)%dopar%{
  require(MCMCpack)
  require(matlab)
  require(mnormt)
  require(Matrix)
  
  Nsam <- 20000
  M = 5
  alpha0 <- 1
  parm0 <-c(0,10,.977) ## sig0  = parm0[2], where sig0 is the variance
  phiC_nw <- rep(0,p)
  phiC_new <- rnorm(p)
  Cii <- 1:p
  #Sig<- diag(rep(1,p))
  delt0 <- rep(0,p)
  X <- rmnorm(n=n1,mean=rep(0,p),varcov = Sig)
  q <- floor(p*(p0[idf]))
  mu2 <- mu2vl[idf] 
  if(q>0){meany <- c(rep(mu2,q),rep(0,p-q))} else{
    meany <- rep(0,p)
  }
  Y <- rmnorm(n=n2,mean=meany,varcov = Sig)
  out1 <- McMUp_alg8KnSig2(X,Y,phiC_nw,Nsam,alpha0,parm0,M,Sig) 
  if(q > 0){c(p0,mu2vl[idf],mean(rowSums(cbind(out1[,1:q]!=0,out1[,-c(1:q)]==0))==ncol(out1)),colMeans(out1))} else {
    c(p0,mu2vl[idf],mean(rowSums(out1==0)==ncol(out1)),colMeans(out1))
  } 
}

stopCluster(cl)


cat("\n Simulation Donce \n")
colnames(res) <- c("P0","mu2","Prob_true", paste("parm",1:p,sep=""))

path_maj <-"/data/rszoh/"

dir.create(paste(path_maj,"Test_Based_on_Clustering",sep=""),showWarnings=FALSE)
dir.create(paste(path_maj,"Test_Based_on_Clustering/wrt_mu2",sep=""),showWarnings=FALSE)
loc <- paste(path_maj,"Test_Based_on_Clustering/wrt_mu2",sep="") 
write.csv(res,file=paste(loc,"/Out_",round(mu2vl[idf],2),".csv",sep=""))

cat("\n ** Done saving output ** \n")