# Version 8/4/12
# U_{ij}'s can come from different clusters for j=1,2,...,m_i.



rm(list=ls())




	###################
	### Get Seed No ###
	###################

seedno <- as.integer(Sys.getenv("seedno"))



	###########################
	### Set MCMC Parameters ###
	###########################

simsize <- 5000
burnin <- 3000

simsize.mh.u <- 1





	######################################
	### Add Libraries and Source Files ###
	######################################


# These libraries must be available to R
library(mvtnorm,lib.loc="~/RLibs") 
library(msm,lib.loc="~/RLibs") 

source("bayesdeconvfuncts_PSBP_abhra.R")
source("bayesdeconvfuncts_abhra.R")





	#################################
	### Set Simulation Parameters ###
	#################################



#sample.size.choices <- c(250)
#densityno.choices <- c(1)
#error.true.choices <- c("ugly101")
#no.replicates <- c(3)

sample.size.choices <- c(1000)
densityno.choices <- c(2)
error.true.choices <- c("ugly101","ugly102","ugly103")
no.replicates <- c(3)

#sample.size.choices <- c(250,500,1000)
#densityno.choices <- c(1,2)
#error.true.choices <- c("normal","skewnormal","ugly1","ugly2","ugly3","ugly4","ugly5","ugly6","ugly101","ugly102","ugly103")
#no.replicates <- c(2,3,4,5)

for(reps in no.replicates)
{
for(n in sample.size.choices) # obs
{
for(densityno in densityno.choices)
{
if(densityno==1)
	prop.x.true <- c(1/2,1/2)
if(densityno==2)
	prop.x.true <- c(4/5,1/5)
for(error.true in error.true.choices)
{


mis <- rep(reps,length=n)	# round(2+runif(n)*5) # reps per observation




	###################
	### Set Seed No ###
	###################

set.seed(seedno)






	#####################
	### Generate Data ###
	#####################


# Generate x
xs.true <- c(rnorm(prop.x.true[1]*n,0,.75),rnorm(prop.x.true[2]*n,3,.75))



# Generate es.true, us.true, ws
e.grid <- seq(-4,4,length=500)
if(error.true=="normal")
	{
	ws <- rnorm(sum(mis),mean=rep(xs.true,times=mis),sd=rep(abs(1+(xs.true/4)),times=mis))
	density.e.true <- dnorm(e.grid)
	}
if(error.true=="skewnormal")
	{
	skewness.true <- 7
	es.true <- rskewnorm(sum(mis),0,1,skewness.true)
	density.e.true <- dskewnorm(e.grid,mean=rep(0,times=length(e.grid)),rep(1,times=length(e.grid)),skewness.true)
	us.true <- es.true*rep(abs(1+xs.true/4),times=mis)
	ws <- rep(xs.true,times=mis) + us.true
	}
if(error.true=="ugly1")
	{
	pi= c(1)
	p = c(0.40)
	mu_curl = c(2)
	sigmasq1 = c(2)
	sigmasq2 = c(1)
	
	ugerrors <- mixnorm(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,F)
	es.true <- ugerrors$es
	density.e.true <- ugerrors$density
	us.true <- es.true*rep(abs(1+xs.true/4),times=mis)
	ws <- rep(xs.true,times=mis) + us.true
	}
if(error.true=="ugly2")
	{
	pi= c(1)
	p = c(0.5)
	mu_curl = c(2)
	sigmasq1 = c(1)
	sigmasq2 = c(1)
	
	ugerrors <- mixnorm(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,F)
	es.true <- ugerrors$es
	density.e.true <- ugerrors$density
	us.true <- es.true*rep(abs(1+xs.true/4),times=mis)
	ws <- rep(xs.true,times=mis) + us.true
	}
if(error.true=="ugly3")
	{
	pi= c(0.3,0.7)
	p = c(0.60,0.50)
	mu_curl = c(5,0)
	sigmasq1 = c(1,4)
	sigmasq2 = c(2,1)
	
	ugerrors <- mixnorm(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,F)
	es.true <- ugerrors$es
	density.e.true <- ugerrors$density
	us.true <- es.true*rep(abs(1+xs.true/4),times=mis)
	ws <- rep(xs.true,times=mis) + us.true
	}
if(error.true=="ugly4")
	{
	pi= c(0.3,0.7)
	p = c(0.60,0.50)
	mu_curl = c(0,4)
	sigmasq1 = c(.5,4)
	sigmasq2 = c(.5,4)
	
	ugerrors <- mixnorm(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,F)
	es.true <- ugerrors$es
	density.e.true <- ugerrors$density
	us.true <- es.true*rep(abs(1+xs.true/4),times=mis)
	ws <- rep(xs.true,times=mis) + us.true
	}
if(error.true=="ht_ugly1")
	{
	pi= c(0.8,0.2)
	p = c(0.60,0.50)
	mu_curl = c(0,0)
	sigmasq1 = c(0.25,5)
	sigmasq2 = c(0.25,5)
	
	ugerrors <- mixnorm(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,F)
	es.true <- ugerrors$es
	density.e.true <- ugerrors$density
	us.true <- es.true*rep(abs(1+xs.true/4),times=mis)
	ws <- rep(xs.true,times=mis) + us.true
	}
if(error.true=="ht_ugly2")
	{
	pi = c(1)
	mu = c(0)
	b = c(2)
	
	ugerrors <- mixlaplace(n,pi,mu,b,e.grid,F)
	es.true <- ugerrors$es
	density.e.true <- ugerrors$density
	us.true <- es.true*rep(abs(1+xs.true/4),times=mis)
	ws <- rep(xs.true,times=mis) + us.true
	}




if((error.true=="ugly101")||(error.true=="ugly102")||(error.true=="ugly103"))
	{
	density.e.true <- numeric(length(e.grid))

	max.z.u = 8
	alpha = rep(2,(max.z.u-1))
	beta = rep(0.5,(max.z.u-1))
	xstar = c(-1.9,-1,0,1,2.5,4,5.5)

	var.grid= seq(-2,6,len=500)
	var.quad = (1+var.grid/4)^2
	var.NP = numeric(length(var.grid))

	fr_PSBP <- function(sigma_sqs)			# function to be maximized
		{
		var.quad = (1+var.grid/4)^2
		var.NP = numeric(length(var.grid)) 
		for(ii in 1:length(var.grid))
			var.NP[ii] = var.NP[ii] + sum(pi.fn(var.grid[ii],alpha,beta,xstar,max.z.u)*sigma_sqs)
		y=sum((var.quad-var.NP)^2)
		return(y)
		}	
	sigma_sqs = optim(rep(1,max.z.u), fr_PSBP, method = "BFGS")$par
	}
if(error.true=="ugly101")
	{
	us.true <- numeric(sum(mis))
	for(ii in 1:n)
		{
		prob <- pi.fn(xs.true[ii],alpha,beta,xstar,max.z.u)
		z.u <- sample(1:max.z.u,mis[ii],T,prob)
		us.true[(sum(mis[1:ii-1])+1):sum(mis[1:ii])] <- rnorm(mis[ii],0,sqrt(sigma_sqs[z.u]))
		}
	ws <- rep(xs.true,times=mis)+us.true
	}
if(error.true=="ugly102")
	{
	skewness <- 7
	us.true <- numeric(sum(mis))
	for(ii in 1:n)
		{
		prob <- pi.fn(xs.true[ii],alpha,beta,xstar,max.z.u)
		z.u <- sample(1:max.z.u,mis[ii],T,prob)
		us.true[(sum(mis[1:ii-1])+1):sum(mis[1:ii])] <- rskewnorm(mis[ii],0,sqrt(sigma_sqs[z.u]),skewness)
		}
	ws <- rep(xs.true,times=mis)+us.true
	}
if(error.true=="ugly103")
	{
	skewness <- c(-7,-3,-1,0,1,3,7)
	us.true <- numeric(sum(mis))
	for(ii in 1:n)
		{
		prob <- pi.fn(xs.true[ii],alpha,beta,xstar,max.z.u)
		z.u <- sample(1:max.z.u,mis[ii],T,prob)
		us.true[(sum(mis[1:ii-1])+1):sum(mis[1:ii])] <- rskewnorm(mis[ii],0,sqrt(sigma_sqs[z.u]),skewness[z.u])
		}
	ws <- rep(xs.true,times=mis)+us.true
	}









	#################################
	### Priors and Initial Values ###
	#################################

### Initialization and prior of xs and us

inds <- rep(1:n,times=mis)
wbars <- tapply(ws,inds,"mean")  
xs <- current.xs <- start.xs <- as.vector(wbars)
s2is <- as.vector(tapply(ws,inds,var))
us <- ws - rep(xs,times=mis)
x.grid <- seq(min(start.xs),max(start.xs),length=100)

alpha.x=0.1
if(n==250)
	alpha.x=0.1
if(n==500)
	alpha.x=0.01
if(n==1000)
	alpha.x=0.001


# Normal-Inverse-Gamma
#nu0.x = 1/5
#gama0.x = 3                   # Note the spelling, also not that this must be greater than 2
#mu0.x = 2
#sigmasq0.x = var(xs)*(gama0.x-1)/(1+1/nu0.x)
#student.x.grid = sqrt(gama0.x)*(x.grid-mu0.x)/(sqrt(sigmasq0.x)*sqrt(1+1/nu0.x))

# Normal-Inverse-Scaled-Chisquare
nu0.x = 1/5
gama0.x = 5	                  		# Note the spelling, also not that this must be greater than 2
mu0.x = 2
sigmasq0.x = 1*(gama0.x-2)/gama0.x		# var(xs)/((1+1/nu0.x)*(gama0.x/(gama0.x-2)))
student.x.grid = (x.grid-mu0.x)/(sqrt(sigmasq0.x)*sqrt(1+1/nu0.x))

z.x <- 1:n	
mu.x <- xs						# unique values
sigmasq.x <- rep(var(xs)/5,n)			# unique values



var.grid <- seq(-2,5,length=100)



### Initialization and prior of PSBP parameters

z.u <- rep(1,n)
max.z.u <- 10

mu0.alpha <- 1
sigmasq0.alpha <- 2
alpha <- rep(mu0.alpha,max.z.u-1)

mu0.beta <- 2
sigmasq0.beta <- 1
beta <- rep(mu0.beta,max.z.u-1)

xstarspace <- seq(min(wbars),max(wbars),len=50)
xstar <- sample(xstarspace,max.z.u-1,FALSE,rep(1,length(xstarspace)))
xstar[1:5] <- seq(min(wbars)+0.5,max(wbars)-0.5,len=5)

p <- rep(0.5,max.z.u)
mu <- rep(0,max.z.u)
sigmasq1 <- c(seq(0.1,2,len=5),seq(0.1,2,len=max.z.u-5))
sigmasq2 <- c(seq(0.1,2,len=5),seq(0.1,2,len=max.z.u-5))
params.u <- rbind(p,mu,sigmasq1,sigmasq2)





	###############################
	### Storage for MCMC Output ###
	###############################


density.x.est <- numeric(length(x.grid))
var.est <- numeric(length(var.grid))

accepted.params.u <- matrix(0,nrow=simsize,ncol=max.z.u)
accepted.xstar <- matrix(0,nrow=simsize,ncol=max.z.u-1)




	##################
	### Start MCMC ###
	##################



for (iii in 1:simsize)
	{
	if(iii%%1000==0)
		print(iii)


	### Updating z.x
	
	#student.x = sqrt(gama0.x)*(xs-mu0.x)/(sqrt(sigmasq0.x)*sqrt(1+1/nu0.x))
	#marginal.tpdfs = dt(student.x,df=2*gama0.x)
	
	student.x = (xs-mu0.x)/(sqrt(sigmasq0.x)*sqrt(1+1/nu0.x))
	marginal.tpdfs = dt(student.x,df=gama0.x)
	
	for(ii in 1:n)
		{
		prob.x = tabulate(z.x) 
        	k.x = length(unique(z.x[-ii]))
        	if((prob.x[z.x[ii]]==1)&&(z.x[ii]<max(z.x)))    # z.x[ii] appears only once in z.x
			{
			temp.mu.x = mu.x[z.x[ii]]
			temp.sigmasq.x = sigmasq.x[z.x[ii]]
			for(jj in z.x[ii]:k.x)
				{
				mu.x[jj] = mu.x[jj+1]
				sigmasq.x[jj] = sigmasq.x[jj+1]
				}
			z.x[z.x > z.x[ii]] = z.x[z.x > z.x[ii]]-1
			z.x[ii] = k.x+1         # Necessary steps
			mu.x[k.x+1] = temp.mu.x
			sigmasq.x[k.x+1] = temp.sigmasq.x  
			}
        	prob.minus.x = tabulate(z.x[-ii])
		for(kk in 1:k.x)
			{
			likelihood = dnorm(xs[ii],mu.x[kk],sqrt(sigmasq.x[kk]))
			prob.x[kk] = prob.minus.x[kk] * likelihood
			}
		prob.x[k.x+1] = alpha.x * marginal.tpdfs[ii]
		prob.x = prob.x/sum(prob.x)    		# NOT really necessary        
		z.x[ii] = sample(k.x+1,1,TRUE,prob.x)   	# New z.x[ii] drawn
		if(z.x[ii]==(k.x+1))
			{
			#sigmasq.x[k.x+1] = 1/rgamma(1,shape=gama0.x,rate=sigmasq0.x)			# Normal-Inverse-Gamma
			sigmasq.x[k.x+1] = 1/rgamma(1,shape=gama0.x/2,rate=gama0.x*sigmasq0.x/2)	# Normal-Inverse-Scaled-Chisquare
			mu.x[k.x+1] = rnorm(1,mu0.x,sqrt(sigmasq.x[k.x+1]/nu0.x))
			}
		}
    
    
	### Updating mu.x, sigmsq.x
	
	k.x = max(z.x)                # Number of clusters
	for(kk in 1:k.x)
		{
		temp = which(z.x==kk)
		xspool = xs[temp]

		#nutemp.x = nu0.x + length(xspool)
		#mutemp.x = (nu0.x*mu0.x + sum(xspool)) /(nutemp.x)
		#gamatemp.x = gama0.x + length(xspool)/2				# Normal-Inverse-Gamma
		#sigmasqtemp.x = sigmasq0.x + (1/2)*(sum(xspool^2) + nu0.x*mu0.x^2 - nutemp.x*mutemp.x^2)
		#sigmasq.x[kk] = 1/rgamma(1,shape=gamatemp.x,rate=sigmasqtemp.x)        
		#mu.x[kk] = rnorm(1,mutemp.x,sqrt(sigmasq.x[kk]/nutemp.x))

		nutemp.x = nu0.x + length(xspool)
		mutemp.x = (nu0.x*mu0.x + sum(xspool)) / (nutemp.x)
		gamatemp.x = gama0.x + length(xspool)				# Normal-Inverse-Scaled-Chisquare
		sigmasqtemp.x = (gama0.x*sigmasq0.x + sum(xspool^2) + nu0.x*mu0.x^2 - nutemp.x*mutemp.x^2) / gamatemp.x	
		sigmasq.x[kk] = 1/rgamma(1,shape=gamatemp.x/2,rate=(gamatemp.x*sigmasqtemp.x)/2)
		mu.x[kk] = rnorm(1,mutemp.x,sqrt(sigmasq.x[kk]/nutemp.x))
		}




	### Updating xs (and us)
	prob.x = tabulate(z.x) 
      if(T)
	{
	proposed.xs = rnorm(n,current.xs,sd=diff(range(start.xs))/6)	
	proposed.us = ws - rep(proposed.xs,times=mis)
	proposed.prior = current.prior = numeric(n)
	for(kk in 1:k.x)
		{
		proposed.prior = proposed.prior + prob.x[kk]*dnorm(proposed.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
		current.prior = current.prior + prob.x[kk]*dnorm(current.xs,mu.x[kk],sqrt(sigmasq.x[kk]))
		}
	mh.ratio = proposed.prior/current.prior
	current.likelihood = proposed.likelihood = numeric(n)
	for(ii in 1:n)
		{
		temp = (sum(mis[1:ii-1])+1):sum(mis[1:ii])

		current.likelihood[ii] = prod(likelihood.fn2(us[temp],params.u[,z.u[ii]]))
		proposed.likelihood[ii] = prod(likelihood.fn2(proposed.us[temp],params.u[,z.u[ii]]))
		}
	mh.ratio = mh.ratio * proposed.likelihood/current.likelihood	

	u = runif(n)
	inds.to.replace =(1:n)[u<mh.ratio]
	xs[inds.to.replace] = current.xs[inds.to.replace] = proposed.xs[inds.to.replace]
	
	us = ws - rep(xs,times=mis)
	}






	### Updating z.u
	if(T)
	{
	for(jj in 1:sum(mis))
		{
		likelihoods.u <- likelihood.fn1(us[jj],params.u)
		prob.u <- pi.fn(xs[inds[jj]],alpha,beta,xstar,max.z.u) * likelihoods.u
		z.u[jj] <- sample(1:max.z.u,1,TRUE,prob.u)
		}



	### Updating params.u


	for(rr in 1:simsize.mh.u)
		{
		for(kk in 1:max.z.u)
			{
			temp = which(z.u==kk)
			uspool = us[temp]
			
			proposed.params.u <- r.proposal.params.u(params.u[,kk])
			proposed.log.likelihood = sum(log(likelihood.fn2(uspool,proposed.params.u)))
			current.log.likelihood = sum(log(likelihood.fn2(uspool,params.u[,kk])))

	
			log.acc.prob <- proposed.log.likelihood-current.log.likelihood
			if(log(runif(1))<log.acc.prob)
				{
				params.u[,kk] <- proposed.params.u
				accepted.params.u[iii,kk] <- 1
				}
			} 
		}


	}




	
	### Updating auxiliary Z.u	# CAUTION: Z.u and z.u
	if(T)
	{
	Z.u <- matrix(0,nrow=sum(mis),ncol=max.z.u-1)
	for(jj in 1:sum(mis))
		{
		#print(jj)
		KK <- z.u[jj]
		if(KK>1){
			for(kk in 1:(KK-1))
				Z.u[jj,kk] <- rtnorm(1,mean=alpha[kk]-beta[kk]*(xs[inds[jj]]-xstar[kk])^2,sd=1,lower=-Inf,upper=0)
			}
		if(KK<max.z.u){
			Z.u[jj,KK] <- rtnorm(1,mean=alpha[KK]-beta[KK]*(xs[inds[jj]]-xstar[KK])^2,sd=1,lower=0,upper=Inf)
			}
		}
	 

	### Updating alpha, beta, xstar
	for(kk in 1:(max.z.u-1))
		{
		#print(kk)
	
		temp = which(z.u>=kk)

		post.sigmasq.alpha = 1/(1/sigmasq0.alpha + length(temp)) 
		post.mu.alpha = mu0.alpha/sigmasq0.alpha + sum(Z.u[temp,kk]+beta[kk]*(xs[inds[temp]]-xstar[kk])^2)
		post.mu.alpha = post.mu.alpha * post.sigmasq.alpha
		alpha[kk] = rtnorm(1,mean=post.mu.alpha,sd=sqrt(post.sigmasq.alpha),lower=0,upper=Inf)

		post.sigmasq.beta = 1/(1/sigmasq0.beta + sum((xs[inds[temp]]-xstar[kk])^4))
		post.mu.beta = mu0.beta/sigmasq0.beta - sum((Z.u[temp,kk]-alpha[kk])*(xs[inds[temp]]-xstar[kk])^2)
		post.mu.beta = post.mu.beta * post.sigmasq.beta
		beta[kk] = rtnorm(1,mean=post.mu.beta,sd=sqrt(post.sigmasq.beta),lower=0,upper=Inf)

		log.prob.xstar = numeric(length(xstarspace))
		for(jj in 1:length(xstarspace))
			log.prob.xstar[jj] = sum(log(dnorm(Z.u[temp,kk],alpha[kk]-beta[kk]*(xs[inds[temp]]-xstarspace[jj])^2,sd=1)))
		if(max(log.prob.xstar)>-Inf)
			log.prob.xstar = log.prob.xstar - max(log.prob.xstar)
		#print(log.prob.xstar)
		if(sum(exp(log.prob.xstar))>0)
			{
			xstar[kk] = sample(xstarspace,1,TRUE,exp(log.prob.xstar))
			accepted.xstar[iii,kk] = 1
			#print("xstar accepted")
			}
		}
	}




	if(iii>burnin)
		{
		prob.x = tabulate(z.x)/n
        	for(kk in 1:max(z.x))
			density.x.est = density.x.est + prob.x[kk]*dnorm(x.grid,mu.x[kk],sqrt(sigmasq.x[kk]))

		for(ii in 1:length(var.grid))
			var.est[ii] = var.est[ii] + sum(pi.fn(var.grid[ii],alpha,beta,xstar,max.z.u)*var.fn1(params.u))
		}

	}





density.x.true <- prop.x.true[1]*dnorm(x.grid,0,.75)+prop.x.true[2]*dnorm(x.grid,3,.75)
density.x.est <- density.x.est/(simsize-burnin)
var.est <- var.est/(simsize-burnin)

var.true <- (1+var.grid/4)^2



if(F)
	{
	par(mfrow=c(3,1))

	plot(x.grid,density.x.true,ylim=c(0,max(density.x.true,density.x.est)),type="l")
	points(x.grid,density.x.est,type="l",lty=2)


	plot(xs.true,s2is,pch="*")
	points(var.grid,var.true,ylim=c(0,max(var.true,var.est)),type="l")
	points(var.grid,var.est,type="l",lty=2)

	plot(xs.true,xs)
	par(mfrow=c(1,1))
	}



filename <- paste("MISEs",error.true,"PSBP",n,reps,densityno,sep="_")
filename <- paste(filename,".RData",sep="")
load(filename)


VarFns <- rbind(VarFns,c(seedno,var.est)) 
Densities_x <- rbind(Densities_x,c(seedno,density.x.est))


delta.x <- x.grid[2] - x.grid[1]
MISE <- sum((density.x.est-density.x.true)^2)*delta.x


MISEs <-  rbind(MISEs,c(seedno,MISE))
MISEs[,2] <- MISEs[order(MISEs[,1]),2]
MISEs[,1] <- sort(MISEs[,1])


save(MISEs,VarFns,Densities_x,file=filename)


}	# loop for error.true.choices ends here
}	# loop for densityno.choices ends here
}	# loop for sample.size.choices ends here
}	# loop for no.replicates ends here







