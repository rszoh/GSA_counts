
uglyerrors <- function(n,p,mu,sigmasq,e.grid,plot=TRUE)
	{
	es <- c(rnorm(p[1]*n,mu[1],sqrt(sigmasq[1])),rnorm(p[2]*n,mu[2],sqrt(sigmasq[2])),rnorm(p[3]*n,mu[3],sqrt(sigmasq[3])))

	sigmasq.e <- sigmasq/sum(p*(mu^2+sigmasq))
	mu.e <- (mu-sum(p*mu))/sqrt(sum(p*(mu^2+sigmasq)))
	
	es <- (es-sum(p*mu))/sqrt(sum(p*(mu^2+sigmasq)))
	es <- es[sample(1:n,n,TRUE,prob=rep(1,n))]


	density <- ( p[1]*dnorm(e.grid,mu.e[1],sqrt(sigmasq.e[1]))
	 	+ p[2]*dnorm(e.grid,mu.e[2],sqrt(sigmasq.e[2]))
	 	+ p[3]*dnorm(e.grid,mu.e[3],sqrt(sigmasq.e[3]))  )

	if(plot)
		{
		par(mfrow=c(2,1))
		hist(es,xlim=c(min(e.grid),max(e.grid)),breaks=30)
		plot(e.grid,q,type="l")
		par(mfrow=c(1,1))
		}

	return(list(es=es,density=density))
	}



#############################################################


rskewnorm <- function(n,mean,sd,skewness)
	{
	c = 2/pi  
	delta = skewness/sqrt(1+skewness^2)
	xi = - (sqrt(c) * skewness)/sqrt(1 + skewness^2 * (1-c))
	omega = sqrt(1+xi^2)

	y = delta * abs(rnorm(n,0,1)) + (1-delta^2)^(1/2) * rnorm(n,0,1)
	y = xi + omega * y
	return(mean+sd*y)
	}

dskewnorm <- function(x,mean,sd,skewness)
	{
	c = 2/pi
	delta = skewness/sqrt(1+skewness^2)
	zeta1 = delta * sqrt(c)
	zeta2 = sqrt(1 - c*delta^2)
	y = rep(0,length(x))
	xmod = zeta1 + zeta2*(x-mean)/sd
	for(i in 1:length(x))
      	y[i] = (2*zeta2/sd[i]) * dnorm(xmod[i]) * pnorm(skewness*xmod[i])
	return(y)
	}




#############################################################

pi.fn <- function(x,alpha,beta,xstar,K)
	{
	#w <- pnorm(alpha-beta*abs(x-xstar))
	#w <- pnorm(alpha-beta*abs(x-xstar)^1.7)
	w <- pnorm(alpha-beta*abs(x-xstar)^2)
	pi <- numeric(K)
	pi[1] <- w[1]
	for(k in 2:(K-1))
		pi[k] <- w[k]*prod(1-w[1:(k-1)])
	pi[K] <- prod(w[1:(K-1)])
	#pi[K] <- 1-sum(pi[1:(K-1)])
	return(pi)
	}

K <- 8
alpha <- rep(3,(K-1))
beta <- rep(2,(K-1))
xstar <- seq(-3,3,len=(K-1))
round(pi.fn(-1,alpha,beta,xstar,K),3)



#############################################################

likelihood.fn1 <- function(us,params)
	{
	p <- params[1,]
	mu <- params[2,]
	sigma1 <- sqrt(params[3,])
	sigma2 <- sqrt(params[4,])
	
	c1 <- (1-p)/sqrt(p^2+(1-p)^2)
	c2 <- -p/sqrt(p^2+(1-p)^2)

	mu1 <- c1*mu
	mu2 <- c2*mu

	y <- p*dnorm(us,mu1,sigma1) + (1-p)*dnorm(us,mu2,sigma2) 
	return(y)
	}

us <- 2
K <- 10
p <- rbeta(K,1,1)
mu <- rnorm(K,0,3)
sigmasq1 <- 1/rgamma(K,2,2)
sigmasq2 <- 1/rgamma(K,2,2)
params <- rbind(p,mu,sigmasq1,sigmasq2)

likelihood.fn1(us,params)




likelihood.fn2 <- function(us,params)
	{
	p <- params[1]
	mu <- params[2]
	sigma1 <- sqrt(params[3])
	sigma2 <- sqrt(params[4])
	
	c1 <- (1-p)/sqrt(p^2+(1-p)^2)
	c2 <- -p/sqrt(p^2+(1-p)^2)

	mu1 <- c1*mu
	mu2 <- c2*mu

	y <- p*dnorm(us,mu1,sigma1) + (1-p)*dnorm(us,mu2,sigma2) 
	return(y)
	}

us <- c(-2,2)
p <- rbeta(1,1,1)
mu <- rnorm(1,0,3)
sigmasq1 <- 1/rgamma(1,2,2)
sigmasq2 <- 1/rgamma(1,2,2)
params <- rbind(p,mu,sigmasq1,sigmasq2)

likelihood.fn2(us,params)




likelihood.fn3 <- function(ws,xs,params)
	{
	p <- params[1,]
	mu <- params[2,]
	sigma1 <- sqrt(params[3,])
	sigma2 <- sqrt(params[4,])
	
	c1 <- (1-p)/sqrt(p^2+(1-p)^2)
	c2 <- -p/sqrt(p^2+(1-p)^2)

	mu1 <- xs+c1*mu
	mu2 <- xs+c2*mu

	y <- p*dnorm(ws,mu1,sigma1) + (1-p)*dnorm(ws,mu2,sigma2) 
	return(y)
	}

ws <- 2
xs <- 1
K <- 10
p <- rbeta(K,1,1)
mu <- rnorm(K,0,3)
sigmasq1 <- 1/rgamma(K,2,2)
sigmasq2 <- 1/rgamma(K,2,2)
params <- rbind(p,mu,sigmasq1,sigmasq2)

likelihood.fn3(ws,xs,params)



var.fn1 <- function(params)
	{
	p <- params[1,]
	mu <- params[2,]
	sigma1 <- sqrt(params[3,])
	sigma2 <- sqrt(params[4,])
	
	c1 <- (1-p)/sqrt(p^2+(1-p)^2)
	c2 <- -p/sqrt(p^2+(1-p)^2)

	mu1 <- c1*mu
	mu2 <- c2*mu

	y <- p*(mu1^2+sigma1^2) + (1-p)*(mu2^2+sigma2^2) 
	return(y)
	}


#############################################################

plot.likelihood.fn <- function(params)
	{
	us <- seq(-10,10,len=100)
	params <- params

	p <- params[1]
	c1 <- (1-p)/sqrt(p^2+(1-p)^2)
	c2 <- -p/sqrt(p^2+(1-p)^2)

	mu <- params[2]
	mu1 <- c1*mu
	mu2 <- c2*mu

	sigma1 <- sqrt(params[3])
	sigma2 <- sqrt(params[4])

	y <- p*dnorm(us,mu1,sigma1) + (1-p)*dnorm(us,mu2,sigma2) 
	plot(us,y,type="l")
	}

K <- 1
p <- rbeta(K,1,1)
mu <- rnorm(K,0,3)
sigmasq1 <- 1/rgamma(K,2,2)
sigmasq2 <- 1/rgamma(K,2,2)
params <- rbind(p,mu,sigmasq1,sigmasq2)

#plot.likelihood.fn(params)


#############################################################

r.proposal.params.u <- function(params)
	{
	p <- rtnorm(1,mean=params[1],sd=0.1,lower=0,upper=1)
	mu <- rnorm(1,params[2],sd=0.1)
	sigmasq1 <- rtnorm(1,mean=params[3],sd=0.1,lower=0,upper=Inf)
	sigmasq2 <- rtnorm(1,mean=params[4],sd=0.1,lower=0,upper=Inf)
	return(rbind(p,mu,sigmasq1,sigmasq2))
	}

r.proposal.params.u(params[,1])





