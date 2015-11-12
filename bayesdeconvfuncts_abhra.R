
#############################################################


P.mat <- function(K)
	{
	# penalty matrix for density
	D <- diag(rep(1,K))
	D <- diff(diff(D))
	P <- t(D)%*%D 
	return(P)
	}


#############################################################


B.basis <- function(x,knots)
	{
	delta <- knots[2]-knots[1]
	n <- length(x)
	K <- length(knots)
	B <- matrix(0,n,K+1)
	for (jj in 1:(K-1))
      	{
		act.inds <- (1:n)[(x>=knots[jj])&(x<=knots[jj+1])]
		act.x <- x[act.inds]
		resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
        
		B[act.inds,jj] <- (1/2)*(1-resc.x)^2
		B[act.inds,jj+1] <- -(resc.x^2)+resc.x+1/2
		B[act.inds,jj+2] <- (resc.x^2)/2
		}
	return(B)
	}


#############################################################


rskewnorm <- function(n,mean,sd,skewness)
	{
	c = 2/3.1415926
	delta = skewness/sqrt(1+skewness^2)
	xi = - (sqrt(c) * skewness)/sqrt(1 + skewness^2 * (1-c))
	omega = sqrt(1+xi^2)

	y = delta * abs(rnorm(n,0,1)) + (1-delta^2)^(1/2) * rnorm(n,0,1)
	y = xi + omega * y
	return(mean+sd*y)
	}


dskewnorm <- function(x,mean,sd,skewness)
	{
	c = 2/3.1415926
	delta = skewness/sqrt(1+skewness^2)
	zeta1 = delta * sqrt(c)
	zeta2 = sqrt(1 - c*delta^2)
	y = numeric(length(x))
	xmod = zeta1 + zeta2*(x-mean)/sd
	for(i in 1:length(x))
      	y[i] = (2*zeta2/sd[i]) * dnorm(xmod[i]) * pnorm(skewness*xmod[i])
	return(y)
	}


#############################################################


d.restricted.mix.norm <- function(x,mean,sd,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		d = params[,2]
		mu_curl = params[,3]
		}
	else
		{
		p = params[1]
		d = params[2]
		mu_curl = params[3]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)

    
    	mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	sigma1 = sqrt(d/p - mu1^2)
	sigma2 = sqrt((1-d)/(1-p) - mu2^2)

    	y = p*dnorm(x,mean+sd*mu1,sd*sigma1) + (1-p)*dnorm(x,mean+sd*mu2,sd*sigma2)
	return(y)
	}


r.restricted.mix.norm <- function(n,params)
	{
	p = params[1]
	d = params[2]
	mu_curl = params[3]

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)

    
    	mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	sigma1 = sqrt(d/p - mu1^2)
	sigma2 = sqrt((1-d)/(1-p) - mu2^2)

    	inds = sample(0:1,n,TRUE,prob=c(p,(1-p)))
	y = numeric(n)
	temp = which(inds==0)
	y[temp] = rnorm(length(temp),mu1,sigma1) 
	temp = which(inds==1)
	y[temp] = rnorm(length(temp),mu2,sigma2)
	return(y)
	}


mixnorm <- function(n,pi,p,d,mu_curl,e.grid,plot=TRUE)
	{
	m = length(pi)
	y = numeric(n)
	density <- numeric(length(e.grid))
	inds = sample(1:m,n,TRUE,prob=pi)
	for(ii in 1:m)
		{
		temp = which(inds==ii)
		y[temp] = r.restricted.mix.norm(length(temp),c(p[ii],d[ii],mu_curl[ii]))
		density <- density + pi[ii]*d.restricted.mix.norm(e.grid,mean=0,sd=1,c(p[ii],d[ii],mu_curl[ii]))
		}
	
	if(plot)
		{
		hist(y,xlim=c(min(e.grid),max(e.grid)),breaks=30,freq=FALSE)
		points(e.grid,density,type="l")
		}

	return(list(es=y,density=density))
	}



r.beta.proposal.params.restricted.mix.norm <- function(p.a,p.b,d.a,d.b)
	{
	p = rbeta(1,p.a,p.b)
	d = rbeta(1,d.a,d.b)

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
	
	d1 = sqrt(d/p)
	d2 = sqrt((1-d)/(1-p))
    
	mu.curl.max = min(d1/c1,-d2/c2)
	mu.curl.min = max(-d1/c1,d2/c2)
    
	#mu_curl = runif(mu.curl.min,mu_curl.max)
	mu_curl = runif(1,0,mu.curl.max)
	
	y = c(p,d,mu_curl)
	return(y)
	}    


r.tnorm.proposal.params.restricted.mix.norm <- function(params)
	{
	current.p <- params[1]
	current.d <- params[2]
	current.mu_curl <- params[3]

	p = rtnorm(1,current.p,0.01,lower=0,upper=1)
	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
    
	d = rtnorm(1,current.d,0.01,lower=0,upper=1)
	d1 = sqrt(d/p)
	d2 = sqrt((1-d)/(1-p))
    
	mu.curl.max = min(d1/c1,-d2/c2)
	mu.curl.min = max(-d1/c1,d2/c2)
    
	#mu_curl = runif(mu.curl.min,mu_curl.max)
	mu_curl = runif(1,0,mu.curl.max)

	y = c(p,d,mu_curl)
	return(y)
	}    


rawmoments_mixnorm <- function(order,pi,p,d,mu_curl)
	{
	m = length(pi)
	moment = 0
	for(ii in 1:m)
		{
		c1 = (1-p[ii])/(p[ii]^2+(1-p[ii])^2)
		c2 = -p[ii]/(p[ii]^2+(1-p[ii])^2)

    		mu1 = c1*mu_curl[ii]
		mu2 = c2*mu_curl[ii]

		sigma1 = sqrt(d[ii]/p[ii] - mu1^2)
		sigma2 = sqrt((1-d[ii])/(1-p[ii]) - mu2^2)

		if(order==2)
			moment = moment + pi[ii]*(p[ii]*(mu1^2+sigma1^2)+(1-p[ii])*(mu2^2+sigma2^2))
		else if(order==3)
			moment = moment + pi[ii]*(p[ii]*(mu1^3+3*mu1*sigma1^2)+(1-p[ii])*(mu2^3+3*mu2*sigma2^2))
		else if(order==4)
			moment = moment + pi[ii]*(p[ii]*(mu1^4+6*(mu1^2)*sigma1^2+3*sigma1^4)+(1-p[ii])*(mu2^4+6*(mu2^2)*sigma2^2+3*sigma2^4))
		else
			stop("Only second, third and fourth order raw moments are supported.")
		}
	return(moment)
	}


#############################################################


fr <- function(thetas)			# function to be maximized
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)
	y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum(rep(log(vars),times=mis))/2 - sum(us^2/rep(vars,times=mis))/2
	return(-y)
	}


############################################################


gr <- function(thetas)			# gradient function of fr
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)
	y = - P.t %*% thetas / s2t
	B.basis.components = matrix(0,nrow=K.t+1,ncol=n)
	for(kk in 1:(K.t+1))
		{
   		thetas.new = rep(0,K.t+1)
		thetas.new[kk] = 1
		B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
		}
	for(kk in 1:(K.t+1))
		for(ii in 1:n)
			for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
				y[kk] = y[kk] - (1 - (us[jj]^2)/vars[ii]) * B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
	return(-y)
	}


############################################################


prop.sig.thetas.fn <- function(thetas,s2t)
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)
	prop.sig.thetas = matrix(0,nrow=K.t+1,ncol=K.t+1)
	B.basis.components = matrix(0,nrow=K.t+1,ncol=n);
	for(kk in 1:(K.t+1))
		{
   		thetas.new = rep(0,K.t+1)
		thetas.new[kk] = 1
		B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
		}

	for(kk in 1:(K.t+1))
		{
		for(ll in kk:(K.t+1))
			{
			if(kk==ll)
	  			for(ii in 1:n)
	     				for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
	        				prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2) + (1-(us[jj]^2)/vars[ii])*B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
			else
	  			for(ii in 1:n)
	     				for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
	       				prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2)
			} 
		}
	for(kk in 2:(K.t+1))
   		for(ll in 1:(kk-1))
			prop.sig.thetas[kk,ll] = prop.sig.thetas[ll,kk]
	prop.sig.thetas = prop.sig.thetas + P.t/s2t
	prop.sig.thetas = round(solve(prop.sig.thetas),4)

	#prop.sig.thetas = prop.sig.thetas + n*sum(mis)*P.t/s2t
	#if(is.positive.definite(prop.sig.thetas)==FALSE)
	#	prop.sig.thetas = make.positive.definite(prop.sig.thetas,tol=.00001)
	#prop.sig.thetas = round(solve(prop.sig.thetas),4)
	
	return(prop.sig.thetas)
	}



