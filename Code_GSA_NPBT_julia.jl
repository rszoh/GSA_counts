#module SimModule
using Distributions
#using NumericExtensions
using Devectorize
#using Permutations
#using GramSchmidt
#using IterativeSolvers
using Iterators
#using Gadfly
using DataFrames



##function Loglikhodfunc(val,i,X, Y,Sig,delt0)
##    n1 = size(X,1) ## number of rows
# #   n2 = size(Y,1) ## number of rows
# #   p = size(Y,2)  ## number of variables
# #   delt = delt0
#    delt[i] = delt0
#    mud = zeros(p)
#    @devec mud = mean(Y,1) - mean(X,2) - delt'
#    tp = inv(Sig)*(n2.*At_mul_B(mud,mud) + At_mul_B(mean(Y,1),mean(Y,1)) + At_mul_B(mean(X,1),mean(X,1)))
#    out = -.5*(n1+n2-1)*(logdet(Sig)) - .5*sum(diag(tp))
#    return out
#end


function Loglikhodfunc_knSig(val, i, X, Y, Sig, delt0) ### known sigma
    n1 = size(X,2) ## number of col
    n2 = size(Y,2) ## number of col
    p = size(Y,1)  ## number of variables p
    delt1 = copy(delt0)
    delt = zeros(p,1)
    delt = reshape(delt1,p,1)
    delt[i] = val
    mud = zeros(p)
    @devec mud = mean(Y,2) - mean(X,2) - delt ## vector of length p
    n0  = 1/(1/n1 + 1/n2)
   # return -.5*trace(n0.*inv(Sig)*A_mul_Bt(mud,mud))
   return -.5*trace(n0.*At_ldiv_B(Sig,mud)*mud')
end

function UpdatCi_alg8_KnSig2(X, Y, phiC_nw, alpha0, parm0, M, Sig)
   p = size(Y,1) ## number of variables p
   phiC = zeros(p)
   phiC0 = zeros(p)
    phiC = copy(phiC_nw)
    phiC0 = copy(phiC_nw) ## copy of vector phic
   mu0 = parm0[1]
   sig0 = parm0[2]
   p0 = parm0[3]
   cp = 0
   K_m = 0
   h = 0
   for i in 1:p
       deleteat!(phiC,i) ## delete the ith element of phiC
       ival = unique(phiC)
       K_m = length(ival)
       h = K_m + M
       phiC_tp = zeros(h)
       prob0 = zeros(h)
        phiC_tp[1:K_m] = vec(ival)
         
       if sum(phiC .== phiC0[i]) > 0
           phiC_tp[(K_m+1):h] = rand(Binomial(1, 1 - p0),M).*rand(Normal(mu0,sqrt(sig0)),M)
       else
           phiC_tp[(K_m+1):h] = [phiC0[i], rand(Binomial(1, 1 - p0),M-1).*rand(Normal(mu0,sqrt(sig0)),M-1)]
       end
             
       for l0 in 1:K_m
           prob0[l0] = Loglikhodfunc_knSig(phiC_tp[l0],i,X, Y,Sig,phiC0) + log(sum(phiC .== ival[l0])/(p-1+alpha0))
       end
       for l0 in (K_m+1):h
           prob0[l0] = log(alpha0) - log(M) - log(p-1+alpha0) + Loglikhodfunc_knSig(phiC_tp[l0],i,X, Y,Sig,phiC0)
       end
      # @show sum(exp(prob0 - maximum(prob0))./ sum(exp(prob0 - maximum(prob0))))
      #@devec prob0 = exp(prob0 - maximum(prob0))./ sum(exp(prob0 - maximum(prob0)))
      #@show sum(prob0)
       phiC0[i] = phiC_tp[vec(rand(Multinomial(1,exp(prob0 - maximum(prob0))./ sum(exp(prob0 - maximum(prob0)))),1) .== 1)][1,1]
        #append!(phiC,0.0)
         phiC = copy(phiC0)
   end
    idf = collect(1:length(phiC0))[phiC0 .== 0.0]
    if length(idf) > 0
        phiC0[idf] = rep(0.0,length(idf))
    end
   return phiC0
end

function Post_delt_KnSig(delti, X, Y, parm0, Sig0, delt0)
    n1 = size(X,2) ## number of col
    n2 = size(Y,2) ## number of col
    p = size(Y,1)  ## number of variables
    delt = copy(delt0)
    #@show delt length(delt)
    n0 = 1/(1/n1 + 1/n2)
    m0 = copy(parm0[1])
    tau0 = copy(parm0[2])
    p0 = copy(parm0[3])
    Sig = copy(Sig0)
    scale!(Sig,1/n0)
    #id0 = [1:length(delt)][delt .== delti] ## where are the delti identical
    id0 = collect(1:length(delt))[delt .== delti] ## where are the delti identical
    if length(id0) == 0
        error("Something is wrong ")
    end
    D = zeros(p)
     @devec D = mean(Y,2) - mean(X,2) ### matrix of 1 by p 
    #invSig = zeros(p,p)
    if length(id0) < length(delt)
        #id1 = [1:length(delt)][delt .!= delti]
        id1 = collect(1:length(delt))[delt .!= delti]
        #invSig = inv(Sig[id1,id1])
        mu0 = mean(D[id0,1] + A_rdiv_Bt(Sig[id0,id1],Sig[id1,id1])*(delt[id1] - D[id1,1]))
        on00 = ones(length(id0))./length(id0)
        Sig00 = (on00'*(Sig[id0,id0] - A_rdiv_Bt(Sig[id0,id1],Sig[id1,id1])*Sig[id1,id0])*on00)[1,1]
      else
        #D = mean(Y,2)' - mean(X,2)' ### matrix of 1 by p
        mu0 = mean(D)
        on00 = ones(length(id0))./length(id0)
        Sig00 = (on00'*Sig[id0,id0]*on00)
    end
       Sigs = 1.0 ./(1.0/tau0 + 1.0/Sig00[1])
       mus = mu0*Sigs/Sig00[1]
       #mus p0 pdf(Normal(0,sqrt(tau0)),0.0) pdf(Normal(mus,sqrt(Sigs)), 0.0)
       pis = 1.0 ./(1.0 + (1.0-p0)*pdf(Normal(0, sqrt(tau0)),0.0)/(p0*pdf(Normal(mus, sqrt(Sigs)),0.0)))
       out = (1.0 - rand(Binomial(1, pis),1)).*rand(Normal(mus,sqrt(Sigs)),1)
       return out == 0.0 ? abs(out) : abs(out)
end
      
     
function  McMUp_alg8KnSig2(X, Y, phiC, Nsam, alpha0, parm0, M, Sig0)
   out_phiC = zeros(0)
    n1 = size(X,2) ## number of col
    n2 = size(Y,2) ## number of col
    p = size(Y,1)  ## number of variables
    phiC_new0 = copy(phiC) 
    Sig = copy(Sig0)
    phiC_new00 = zeros(p)
    phiC_new01 = zeros(p)
    out_res = zeros(Nsam, p)
    cpt = Nsam/10
    println("Start MCMC")
   for m in 1:Nsam
    phiC_new = UpdatCi_alg8_KnSig2(X,Y,phiC_new0,alpha0,parm0,M,Sig) ## store the new Phics
    @devec phiC_new00 = phiC_new  ## assign the new PhiC to two vectors
    @devec phiC_new01 = phiC_new  ## assign the new PhiC to two vectors
    ival = unique(phiC_new) ## find the unique values of Phic
     for lx in 1:length(ival)
      #id00 = [1:length(phiC_new01)][phiC_new0 .== ival[lx]]
      id00 = collect(1:length(phiC_new))[phiC_new .== ival[lx]]
      phiC_new00[id00] = rep(Post_delt_KnSig(ival[lx],X,Y,parm0,Sig,phiC_new)[1],length(id00))
     end
     @devec phiC_new = phiC_new00
     append!(out_phiC, phiC_new)
     #if m%cpt ==0 
     #println("Done with Iter -- $m")
     #end
   end
    out_res = reshape(out_phiC,p,Nsam)'
    println("Done MCMC")
    return out_res
end

function SIM_data(parm0) ## parm contains important info to simulate data and other things!
  @devec parm = parm0
  const  n1 = floor(Integer, parm[1])
   const n2 = floor(Integer,parm[2])
   const p = floor(Integer,parm[3])
  const  p0 = parm[4]
   const mu2 = parm[5]
  const Nsam = floor(Integer, parm[6])
   const tau0 = parm[7]
   const m0 = parm[8]
   const alpha0 = parm[9]
   const M = floor(Integer,parm[10])
   
   # println("Form true Sig  -- ")
    B = .85*eye(25,25) + .15*ones(25,25)
    Sig0 = kron(eye(8),B)
    Sig0 = eye(p) 
    # println("Simulate data")
     X0 = rand(MvNormal(zeros(p), Sig0), n1)
     q = floor(Integer, p0*p) ## number of non-zeros entries
     if q > 0
      meany = zeros(p)
      meany[1:q] = rep(mu2,q)
      Y0 = rand(MvNormal(meany, Sig0), n2)
  else
      meany = zeros(p)
      Y0 = rand(MvNormal(meany, Sig0), n2)
  end
  phiC_nw0 = zeros(p)
  phiC_new0 = randn(p)
  
  out1  =  McMUp_alg8KnSig2(X0,Y0,phiC_nw0,Nsam,alpha0,[m0, tau0, alpha0],M,Sig0)
  #println("MCMC Done!") 
  if q == 0
   ot = [parm mean(prod(out1 .== 0.0,2)) mean(out1,1)]
  else
    ot = [parm mean(prod(out1[:,1:q] .!= 0.0,2).*prod(out1[:,(q+1):end] .== 0.0, 2))  mean(out1,1)]
  end
  return ot
end
                                                                        
           