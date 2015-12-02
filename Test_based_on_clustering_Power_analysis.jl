 blas_set_num_threads(2)
 #cd("U:\\Roger_Abhra_project")
include("Code_GSA_NPBT_julia.jl")

N = 50 ## number of samples
Nsam = 1000   ## number of samples
mu0 = 0.0
tau0 = 10 ## prior variance for delt_i
p0 = [.01,.05,.5]
mu2vl = sqrt(tau0).*[.1,.5,1.5,3.0,6.0] ### where c(.1,.5,1.5,3,6) are the SNR 
p  = 200
n1 = 50
n2 = 50
M = 5
alpha0 = 1.0
parm0 = hcat(map(x->vcat(x...),product(n1,n2,p,p0,mu2vl,Nsam,tau0,mu0,alpha0,M))...)'  ### all parameters combinasion

parm = repmat(vcat([n1 n2 p 0.0 0.0 Nsam tau0 mu0 alpha0 M],parm0),3,1)

parm = parm[1:3,:]
#@time out = SIM_data(parm[1,:])
error("Stop here")
#path00 = pwd()
#push!(LOAD_PATH,path00)

size(parm)
  println("Start here ... ")
  println("path $(pwd())")
  nprocs0 = min(3,size(parm,1))
  addprocs(nprocs0)
  Val = { parm[i,:] for i in 1:size(parm,1)}
  # Val = Any[parm[i,:] for i in 1:size(parm,1)]
  #@everywhere include("Code_GSA_NPBT_julia.jl")
  require("Code_GSA_NPBT_julia.jl")
  #using SimModule
  ot = pmap(SIM_data, Val)
  rmprocs(workers())
  ot 
 outf=ot[1]
 for i in 2:(size(ot,1))
 outf = vcat(outf,ot[i])
end
println("It is all done")
 data = convert(DataFrame,outf)

 #names!(data,[:n1,:n2,:p,:p0,:mu2,:Nsam,:tau0,:mu0,:alpha0,:m,:Log_bay2,:P_val,:P_val2])
 #mkdir("/data/rszoh/Test_Based_on_Clustering/Julia_output/")
 #writetable("/data/rszoh/Test_Based_on_Clustering/Julia_output/output_11_30_15.csv",data)
 writetable("/data/rszoh/Test_Based_on_Clustering/Julia_output/output_11_30_15.csv",data)
 
 println("Done and Done !!")

error("Stop here")

 
 Nsam = 100
  M = 5
  alpha0 = 1
  p = 30
  n1 = 10 
  n2 = 10
  parm00 = [0,10,.977] ## sig0  = parm0[2], where sig0 is the variance
  phiC_nw0 = rep(0,p)
  phiC_new0 = randn(p)
  Cii0  = 1:p
  Sig0 = eye(p)
  delt00 = zeros(p,1)
  X0 = rand(MvNormal(zeros(p), Sig0), n1)
  q = 2
    mu2 = 5 
  meany = zeros(p)
  meany[1:q] = rep(mu2,q)
  Y0 = rand(MvNormal(meany, Sig0), n2)
  
  out = Loglikhodfunc_knSig(1,1,X0,Y0,Sig0,delt00)
  
  out0 = UpdatCi_alg8_KnSig2(X0,Y0,phiC_new0,alpha0,parm00,M,Sig0)
  
  out0
  Post_delt_KnSig(out0[1],X0,Y0,parm00,Sig0,out0)
  
  out1  =  McMUp_alg8KnSig2(X0,Y0,out0,Nsam,alpha0,parm00,M,Sig0)
  
  mean(rowSums(out1==0)==ncol(out1))
  
function essai(x::Float64)
 if x < 0
  println("you enter x_neg= $x")
  println("this is a negative number")
  y = x+100
  println("this is Y= $y")
 else
  println("you enter x_pos= $x")
  println("this is a positive number")
  y = x+100
  println("this is Y= $y")
 end
 return x
end