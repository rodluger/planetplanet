using PyPlot
using Optim
using JLD

#PyPlot.matplotlib[:rc]("text", usetex=true) # allow tex rendering
#PyPlot.matplotlib[:rcParams]["text.latex.unicode"] = true
#PyPlot.matplotlib[:rc]("text.latex",preamble=L"\usepackage{amsmath}")

@load "MC_circular_1000000sig_rho_sigtheta.jld"

nsig_theta = 20
sig_theta1 = 0.0001
sig_theta2 = 0.05

dsig_theta = sig_theta * log(10.) * log10(sig_theta2/sig_theta1)/nsig_theta

sig_ecc1 = 0.0
sig_ecc2 = 0.0
nsig_ecc = 1

dsig_ecc = 1.0

rho_star1 = 48.0
rho_star2 = 54.0
nrho = 30

drho = (rho_star2-rho_star1)/nrho

prob_theta_sum = zeros(nsig_theta)
prob_rho_sum = zeros(nrho)
  for j=1:nsig_theta
  for k=1:nsig_ecc
    prob_theta_sum[j] += drho*sum(prob[j,:,k])
    for m=1:nrho
      prob_rho_sum[m] +=  dsig_theta[j]*sum(prob[j,m,k])
    end
  end
  end

axis([1e-3,.2,0,.5])

#function generalized_logistic(t,param)
#A = param[1]
#B = param[2]
#C = param[3]
#nu = param[4]
#logecc = B.*log10(t./abs(C))
#return A.*(1.-logecc./(1.+abs(logecc).^nu).^(1./nu))
#end


#function chisquare(param)
#return sum((generalized_logistic(sig_ecc,param)-p).^2)
#end

#param = [.24,3.7,.019,2.25]

#result = optimize(chisquare,param)

#param = result.minimizer

using PyPlot

clf()
prob_theta_sum /= maximum(prob_theta_sum)
sig_theta *= 180./pi
semilogx(sig_theta,prob_theta_sum,linewidth=3,alpha=0.5)
cum_prob = zeros(nsig_theta)
cum_prob[1] = dsig_theta[1]*.5*prob_theta_sum[1]
for i=2:nsig_theta
  cum_prob[i] = cum_prob[i-1]+dsig_theta[i]*prob_theta_sum[i]
end
scatter(sig_theta,prob_theta_sum,s=64)
plot(sig_theta,cum_prob/maximum(cum_prob),linewidth=3,alpha=0.5)

axis([6e-3,5,0,maximum(prob_theta_sum)*1.05])
#xlabel(L"\boldsymbol{\sigma_\theta {[\rm deg]}}}",fontsize=20)
#xlabel(L"\boldsymbol{\sigma}_\theta {[\rm deg]}}",fontsize=20)
#xlabel(L"\sigma_\theta {[\rm deg]}}",fontsize=20)
xlabel("Polar angle scatter [deg]",fontsize=20,fontweight="bold")
ylabel("Relative probability",fontsize=20,fontweight="bold")
tick_params("both",labelsize=20)

tight_layout()
savefig("sigma_theta_circular.pdf", bbox_inches="tight")

read(STDIN,Char)
clf()
scatter(rho_star,prob_rho_sum,s=64)
rho_grid = linspace(45,56,1000)

#plot(rho_grid,maximum(prob_rho_sum)*(exp(-(rho_grid-51.15).^2*.5/.6^2)+.03*exp(-(rho_grid-51.3).^2*.5/2.^2)))
#plot(rho_grid,maximum(prob_rho_sum)*exp(-(rho_grid-51.25).^2*.5/.65^2))

#xlabel(L"\boldsymbol{\rho_* [\rho_\odot]}",fontsize=20)
#xlabel(L"\mathbf{\rho_* [\rho_\odot]}",fontsize=20)
xlabel("Stellar density [Solar density]",fontsize=20,fontweight="bold")
ylabel("Relative probability",fontsize=20,fontweight="bold")
axis([45,56,0,1.0])
tick_params("both",labelsize=20)

param2=zeros(6)
param2[1]=.65*maximum(prob_rho_sum)
param2[2]=51.0
param2[3]=0.67
param2[4]=.35*maximum(prob_rho_sum)
param2[5]=51.0
param2[6]=1.5

function chi_normal2(param)
chi = 0.0
for i=1:nrho
  chi += (param[1]*exp(-(rho_star[i]-param[2])^2*.5/param[3]^2)
         +param[4]*exp(-(rho_star[i]-param[5])^2*.5/param[6]^2)-prob_rho_sum[i])^2
end
return chi 
end

function chi_normal(param)
chi = 0.0
for i=1:nrho
  chi += (param[1]*exp(-(rho_star[i]-param[2])^2*.5/param[3]^2)-prob_rho_sum[i])^2
end
return chi 
end

#result2 = optimize(chi_normal2,param2)
param=zeros(3)
param[1]=.65*maximum(prob_rho_sum)
param[2]=51.0
param[3]=0.67
result = optimize(chi_normal,param)
#println(chi_normal(param))

param= result.minimizer
#plot(rho_grid,param2[1]*exp(-(rho_grid-param2[2]).^2*.5/param2[3]^2)+param2[4]*exp(-(rho_grid-param2[5]).^2*.5/param2[6]^2))
plot(rho_grid,param[1]*exp(-(rho_grid-param[2]).^2*.5/param[3]^2),linewidth=3,alpha=0.5)
prob_rho_gillon = zeros(1000)
for i=1:1000
  if rho_grid[i] > 50.7
#    prob_rho_gillon[i] = (param2[1]+param2[4])*exp(-0.5*(rho_grid[i]-50.7)^2/1.2^2)
    prob_rho_gillon[i] = param[1]*exp(-0.5*(rho_grid[i]-50.7)^2/1.2^2)
  else
#    prob_rho_gillon[i] = (param2[1]+param2[4])*exp(-0.5*(rho_grid[i]-50.7)^2/2.2^2)
    prob_rho_gillon[i] = param[1]*exp(-0.5*(rho_grid[i]-50.7)^2/2.2^2)
  end
end
plot(rho_grid,prob_rho_gillon,linewidth=3,alpha=0.5)
println("Density: ",param[2],"+-",param[3])
tight_layout()
savefig("rho_star_circular.pdf", bbox_inches="tight")

#read(STDIN,Char)
#clf()
#semilogx(sig_theta,prob_theta_sum,alpha=0.5)
#scatter(sig_theta,prob_theta_sum,linewidth=3,alpha=0.5)
#plot(sig_theta,cum_prob/maximum(cum_prob),linewidth=3,alpha=0.5)

#axis([6e-3,5,0,maximum(prob_theta_sum)*1.05])
#xlabel(L"\sigma_\theta [deg]")
#ylabel("Relative probability")
#
#function quartic(t,param)
#return exp(-(param[1]+param[2]*t.^2+param[3]*t.^4))
#end
#
#p = prob_theta_sum
#param = [1.,1.,1.]
#
#function chisquare(param)
#return sum((quartic(sig_theta,param)-p).^2)
#end
#
#
#result = optimize(chisquare,param)
#
#param = result.minimizer
#
#plot(sig_theta,quartic(sig_theta,param))
