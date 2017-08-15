# Monte carlo transit probabilities as a function of viewer angle and vector distribution.
# Add in the sizes of the planets [x]
#using PyPlot
using CGS
using JLD

function monte_carlo_circular(nsim::Int)

#x=0.0; y=0.0; length=2.0
length=2.0
#function rand2()
# Generates a random unit vector
#while length > 1.0
#  x=2.0*rand()-1.0
#  y=2.0*rand()-1.0
#  length = x*x+y*y
#end
#length= sqrt(length)
#x /=length
#y /=length
#return x/length,y/length
#end

nsig_theta = 20
sig_theta1 = 0.0001
sig_theta2 = 0.05
#nsig_theta = 10
#sig_theta1 = 0.00003
#sig_theta2 = 0.0600
#rho_star1 = 50.7
#rho_star2 = 50.7
rho_star1 = 48.0
rho_star2 = 54.0
#rho_star1 = 38.0
#rho_star2 = 43.0
nrho = 30
# eccentricity:
nsig_ecc = 1
sig_theta= zeros(nsig_theta)
sig_ecc= zeros(nsig_ecc)
rho_star = zeros(nrho) # In units of the Solar density
Grho_sun = GRAV*MSUN/(4pi/3*RSUN^3)
println("G*rho_sun: ",Grho_sun)
prob = zeros(nsig_theta,nrho,nsig_ecc)
prob_rho = zeros(nrho)
prob_sigma = zeros(nsig_theta)
prob_ecc= zeros(nsig_ecc)
#nsim = 100000000
#nsim = 400000000
nplanet = 7
#imu = [89.65, 89.670, 89.750, 89.86, 89.680, 89.710, 89.80] * pi / 180
#isig = [0.245, 0.170, 0.160, 0.110, 0.034, 0.025, 0.075] * pi / 180
#imin = 1./ [20.50, 28.08, 39.55, 51.97, 68.4, 83.2, 109.]
period = [1.51087081,2.4218233,4.049610,6.099615,9.206690,12.35294,18.767]
dur_obs =[0.025284540348065355, 0.029424593451643284, 0.03408749899000376, 0.03972169981981801, 0.04348163437443351, 0.04754588636760185, 0.052297371090347336]
sd_obs = [0.00012137438458874883, 0.00014791019776496287, 0.0004472829006062137, 0.000513238038271098, 0.0003985200459351564, 0.0004535618470311693, 0.0014067902207811922]
rr_obs = [0.08524846380071895,0.08285759783610132, 0.06054351753969651, 0.07199543643724797, 0.08204346713215337, 0.08841765826354074, 0.059324407529117607]
sr_obs = [0.0004978176874500159,0.0006002276533876536,0.0014080585401894483,0.0017882099896087914,0.001324749805415716,0.0014992826946533435,0.0024903810887425457]

aonr_obs = [20.5,28.08,39.55,51.97,68.4,83.2,114.]
b_obs = [.126,.161,.17,.12,.382,.421,.45]

for ip=1:nplanet
  println("Duration of planet ",ip," ",dur_obs[ip]*24.*60.,"+-",sd_obs[ip]*24.*60.)
  println("Inferred density: ",(period[ip]*24.*3600.)/(dur_obs[ip]*24.*3600.)^3/pi^2*3.*((1.+rr_obs[ip])^2-b_obs[ip]^2)^1.5/Grho_sun)
  println("Inferred density: ",(period[ip]/365.25)*(13./dur_obs[ip]/24.)^3*((1.+rr_obs[ip])^2-b_obs[ip]^2)^1.5)
  println("Inferred density: ",aonr_obs[ip]^3/(period[ip]*24.*3600.)^2*3pi/Grho_sun)
end
#read(STDIN,Char)

# Fraction of simulations that have non-zero probability for each sigma_theta:
frac = zeros(Int,nsig_theta,nrho,nsig_ecc)
lnlike_max = -Inf
dur_max = zeros(nplanet)
dur = zeros(nplanet)
b_max = zeros(nplanet)
b = zeros(nplanet)
time_tot = 0.0
aonr = zeros(nplanet)
imin = zeros(nplanet)
iobs = 0.0; sinphiobs = 0.0; cosphiobs=0.0; lnlike=0.0; ip=0 ; theta = 0.0;  rr=0.0; phi=0.0; sinphi=0.0; cosphi=0.0; cosinc=0.0; b0=0.0; dur0=0.0; sinom = 0.0; ecc=0.0
thetap=zeros(nplanet)
imax = 0.0; tmp=0.0
for i=1:nsig_theta
  tic()
  sig_theta[i] = sig_theta1*10.^(log10(sig_theta2/sig_theta1)*(float(i)-0.5)/nsig_theta)
  for j=1:nrho
  rho_star[j] = rho_star1+(rho_star2-rho_star1)*(float(j)-0.5)/nrho
  # Compute semi-major axes of the planets in units of stellar radius:
#  println("aonr: ",aonr)
  fac = (24.*3600.).^(2./3.)*(Grho_sun*rho_star[j]/3./pi)^(1./3.)
  for ip=1:nplanet
    aonr[ip] = fac*period[ip]^(2./3.)
    imin[ip] = 1./aonr[ip]
  end
  for k=1:nsig_ecc
# Assume zero eccentricity:
  sig_ecc[k] = 0.0
  for isim=1:nsim
  # Draw planets:
  # Randomly choose an observer within ~1 deg of edge-on:
#    iobs = .5*pi + 0.07*(rand()-.5)
    iobs = 0.07*(rand()-.5)
#    iobs = 0.15*(rand()-.5)
#    phiobs = 2pi*rand()
#    sinphiobs=sin(phiobs)
#    cosphiobs=sin(phiobs)
#    sinphiobs = randn()
#    cosphiobs = randn()
    length=2.0
#    sinphiobs, cosphiobs = rand2()
    while length > 1.0
      cosphiobs=2.0*rand()-1.0
      sinphiobs=2.0*rand()-1.0
      length = cosphiobs*cosphiobs+sinphiobs*sinphiobs
    end
    length= sqrt(length)
    cosphiobs /=length
    sinphiobs /=length
#    norm = sqrt(sinphiobs^2+cosphiobs^2)
#    sinphiobs /=norm
#    cosphiobs /=norm
#    uobs=randn(2)
#    uobs = uobs./norm(uobs)
  # Now, compute angles:
    lnlike = 0.0
#    ip = 1
#    while lnlike > -Inf && ip <= nplanet
    ip = 7
    while lnlike > -Inf && ip > 0
      theta = randn()*sig_theta[i]
      length=2.0
#      sinphi, cosphi = rand2()
      while length > 1.0
        cosphi=2.0*rand()-1.0
        sinphi=2.0*rand()-1.0
        length = cosphi*cosphi+sinphi*sinphi
      end
      length= sqrt(length)
      cosphi /=length
      sinphi /=length
      sinom=sin(rand()*2pi)
      tmp = cosphiobs*cosphi+sinphiobs*sinphi
      cosinc = .5*abs(sin(theta-iobs)*(tmp+1.)+(tmp-1.)*sin(theta+iobs))
      b0 = aonr[ip]*cosinc
      if abs(b0) > 1.0
        # Discard any simulation with |b| > 1
        lnlike -= Inf
      else
        # Draw from planet radius ratio distribution:
        rr = -1.0
        while rr < 0
          rr = rr_obs[ip]+randn()*sr_obs[ip]
        end
        # Compute impact parameter:
        # Compute the duration:
        dur0 = period[ip]/pi/aonr[ip]*sqrt(((1.+rr)^2-b0^2))
        # Compare to the observed duration:
        lnlike -= .5*((dur0-dur_obs[ip])/sd_obs[ip])^2
        b[ip]=b0
        dur[ip]=dur0
        thetap[ip]=theta
      end
#      ip +=1
      ip -=1
    end
    if lnlike != -Inf
      frac[i,j,k] += 1
      if abs(iobs) > imax && lnlike > log(1e-5)
        imax = abs(iobs)
#        println("imax: ",imax*180./pi," lnlike: ",lnlike," b: ",b," theta: ",thetap*180./pi)
      end
    end
    prob[i,j,k] += exp(lnlike)
    if lnlike > lnlike_max
      lnlike_max = lnlike
      dur_max = dur
      b_max = b
#      println("Maximum likelihood: ",lnlike_max)
#      println("Durations: ",dur_max)
#      println("Impact pm: ",b_max)
    end
  end
  end
#  println(j," ",rho_star[j]," ",prob[i,j,k])
  end
  prob_sigma[i] = maximum(prob[i,:,:])
  time_tot += toc()
#  println(i," ",sig_theta[i]," ",prob_sigma[i])
  println("Time remaining: ",(nsig_theta-i)*time_tot/i)
end
for j=1:nrho
  prob_rho[j] = maximum(prob[:,j,:])
end
#semilogx(sig_theta*180./pi,exp(-abs(sig_theta.*180./pi/.15)*.5)*2750.)
#fig,axes = subplots(1,2)
#ax = axes[1]
#ax[:semilogx](sig_theta*180./pi,prob_sigma)
#ax[:scatter](sig_theta*180/pi,prob_sigma)
#ax[:axis]([2e-3,3,0,maximum(prob_sigma)*1.05])
#ax[:set_xlabel](L"\sigma_\theta [deg]")
#ax[:set_ylabel]("Relative probability")
#ax = axes[2]
#ax[:plot](rho_star,prob_rho)
#ax[:scatter](rho_star,prob_rho)
#ax[:axis]([rho_star1,rho_star2,0,maximum(prob_rho)*1.05])
#ax[:set_xlabel](L"\rho_* [\rho_\odot]")
#ax[:set_ylabel]("Relative probability")
#
#
save(string("MC_circular_",string(nsim),"sig_rho_sigtheta.jld"),"nsim",nsim,"nsig_theta",nsig_theta,"nrho",nrho,"sig_theta",sig_theta,"rho_star",rho_star,"prob_sigma",prob_sigma,"prob_rho",prob_rho,"prob",prob,"frac",frac,"lnlike_max",lnlike_max,"dur_max",dur_max,"b_max",b_max,"imax",imax)
return sig_theta,rho_star,prob
end
